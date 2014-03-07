/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2013 Lawrence Livermore National Security, LLC
 * Description:   Scalable load balancer using tree algorithm.
 *
 ************************************************************************/

#ifndef included_mesh_CascadePartitioner_C
#define included_mesh_CascadePartitioner_C

#include "SAMRAI/mesh/CascadePartitioner.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace mesh {

/*!
 * @brief Make a cycle-zero, single-process group.
 */
void CascadePartitionerGroup::makeSingleProcessGroup(
   tbox::SAMRAI_MPI &mpi,
   const PartitioningParams &pparams,
   TransitLoad &local_load,
   double global_load )
{
   d_mpi = mpi;
   d_cycle_num = 0;
   d_first_lower_rank = d_mpi.getRank();
   d_first_upper_rank = d_first_lower_rank+1;
   d_end_rank = d_first_lower_rank+1;
   d_contact = -1;
   d_our_half = 0;
   d_our_position = Lower;
   d_lower_weight = local_load.getSumLoad();
   d_upper_weight = 0.0;
   d_our_weight = &d_lower_weight;
   d_far_weight = 0;
   d_local_load = &local_load;
   d_global_load_avg = global_load/d_mpi.getSize();
   d_pparams = &pparams;
}

/*!
 * @brief Make a combo group consistng of the given half-group and
 * its corresponding partner, which this method will determine.
 *
 * Local process provides the half containing the local rank
 * (our_half).  Data for other half (the far half) is obtained by
 * communication.
 */
void CascadePartitionerGroup::makeComboGroup( CascadePartitionerGroup &our_half )
{
   /*
    * Set up the groups
    */

   d_our_half = &our_half;
   d_cycle_num = our_half.d_cycle_num+1;
   d_mpi = our_half.d_mpi;
   d_local_load = 0; // For single-process groups only.

   int group_size = 1 << d_cycle_num;
   int group_num = d_mpi.getRank()/group_size;

   d_first_lower_rank = group_size*group_num;
   d_first_upper_rank = tbox::MathUtilities<int>::Min(
      d_first_lower_rank + group_size/2, d_mpi.getSize());
   d_end_rank = tbox::MathUtilities<int>::Min(
      d_first_lower_rank + group_size, d_mpi.getSize());

   int relative_rank = d_mpi.getRank() - d_first_lower_rank;

   d_our_position = relative_rank >= group_size/2 ? Upper : Lower;
   d_contact = d_our_position == Lower ?
      d_mpi.getRank() + group_size/2 :
      d_mpi.getRank() - group_size/2;
   if ( d_contact >= d_mpi.getSize() ) {
      d_contact = -1;
   }
   d_our_weight = d_our_position == Lower ? &d_lower_weight : &d_upper_weight;
   d_far_weight = d_our_position == Lower ? &d_upper_weight : &d_lower_weight;

   /*
    * Record weights of the two halves.  Needs communication for far
    * half.
    */
   double our_weight = d_our_half->getComboWeight();
   double far_weight = 0.0;
   if ( d_contact >= 0 ) {
      tbox::MessageStream send_msg;
      send_msg << our_weight << d_our_half_may_supply;

      std::vector<char> tmp_buffer(send_msg.getCurrentSize());

      tbox::SAMRAI_MPI::Status status;
      d_mpi.Sendrecv( (void*)(send_msg.getBufferStart()), 1, MPI_CHAR, d_contact, 0,
                      &tmp_buffer[0], 1, MPI_CHAR, d_contact, 1,
                      &status );

      tbox::MessageStream recv_msg( tmp_buffer.size(),
                                    tbox::MessageStream::Read,
                                    &tmp_buffer[0], false );
      recv_msg >> far_weight >> d_contact_may_supply;
   }
   d_lower_weight = d_our_position == Lower ? our_weight : far_weight;
   d_upper_weight = d_our_position == Upper ? our_weight : far_weight;

   d_global_load_avg = our_half.d_global_load_avg;
   d_pparams = our_half.d_pparams;
}



/*
 *************************************************************************
 * If one half has a positive surplus and the other half has a
 * negative surplus, the former supplies work to the latter.  Amount
 * supplied is ideally the minimum of the supplier's surplus and the
 * demander's deficit.  (Actual ammounts are limited by load cutting
 * restrictions.)
 *
 * This method records estimates of the weight changes to the groups
 * it knows about, It doesn't record the actual weight changes because
 * that happens remotely.  Each process on the supply half send a
 * message received by its contact on the demand half.  The messages
 * has the actual work to be transfered.
 *************************************************************************
 */
void
CascadePartitionerGroup::balanceConstituentHalves()
{

   if ( ourSurplus() >  d_pparams->getLoadComparisonTol() &&
        farSurplus() < -d_pparams->getLoadComparisonTol() ) {

      double supplied_weight = supplyLoadFromOurHalf(
         -farSurplus(), d_contact );

      recordDemandReceivedByFarHalf(supplied_weight);
   }

   else if ( farSurplus() >  d_pparams->getLoadComparisonTol() &&
             ourSurplus() < -d_pparams->getLoadComparisonTol() ) {

      if ( d_contact_may_supply ) {
         d_comm.setPeerRank(d_contact);
         d_comm.beginRecv();
      }

      double supplied_weight = supplyLoadFromFarHalf(
         -ourSurplus() );

      recordDemandReceivedByOurHalf(supplied_weight);

      unpackSuppliedLoad();
   }

   return;
}



/*
 *************************************************************************
 * Supply specified amount of load to another group with a work
 * demand.  Any load supplied by local process is to be sent to the
 * designated taker, a process in the demand group.  Give priority
 * to supplies closest to the taker in rank space.
 *
 * 1. If group is single-process, try to remove load from d_local_load
 *    and put it in d_shipment.
 * 2. Else:
 *    a: Supply load from the half closest to the taker.
 *    b: If step a didn't supply enough, remove some from the other half.
 *
 * This method is recursive by the sequence
 * supplyLoad-supplyLoadFromOurHalf-supplyLoad.
 *************************************************************************
 */
double
CascadePartitionerGroup::supplyLoad( double amount, int taker )
{
   TBOX_ASSERT( amount > 0.0 );
   double removed = 0.0;

   if ( d_cycle_num == 0 ) {
      // Group is single-process, not two halves.
      d_shipment = boost::shared_ptr<TransitLoad>(d_local_load->clone());
      d_shipment->adjustLoad( *d_local_load, amount, amount, amount );
      removed = d_shipment->getSumLoad();
      d_lower_weight -= removed;
      sendMyShipment(taker);
   }

   /*
    * Group contains two halves.  First, supply from the half closer
    * to taker.  If more is needed, supply from other half.
    */
   else if ( ( d_our_position == Lower && taker <  d_first_lower_rank ) ||
             ( d_our_position == Upper && taker >= d_first_upper_rank ) ) {
      removed = supplyLoadFromOurHalf( amount, taker );
      if ( removed < amount ) {
         removed += supplyLoadFromFarHalf( amount-removed );
      }
   }
   else {
      removed = supplyLoadFromFarHalf( amount );
      if ( removed < amount ) {
         removed += supplyLoadFromOurHalf( amount-removed, taker );
      }
   }
   return removed;
}



/*
 *************************************************************************
 *************************************************************************
 */
void
CascadePartitionerGroup::sendMyShipment( int taker )
{
   tbox::MessageStream msg;
   msg << *d_shipment;
   d_comm.setPeerRank(taker);
   d_comm.beginSend( (const char*)(msg.getBufferStart()), msg.getCurrentSize() );
   return;
}



/*
 *************************************************************************
 *************************************************************************
 */
void
CascadePartitionerGroup::unpackSuppliedLoad()
{
   d_comm.completeCurrentOperation();
   tbox::MessageStream recv_msg( d_comm.getRecvSize(), tbox::MessageStream::Read,
                                 d_comm.getRecvData(), true );
   recv_msg >> d_local_load;
   return;
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif

#endif
