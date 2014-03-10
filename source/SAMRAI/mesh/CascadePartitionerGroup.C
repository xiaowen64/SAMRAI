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
   const CascadePartitioner *common_data,
   TransitLoad &local_load )
{
   d_common = common_data;
   d_cycle_num = 0;
   d_first_lower_rank = d_common->d_mpi.getRank();
   d_first_upper_rank = d_first_lower_rank+1;
   d_end_rank = d_first_lower_rank+1;
   d_contact = -1;
   d_our_half = 0;
   d_our_position = Lower;
   d_lower_work = local_load.getSumLoad();
   d_upper_work = 0.0;
   d_our_work = &d_lower_work;
   d_far_work = 0;
   d_local_load = &local_load;
   d_lower_capacity = d_common->d_global_load_avg;
   d_upper_capacity = 0.0;
   if ( d_common->d_print_steps ) {
      tbox::plog << "CascadePartitionerGroup::makeSingleProcessGroup:\n";
      printClassData( tbox::plog, "\t" );
   }
}

/*!
 * @brief Make a combined group consistng of the given half-group and
 * its corresponding partner, which this method will determine.
 *
 * Local process provides the half containing the local rank
 * (our_half).  Data for other half (the far half) is obtained by
 * communication.
 */
void CascadePartitionerGroup::makeCombinedGroup( CascadePartitionerGroup &our_half )
{
   /*
    * Set up the groups
    */

   d_common = our_half.d_common;
   d_our_half = &our_half;
   d_cycle_num = our_half.d_cycle_num+1;
   d_local_load = 0; // For single-process groups only.

   int group_size = 1 << d_cycle_num;
   int group_num = d_common->d_mpi.getRank()/group_size;

   d_first_lower_rank = group_size*group_num;
   d_first_upper_rank = tbox::MathUtilities<int>::Min(
      d_first_lower_rank + group_size/2, d_common->d_mpi.getSize());
   d_end_rank = tbox::MathUtilities<int>::Min(
      d_first_lower_rank + group_size, d_common->d_mpi.getSize());

   int relative_rank = d_common->d_mpi.getRank() - d_first_lower_rank;

   d_our_position = relative_rank >= group_size/2 ? Upper : Lower;
   d_contact = d_our_position == Lower ?
      d_common->d_mpi.getRank() + group_size/2 :
      d_common->d_mpi.getRank() - group_size/2;
   if ( d_contact >= d_common->d_mpi.getSize() ) {
      d_contact = -1;
   }
   d_our_work = d_our_position == Lower ? &d_lower_work : &d_upper_work;
   d_far_work = d_our_position == Lower ? &d_upper_work : &d_lower_work;

   /*
    * Record works of the two halves.  Needs communication to get
    * data about the far half of the group.
    */
   double our_work = d_our_half->getCombinedWork();
   double far_work = 0.0;
   if ( d_contact >= 0 ) {
      tbox::MessageStream send_msg;
      send_msg << our_work << d_our_half_may_supply;

      std::vector<char> tmp_buffer(send_msg.getCurrentSize());

      tbox::SAMRAI_MPI::Status status;
      d_common->d_mpi.Sendrecv( (void*)(send_msg.getBufferStart()), 1, MPI_CHAR, d_contact, 0,
                                &tmp_buffer[0], 1, MPI_CHAR, d_contact, 1,
                                &status );

      tbox::MessageStream recv_msg( tmp_buffer.size(),
                                    tbox::MessageStream::Read,
                                    &tmp_buffer[0], false );
      recv_msg >> far_work >> d_contact_may_supply;
   }
   d_lower_work = d_our_position == Lower ? our_work : far_work;
   d_upper_work = d_our_position == Upper ? our_work : far_work;

   d_lower_capacity = d_common->d_global_load_avg*(d_first_upper_rank-d_first_lower_rank);
   d_upper_capacity = d_common->d_global_load_avg*(d_end_rank-d_first_upper_rank);

   if ( d_common->d_print_steps ) {
      tbox::plog << "CascadePartitionerGroup::makeCombinedGroup:\n";
      printClassData( tbox::plog, "\t" );
   }
}



/*
 *************************************************************************
 * If one half has a positive surplus and the other half has a
 * negative surplus, the former supplies work to the latter.  Amount
 * supplied is ideally the minimum of the supplier's surplus and the
 * requestor's deficit.  (Actual ammounts are affected by load cutting
 * restrictions.)
 *
 * This method records estimates of the work changes to the groups
 * it knows about, It doesn't record the actual work changes because
 * that happens remotely.  Each process on the supply half send a
 * message received by its contact on the requesting half.  The messages
 * has the actual work to be transfered.
 *************************************************************************
 */
void
CascadePartitionerGroup::balanceConstituentHalves()
{

   if ( ourSurplus() >  d_common->d_pparams->getLoadComparisonTol() &&
        farSurplus() < -d_common->d_pparams->getLoadComparisonTol() ) {

      double work_supplied = supplyWorkFromOurHalf(
         -farSurplus(), d_contact );

      recordWorkTakenByFarHalf(work_supplied);
   }

   else if ( farSurplus() >  d_common->d_pparams->getLoadComparisonTol() &&
             ourSurplus() < -d_common->d_pparams->getLoadComparisonTol() ) {

      if ( d_contact_may_supply ) {
         d_comm.setPeerRank(d_contact);
         d_comm.beginRecv();
      }

      double work_supplied = supplyWorkFromFarHalf(
         -ourSurplus() );

      recordWorkTakenByOurHalf(work_supplied);

      unpackSuppliedLoad();
   }

   return;
}



/*
 *************************************************************************
 * Supply work_requested of load to another group requesting work.
 * Any load supplied by local process is to be sent to the designated
 * taker, a process in the requesting group.  Give priority to
 * supplies closest to the taker in rank space.
 *
 * 1. If group is single-process, remove load from d_local_load
 *    and put it in d_shipment.
 * 2. Else:
 *    A: Remove load from the half closest to the taker.
 *    B: If step A didn't supply enough, remove some from the other half.
 *
 * This method is recursive by the sequence
 * supplyWork-supplyWorkFromOurHalf-supplyWork.
 *************************************************************************
 */
double
CascadePartitionerGroup::supplyWork( double work_requested, int taker )
{
   TBOX_ASSERT( work_requested > 0.0 );
   TBOX_ASSERT( d_our_half_may_supply );  // Should be checked by calling method.
   double work_supplied = 0.0;

   if ( d_cycle_num == 0 ) {
      // Group is single-process, not two halves.
      d_shipment = boost::shared_ptr<TransitLoad>(d_local_load->clone());
      d_shipment->adjustLoad( *d_local_load, work_requested, work_requested, work_requested );
      work_supplied = d_shipment->getSumLoad();
      d_lower_work -= work_supplied;
      sendMyShipment(taker);
   }

   else if ( ( d_our_position == Lower && taker <  d_first_lower_rank ) ||
             ( d_our_position == Upper && taker >= d_first_upper_rank ) ) {
      work_supplied = supplyWorkFromOurHalf( work_requested, taker );
      if ( work_supplied < work_requested ) {
         work_supplied += supplyWorkFromFarHalf( work_requested-work_supplied );
      }
   }
   else {
      work_supplied = supplyWorkFromFarHalf( work_requested );
      if ( work_supplied < work_requested ) {
         work_supplied += supplyWorkFromOurHalf( work_requested-work_supplied, taker );
      }
   }
   return work_supplied;
}



/*
 *************************************************************************
 * Try to supply an amount of work from our half of the group.
 *
 * This method is recursive by the sequence
 * supplyWorkFromOurHalf-supplyWork-supplyWorkFromOurHalf.
 *************************************************************************
 */
double CascadePartitionerGroup::supplyWorkFromOurHalf( double work_requested, int taker ) {
   TBOX_ASSERT( work_requested > 0.0 );
   double work_supplied = 0.0;
   if ( d_our_half_may_supply &&
        ourSurplus() >= d_common->d_pparams->getLoadComparisonTol() ) {
      work_supplied = d_our_half->supplyWork( work_requested, taker );
      *d_our_work -= work_supplied;
   }
   return work_supplied;
}



/*
 *************************************************************************
 * Symbolically try to supply an amount of work from the
 * half of the group not containing the local process.
 *
 * No real work is exchanged because the local process is not in
 * the far half.  This method just estimates what the far half
 * would give away.
 *************************************************************************
 */
double CascadePartitionerGroup::supplyWorkFromFarHalf( double work_requested ) {
   TBOX_ASSERT( work_requested > 0.0 );
   double work_supplied = 0.0;
   if ( d_far_half_may_supply &&
        farSurplus() >= d_common->d_pparams->getLoadComparisonTol() ) {
      work_supplied = tbox::MathUtilities<double>::Min( work_requested, *d_far_work );
      *d_far_work -= work_supplied;
   }
   return work_supplied;
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



/*
 *************************************************************************
 *************************************************************************
 */
void
CascadePartitionerGroup::printClassData( std::ostream &co, const std::string &border ) const
{
   co << border << "cycle " << d_cycle_num
      << "  [" << d_first_lower_rank << ',' << d_first_upper_rank << ',' << d_end_rank
      << ")  group_size=" << d_end_rank-d_first_lower_rank << '='
      << d_first_upper_rank-d_first_lower_rank << '+' << d_end_rank-d_first_upper_rank << "\n"
      << border << " lower_work=" << d_lower_work << '/' << d_lower_capacity
      << "  " << " upper_work=" << d_upper_work << '/' << d_upper_capacity << '\n'
      << border << " our_position=" << d_our_position
      << "  our_half_may_supply: " << d_our_half_may_supply
      << "  far_half_may_supply: " << d_far_half_may_supply
      << "  contact=" << d_contact << "  contact_may_supply: " << d_contact_may_supply
      << '\n';
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
