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
   d_lower_begin = d_common->d_mpi.getRank();
   d_upper_begin = d_lower_begin+1;
   d_upper_end = d_lower_begin+1;
   d_contact = -1;
   d_our_half = 0;
   d_our_position = Lower;

   d_lower_work = local_load.getSumLoad();
   d_upper_work = 0.0;
   d_our_work = &d_lower_work;
   d_far_work = 0;

   d_lower_capacity = d_common->d_global_load_avg;
   d_upper_capacity = 0.0;

   d_our_half_may_supply = d_local_may_supply =
      d_lower_work > ( d_lower_capacity + d_common->d_pparams->getLoadComparisonTol() );

   d_far_half_may_supply = d_contact_may_supply = false; // No far half.

   if ( d_common->d_print_steps ) {
      tbox::plog << "CascadePartitionerGroup::makeSingleProcessGroup: leaving\n";
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
   if ( our_half.d_common->d_print_steps ) {
      tbox::plog << "CascadePartitionerGroup::makeCombinedGroup: entered\n";
   }

   /*
    * Set up the groups
    */

   d_common = our_half.d_common;
   d_our_half = &our_half;
   d_cycle_num = our_half.d_cycle_num+1;

   int group_size = 1 << d_cycle_num;
   int group_num = d_common->d_mpi.getRank()/group_size;

   d_lower_begin = group_size*group_num;
   d_upper_begin = tbox::MathUtilities<int>::Min(
      d_lower_begin + group_size/2, d_common->d_mpi.getSize());
   d_upper_end = tbox::MathUtilities<int>::Min(
      d_lower_begin + group_size, d_common->d_mpi.getSize());

   d_lower_capacity = d_common->d_global_load_avg*(d_upper_begin-d_lower_begin);
   d_upper_capacity = d_common->d_global_load_avg*(d_upper_end-d_upper_begin);

   int relative_rank = d_common->d_mpi.getRank() - d_lower_begin;

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
    * Determine and record work of the two halves.  Needs
    * communication to get data about the far half of the group.
    */

   *d_our_work = d_our_half->getCombinedWork();
   d_our_half_may_supply = *d_our_work > d_common->d_pparams->getLoadComparisonTol() +
      ( d_our_position == Lower ? d_lower_capacity : d_upper_capacity ) ;
   d_local_may_supply = d_our_half->d_local_may_supply && d_our_half_may_supply;

   double our_work = d_our_half->getCombinedWork();
   double far_work = 0.0;
   if ( d_contact >= 0 ) {
      tbox::MessageStream send_msg;
      send_msg << our_work << d_our_half_may_supply << d_local_may_supply;

      std::vector<char> tmp_buffer(send_msg.getCurrentSize());

      tbox::SAMRAI_MPI::Status status;
      d_common->d_mpi.Sendrecv(
         (void*)(send_msg.getBufferStart()),
         send_msg.getCurrentSize(),
         MPI_CHAR,
         d_contact,
         CascadePartitionerGroup_TAG_InfoExchange,
         &tmp_buffer[0],
         send_msg.getCurrentSize(),
         MPI_CHAR,
         d_contact,
         CascadePartitionerGroup_TAG_InfoExchange,
         &status );

      tbox::MessageStream recv_msg( tmp_buffer.size(),
                                    tbox::MessageStream::Read,
                                    &tmp_buffer[0], false );
      recv_msg >> far_work >> d_far_half_may_supply >> d_contact_may_supply;
   }

   d_lower_work = d_our_position == Lower ? our_work : far_work;
   d_upper_work = d_our_position == Upper ? our_work : far_work;

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

      double work_supplied = supplyWorkFromOurHalf( -farSurplus(), d_contact );

      // Record work taken by the far half.
      *d_far_work += work_supplied;
      d_far_half_may_supply = d_contact_may_supply = false;

      if ( d_local_may_supply ) {
         sendMyShipment(d_contact);
      }

      if ( d_common->d_print_steps ) {
         tbox::plog << "CascadePartitionerGroup::balanceConstituentHalves:"
                    << "  supplied " << work_supplied
                    << " from our half to far half.  d_local_may_supply="
                    << d_local_may_supply << "\n";
      }
   }

   else if ( farSurplus() >  d_common->d_pparams->getLoadComparisonTol() &&
             ourSurplus() < -d_common->d_pparams->getLoadComparisonTol() ) {

      if ( d_contact_may_supply ) {
         d_common->d_comm_peer.setPeerRank(d_contact);
         d_common->d_comm_peer.setMPITag( CascadePartitionerGroup_TAG_LoadTransfer0,
                                          CascadePartitionerGroup_TAG_LoadTransfer1 );
         d_common->d_comm_peer.beginRecv();
      }

      double work_supplied = supplyWorkFromFarHalf( -ourSurplus() );

      // Record work taken by our half.
      *d_our_work += work_supplied;
      d_our_half_may_supply = d_local_may_supply = false;

      if ( d_common->d_print_steps ) {
         tbox::plog << "CascadePartitionerGroup::balanceConstituentHalves:"
                    << "  recorded supply of " << work_supplied
                    << " from far half to our half.  d_contact_may_supply="
                    << d_contact_may_supply << "\n";
      }

      if ( d_contact_may_supply ) {
         receiveAndUnpackSuppliedLoad();
      }
   }

   // Complete the load send, if there was any.
   d_common->d_comm_stage.advanceAll();
   while ( d_common->d_comm_stage.numberOfCompletedMembers() > 0 ) {
      d_common->d_comm_stage.popCompletionQueue();
   }

   if ( d_common->d_print_steps ) {
      tbox::plog << "CascadePartitionerGroup::balanceConstituentHalves: leaving with state:\n";
      printClassData( tbox::plog, "\t" );
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
 * 1. If group is single-process, remove load from d_common->d_local_load
 *    and put it in d_common->d_shipment.
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
   double work_supplied = 0.0;

   if ( d_cycle_num == 0 ) {
      if ( d_local_may_supply ) {
         /*
          * Estimate the work that can be supplied.  Even though we can
          * compute the work_suppliled exactly, we must use the estimate
          * based on perfect work division, not actual work division, so
          * that our record matches the record remote processes keep on
          * us.
          */
         work_supplied = tbox::MathUtilities<double>::Min(
            work_requested, d_common->d_local_load->getSumLoad() );
         d_lower_work -= work_supplied;

         d_common->d_shipment->adjustLoad(
            *d_common->d_local_load,
            work_requested,
            work_requested,
            work_requested );
         if ( d_common->d_print_steps ) {
            tbox::plog << "CascadePartitionerGroup::supplyWork giving to " << taker << ": ";
            d_common->d_shipment->recursivePrint();
            tbox::plog << "CascadePartitionerGroup::supplyWork keeping: ";
            d_common->d_local_load->recursivePrint();
         }

         d_local_may_supply =
            d_common->d_local_load->getSumLoad() > ( d_lower_capacity + d_common->d_pparams->getLoadComparisonTol() );
         d_our_half_may_supply = ourSurplus() > d_common->d_pparams->getLoadComparisonTol();
      }
   }

   else if ( ( d_our_position == Lower && taker <  d_upper_begin ) ||
             ( d_our_position == Upper && taker >= d_upper_begin ) ) {
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
 * Supply an amount of work from our half of the group.
 *
 * This method is recursive by the sequence
 * supplyWorkFromOurHalf-supplyWork-supplyWorkFromOurHalf.
 *
 * Note: Don't assume this work will be added to the far half.
 * It may be combined with other works and added to a bigger group.
 *************************************************************************
 */
double CascadePartitionerGroup::supplyWorkFromOurHalf( double work_requested, int taker ) {
   TBOX_ASSERT( work_requested > 0.0 );
   double work_supplied = 0.0;
   if ( d_our_half_may_supply ) {
      work_supplied = d_our_half->supplyWork( work_requested, taker );
      *d_our_work -= work_supplied;
      d_our_half_may_supply = ourSurplus() > d_common->d_pparams->getLoadComparisonTol();
   }
   return work_supplied;
}



/*
 *************************************************************************
 * Symbolically supply an amount of work from the far half of the
 * group.
 *
 * No real work is exchanged because the local process is not in
 * the far half.  This method just estimates what the far half
 * would give away.
 *
 * Note: Don't assume this work will be added to our half.
 * It may be combined with other works and added to a bigger group.
 *************************************************************************
 */
double CascadePartitionerGroup::supplyWorkFromFarHalf( double work_requested ) {
   TBOX_ASSERT( work_requested > 0.0 );
   double work_supplied = 0.0;
   if ( d_far_half_may_supply ) {
      work_supplied = tbox::MathUtilities<double>::Min( work_requested, *d_far_work );
      *d_far_work -= work_supplied;
      d_far_half_may_supply = farSurplus() > d_common->d_pparams->getLoadComparisonTol();
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
   if ( d_common->d_print_steps ) {
      tbox::plog << "CascadePartitionerGroup::sendMyShipment: sending to " << taker << ":\n";
      d_common->d_shipment->recursivePrint(tbox::plog, "Send: ", 2);
   }
   tbox::MessageStream msg;
   msg << *d_common->d_shipment;
   d_common->d_comm_peer.setPeerRank(taker);
   d_common->d_comm_peer.setMPITag( CascadePartitionerGroup_TAG_LoadTransfer0,
                                    CascadePartitionerGroup_TAG_LoadTransfer1 );
   d_common->d_comm_peer.beginSend( (const char*)(msg.getBufferStart()),
                                    msg.getCurrentSize() );
   d_common->d_shipment->clear();
   return;
}



/*
 *************************************************************************
 *************************************************************************
 */
void
CascadePartitionerGroup::receiveAndUnpackSuppliedLoad()
{
   d_common->d_comm_peer.completeCurrentOperation();
   tbox::MessageStream recv_msg( d_common->d_comm_peer.getRecvSize(),
                                 tbox::MessageStream::Read,
                                 d_common->d_comm_peer.getRecvData(),
                                 true );
   recv_msg >> *d_common->d_local_load;
   if ( d_common->d_print_steps ) {
      tbox::plog << "CascadePartitionerGroup::receiveAndUnpackSuppliedLoad: updated d_local_load:\n";
     d_common->d_local_load->recursivePrint(tbox::plog, "Curr: ", 2);
   }
   return;
}



/*
 *************************************************************************
 *************************************************************************
 */
void
CascadePartitionerGroup::printClassData( std::ostream &co, const std::string &border ) const
{
   const std::string indent( border + std::string(d_cycle_num,' ') + std::string(d_cycle_num,' ') );
   co << indent << "cycle " << d_cycle_num
      << "  [" << d_lower_begin << ',' << d_upper_begin << ',' << d_upper_end
      << ")  group_size=" << d_upper_end-d_lower_begin << '='
      << d_upper_begin-d_lower_begin << '+' << d_upper_end-d_upper_begin
      << "  this=" << this << "  our_half=" << d_our_half
      << '\n' << indent
      << "our_position=" << d_our_position << "  contact=" << d_contact
      << "  lower_work=" << d_lower_work << '/' << d_lower_capacity
      << "  " << " upper_work=" << d_upper_work << '/' << d_upper_capacity
      << '\n' << indent
      << "our_half_may_supply=" << d_our_half_may_supply
      << "  local_may_supply=" << d_local_may_supply
      << "  far_half_may_supply=" << d_far_half_may_supply
      << "  contact_may_supply=" << d_contact_may_supply
      << '\n' << indent
      << "local_load:";
   d_common->d_local_load->recursivePrint(tbox::plog, indent,0);
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
