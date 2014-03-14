/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2013 Lawrence Livermore National Security, LLC
 * Description:   Scalable load balancer using tree algorithm.
 *
 ************************************************************************/

#ifndef included_mesh_CascadePartitionerTree_C
#define included_mesh_CascadePartitionerTree_C

#include "SAMRAI/mesh/CascadePartitionerTree.h"
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

/*
 * Construct the tree for the given CascadePartitioner by constructing
 * the root and recursively constructing its children.
 *
 * Only children relevant to local process is actually constructed.
 * These are the groups containing the local processes and their
 * sibling groups.
 */
CascadePartitionerTree::CascadePartitionerTree(
   const CascadePartitioner &partitioner ) :
   d_common(&partitioner),
   d_generation(0),

   d_begin(0),
   d_end(partitioner.d_mpi.getSize()),
   d_position(-1),

   d_parent(0),
   d_near(0),
   d_far(0),
   d_leaf(0),

   d_work(0),
   d_capacity(d_common->d_global_load_avg*d_common->d_mpi.getSize()),
   d_group_may_supply(false)
{
   d_children[0] = d_children[1] = 0;
   d_contact[0] = d_contact[1] = -1;
   d_process_may_supply[0] = d_process_may_supply[1] = false;

   makeChildren();

   if ( d_common->d_print_steps ) {
      tbox::plog << "CascadePartitionerTree::makeSingleProcessGroup: leaving\n";
      printClassData( tbox::plog, "\t" );
   }
}



/*
 * Construct a child group.
 */
CascadePartitionerTree::CascadePartitionerTree(
   CascadePartitionerTree &parent,
   Position my_position ) :
   d_common(parent.d_common),
   d_generation(1 + parent.d_generation),

   d_begin(parent.d_begin),
   d_end(parent.d_end),
   d_position(my_position),

   d_parent(&parent),
   d_near(0),
   d_far(0),
   d_leaf(0),

   d_work(0.0),
   d_capacity(0),
   d_group_may_supply(false)
{
   d_children[0] = d_children[1] = 0;
   d_contact[0] = d_contact[1] = -1;
   d_process_may_supply[0] = d_process_may_supply[1] = false;

   const int upper_begin = (d_parent->d_begin + d_parent->d_end)/2;
   if ( my_position == Lower ) { d_end = upper_begin; }
   else { d_begin = upper_begin; }

   d_capacity = d_common->d_global_load_avg*(d_end-d_begin);

   /*
    * Assign contacts by pairing this process with the sibling group's
    * process having the same relative rank.  If upper group has an
    * odd rank, pair it with the last rank in lower group (giving last
    * rank in lower group two contacts).
    */
   const int relative_rank = d_common->d_mpi.getRank() - d_begin;
   d_contact[1] = -1;
   if ( my_position == Lower ) {
      d_contact[0] = relative_rank + upper_begin;
      if ( d_common->d_mpi.getRank() == upper_begin+1 &&
           (d_parent->d_end - d_parent->d_begin)%2 ) {
         d_contact[1] = 1 + d_contact[0];
      }
   }
   else {
      d_contact[0] = tbox::MathUtilities<int>::Max(
         relative_rank + d_parent->d_begin, d_begin-1 );
   }

   d_process_may_supply[0] = d_process_may_supply[1] = false;

   // Local process only needs child details of groups containing it.
   if ( d_common->d_mpi.getRank() >= d_begin &&
        d_common->d_mpi.getRank() < d_end ) {
      makeChildren();
   }

   return;
}



/*
 *************************************************************************
 *************************************************************************
 */
void
CascadePartitionerTree::makeChildren()
{
   if ( (d_end-d_begin) > 1 ) {

      d_children[0] = new CascadePartitionerTree( *this, Lower );
      d_children[1] = new CascadePartitionerTree( *this, Upper );

      const bool in_upper = d_common->d_mpi.getRank() >= d_children[1]->d_begin;
      d_near = d_children[in_upper];
      d_far = d_children[!in_upper];

      d_leaf = d_near->d_leaf;
   }

   else {
      // This is a leaf and has no children.
      if ( d_begin == d_common->d_mpi.getRank() ) {
         d_leaf = this;
      }
   }

   return;
}



/*
 *************************************************************************
 *************************************************************************
 */
CascadePartitionerTree::~CascadePartitionerTree()
{
   if ( d_children[0] ) delete d_children[0];
   if ( d_children[1] ) delete d_children[1];
   d_children[0] = d_children[1] = d_near = d_far = d_leaf = 0;
}



/*
 *************************************************************************
 *************************************************************************
 */
void CascadePartitionerTree::balanceAll()
{
   const int lg_size = CascadePartitioner::lgInt(d_common->d_mpi.getSize());

   /*
    * Calling balanceChildren balances current_group's two branches,
    * but within each branch, there may be imbalance.
    * We need lg_size outer cycles, each balancing a generation.
    */
   CascadePartitionerTree *current_group = this;
   for ( int outer_cycle=0; outer_cycle<lg_size; ++outer_cycle ) {
      current_group->balanceChildren();
      current_group = current_group->d_near;
   }

}



/*
 *************************************************************************
 * If one child has a positive surplus and the other has a negative
 * surplus, the former supplies work to the latter.  Amount supplied
 * is ideally the minimum of the supplier's surplus and the
 * requestor's deficit.  (Actual ammounts are affected by load cutting
 * restrictions.)  Note that only children groups will be balanced
 * (not all descendents).
 *
 * This method records estimates of the work changes to the groups it
 * knows about.  It doesn't record the actual work changes because
 * that happens remotely.  Each process in the supply group send a
 * message received by its contact(s) on the requesting group.  The
 * messages has the actual work to be transfered.
 *************************************************************************
 */
void
CascadePartitionerTree::balanceChildren()
{

   if ( d_near->surplus() > d_common->d_pparams->getLoadComparisonTol() &&
        d_far->surplus() < -d_common->d_pparams->getLoadComparisonTol() ) {

      double work_supplied = d_near->supplyWork( -d_far->surplus(), d_contact[0] );

      // Record work taken by the far child.
      d_far->d_work += work_supplied;
      d_far->d_group_may_supply = d_far->d_process_may_supply[0] = d_far->d_process_may_supply[1] = false;

      if ( d_process_may_supply[0] ) {
         sendShipment(d_contact[0]);
      }

      if ( d_common->d_print_steps ) {
         tbox::plog << "CascadePartitionerTree::balanceConstituentHalves:"
                    << "  supplied " << work_supplied
                    << " from our half to far half.  d_process_may_supply="
                    << d_process_may_supply[0] << ',' << d_process_may_supply[1] << "\n";
      }
   }

   else if ( d_far->surplus()  >  d_common->d_pparams->getLoadComparisonTol() &&
             d_near->surplus() < -d_common->d_pparams->getLoadComparisonTol() ) {

      if ( d_far->d_group_may_supply ) {
         if ( d_process_may_supply[0] ) {
            d_common->d_comm_peer[0].setPeerRank(d_contact[0]);
            d_common->d_comm_peer[0].setMPITag( CascadePartitionerTree_TAG_LoadTransfer0,
                                                CascadePartitionerTree_TAG_LoadTransfer1 );
            d_common->d_comm_peer[0].beginRecv();
         }
         if ( d_process_may_supply[1] ) {
            d_common->d_comm_peer[1].setPeerRank(d_contact[1]);
            d_common->d_comm_peer[1].setMPITag( CascadePartitionerTree_TAG_LoadTransfer0,
                                                CascadePartitionerTree_TAG_LoadTransfer1 );
            d_common->d_comm_peer[1].beginRecv();
         }
      }

      double work_supplied = d_far->supplyWork( -d_near->surplus(), d_common->d_mpi.getRank() );

      // Record work taken by near child.
      d_near->d_work += work_supplied;
      d_near->d_group_may_supply = d_near->d_process_may_supply[0] = false;

      if ( d_common->d_print_steps ) {
         tbox::plog << "CascadePartitionerTree::balanceConstituentHalves:"
                    << "  recorded supply of " << work_supplied
                    << " from far half to our half.  d_process_may_supply="
                    << d_process_may_supply[0] << ',' << d_process_may_supply[1] << "\n";
      }

      if ( d_process_may_supply[0] || d_process_may_supply[1] ) {
         receiveAndUnpackSuppliedLoad();
      }
   }

   // Complete the load send, if there was any.
   d_common->d_comm_stage.advanceAll();
   while ( d_common->d_comm_stage.numberOfCompletedMembers() > 0 ) {
      d_common->d_comm_stage.popCompletionQueue();
   }

   if ( d_common->d_print_steps ) {
      tbox::plog << "CascadePartitionerTree::balanceConstituentHalves: leaving with state:\n";
      printClassData( tbox::plog, "\t" );
   }
   return;
}



/*
 *************************************************************************
 * Supply work_requested to another group requesting work.  Any load
 * supplied by local process is will be sent to the designated taker,
 * a process in the requesting group.  Give priority to supplies
 * closest to the taker in rank space.
 *
 * 1. If this is a leaf, remove load from d_common->d_local_load
 *    and put it in d_common->d_shipment.
 * 2. Else:
 *    A: Remove load from the child group closer to the taker.
 *    B: If step A didn't supply enough, remove some from the other child.
 *
 * This method recurses into descendent groups, estimating changes
 * needed to supply the work_requested.  Estimation is necessary
 * because most of the changes take place remotely.  If the recursion
 * reaches the leaf group of the local process, it sets aside the work
 * that the local process supplies.
 *************************************************************************
 */
double
CascadePartitionerTree::supplyWork( double work_requested, int taker )
{
   TBOX_ASSERT( work_requested > 0.0 );
   TBOX_ASSERT( taker < d_begin && taker >= d_end ); // Taker must not be member.
   double est_work_supplied = 0.0; // Estimate of work supplied by this group.

   if ( d_children[0] != 0 ) {
      // Not a leaf:  Supply load from children.
      if ( taker < d_begin ) {
         est_work_supplied = d_children[0]->supplyWork( work_requested, taker );
         if ( est_work_supplied < work_requested ) {
            est_work_supplied += d_children[1]->supplyWork( work_requested-est_work_supplied, taker );
         }
      }
      else {
         est_work_supplied = d_children[1]->supplyWork( work_requested, taker );
         if ( est_work_supplied < work_requested ) {
            est_work_supplied = d_children[0]->supplyWork( work_requested-est_work_supplied, taker );
         }
      }
   }

   else {
      // Is a leaf:  Supply load from self.
      if ( d_group_may_supply ) {
         est_work_supplied = tbox::MathUtilities<double>::Min( work_requested, d_work );
         d_work -= est_work_supplied;

         if ( d_begin == d_common->d_mpi.getRank() ) {
            d_common->d_shipment->adjustLoad(
               *d_common->d_local_load,
               work_requested,
               work_requested,
               work_requested );
            if ( d_common->d_print_steps ) {
               tbox::plog << "CascadePartitionerTree::supplyWork giving to " << taker << ": ";
               d_common->d_shipment->recursivePrint();
               tbox::plog << "CascadePartitionerTree::supplyWork keeping: ";
               d_common->d_local_load->recursivePrint();
            }

            d_process_may_supply[0] =
               d_common->d_local_load->getSumLoad() - d_capacity > d_common->d_pparams->getLoadComparisonTol();

            d_group_may_supply = surplus() > d_common->d_pparams->getLoadComparisonTol();
         }

      }
   }

   return est_work_supplied;
}



/*
 *************************************************************************
 *************************************************************************
 */
void
CascadePartitionerTree::sendShipment( int taker )
{
   if ( d_common->d_print_steps ) {
      tbox::plog << "CascadePartitionerTree::sendMyShipment: sending to " << taker << ":\n";
      d_common->d_shipment->recursivePrint(tbox::plog, "Send: ", 2);
   }
   tbox::MessageStream msg;
   msg << *d_common->d_shipment;
   d_common->d_comm_peer[0].setPeerRank(taker);
   d_common->d_comm_peer[0].setMPITag( CascadePartitionerTree_TAG_LoadTransfer0,
                                       CascadePartitionerTree_TAG_LoadTransfer1 );
   d_common->d_comm_peer[0].beginSend( (const char*)(msg.getBufferStart()),
                                       msg.getCurrentSize() );
   d_common->d_shipment->clear();
   return;
}



/*
 *************************************************************************
 *************************************************************************
 */
void
CascadePartitionerTree::receiveAndUnpackSuppliedLoad()
{
   while ( d_common->d_comm_stage.advanceAny() ) {
      tbox::AsyncCommStage::Member *completed = d_common->d_comm_stage.popCompletionQueue();
      tbox::AsyncCommPeer<char> *comm_peer = static_cast<tbox::AsyncCommPeer<char>*>(completed);
      tbox::MessageStream recv_msg( comm_peer->getRecvSize(),
                                    tbox::MessageStream::Read,
                                    comm_peer->getRecvData(),
                                    true );
      recv_msg >> *d_common->d_local_load;
      if ( d_common->d_print_steps ) {
         tbox::plog << "CascadePartitionerTree::receiveAndUnpackSuppliedLoad: updated d_local_load:\n";
         d_common->d_local_load->recursivePrint(tbox::plog, "Curr: ", 2);
      }
   }
   return;
}



/*
 *************************************************************************
 *************************************************************************
 */
int
CascadePartitionerTree::cycleNum() const
{
   return CascadePartitioner::lgInt(d_common->d_mpi.getSize()) - d_generation;
}



/*
 *************************************************************************
 *************************************************************************
 */
void
CascadePartitionerTree::printClassData( std::ostream &co, const std::string &border ) const
{
   const std::string indent( border + std::string(cycleNum(),' ') + std::string(cycleNum(),' ') );
   co << indent << "generation=" << d_generation << "  cycle=" << cycleNum()
      << "  [" << d_begin << ',' << d_end << ")  group_size="
      << d_end-d_begin << '='
      << "  this=" << this
      << '\n' << indent
      << "position=" << d_position << "  contact=" << d_contact[0] << ',' << d_contact[1]
      << "  work=" << d_work << '/' << d_capacity
      << '\n' << indent
      << "group_may_supply=" << d_group_may_supply
      << "  process_may_supply=" << d_process_may_supply[0] << ',' << d_process_may_supply[1]
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
