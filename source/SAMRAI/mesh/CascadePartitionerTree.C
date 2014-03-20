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
   d_gen_num(0),

   d_begin(0),
   d_end(partitioner.d_mpi.getSize()),
   d_position(-1),

   d_parent(0),
   d_near(0),
   d_far(0),
   d_leaf(0),

   d_work(partitioner.d_local_load->getSumLoad()),
   d_capacity(d_common->d_global_load_avg*d_common->d_mpi.getSize()),
   d_group_may_supply(false)
{
   d_children[0] = d_children[1] = 0;
   d_contact[0] = d_contact[1] = -1;
   d_process_may_supply[0] = d_process_may_supply[1] = false;

   makeChildren();

   if ( d_common->d_print_steps ) {
      tbox::plog << "CascadePartitionerTree::root constructor: leaving\n";
      printClassData( tbox::plog, "\t" );
   }
}



/*
 * Construct a child group.  Child group has either the first half or
 * the second half of the parent group.  If parent group has odd
 * number of processes, the extra process is placed in the Upper child
 * group.
 */
CascadePartitionerTree::CascadePartitionerTree(
   CascadePartitionerTree &parent,
   Position group_position ) :
   d_common(parent.d_common),
   d_gen_num(1 + parent.d_gen_num),

   d_begin(parent.d_begin),
   d_end(parent.d_end),
   d_position(group_position),

   d_parent(&parent),
   d_near(0),
   d_far(0),
   d_leaf(0),

   d_work(-1.0),
   d_capacity(0),
   d_group_may_supply(false)
{
   d_children[0] = d_children[1] = 0;
   d_contact[0] = d_contact[1] = -1;
   d_process_may_supply[0] = d_process_may_supply[1] = false;

   const int upper_begin = (d_parent->d_begin + d_parent->d_end)/2;
   if ( group_position == Lower ) { d_end = upper_begin; }
   else { d_begin = upper_begin; }

   d_capacity = d_common->d_global_load_avg*(d_end-d_begin);

   /*
    * Assign contacts by pairing this process with the sibling group's
    * process having the same relative rank.  If upper group has an
    * odd rank, pair it with the last rank in lower group (giving last
    * rank in lower group two contacts).
    */
   const int relative_rank = d_common->d_mpi.getRank() < upper_begin ?
      d_common->d_mpi.getRank() - d_parent->d_begin :
      d_common->d_mpi.getRank() - upper_begin ;
   if ( group_position == Lower ) {
      d_contact[0] = relative_rank + upper_begin;
      if ( d_common->d_mpi.getRank() == upper_begin-1 &&
           (d_parent->d_end - d_parent->d_begin)%2 ) {
         d_contact[1] = 1 + d_contact[0];
      }
   }
   else {
      d_contact[0] = tbox::MathUtilities<int>::Min(
         relative_rank + d_parent->d_begin, upper_begin-1 );
   }

   d_process_may_supply[0] = d_process_may_supply[1] = false;

   // If local process is in this group, make its children.  Otherwise skip.
   if ( d_common->d_mpi.getRank() >= d_begin &&
        d_common->d_mpi.getRank() < d_end ) {
      makeChildren();
   }


   if ( d_common->d_print_steps ) {
      tbox::plog << "CascadePartitionerTree::non-root constructor: leaving\n";
      printClassData( tbox::plog, "\t" );
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

      const bool in_upper_branch = d_common->d_mpi.getRank() >= d_children[1]->d_begin;
      d_near = d_children[in_upper_branch];
      d_far = d_children[!in_upper_branch];

      d_leaf = d_near->d_leaf;
   }

   else {
      // This is a local leaf.
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
 * This method contains the looping structure required to balance all
 * groups.  Starting at the leaves, combine sibling data to balance
 * them, then continue toward the root.
 *************************************************************************
 */
void CascadePartitionerTree::balanceAll()
{
   int top_group_num = 0; // For debugging only.

   /*
    * Balancing a group doesn't guarantee all its descendent groups
    * are balanced.  This loop goes through all generations to balance
    * each group (top_group) in turn.
    */
   for ( CascadePartitionerTree *top_group=this; top_group!=d_leaf;
         top_group = top_group->d_near ) {

      tbox::plog << "\nCascadePartitionerTree::balanceAll balancing outer top_group "
                 << top_group->d_gen_num
                 << "  with exact local_load=" << d_common->d_local_load->getSumLoad()
                 << std::endl;

      d_leaf->recomputeLeafData();

      /*
       * Balance top_group by sweeping from leaf to top_group, shifting
       * load from surplus child to deficit child as we go.
       */
      for ( CascadePartitionerTree *current_group = d_leaf->d_parent;
            current_group != 0 && current_group->d_near != top_group;
            current_group = current_group->d_parent ) {

         current_group->combineChildren();
         if ( d_common->d_print_steps ) {
            tbox::plog << "\nCascadePartitionerTree::balanceAll outer top_group "
                       << top_group->d_gen_num << "  shuffling generation "
                       << current_group->d_gen_num << "\n";
            current_group->printClassData( tbox::plog, "\t" );
            tbox::plog << "\tchild 0:\n";
            current_group->d_children[0]->printClassData( tbox::plog, "\t" );
            tbox::plog << "\tchild 1:\n";
            current_group->d_children[1]->printClassData( tbox::plog, "\t" );
         }

         current_group->balanceChildren();
         if ( d_common->d_print_steps ) {
            tbox::plog << "\nCascadePartitionerTree::balanceAll outer top_group "
                       << top_group->d_gen_num << "  shuffled generation "
                       << current_group->d_gen_num << "\n";
            current_group->printClassData( tbox::plog, "\t" );
            tbox::plog << "\tchild 0:\n";
            current_group->d_children[0]->printClassData( tbox::plog, "\t" );
            tbox::plog << "\tchild 1:\n";
            current_group->d_children[1]->printClassData( tbox::plog, "\t" );
         }

      } // Inner loop.

      ++top_group_num; // For debugging only.

   } // Outer loop.

}



/*
 *************************************************************************
 * Combine near and far children data (using communication) to compute
 * work-related data for this group.
 *
 * NOTE: This method can probably be re-organized as if there's only
 * one contact.  We only need to exchange data with the second contact
 * if the the far group may supply and the near group has a deficit.
 *************************************************************************
 */
void
CascadePartitionerTree::combineChildren()
{
   if ( false && d_common->d_print_steps ) {
      tbox::plog << "CascadePartitionerTree::combineChildren: entered with state:\n";
      printClassData( tbox::plog, "\t" );
      tbox::plog << "\tchild 0:\n";
      d_children[0]->printClassData( tbox::plog, "\t" );
      tbox::plog << "\tchild 1:\n";
      d_children[1]->printClassData( tbox::plog, "\t" );
   }

   /*
    * Determine and record work of the two halves.  Needs
    * communication to get data about the far half of the group.
    */

   tbox::MessageStream send_msg;
   send_msg << d_near->d_work << d_near->d_group_may_supply << d_near->d_process_may_supply[0];

   for ( int i=0; i<2; ++i ) {
      if ( d_near->d_contact[i] >= 0 ) {

         d_common->d_comm_peer[i].setPeerRank(d_near->d_contact[i]);
         d_common->d_comm_peer[i].setMPITag( CascadePartitionerTree_TAG_InfoExchange0,
                                             CascadePartitionerTree_TAG_InfoExchange1 );
         d_common->d_comm_peer[i].limitFirstDataLength( send_msg.getCurrentSize() );
         d_common->d_comm_peer[i].beginRecv(true);

         d_common->d_comm_peer[2+i].setPeerRank(d_near->d_contact[i]);
         d_common->d_comm_peer[2+i].setMPITag( CascadePartitionerTree_TAG_InfoExchange0,
                                               CascadePartitionerTree_TAG_InfoExchange1 );
         d_common->d_comm_peer[2+i].limitFirstDataLength( send_msg.getCurrentSize() );
         d_common->d_comm_peer[2+i].beginSend( static_cast<const char*>(send_msg.getBufferStart()),
                                               send_msg.getCurrentSize(), true );

      }
   }

   d_far->d_work = 0.0;
   d_far->d_group_may_supply = true;
   while ( d_common->d_comm_stage.numberOfCompletedMembers() > 0 ||
           d_common->d_comm_stage.advanceAny() ) {

      tbox::AsyncCommPeer<char> *completed = static_cast<tbox::AsyncCommPeer<char>*>(
         d_common->d_comm_stage.popCompletionQueue() );

      const int i = completed - d_common->d_comm_peer;
      if ( i < 2 ) {
         // This was a receive.
         tbox::MessageStream recv_msg(completed->getRecvSize(),
                                      tbox::MessageStream::Read,
                                      completed->getRecvData(),
                                      false);
         recv_msg >> d_far->d_work >> d_far->d_group_may_supply >> d_far->d_process_may_supply[i];
      }
   }
   TBOX_ASSERT( d_common->d_comm_stage.numberOfPendingMembers() == 0 );

   // Work data is function of children data.
   d_work = d_children[0]->d_work + d_children[1]->d_work;
   d_group_may_supply = d_children[0]->d_group_may_supply || d_children[1]->d_group_may_supply;

   // If process still may supply for near child, it may supply for this group.
   d_process_may_supply[0] = d_near->d_process_may_supply[0];

   if ( false && d_common->d_print_steps ) {
      tbox::plog << "CascadePartitionerTree::combineChildren: leaving with state:\n";
      printClassData( tbox::plog, "\t" );
      tbox::plog << "\tchild 0:\n";
      d_children[0]->printClassData( tbox::plog, "\t" );
      tbox::plog << "\tchild 1:\n";
      d_children[1]->printClassData( tbox::plog, "\t" );
      tbox::plog << "\tlocal_load:";
      d_common->d_local_load->recursivePrint(tbox::plog, "\t",0);
   }
   return;
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
 *
 * Note: This method balances only the children groups (the lower
 * child and upper), not all descendents.  To balance the near grand
 * child, call d_near->balanceChildren().
 *************************************************************************
 */
void
CascadePartitionerTree::balanceChildren()
{
   if ( d_common->d_print_steps ) {
      tbox::plog << "CascadePartitionerTree::balanceChildren: entered\n";
   }

   TBOX_ASSERT( d_common->d_shipment->empty() );

   if ( d_near->surplus() > d_common->d_pparams->getLoadComparisonTol() &&
        d_far->surplus() < -d_common->d_pparams->getLoadComparisonTol() ) {
      // Outgoing work, from near child to far child.

      if ( d_near->d_process_may_supply[0] ) {

         double work_supplied = d_near->supplyWork( -d_far->surplus(), d_near->d_contact[0] );

         // Record work taken by the far child.
         d_far->d_work += work_supplied;
         d_far->d_group_may_supply = d_far->d_process_may_supply[0] = d_far->d_process_may_supply[1] = false;

         if ( d_common->d_print_steps ) {
            tbox::plog << "CascadePartitionerTree::balanceChildren:"
                       << "  record outgoing shipment of " << work_supplied
                       << " from our half to far half.  Send to " << d_near->d_contact[0] << "\n";
         }

         TBOX_ASSERT( d_near->d_contact[0] >= 0 );
         sendShipment(d_near->d_contact[0]); // If 2 contacts in far group, send to the first one only.
      }
   }

   else if ( d_far->surplus()  >  d_common->d_pparams->getLoadComparisonTol() &&
             d_near->surplus() < -d_common->d_pparams->getLoadComparisonTol() ) {
      // Incoming work, from far child to near child.

      /*
       * Even when a group may supply, it will not supply to a
       * redundant_demand.  It will supply to its first contact
       * everything it wants to send to the group.  A redundant_demand
       * is one from the second contact in the same group.  Local
       * demand is redundant if the near group is bigger than the far
       * group and local process is the last rank in its group.
       */
      const bool redundant_demand = d_near->size() > d_far->size() &&
         d_common->d_mpi.getRank() == d_end-1;

      if ( d_far->d_group_may_supply && !redundant_demand ) {
         if ( d_far->d_process_may_supply[0] ) {
            d_common->d_comm_peer[0].setPeerRank(d_near->d_contact[0]);
            d_common->d_comm_peer[0].setMPITag( CascadePartitionerTree_TAG_LoadTransfer0,
                                                CascadePartitionerTree_TAG_LoadTransfer1 );
            d_common->d_comm_peer[0].beginRecv(true);
            if ( d_common->d_print_steps ) {
               tbox::plog << "CascadePartitionerTree::balanceChildren:"
                          << "  expecting shipment from first contact " << d_near->d_contact[0]
                          << "  redundant_demand=" << redundant_demand
                          << "  d_near->size()=" << d_near->size()
                          << "  d_far->size()=" << d_far->size()
                          << "  d_end=" << d_end
                          << "  d_common->d_mpi.getRank()=" << d_common->d_mpi.getRank()
                          << "\n";
            }
         }
         if ( d_far->d_process_may_supply[1] ) {
            d_common->d_comm_peer[1].setPeerRank(d_near->d_contact[1]);
            d_common->d_comm_peer[1].setMPITag( CascadePartitionerTree_TAG_LoadTransfer0,
                                                CascadePartitionerTree_TAG_LoadTransfer1 );
            d_common->d_comm_peer[1].beginRecv(true);
            if ( d_common->d_print_steps ) {
               tbox::plog << "CascadePartitionerTree::balanceChildren:"
                          << "  expecting shipment from second contact " << d_near->d_contact[1]
                          << "  redundant_demand=" << redundant_demand
                          << "  d_near->size()=" << d_near->size()
                          << "  d_far->size()=" << d_far->size()
                          << "  d_end=" << d_end
                          << "  d_common->d_mpi.getRank()=" << d_common->d_mpi.getRank()
                          << "\n";
            }
         }
      }

      double work_supplied = d_far->supplyWork( -d_near->surplus(), d_common->d_mpi.getRank() );

      // Record work taken by near child group.
      d_near->d_work += work_supplied;
      d_near->d_group_may_supply = d_near->d_process_may_supply[0] = false;

      if ( d_common->d_print_steps ) {
         tbox::plog << "CascadePartitionerTree::balanceChildren:"
                    << "  recorded incoming shipment of " << work_supplied
                    << " from far half to our half.\n";
      }

      if ( d_far->d_process_may_supply[0] || d_far->d_process_may_supply[1] ) {
         receiveAndUnpackSuppliedLoad();
      }
   }

   // Complete the load send, if there was any.
   d_common->d_comm_stage.advanceAll();
   while ( d_common->d_comm_stage.numberOfCompletedMembers() > 0 ) {
      d_common->d_comm_stage.popCompletionQueue();
   }

   if ( d_common->d_print_steps ) {
      tbox::plog << "CascadePartitionerTree::balanceChildren: leaving with state:\n";
      printClassData( tbox::plog, "\t" );
      tbox::plog << "\tchild 0:\n";
      d_children[0]->printClassData( tbox::plog, "\t" );
      tbox::plog << "\tchild 1:\n";
      d_children[1]->printClassData( tbox::plog, "\t" );
      tbox::plog << "\tlocal_load:";
      d_common->d_local_load->recursivePrint(tbox::plog, "\t",0);
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
   TBOX_ASSERT( !containsRank(taker) );
   TBOX_ASSERT( containsRank(d_common->d_mpi.getRank()) || d_children[0] == 0 ); // Only near groups should store children.

   double est_work_supplied = 0.0; // Estimate of work supplied by this group.

   if ( d_children[0] != 0 ) {
      // Near group and not a leaf:  Supply load from children.
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

   else if ( !containsRank(d_common->d_mpi.getRank()) ) {
      // Is a far group: Doesn't matter if it's a leaf.
      est_work_supplied = tbox::MathUtilities<double>::Min( work_requested, surplus() );
      d_work -= est_work_supplied;
   }

   else {
      /*
       * Is a near leaf group: Compute how much work we can supply,
       * update group's work estimate and apportion the load to
       * supply.
       */
      if ( d_group_may_supply ) {
         est_work_supplied = tbox::MathUtilities<double>::Min( work_requested, surplus() );
         d_work -= est_work_supplied;

         const double tolerance = d_common->d_flexible_load_tol*d_common->d_global_load_avg;
         d_common->d_shipment->adjustLoad(
            *d_common->d_local_load,
            est_work_supplied,
            est_work_supplied-tolerance,
            est_work_supplied+tolerance );
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
      tbox::plog << "CascadePartitionerTree::sendMyShipment: sending to " << taker << ' ';
      d_common->d_shipment->recursivePrint(tbox::plog, "", 0);
      tbox::plog << " leaving d_local_load with ";
      d_common->d_local_load->recursivePrint(tbox::plog, "", 0);
      tbox::plog << '\n';
   }
   tbox::MessageStream msg;
   msg << *d_common->d_shipment;
   d_common->d_comm_peer[0].setPeerRank(taker);
   d_common->d_comm_peer[0].setMPITag( CascadePartitionerTree_TAG_LoadTransfer0,
                                       CascadePartitionerTree_TAG_LoadTransfer1 );
   d_common->d_comm_peer[0].beginSend( (const char*)(msg.getBufferStart()),
                                       msg.getCurrentSize(), true );
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
   while ( d_common->d_comm_stage.numberOfCompletedMembers() > 0 ||
           d_common->d_comm_stage.advanceAny() ) {
      tbox::AsyncCommStage::Member *completed = d_common->d_comm_stage.popCompletionQueue();
      tbox::AsyncCommPeer<char> *comm_peer = static_cast<tbox::AsyncCommPeer<char>*>(completed);
      tbox::MessageStream recv_msg( comm_peer->getRecvSize(),
                                    tbox::MessageStream::Read,
                                    comm_peer->getRecvData(),
                                    true );
      recv_msg >> *d_common->d_shipment;
      d_common->d_local_load->insertAll(*d_common->d_shipment);
      if ( d_common->d_print_steps ) {
         tbox::plog << "CascadePartitionerTree::receiveAndUnpackSuppliedLoad: received ";
         d_common->d_shipment->recursivePrint(tbox::plog, "", 0);
         tbox::plog << " from process " << comm_peer->getPeerRank()
                    << " and updated d_local_load to ";
         d_common->d_local_load->recursivePrint(tbox::plog, "", 0);
         tbox::plog << '\n';
      }
      d_common->d_shipment->clear();
   }
   return;
}



/*
 *************************************************************************
 * Should only be called for leaves.  For bigger groups, communication
 * is required and is determined as part of work exchange.
 *************************************************************************
 */
void
CascadePartitionerTree::recomputeLeafData()
{
   TBOX_ASSERT( this == d_leaf ); // Should only be called for leaves.
   d_work = d_common->d_local_load->getSumLoad();
   d_group_may_supply = d_process_may_supply[0] =
      d_work > d_capacity + d_common->d_pparams->getLoadComparisonTol() ;
   /*
    * TODO: Try not resetting d_group_may_supply.  Seems bad that a
    * leaf that received slightly more than its demand would now try
    * to send that slight bit elsewhere.
    */
}



/*
 *************************************************************************
 *************************************************************************
 */
int
CascadePartitionerTree::cycleNum() const
{
   return CascadePartitioner::lgInt(d_common->d_mpi.getSize()) - d_gen_num;
}



/*
 *************************************************************************
 *************************************************************************
 */
void
CascadePartitionerTree::printClassData( std::ostream &co, const std::string &border ) const
{
   const std::string indent( border + std::string(d_gen_num,' ') + std::string(d_gen_num,' ') );
   co << indent << "gen_num=" << d_gen_num << "  cycle=" << cycleNum()
      << "  [" << d_begin << ',' << d_end << ")  group_size=" << d_end-d_begin
      << "  local leaf=" << (this==d_leaf) << "  this=" << this
      << "  near=" << d_near << "  far=" << d_far
      << '\n' << indent
      << "position=" << d_position << "  contact=" << d_contact[0] << ',' << d_contact[1]
      << "  work=" << d_work << '/' << d_capacity
      << "  local_load=" << d_common->d_local_load->getSumLoad()
      << '\n' << indent
      << "group_may_supply=" << d_group_may_supply
      << "  process_may_supply=" << d_process_may_supply[0] << ',' << d_process_may_supply[1]
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
