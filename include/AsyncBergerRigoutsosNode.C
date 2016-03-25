/*
 * File:         $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-1/source/mesh/clustering/AsyncBergerRigoutsosNode.C $
 * Copyright:    (c) 1997-2007 Lawrence Livermore National Security, LLC
 * Revision:     $LastChangedRevision: 1889 $
 * Modified:     $LastChangedDate: 2008-01-22 16:46:52 -0800 (Tue, 22 Jan 2008) $
 * Description:  Node in asynchronous Berger-Rigoutsos dendogram
 */

#ifndef included_mesh_AsyncBergerRigoutsosNode_C
#define included_mesh_AsyncBergerRigoutsosNode_C

#include "AsyncBergerRigoutsosNode.h"
#include "CellData.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"

#include <set>

namespace SAMRAI {
namespace mesh {

/*
********************************************************************
*                                                                  *
* Constructor.                                                     *
*                                                                  *
* Set the node id and generation number to their proper values     *
* (see documentation for those parameters).                        *
*                                                                  *
* If constructing a root node,                                     *
*   - set MPI tag                                                  *
*   - set generation number to 1                                   *
*   - set the wait phase for starting                              *
*                                                                  *
* If possible, give the node id d_pos the position number in the   *
* dendogram.  If that turns out to be too big, set d_pos to -1 if  *
* it is the left child and -2 if it is the right child.            *
*                                                                  *
********************************************************************
*/
template<int DIM>
AsyncBergerRigoutsosNode<DIM>::AsyncBergerRigoutsosNode(
   CommonParams *common_params,
   const hier::Box<DIM> *bound_box,
   AsyncBergerRigoutsosNode *parent,
   const int child_number )
   : d_pos( parent == NULL ? 1 :
           ( parent->d_pos > 0 && 
             parent->d_pos < tbox::MathUtilities<int>::getMax()/2 ) ?
           2*parent->d_pos+child_number :
           ( child_number==0 ? -1 : -2 ) ),
     d_common(common_params),
     d_parent(parent),
     d_lft_child(NULL),
     d_rht_child(NULL),
     d_box( bound_box ? *bound_box : hier::Box<DIM>() ),
     d_owner( parent == NULL ? 0 : -1 ),
     d_group(0),
     d_mpi_tag(-1),
     d_overlap(-1),
     d_box_acceptance(undetermined),
     d_node(),
     d_node_iterator(),
     d_wait_phase( parent == NULL ? to_be_launched : for_data_only),
     d_send_msg(),
     d_recv_msg(),
     d_comm_group(NULL),
     d_generation( parent == NULL ? 1 : d_parent->d_generation+1),
     d_n_cont(0)
{

#ifdef DEBUG_CHECK_ASSERTIONS
   if ( d_parent != NULL && d_parent->d_pos >= 0 && d_pos < 0 ) {
      TBOX_WARNING( "Too many generations for node identification.\n"
                    <<"The node id cannot be increased any further.\n"
                    <<"This affects only the node id, which is only\n"
                    <<"used for analysis and debugging and does not\n"
                    <<"affect the algorithm.\n"
                    <<"Last valid node id is " << d_parent->d_pos << '\n' );
   }
#endif


#ifdef DEBUG_CHECK_ASSERTIONS
   d_node_iterator = GraphNodeContainer().end();
#endif

   if ( d_parent == NULL ) {
#if defined(HAVE_MPI)
      int group_size = 1;
      MPI_Comm_size(d_common->mpi_communicator, &group_size);
      d_group.resizeArray(group_size);
      for ( int i=0; i<group_size; ++i ) {
         d_group[i] = i;
      }
#else
      d_group.resizeArray(1);
      d_group[0] = 0;
#endif
   }

   if ( d_parent == NULL ) {
      if ( d_common->rank == 0 ) {
         claimMPITag();
      }
      else {
         d_mpi_tag = 0;
      }
#ifdef DEBUG_CHECK_ASSERTIONS
      /*
       * Even though owner has different way of getting MPI tag,
       * all processes should start with the same mpi tag.
       */
      TBOX_ASSERT( d_mpi_tag == 0 );
#endif
   }


   if ( d_parent == NULL && d_owner == d_common->rank ) {
      ++(d_common->num_nodes_owned);
      ++(d_common->max_nodes_owned);
   }
   ++(d_common->num_nodes_allocated);
   ++(d_common->max_nodes_allocated);
   if ( d_common->max_generation < d_generation ) {
      d_common->max_generation = d_generation;
   }

   if ( d_common->log_node_history ) {
      tbox::plog << d_common->num_nodes_allocated << "-allocated  "
                 << d_common->num_nodes_active << "-active  "
                 << d_common->num_nodes_owned << "-owned  "
                 << d_common->num_nodes_completed << "-completed  "
                 << "Construct " << d_generation << ':' << d_pos
                 << ".\n";
   }

   return;
}


template<int DIM>
AsyncBergerRigoutsosNode<DIM>::~AsyncBergerRigoutsosNode(void)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   /*
    * Forbid deleting a node that is running because there may
    * be pending communication (by the node or its children).
    * Note that this is NOT an extra restriction over the
    * recursive implementation.
    */
   if ( d_wait_phase != for_data_only &&
        d_wait_phase != to_be_launched &&
        d_wait_phase != completed ) {
      TBOX_ERROR("Should not delete a node that is currently running\n"
                 <<"the Berger-Rigoutsos algorithm because there\n"
                 <<"may be pending communications.");
   }
#endif

   if ( d_comm_group != NULL ) {
      if ( ! d_comm_group->isDone() ) {
         TBOX_ERROR("Library error: Destructing a node with an unfinished\n"
                    <<"communication tree is bad because it leaves\n"
                    <<"outstanding MPI messages.");
      }
      // d_common->job_relauncher.getCommStage().deallocateCommGroup(d_comm_group);
      delete d_comm_group;
      d_comm_group = NULL;
   }

   --(d_common->num_nodes_allocated);

   if ( d_parent != NULL && d_common->log_node_history ) {
      tbox::plog << d_common->num_nodes_allocated << "-allocated  "
                 << d_common->num_nodes_active << "-active  "
                 << d_common->num_nodes_owned << "-owned  "
                 << d_common->num_nodes_completed << "-completed  "
                 << "Destruct " << d_generation << ':' << d_pos
                 << "  " << d_node
                 << "  " << d_box
                 << ".\n";
   }

   d_wait_phase = deallocated;

   return;
}



/*
********************************************************************
* Preprocess, run the asynchronous JobRelauncher and postprocess.  *
********************************************************************
*/
template<int DIM>
void AsyncBergerRigoutsosNode<DIM>::runAlgorithm()
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( d_parent == NULL );
#endif

   d_common->t_cluster->start();

   /*
    * If compute_edges == 1:
    *   - Compute edges from tagged level to new levels.
    *     These edges are organized around the tagged nodes.
    *     They do not need to be shared with the owners of the
    *     new nodes.
    *
    * If compute_edges == 2:
    *   - Compute edges as in compute_edges == 1 case.
    *   - Owners of new edges send new edge data to owners
    *     of new nodes.  This creates the neighbor data
    *     organized around the new nodes.
    */

   if ( d_common->compute_edges > 0 ) {

      d_common->tag_node_set.setTo(*d_common->level);

      /*
       * Create empty neighbor lists for nodes on tagged layer.
       * As new nodes are finalized, they will be added to
       * these lists.
       */
      const GraphNodeContainer &tag_node_container =
         d_common->tag_node_set.getNodeContainer();
      for ( typename GraphNodeContainer::const_iterator
               ni=tag_node_container.begin();
            ni != tag_node_container.end(); ++ni ) {
         d_common->tag_cnect_new[(*ni).getLocalIndex()];
      }
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT( tag_node_container.size() == d_common->tag_cnect_new.size() );
#endif

   }

   if ( d_common->compute_edges > 1 ) {
      d_common->new_node_set.setRefinementRatio(d_common->level->getRatio());
   }

   d_common->job_relauncher.runAlgorithm( this );

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( d_wait_phase == completed );
   if ( d_common->compute_edges > 2 ) {
      // Each new node should have its own neighbor list.
      TBOX_ASSERT( d_common->new_node_set.getNodeContainer().size() ==
                   d_common->new_cnect_tag.size() );
   }
#endif

   d_common->t_cluster->stop();

   /*
    * Share edges with owners, if requested.
    * This is a one-time operation that is not considered a part
    * of continueAlgorithm(), so it lies outside that timimg.
    */
   if ( d_common->compute_edges > 1 ) {
      shareNewEdgesWithOwners();
   }

   return;
}



/*
********************************************************************
* This method looks messy, but it is just the BR agorithm,         *
* with multiple pause and continue points implemented by           *
* the goto and labels.  Each pause point is accompanied by         *
* a line setting d_wait_phase so that the algorithm can            *
* continue where it left off when this method is called again.     *
* The BR algorithm is not completed until this method returns      *
* the WaitPhase value "completed".                                 *
********************************************************************
*/
template<int DIM>
typename AsyncBergerRigoutsosNode<DIM>::WaitPhase
AsyncBergerRigoutsosNode<DIM>::continueAlgorithm()
{
   d_common->t_continue_algorithm->start();
   ++d_n_cont;


#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( d_parent == NULL || d_parent->d_wait_phase != completed );
   TBOX_ASSERT( ! inRelaunchQueue(this) );
#endif

   /*
    * Skip right to where we left off,
    * which is specified by the wait phase variable.
    */
   switch ( d_wait_phase ) {
   case for_data_only:
      TBOX_ERROR("Library error: Attempt to execute data-only node.");
   case to_be_launched:
      goto TO_BE_LAUNCHED;
   case reduce_histogram:
      goto REDUCE_HISTOGRAM;
   case bcast_acceptability:
      goto BCAST_ACCEPTABILITY;
   case gather_grouping_criteria:
      goto GATHER_GROUPING_CRITERIA;
   case bcast_child_groups:
      goto BCAST_CHILD_GROUPS;
   case run_children:
      goto RUN_CHILDREN;
   case bcast_to_dropouts:
      goto BCAST_TO_DROPOUTS;
   case completed:
      TBOX_ERROR("Library error: Senseless continuation of completed node.");
   default:
      TBOX_ERROR("Library error: Nonexistent phase.");
   }

   bool sub_completed;

   /*
    * Delegated tasks: Major tasks are delegated to private methods.
    * These methods may check whether the process is the owner or
    * just a contributor and fork appropriately.  The communication
    * checking tasks return whether communication is completed, but
    * they do NOT change the d_wait_phase variable, which is done
    * in this function.
    */


  TO_BE_LAUNCHED:

   ++(d_common->num_nodes_active);
   if ( d_common->max_nodes_active < d_common->num_nodes_active ) {
      d_common->max_nodes_active = d_common->num_nodes_active;
   }

   if ( d_common->log_node_history ) {
      tbox::plog << d_common->num_nodes_allocated << "-allocated  "
                 << d_common->num_nodes_active << "-active  "
                 << d_common->num_nodes_owned << "-owned  "
                 << d_common->num_nodes_completed << "-completed  "
                 << "Commence " << d_generation << ':' << d_pos
                 << "  " << d_box
                 << "  accept=" << d_box_acceptance
                 << "  ovlap=" << d_overlap
                 << "  owner=" << d_owner
                 << "  gsize=" << d_group.size()
                 << ".\n";
   }

   if ( d_parent == NULL || d_overlap > 0 || d_common->rank == d_owner ) {

#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT( inGroup(d_group) );
#endif

      // Set up communication group for operations in participating group.
      d_comm_group = d_common->job_relauncher.getCommStage().
         allocateCommGroup(computeCommunicationTreeDegree(d_group.size()),
                           this);
      d_comm_group->setMPICommunicator(d_common->mpi_communicator);
      d_comm_group->setGroupAndRootRank(d_group, d_owner);


      makeLocalTagHistogram();


      if ( d_group.size() > 1 ) {
         d_common->t_reduce_histogram->start();
         reduceHistogram_start();
        REDUCE_HISTOGRAM:
         if ( ! d_common->t_reduce_histogram->isRunning() )
            d_common->t_reduce_histogram->start();
         sub_completed = reduceHistogram_check();
         d_common->t_reduce_histogram->stop();
         if ( !sub_completed ) {
            d_wait_phase = reduce_histogram;
            goto RETURN;
         }
      }


      if ( d_common->rank == d_owner ) {
         computeMinimalBoundingBoxForTags();
         acceptOrSplitBox();
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT( boxAccepted() || boxRejected() ||
                      ( boxHasNoTag() && d_parent == 0 ) );
#endif
         if ( ! boxHasNoTag() ) {
            /*
             * A layer node is created even if box is not acceptable,
             * so that the children can reference its local index in case
             * the box is later accepted based on the combined tolerance
             * of the children.  The node would be erased later if
             * it is not finally accepted.
             */
            createLayerNode();
         }
      }


      if ( d_group.size() > 1 ) {
         d_common->t_bcast_acceptability->start();
         broadcastAcceptability_start();
        BCAST_ACCEPTABILITY:
         if ( ! d_common->t_bcast_acceptability->isRunning() )
            d_common->t_bcast_acceptability->start();
         sub_completed = broadcastAcceptability_check();
         d_common->t_bcast_acceptability->stop();
         if ( !sub_completed ) {
            d_wait_phase = bcast_acceptability;
            goto RETURN;
         }
      }
#ifdef DEBUG_CHECK_ASSERTIONS
      if ( d_common->rank == d_owner ) {
         TBOX_ASSERT( d_box_acceptance == accepted_by_calculation ||
                      d_box_acceptance == rejected_by_calculation ||
                      d_box_acceptance == hasnotag_by_owner );
      }
      else {
         TBOX_ASSERT( d_box_acceptance == accepted_by_owner ||
                      d_box_acceptance == rejected_by_owner ||
                      d_box_acceptance == hasnotag_by_owner );
      }
#endif


      if ( boxRejected() ) {

         if ( ! d_common->use_level_boxes ) {

            /*
             * Compute children groups and owners without assuming
             * entire mesh structure is known locally.
             */
            countOverlapWithLocalPatches();

            if ( d_group.size() > 1 ) {
               d_common->t_gather_grouping_criteria->start();
               gatherGroupingCriteria_start();
              GATHER_GROUPING_CRITERIA:
               if ( ! d_common->t_gather_grouping_criteria->isRunning() )
                  d_common->t_gather_grouping_criteria->start();
               sub_completed = gatherGroupingCriteria_check();
               d_common->t_gather_grouping_criteria->stop();
               if ( !sub_completed ) {
                  d_wait_phase = gather_grouping_criteria;
                  goto RETURN;
               }
            }

            if ( d_common->rank == d_owner ) {
               formChildGroups();
            }

            if ( d_group.size() > 1 ) {
               d_common->t_bcast_child_groups->start();
               broadcastChildGroups_start();
              BCAST_CHILD_GROUPS:
               if ( ! d_common->t_bcast_child_groups->isRunning() )
                  d_common->t_bcast_child_groups->start();
               sub_completed = broadcastChildGroups_check();
               d_common->t_bcast_child_groups->stop();
               if ( !sub_completed ) {
                  d_wait_phase = bcast_child_groups;
                  goto RETURN;
               }
            }

         }
         else {
            /*
             * Compute children groups and owners using local copy
             * of entire mesh structure.
             */
            formChildGroupsUsingLevelBoxes();
         }

         if ( d_lft_child->d_owner == d_common->rank ) {
            ++(d_common->num_nodes_owned);
            ++(d_common->max_nodes_owned);
         }
         if ( d_rht_child->d_owner == d_common->rank ) {
            ++(d_common->num_nodes_owned);
            ++(d_common->max_nodes_owned);
         }

         runChildren_start();
        RUN_CHILDREN:
         sub_completed = runChildren_check();
         if ( !sub_completed ) {
            d_wait_phase = run_children;
            goto RETURN;
         }
      }

      else if ( boxAccepted() ) {
         if ( d_common->rank == d_owner ) {
            ++(d_common->num_boxes_generated);
         }
      }

      else {
         // Box has no tag.
      }

      // All done with communication within participating group.
      delete d_comm_group;
      d_comm_group = NULL;

   }

   else {
#ifdef DEBUG_CHECK_ASSERTIONS
      /*
       * This process is not in the group that decides on the box for
       * this dendogram node.
       */
      TBOX_ASSERT( ! inGroup(d_group) );
#endif
   }


   if ( d_parent == NULL ) {
      /*
       * Compute edges and set up edge sharing data.
       * This is usually done by a node's parent in the
       * runChildren_check() method because only the
       * parent can know if the node's box will be
       * kept or recombined with the sibling.
       * But the root node must do this itself because it has no parent.
       */
      if ( d_common->compute_edges > 0 && boxAccepted() ) {
         computeNewGraphEdges();
      }
   }


#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( d_lft_child == NULL );
   TBOX_ASSERT( d_rht_child == NULL );
   TBOX_ASSERT( ! inRelaunchQueue(this) );
#endif


   /*
    * Broadcast the result to dropouts.
    * Dropout processes are those that participated in the
    * parent but not in this dendogram node.  They need the
    * result to perform combined efficiency check for the
    * parent.
    *
    * Processes that should participate in the dropout broadcast
    * are the dropouts (processes with zero overlap) and the owner.
    *
    * Broadcast to dropouts is only needed if:
    *
    *    - In multi-owner mode and edge-computing mode.
    *      In single-owner mode, only the original owner needs
    *      the final result, and it participates everywhere,
    *      so there is no need for this phase.
    *      When computing edges, participant processors must
    *      know results to do recombination check, to determine
    *      if parent box is preferred.
    *
    *    - This is NOT the root dendogram node.  The root node
    *      has no parent and no corresponding dropout group.
    *
    *    - Dropout group is not empty.  Number of dropouts
    *      is the difference between parent group size and this
    *      group size.
    */
   if ( d_overlap == 0 || d_common->rank == d_owner ) {

      if ( ( d_common->owner_mode != SINGLE_OWNER ||
             d_common->compute_edges > 0 ) &&
           d_parent != NULL &&
           d_parent->d_group.size() > d_group.size() ) {

         d_common->t_bcast_to_dropouts->start();
         {
            // Create the communication group for the dropouts.
            tbox::Array<int> dropouts;
            computeDropoutGroup( d_parent->d_group,
                                 d_group,
                                 dropouts,
                                 d_owner );
            d_comm_group = d_common->job_relauncher.getCommStage().
               allocateCommGroup(
                  computeCommunicationTreeDegree(dropouts.size()),
                  this);
            d_comm_group->setMPICommunicator(d_common->mpi_communicator);
            d_comm_group->setGroupAndRootIndex(dropouts, 0);
         }

         broadcastToDropouts_start();
        BCAST_TO_DROPOUTS:
         if ( ! d_common->t_bcast_to_dropouts->isRunning() )
            d_common->t_bcast_to_dropouts->start();
         sub_completed = broadcastToDropouts_check();
         d_common->t_bcast_to_dropouts->stop();

         if ( !sub_completed ) {
            d_wait_phase = bcast_to_dropouts;
            goto RETURN;
         }

         if ( d_common->log_node_history && d_common->rank != d_owner ) {
            tbox::plog << d_common->num_nodes_allocated << "-allocated  "
                       << d_common->num_nodes_active << "-active  "
                       << d_common->num_nodes_owned << "-owned  "
                       << d_common->num_nodes_completed << "-completed  "
                       << "DO Recv " << d_generation << ':' << d_pos
                       << "  " << d_node
                       << "  accept=" << d_box_acceptance
                       << ".\n";
         }

         delete d_comm_group;
         d_comm_group = NULL;
      }
   }


   d_wait_phase = completed;

   if ( d_comm_group != NULL ) {
      // No further communication.  Deallocate the communication group.
      // d_common->job_relauncher.getCommStage().deallocateCommGroup(d_comm_group);
      delete d_comm_group;
      d_comm_group = NULL;
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( d_common->num_nodes_owned >= 0 );
#endif

   // Adjust counters.
   --(d_common->num_nodes_active);
   ++(d_common->num_nodes_completed);
   if ( d_owner == d_common->rank ) {
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT( d_common->num_nodes_owned > 0 );
#endif
      --(d_common->num_nodes_owned);
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT( d_common->num_nodes_owned >= 0 );
#endif
   }
   d_common->num_conts_to_complete += d_n_cont;
   if ( d_common->max_conts_to_complete < d_n_cont ) {
      d_common->max_conts_to_complete = d_n_cont;
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( d_common->num_nodes_owned >= 0 );
#endif

   if ( d_common->log_node_history ) {
      tbox::plog << d_common->num_nodes_allocated << "-allocated  "
                 << d_common->num_nodes_active << "-active  "
                 << d_common->num_nodes_owned << "-owned  "
                 << d_common->num_nodes_completed << "-completed  "
                 << "Complete " << d_generation << ':' << d_pos
                 << "  " << d_node
                 << "  accept=" << d_box_acceptance
                 << ".\n";
   }


   /*
    * Recall that an dendogram node waiting for its children
    * is not placed in the relaunch queue (because it is
    * pointless to relaunch it until the children are completed).
    * Therefore, to eventually continue that node, its last
    * child to complete must put it on the queue.  If this node
    * and its sibling are completed, put the parent on the FRONT
    * queue to be checked immediately.
    */
   if ( d_parent != NULL &&
        d_parent->d_lft_child->d_wait_phase == completed &&
        d_parent->d_rht_child->d_wait_phase == completed ) {
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT( d_parent->d_wait_phase == run_children );
      TBOX_ASSERT( ! inRelaunchQueue(d_parent) );
#endif
      d_common->job_relauncher.getRelaunchQueue().appendItem(d_parent);
   }


  RETURN:


#ifdef DEBUG_CHECK_ASSERTIONS
   tbox::List<tbox::RelaunchableJob*>::Iterator li =
      inRelaunchQueue(this);
   if ( li ) {
      TBOX_ERROR("Library error: adding a duplicate to the leaf queue.");
   }
#endif

#ifdef DEBUG_CHECK_ASSERTIONS
   if ( d_wait_phase != completed && d_wait_phase != run_children ) {
      TBOX_ASSERT( ! d_comm_group->isDone() );
      TBOX_ASSERT( ! d_common->job_relauncher.getCommStage().isDone() );
   }
   if ( d_wait_phase == run_children ) {
      TBOX_ASSERT( ! d_common->job_relauncher.getRelaunchQueue().isEmpty() );
   }
#endif

   d_common->t_continue_algorithm->stop();

   return d_wait_phase;
}




template<int DIM>
void AsyncBergerRigoutsosNode<DIM>::runChildren_start()
{
   /*
    * Children were created to store temporary data
    * and determine participation. Now, run them.
    */

#ifdef DEBUG_CHECK_ASSERTIONS
   /*
    * Should only be here if box is rejected based on calculation.
    */
   TBOX_ASSERT( d_box_acceptance == rejected_by_calculation ||
                d_box_acceptance == rejected_by_owner );
#endif

   d_lft_child->d_wait_phase = to_be_launched;
   d_rht_child->d_wait_phase = to_be_launched;

   /*
    * Queue the children so they get executed.
    * Put them at the front because they have
    * immediate computation (compute histogram)
    * to perform.  Put the left child in front
    * of the right to more closely match the
    * progression of the recursive BR (not essential).
    */
   d_common->job_relauncher.getRelaunchQueue().appendItem(d_rht_child);
   d_common->job_relauncher.getRelaunchQueue().appendItem(d_lft_child);
   return;
}





/*
********************************************************************
* Check for combined tolerance.                                    *
* If both children accepted their boxes without further splitting  *
* but their combined efficiency is not good enough to make         *
* the splitting worth accepting, use the current box instead       *
* of the children boxes.  Otherwise, use the children boxes.       *
********************************************************************
*/
template<int DIM>
bool AsyncBergerRigoutsosNode<DIM>::runChildren_check()
{
   if ( d_lft_child->d_wait_phase != completed ||
        d_rht_child->d_wait_phase != completed ) {
      return false;
   }

   if ( d_lft_child->boxAccepted() &&
        d_rht_child->boxAccepted() &&
        ( d_lft_child->d_box.size() + d_rht_child->d_box.size() >=
          d_common->combine_tol*d_box.size() ) ) {

      // Discard childrens' graph nodes in favor of recombination.

      d_box_acceptance = accepted_by_recombination;

      if ( d_common->log_node_history ) {
         tbox::plog << d_common->num_nodes_allocated << "-allocated  "
                    << d_common->num_nodes_active << "-active  "
                    << d_common->num_nodes_owned << "-owned  "
                    << d_common->num_nodes_completed << "-completed  "
                    << "Recombine " << d_generation << ':' << d_pos
                    << "  " << d_node
                    << " <= " << d_lft_child->d_node
                    << " + " << d_rht_child->d_node
                    << ".\n";
      }

      if ( d_lft_child->d_owner == d_common->rank ) {
         d_lft_child->eraseLayerNode();
         d_lft_child->d_box_acceptance = rejected_by_recombination;
         --(d_common->num_boxes_generated);
      }

      if ( d_rht_child->d_owner == d_common->rank ) {
         d_rht_child->eraseLayerNode();
         d_rht_child->d_box_acceptance = rejected_by_recombination;
         --(d_common->num_boxes_generated);
      }

      if ( d_owner == d_common->rank ) {
         ++(d_common->num_boxes_generated);
      }

   }
   else {

      // Accept childrens' results, discarding graph node.

      if ( d_owner == d_common->rank ) {
         eraseLayerNode();
      }
      if ( d_common->compute_edges > 0 ) {
         if ( d_lft_child->boxAccepted() &&
              d_lft_child->d_box_acceptance != accepted_by_dropout_bcast ) {
            d_lft_child->computeNewGraphEdges();
         }
         if ( d_rht_child->boxAccepted() &&
              d_rht_child->d_box_acceptance != accepted_by_dropout_bcast ) {
            d_rht_child->computeNewGraphEdges();
         }
      }

   }

   /*
    * No longer need children nodes after this point.
    */
   delete d_lft_child;
   delete d_rht_child;
   d_lft_child = d_rht_child = NULL;

   return true;
}





template<int DIM>
AsyncBergerRigoutsosNode<DIM>::CommonParams::CommonParams(
   const tbox::Pointer<hier::PatchLevel<DIM> > level_,
   const int tag_data_index_,
   const int tag_val_ ,
   const hier::IntVector<DIM> min_box_,
   const double efficiency_tol_,
   const double combine_tol_,
   const tbox::SAMRAI_MPI::comm mpi_communicator_)
   :
     new_node_set(),
     tag_cnect_new(),
     edge_senders(),
     edge_messages(),
     level(level_),
     // Box<DIM> clustering parameters ...
     tag_data_index(tag_data_index_),
     tag_val(tag_val_),
     min_box(min_box_),
     efficiency_tol(efficiency_tol_),
     combine_tol(combine_tol_),
     // Implementation flags ...
     compute_edges(2),
     max_gcw(hier::IntVector<DIM>(1)),
     owner_mode(MOST_OVERLAP),
     use_level_boxes(false),
     // Communication parameters ...
     mpi_communicator(mpi_communicator_),
     rank(getRank(mpi_communicator)),
     nproc(getProcCount(mpi_communicator)),
     available_mpi_tag(-1),
     // Analysis support ...
     log_node_history(false),
     num_nodes_allocated(0),
     max_nodes_allocated(0),
     num_nodes_active(0),
     max_nodes_active(0),
     num_nodes_owned(0),
     max_nodes_owned(0),
     num_nodes_completed(0),
     max_generation(0),
     num_boxes_generated(0),
     num_conts_to_complete(0),
     max_conts_to_complete(0)
{
   /*
    * Reserve the tag upper bound for the edge-sharing phase.
    * Divide the rest into tag pools divided among all processes.
    */
#if defined(HAVE_MPI)
   int *tag_upper_bound_ptr, flag;
   MPI_Attr_get( MPI_COMM_WORLD,
                 MPI_TAG_UB,
                 &tag_upper_bound_ptr,
                 &flag );
   tag_upper_bound = *tag_upper_bound_ptr;
#else
   tag_upper_bound = 1000000; // Parameter only used with MPI.
#endif

   available_mpi_tag = tag_upper_bound/nproc*rank;

   // Timers

   t_cluster = tbox::TimerManager::getManager()->
      getTimer("mesh::AsyncBergerRigoutsosNode::cluster");
   t_continue_algorithm = tbox::TimerManager::getManager()->
      getTimer("mesh::AsyncBergerRigoutsosNode::continueAlgorithm()");

   t_compute_new_graph_edges = tbox::TimerManager::getManager()->
      getTimer("mesh::AsyncBergerRigoutsosNode::computeNewGraphEdges()");
   t_share_new_edges = tbox::TimerManager::getManager()->
      getTimer("mesh::AsyncBergerRigoutsosNode::shareNewEdgesWithOwners()");
   t_share_new_edges_send = tbox::TimerManager::getManager()->
      getTimer("mesh::AsyncBergerRigoutsosNode::shareNewEdgesWithOwners()_send");
   t_share_new_edges_recv = tbox::TimerManager::getManager()->
      getTimer("mesh::AsyncBergerRigoutsosNode::shareNewEdgesWithOwners()_recv");
   t_share_new_edges_unpack = tbox::TimerManager::getManager()->
      getTimer("mesh::AsyncBergerRigoutsosNode::shareNewEdgesWithOwners()_unpack");

   // Multi-stage timers.
   t_reduce_histogram = tbox::TimerManager::getManager()->
      getTimer("mesh::AsyncBergerRigoutsosNode::reduce_histogram");
   t_bcast_acceptability = tbox::TimerManager::getManager()->
      getTimer("mesh::AsyncBergerRigoutsosNode::bcast_acceptability");
   t_gather_grouping_criteria = tbox::TimerManager::getManager()->
      getTimer("mesh::AsyncBergerRigoutsosNode::gather_grouping_criteria");
   t_bcast_child_groups = tbox::TimerManager::getManager()->
      getTimer("mesh::AsyncBergerRigoutsosNode::bcast_child_groups");
   t_bcast_to_dropouts = tbox::TimerManager::getManager()->
      getTimer("mesh::AsyncBergerRigoutsosNode::bcast_to_dropouts");

   return;
}




template<int DIM>
int AsyncBergerRigoutsosNode<DIM>::CommonParams::getRank(
   const tbox::SAMRAI_MPI::comm &mpi_communicator_ )
{
#ifdef HAVE_MPI
   int rank = 0;
   MPI_Comm_rank(mpi_communicator_, &rank);
   return rank;
#else
   return 0;
#endif
}




template<int DIM>
int AsyncBergerRigoutsosNode<DIM>::CommonParams::getProcCount(
   const tbox::SAMRAI_MPI::comm &mpi_communicator_ )
{
#ifdef HAVE_MPI
   int proc_count = 0;
   MPI_Comm_size(mpi_communicator_, &proc_count);
   return proc_count;
#else
   return 1;
#endif
}



/*
********************************************************************
*                                                                  *
* Asynchronous methods: these methods have _start and _check       *
* suffices.  They involve initiating some task and checking        *
* whether that task is completed.                                  *
*                                                                  *
********************************************************************
*/


template<int DIM>
void AsyncBergerRigoutsosNode<DIM>::reduceHistogram_start()
{
   if ( d_group.size() == 1 ) {
      return;
   }
   d_comm_group->setMPITag(d_mpi_tag + reduce_histogram_tag);
   const int hist_size = getHistogramBufferSize(d_box);
   if ( d_common->rank == d_owner ) {
      d_recv_msg.setNull();
      d_recv_msg.resizeArray( hist_size );
      putHistogramToBuffer( d_recv_msg.getPointer() );
      d_comm_group->beginSumReduce( d_recv_msg.getPointer(), hist_size );
   }
   else {
      d_send_msg.setNull();
      d_send_msg.resizeArray(hist_size);
      putHistogramToBuffer( d_send_msg.getPointer() );
      d_comm_group->beginSumReduce( d_send_msg.getPointer(), hist_size );
   }
   return;
}


template<int DIM>
bool AsyncBergerRigoutsosNode<DIM>::reduceHistogram_check()
{
   if ( d_group.size() == 1 ) {
      return true;
   }
   if ( d_comm_group->checkOperation() && d_common->rank == d_owner ) {
      getHistogramFromBuffer( d_recv_msg.getPointer() );
   }
   return d_comm_group->isDone();
}


template<int DIM>
void AsyncBergerRigoutsosNode<DIM>::broadcastAcceptability_start()
{
   if ( d_group.size() == 1 ) {
      return;
   }
   d_comm_group->setMPITag( d_mpi_tag + bcast_acceptability_tag );
   /*
    * Items communicated:
    * - local index of node
    * - whether box is accepted
    * - in case box is accepted:
    *   . box (which may have been trimmed to minimal tag bounding box)
    * - in case box is rejected:
    *   . left/right child boxes
    *   . left/right child MPI tags
    */

   const int buffer_size = 1           // Acceptability flag.
                         + 1           // Local index of node.
                         + DIM*2      // Box<DIM>.
                         + DIM*4      // Children boxes.
                         + 2           // Children MPI tags
                         ;

   if ( d_common->rank == d_owner ) {
#ifdef DEBUG_CHECK_ASSERTION
      TBOX_ASSERT( d_box_acceptance == rejected_by_calculation ||
                   d_box_acceptance == accepted_by_calculation ||
                   ( d_parent == NULL && d_box_acceptance == hasnotag_by_owner )
                   );
#endif
      d_send_msg.setNull();
      d_send_msg.resizeArray( buffer_size );
      int *ptr = d_send_msg.getPointer();
      *(ptr++) = d_box_acceptance >= 0 ?
         d_box_acceptance + 2 /* indicate remote decision */ :
         d_box_acceptance;
      if ( ! boxHasNoTag() ) {
         *(ptr++) = d_node.getLocalIndex();
         ptr = putBoxToBuffer(d_box, ptr);
         if ( boxRejected() ) {
            ptr = putBoxToBuffer(d_lft_child->d_box, ptr);
            ptr = putBoxToBuffer(d_rht_child->d_box, ptr);
            *(ptr++) = d_lft_child->d_mpi_tag;
            *(ptr++) = d_rht_child->d_mpi_tag;
         }
      }
      d_comm_group->beginBcast( d_send_msg.getPointer(), buffer_size );
   }
   else {
      d_recv_msg.setNull();
      d_recv_msg.resizeArray( buffer_size );
      d_comm_group->beginBcast( d_recv_msg.getPointer(), buffer_size );
   }
   return;
}


template<int DIM>
bool AsyncBergerRigoutsosNode<DIM>::broadcastAcceptability_check()
{
   if ( d_group.size() == 1 ) {
      return true;
   }
   if ( d_comm_group->checkBcast() && d_common->rank != d_owner ) {

      int *ptr = d_recv_msg.getPointer();
      d_box_acceptance = intToBoxAcceptance (*(ptr++) );
      if ( ! boxHasNoTag() ) {
         const typename GraphNode::LocalIndex node_local_index = *(ptr++);
         ptr = getBoxFromBuffer(d_box, ptr);
         d_node = GraphNode( d_box, node_local_index, d_owner );
      }

#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT( d_node.getLocalIndex() >= 0 );
      TBOX_ASSERT( boxAccepted() || boxRejected() ||
                   ( boxHasNoTag() && d_parent == NULL ) );
      TBOX_ASSERT( d_box.numberCells() >= d_common->min_box );
#endif

      if ( boxRejected() ) {
           
         /*
          * The owner formed its children earlier so it can
          * use their parameters while determining which to run.
          * Contributors create the children when the receive
          * the d_box_acceptance flag indicates tha further
          * branching is required.
          */
         d_lft_child = new AsyncBergerRigoutsosNode( d_common,
                                                     NULL,
                                                     this,
                                                     0);
         d_rht_child = new AsyncBergerRigoutsosNode( d_common,
                                                     NULL,
                                                     this,
                                                     1);

         ptr = getBoxFromBuffer(d_lft_child->d_box, ptr);
         ptr = getBoxFromBuffer(d_rht_child->d_box, ptr);

         d_lft_child->d_mpi_tag = *(ptr++);
         d_rht_child->d_mpi_tag = *(ptr++);

#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT( d_lft_child->d_box.numberCells() >= d_common->min_box );
         TBOX_ASSERT( d_rht_child->d_box.numberCells() >= d_common->min_box );
         TBOX_ASSERT( d_lft_child->d_mpi_tag > -1 );
         TBOX_ASSERT( d_rht_child->d_mpi_tag > -1 );
#endif
         if ( d_common->log_node_history ) {
            tbox::plog << d_common->num_nodes_allocated << "-allocated  "
                       << d_common->num_nodes_active << "-active  "
                       << d_common->num_nodes_owned << "-owned  "
                       << d_common->num_nodes_completed << "-completed  "
                       << "Rm Split " << d_generation << ':' << d_pos
                       << "  " << d_node
                       << " => " << d_lft_child->d_box
                       << " + " << d_rht_child->d_box
                       << ".\n";
         }

      }
   }
   return d_comm_group->isDone();
}


template<int DIM>
void AsyncBergerRigoutsosNode<DIM>::gatherGroupingCriteria_start()
{
   if ( d_group.size() == 1 ) {
      return;
   }
   d_comm_group->setMPITag( d_mpi_tag + gather_grouping_criteria_tag );

   if ( d_common->rank == d_owner ) {
      d_recv_msg.setNull();
      d_recv_msg.resizeArray( 4*d_group.size() );
      d_comm_group->beginGather( d_recv_msg.getPointer(), 4 );
   }
   else {
      d_send_msg.setNull();
      d_send_msg.resizeArray(4);
      d_send_msg[0] = d_lft_child->d_overlap;
      d_send_msg[1] = d_rht_child->d_overlap;
      // Use negative burden measures for uniformity of criteria comparison.
      d_send_msg[2] = -d_common->num_nodes_owned;
      d_send_msg[3] = -d_common->num_nodes_active;
      d_comm_group->beginGather( d_send_msg.getPointer(), 4 );
   }

   return;
}


template<int DIM>
bool AsyncBergerRigoutsosNode<DIM>::gatherGroupingCriteria_check()
{
   if ( d_group.size() == 1 ) {
      return true;
   }
   d_comm_group->checkGather();
   /*
    * Do nothing yet with the overlap data d_recv_msg.
    * We extract it in formChildGroups().
    */
   return d_comm_group->isDone();
}


template<int DIM>
void AsyncBergerRigoutsosNode<DIM>::broadcastChildGroups_start()
{
   if ( d_group.size() == 1 ) {
      return;
   }
   /*
    * Items communicated:
    * - left/right owner
    * - left/right group
    */
   d_comm_group->setMPITag(d_mpi_tag + bcast_child_groups_tag);

   if ( d_common->rank == d_owner ) {

      const int buffer_size = 2                // Left/right owners.
                            + 2                // Left/right group sizes.
                            + d_lft_child->d_group.size() // Left group.
                            + d_rht_child->d_group.size() // Right groups.
                            ;

      int i;

      d_send_msg.setNull();
      d_send_msg.resizeArray( buffer_size );
      int *ptr = d_send_msg.getPointer();

      *(ptr++) = d_lft_child->d_owner;
      *(ptr++) = d_lft_child->d_group.size();
      for ( i=0; i<d_lft_child->d_group.size(); ++i ) {
         *(ptr++) = d_lft_child->d_group[i];
      }
      *(ptr++) = d_rht_child->d_owner;
      *(ptr++) = d_rht_child->d_group.size();
      for ( i=0; i<d_rht_child->d_group.size(); ++i ) {
         *(ptr++) = d_rht_child->d_group[i];
      }

      d_comm_group->beginBcast( d_send_msg.getPointer(), buffer_size );
   }
   else {
      const int buffer_size = 2                // Left/right owners.
                            + 2                // Left/right group sizes.
                            + 2*d_group.size() // Left/right groups.
                            ;
      d_recv_msg.setNull();
      d_recv_msg.resizeArray( buffer_size );

      d_comm_group->beginBcast( d_recv_msg.getPointer(), buffer_size );
   }

   return;
}


template<int DIM>
bool AsyncBergerRigoutsosNode<DIM>::broadcastChildGroups_check()
{
   if ( d_group.size() == 1 ) {
      return true;
   }
   if ( d_comm_group->checkBcast() && d_common->rank != d_owner ) {

      int i;
      int *ptr = d_recv_msg.getPointer();

      d_lft_child->d_owner = *(ptr++);
      d_lft_child->d_group.resizeArray( *(ptr++) );
      for ( i=0; i<d_lft_child->d_group.size(); ++i ) {
         d_lft_child->d_group[i] = *(ptr++);
      }
      d_rht_child->d_owner = *(ptr++);
      d_rht_child->d_group.resizeArray( *(ptr++) );
      for ( i=0; i<d_rht_child->d_group.size(); ++i ) {
         d_rht_child->d_group[i] = *(ptr++);
      }

#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT( d_lft_child->d_owner >= 0 );
      TBOX_ASSERT( d_lft_child->d_group.size() > 0 );
      TBOX_ASSERT( ( d_lft_child->d_overlap > 0 ) ==
                   inGroup(d_lft_child->d_group) );
      TBOX_ASSERT( d_rht_child->d_owner >= 0 );
      TBOX_ASSERT( d_rht_child->d_group.size() > 0 );
      TBOX_ASSERT( ( d_rht_child->d_overlap > 0 ) ==
                   inGroup(d_rht_child->d_group) );
#endif

   }

   return d_comm_group->isDone();
}




template<int DIM>
void AsyncBergerRigoutsosNode<DIM>::broadcastToDropouts_start()
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( d_common->rank == d_owner || d_overlap == 0 );
#endif
   d_comm_group->setMPITag( d_mpi_tag + bcast_to_dropouts_tag );

   const int buffer_size = 1      // d_box_acceptance
                         + 1      // local index of graph node
                         + DIM*2  // d_box (in case it got reduced)
                         ;
   d_send_msg.setNull();
   d_recv_msg.setNull();
   if ( d_common->rank == d_owner ) {
      d_send_msg.resizeArray( buffer_size );
      d_send_msg[0] = d_box_acceptance;
      d_send_msg[1] = d_node.getLocalIndex();
      putBoxToBuffer( d_box, &d_send_msg[2] );
      d_comm_group->beginBcast( d_send_msg.getPointer(),
                               buffer_size );
   }
   else {
      d_recv_msg.resizeArray( buffer_size );
      d_comm_group->beginBcast( d_recv_msg.getPointer(),
                               buffer_size );
   }
   return;
}




template<int DIM>
bool AsyncBergerRigoutsosNode<DIM>::broadcastToDropouts_check()
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( d_common->rank == d_owner || d_overlap == 0 );
#endif
   if ( d_comm_group->checkBcast() ) {
      if ( d_common->rank != d_owner ) {
#ifdef DEBUG_CHECK_ASSERTIONS
         /*
          * We check for the case of the box having no tags,
          * to keeps things explicit and help detect bugs.
          * But in fact, having no tags is impossible
          * in the broadcastToDropout step, because it is
          * only possible for the root dendogram node,
          * which has no dropout group.
          */
         TBOX_ASSERT( d_recv_msg[0] >= 0 );
#endif
         d_box_acceptance = intToBoxAcceptance( ( d_recv_msg[0] % 2 ) +
                                                rejected_by_dropout_bcast );
         const int local_index = d_recv_msg[1];
         getBoxFromBuffer( d_box, &d_recv_msg[2] );
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT( d_box.numberCells() >= d_common->min_box );
#endif
         d_node = GraphNode( d_box, local_index, d_owner );
      }
   }
   return d_comm_group->isDone();
}



/*
********************************************************************
* Utility computations using local data.                           *
********************************************************************
*/



template<int DIM>
void AsyncBergerRigoutsosNode<DIM>::makeLocalTagHistogram()
{

   int d;

   /*
    * Compute the histogram size and allocate space for it.
    */
   for ( d=0; d<DIM; ++d ) {
      d_histogram[d].setNull();
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT( d_box.numberCells(d) > 0 );
#endif
      d_histogram[d].resizeArray( d_box.numberCells(d) );
      int i;
      for ( i=0; i<d_histogram[d].size(); ++i ) {
         d_histogram[d][i] = 0;
      }
   }

   /*
    * Accumulate tag counts in the histogram variable.
    */
   hier::PatchLevel<DIM> &level = *d_common->level;
   for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
      hier::Patch<DIM> &patch = *level.getPatch(ip());

      hier::Box<DIM> intersection = patch.getBox() * d_box;

      if ( !(intersection.empty()) ) {

         tbox::Pointer< pdat::CellData<DIM,int> >
            tag_data_ = patch.getPatchData(d_common->tag_data_index);

         pdat::CellData<DIM,int> &tag_data = *tag_data_;

         for (pdat::CellIterator<DIM> ci(intersection); ci; ci++) {
            if ( tag_data(*ci) == d_common->tag_val ) {
               const int *hi = (int*)( *ci - d_box.lower() );
               for ( d=0; d<DIM; ++d ) {
                  ++(d_histogram[d][hi[d]]);
               }
            }
         }
      }
   }
   return;
}


/*
********************************************************************
* Change d_box to that of the minimal bounding box for tags.       *
* If d_box is changed, reduce d_histogram to new d_box.            *
********************************************************************
*/
template<int DIM>
void AsyncBergerRigoutsosNode<DIM>::computeMinimalBoundingBoxForTags()
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( ! d_box.empty() );
#endif

   int d;

   hier::Index<DIM> new_lower = d_box.lower();
   hier::Index<DIM> new_upper = d_box.upper();

   const hier::IntVector<DIM> &min_box = d_common->min_box;
   hier::IntVector<DIM> box_size = d_box.numberCells();

   /*
    * Bring the lower side of the box up past untagged index planes.
    * Bring the upper side of the box down past untagged index planes.
    * Do not make the box smaller than the min_box requirement.
    */
   for ( d=0; d<DIM; ++d ) {
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT( d_histogram[d].size() != 0 );
#endif
      int *histogram_beg = d_histogram[d].getPointer();
      int *histogram_end = histogram_beg + d_box.numberCells(d) - 1;
      while ( *histogram_beg == 0 &&
              box_size(d) > min_box(d) ) {
         ++new_lower(d);
         ++histogram_beg;
         --box_size(d);
      }
      while ( *histogram_end == 0 &&
              box_size(d) > min_box(d) ) {
         --new_upper(d);
         --histogram_end;
         --box_size(d);
      }
   }

   const hier::Box<DIM> new_box(new_lower,new_upper);
   const hier::IntVector<DIM> new_size = new_box.numberCells();

   if ( new_box != d_box ) {
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT( new_box.numberCells() >= d_common->min_box );
#endif
      /*
       * Save tagged part of the current histogram and reset the box.
       * Is this step really required?  No, we can just keep the
       * shift in a hier::IntVector<DIM> and adjust.
       */
      for ( d=0; d<DIM; ++d ) {
         tbox::Array<int> &h = d_histogram[d];
         const int shift = new_lower(d) - d_box.lower()(d);
         if ( shift > 0 ) {
            int i;
            for ( i=0; i<new_size(d); ++i ) {
               h[i] = h[i+shift];
            }
         }
         h.resizeArray(new_size(d));
      }
      d_box = new_box;
   }

   return;
}



/*
*********************************************************************
* Accept the box or split it, setting d_box_acceptance accordingly. *
*********************************************************************
*/
template<int DIM>
void AsyncBergerRigoutsosNode<DIM>::acceptOrSplitBox()
{

#ifdef DEBUG_CHECK_ASSERTIONS
   if ( d_owner != d_common->rank ) {
      TBOX_ERROR("Only the owner can determine\n"
                 "whether to accept or split a box.\n");
   }
   TBOX_ASSERT( d_box_acceptance == undetermined );
#endif

   /*
    * Box<DIM> d_box is acceptable if
    * - it has a high enough fraction of tagged cells, or
    * - it cannot be split without breaking the minimum
    *   box requirement, or
    *
    * If d_box has no tags:
    * - set d_box_acceptance = hasnotag_by_owner
    * If accepting d_box:
    * - set d_box_acceptance = accepted_by_calculation
    * If rejecting d_box:
    * - set d_box_acceptance = rejected_by_calculation
    * - create left and right children
    * - set children boxes
    * - claim MPI tags for communication at children nodes
    *
    * Instead of writing from scratch,
    * the code to find the split plane was copied
    * from mesh::BergerRigoutsos<DIM>::splitTagBoundBox()
    * and modified.
    */


   if ( d_box_acceptance == undetermined ) {
      /*
       * See if d_box should be accepted based on efficiency.
       */
      int i, num_tagged=0;
      for ( i=0; i<d_histogram[0].size(); ++i ) {
         num_tagged += d_histogram[0][i];
      }
      int num_cells = d_box.size();
      double efficiency = ( num_cells == 0 ? 1.e0 :
                            ((double)num_tagged)/num_cells );
      if ( efficiency >= d_common->efficiency_tol ) {
         d_box_acceptance = accepted_by_calculation;
      }
// Disable Intel warning on real comparison
#ifdef __INTEL_COMPILER
#pragma warning (disable:1572)
#endif
      else if ( efficiency == 0 ) {
         // No tags!  This should be caught at the dendogram root.
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT( d_parent == NULL );
#endif
         d_box_acceptance = hasnotag_by_owner;
         if ( d_common->log_node_history ) {
            tbox::plog << d_common->num_nodes_allocated << "-allocated  "
                       << d_common->num_nodes_active << "-active  "
                       << d_common->num_nodes_owned << "-owned  "
                       << d_common->num_nodes_completed << "-completed  "
                       << "HasNoTag " << d_generation << ':' << d_pos
                       << "  " << d_box
                       << ".\n";
         }
      }
   }


   hier::IntVector<DIM> sorted;

   if ( d_box_acceptance == undetermined ) {
      /*
       * Sort the bound box dimensions from largest to smallest.
       */
      const hier::IntVector<DIM> num_cells = d_box.numberCells();
      int dim;
      for (dim = 0; dim < DIM; dim++) {
         sorted(dim) = dim;
      }
      for (int dim0 = 0; dim0 < DIM-1; dim0++) {
         for (int dim1 = dim0+1; dim1 < DIM; dim1++) {
            if ( num_cells(sorted(dim0)) < num_cells(sorted(dim1)) ) {
               int tmp_dim = sorted(dim0);
               sorted(dim0) = sorted(dim1);
               sorted(dim1) = tmp_dim; 
            }
         }
      }
#ifdef DEBUG_CHECK_ASSERTIONS
      for (dim = 0; dim < DIM-1; dim++) {
         TBOX_ASSERT( num_cells(sorted(dim)) >= num_cells(sorted(dim+1)) );
      }
#endif
   }


   int num_splitable_dim = 0;

   if ( d_box_acceptance == undetermined ) {
      /*
       * Determine number of coordinate directions in bounding box
       * that are splittable according to the minimum box size
       * restriction.  Accept d_box if no such directions exist.
       */
      hier::IntVector<DIM> num_cells = d_box.numberCells();
      for (num_splitable_dim=0; num_splitable_dim<DIM;
           num_splitable_dim++) {
         int tmp_dim = sorted(num_splitable_dim);
         if ( num_cells(tmp_dim) < 2*d_common->min_box(tmp_dim) ) {
            break;
         }
      }
      if (num_splitable_dim == 0) {
         d_box_acceptance = accepted_by_calculation;
      }
   }


   if ( d_box_acceptance == undetermined ) {

      /*
       * Attempt to split box at a zero interior point in the
       * histogram.  Check each splittable direction, from
       * largest to smallest, until zero point found.
       */

      int cut_lo, cut_hi;
      int cut_pt = -(tbox::MathUtilities<int>::getMax());
      int cut_dim = -1; 
      int dim = -1;
      const hier::Index<DIM> box_lo(d_box.lower()); 
      const hier::Index<DIM> box_hi(d_box.upper());
      hier::Index<DIM> lft_hi(box_hi); 
      hier::Index<DIM> rht_lo(box_lo); 

      for (dim = 0; dim < num_splitable_dim; dim++) {
         cut_dim = sorted(dim);
#if 1
         if ( findZeroCutSwath( cut_lo, cut_hi, cut_dim ) ) {
            // Split bound box at cut_pt; cut_dim is splitting dimension.
            TBOX_ASSERT( cut_hi-cut_lo >= 0 );
            lft_hi(cut_dim) = cut_lo - 1;
            rht_lo(cut_dim) = cut_hi + 1;
            break;
         }
#else
         if ( findZeroCutPoint( cut_pt, cut_dim ) ) {
            // Split bound box at cut_pt; cut_dim is splitting dimension.
            lft_hi(cut_dim) = cut_pt;
            rht_lo(cut_dim) = cut_pt + 1;
            break;
         }
#endif
      }

      /*
       * If no zero point found,
       * try Laplacian on longest side of bound box.
       */

      if (dim == num_splitable_dim) {
         cut_dim = sorted(0);
         cutAtLaplacian(cut_pt, 
                        cut_dim, box_lo(cut_dim), box_hi(cut_dim), 
                        d_common->min_box(cut_dim));
         // Split bound box at cut_pt; cut_dim is splitting dimension.
         lft_hi(cut_dim) = cut_pt;
         rht_lo(cut_dim) = cut_pt + 1;
      }


      /*
       * The owner forms its children now so it can use their
       * parameters while determining which to run.
       * Contributors create the children when the receive
       * the d_box_acceptance flag from the owner.
       */
      d_lft_child = new AsyncBergerRigoutsosNode( d_common,
                                                  NULL,
                                                  this,
                                                  0);
      d_rht_child = new AsyncBergerRigoutsosNode( d_common,
                                                  NULL,
                                                  this,
                                                  1);

      d_lft_child->d_box = hier::Box<DIM>(box_lo, lft_hi);
      d_rht_child->d_box = hier::Box<DIM>(rht_lo, box_hi);
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT( d_lft_child->d_box.numberCells() >= d_common->min_box );
      TBOX_ASSERT( d_rht_child->d_box.numberCells() >= d_common->min_box );
#endif

      d_lft_child->claimMPITag();
      d_rht_child->claimMPITag();

      d_box_acceptance = rejected_by_calculation;

      if ( d_common->log_node_history ) {
         tbox::plog << d_common->num_nodes_allocated << "-allocated  "
                    << d_common->num_nodes_active << "-active  "
                    << d_common->num_nodes_owned << "-owned  "
                    << d_common->num_nodes_completed << "-completed  "
                    << "Lc Split " << d_generation << ':' << d_pos
                    << "  " << d_box
                    << " => " << d_lft_child->d_box
                    << " + " << d_rht_child->d_box
                    << ".\n";
      }

   }

   return;
}



/*
********************************************************************
*                                                                  *
* Attempt to find a range with zero histogram value near the       *
* middle of d_box in the given coordinate direction.               *
* Note that the hole is kept more than a minimium distance from    *
* the endpoints of of the index interval.                          *
*                                                                  *
* Note that it is assumed that box indices are cell indices.       *
*                                                                  *
* If a hole is found, cut_lo and cut_hi are set to the             *
* high index planes in the hole.                                   *
*                                                                  *
********************************************************************
*/

template<int DIM>
bool AsyncBergerRigoutsosNode<DIM>::findZeroCutSwath(
   int &cut_lo,
   int &cut_hi,
   const int dim )
{
   const int lo = d_box.lower(dim);
   const int hi = d_box.upper(dim);
   // Compute the limit for the swath.
   const int cut_lo_lim = lo + d_common->min_box(dim);
   const int cut_hi_lim = hi - d_common->min_box(dim);

   /*
    * Start in the middle of the box.
    * Move cut_lo down and cut_hi up until a hole is found.
    * Keep moving in same direction of the hole until the
    * other side of the hole is found.  The two planes form
    * the widest cut possible at the hole.
    */
   cut_lo = cut_hi = (lo + hi)/2;
   while ( (cut_lo>=cut_lo_lim) && (cut_hi<=cut_hi_lim) ) {
      if ( d_histogram[dim][cut_lo-lo] == 0 ) {
         /* The narrow cut is at cut_lo.  Initialize the cut swath here
          * and move cut_lo down until the far side the hole is found.
          */
         cut_hi = cut_lo;
         while ( ( (cut_lo>cut_lo_lim) ) &&
                 ( d_histogram[dim][cut_lo-lo-1] == 0 ) ) {
            --cut_lo;
         }
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT( cut_hi >= cut_lo );
         TBOX_ASSERT( cut_lo - lo >= d_common->min_box(dim) );
         TBOX_ASSERT( hi - cut_hi >= d_common->min_box(dim) );
         for ( int i=cut_lo; i<=cut_hi; ++i ) {
            TBOX_ASSERT( d_histogram[dim][i-lo] == 0 );
         }
#endif
         return(true);
      }
      if ( d_histogram[dim][cut_hi-lo] == 0) {
         /* The narrow cut is at cut_hi.  Initialize the cut swath here
          * and move cut_hi up until the far side the hole is found.
          */
         cut_lo = cut_hi;
         while ( ( (cut_hi<cut_hi_lim) ) &&
                 ( d_histogram[dim][cut_hi-lo+1] == 0 ) ) {
            ++cut_hi;
         }
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT( cut_hi >= cut_lo );
         TBOX_ASSERT( cut_lo - lo >= d_common->min_box(dim) );
         TBOX_ASSERT( hi - cut_hi >= d_common->min_box(dim) );
         for ( int i=cut_lo; i<=cut_hi; ++i ) {
            TBOX_ASSERT( d_histogram[dim][i-lo] == 0 );
         }
#endif
         return(true);
      }
      --cut_lo;
      ++cut_hi;
   }
   
   return(false);
}



/*
********************************************************************
*                                                                  *
* Attempt to find a zero histogram value near the middle of the    *
* index interval (lo, hi) in the given coordinate direction.       *
* Note that the cut_pt is kept more than a minimium distance from  *
* the endpoints of of the index interval. Since box indices are    *
* cell-centered, the cut point value corresponds to the right      *
* edge of the cell whose index is equal to the cut point.          *
*                                                                  *
* Note that it is assumed that box indices are cell indices.       *
*                                                                  *
********************************************************************
*/

template<int DIM>
bool AsyncBergerRigoutsosNode<DIM>::findZeroCutPoint(
   int &cut_pt,
   const int dim )
{
   const int lo = d_box.lower(dim);
   const int hi = d_box.upper(dim);
   const int cut_lo_lim = lo + d_common->min_box(dim) - 1;
   const int cut_hi_lim = hi - d_common->min_box(dim);
   const int cut_mid = (lo + hi)/2;
   int cut_lo = (lo + hi)/2;
   int cut_hi = (lo + hi)/2;

   int ic;
   for ( ic=0;
        ( (cut_lo>=cut_lo_lim) && (cut_hi<=cut_hi_lim) );
        ++ic, --cut_lo, ++cut_hi ) {
      if ( d_histogram[dim][cut_lo-lo] == 0 ) {
         cut_pt = cut_mid-ic;
         return(true);
      }
      if ( d_histogram[dim][cut_hi-lo] == 0) {
         cut_pt = cut_mid+ic;
         return(true);
      }
   }
   
   return(false);
}

/*
***********************************************************************
*                                                                     *
* Attempt to find a point in the given coordinate direction near an   *
* inflection point in the histogram for that direction. Note that the *
* cut point is kept more than a minimium distance from the endpoints  *
* of the index interval (lo, hi).  Also, the box must have at least   *
* three cells along a side to apply the Laplacian test.  If no        *
* inflection point is found, the mid-point of the interval is         *
* returned as the cut point.                                          *
*                                                                     *
* Note that it is assumed that box indices are cell indices.          *
*                                                                     *
***********************************************************************
*/

template<int DIM>
void AsyncBergerRigoutsosNode<DIM>::cutAtLaplacian(
   int &cut_pt, 
   const int dim,
   const int lo, 
   const int hi, 
   const int min_size)
{
   int loc_zero = (lo + hi)/2;

   if ( (hi - lo + 1) >= 3 ) {

      int max_zero = 0;
      const int infpt_lo_lim = tbox::MathUtilities<int>::Max( lo+min_size-1, 
                                                              lo+1 );  
      const int infpt_hi_lim = tbox::MathUtilities<int>::Min( hi-min_size, 
                                                              hi-2 );  

      int last_lap = d_histogram[dim][infpt_lo_lim-1-d_box.lower()(dim)]
                     - 2 * d_histogram[dim][infpt_lo_lim-d_box.lower()(dim)]
                     + d_histogram[dim][infpt_lo_lim + 1-d_box.lower()(dim)];
  
      for (int ic = infpt_lo_lim+1; ic <= infpt_hi_lim+1; ic++) {
         int new_lap = d_histogram[dim][ic - 1-d_box.lower()(dim)]
                       - 2 * d_histogram[dim][ic-d_box.lower()(dim)]
                       + d_histogram[dim][ic + 1-d_box.lower()(dim)];

         if ( ((new_lap < 0) && (last_lap >= 0)) ||
              ((new_lap >= 0) && (last_lap < 0)) ) {

            int delta = new_lap - last_lap;

            if ( delta < 0 ) {
               delta = -delta; 
            }

            if ( delta > max_zero ) {
               loc_zero = ic - 1;
               max_zero = delta;
            }
         }
     
         last_lap = new_lap;
      }

   }

   cut_pt = loc_zero;
   return;
}




/*
********************************************************************
* Create a DNBG layer node in d_common->new_node_set,              *
* where the output boxes of the algorithm is saved.                *
*                                                                  *
* Only the owner should create the layer node this way.            *
* Other processes build layer node using data from owner.          *
********************************************************************
*/
template<int DIM>
void AsyncBergerRigoutsosNode<DIM>::createLayerNode()
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( d_common->rank == d_owner );
#endif
   /*
    * The second argument to addBox states whether
    * a vacant index should be used (if available).
    * If false, the layer nodes are created with strictly
    * ascending indices, which will order the output BoxList<DIM>
    * in the same order as an recursive implementation would
    * (given that this algorithm is run in synchronous,
    * single-owner mode using level boxes).
    */
   d_node_iterator = d_common->new_node_set.addBox(d_box, false);
   d_node = *d_node_iterator;
   return;
}




/*
********************************************************************
* Discard the DNBG layer node.  On the owner, this node is a part  *
* of d_common->new_node_set where it must be removed.              *
* On contributors the layer node can just be ignored.              *
* To prevent busg, the node and its iterator are                   *
* set to unusable values.                                          *
********************************************************************
*/
template<int DIM>
void AsyncBergerRigoutsosNode<DIM>::eraseLayerNode()
{
   if ( d_common->rank == d_owner ) {
      d_common->new_node_set.eraseNode(d_node_iterator);
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   d_node_iterator = GraphNodeContainer().end();
   d_node = GraphNode();
#endif
   return;
}


template<int DIM>
void AsyncBergerRigoutsosNode<DIM>::countOverlapWithLocalPatches()
{
   /*
    * Count overlaps for the left and right sides.
    *
    * Remove the child if it has zero overlap.
    */
   hier::Box<DIM> lft_grown_box = d_lft_child->d_box;
   lft_grown_box.grow( d_common->max_gcw );
   hier::Box<DIM> rht_grown_box = d_rht_child->d_box;
   rht_grown_box.grow( d_common->max_gcw );
   int &lft_overlap = d_lft_child->d_overlap;
   int &rht_overlap = d_rht_child->d_overlap;
   lft_overlap = rht_overlap = 0;

   hier::PatchLevel<DIM> &level = *d_common->level;
   for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {

      const hier::Box<DIM> &patch_box = level.getPatch(ip())->getBox();

      hier::Box<DIM> lft_intersection = patch_box * lft_grown_box;
      lft_overlap += lft_intersection.size();

      hier::Box<DIM> rht_intersection = patch_box * rht_grown_box;
      rht_overlap += rht_intersection.size();

   }
   return;
}





/*
*************************************************************************
* Child groups are subsets of current group.  Each child group          *
* includes processes owning patches that overlap the box of that child. *
* The overlap data has been gathered in d_recv_msg.                     *
* See gatherGroupingCriteria_start() for the format of the message.     *
*************************************************************************
*/
template<int DIM>
void AsyncBergerRigoutsosNode<DIM>::formChildGroups()
{
   /*
    * Form child groups and determine owners from data gathered
    * in the gather_overlap_counts phase.
    */
   if ( d_group.size() == 1 ) {
      // Short cut for trivial groups.
      d_lft_child->d_group.setNull();
      d_rht_child->d_group.setNull();
      d_lft_child->d_group.resizeArray(1);
      d_rht_child->d_group.resizeArray(1);
      d_lft_child->d_group[0] = d_group[0];
      d_rht_child->d_group[0] = d_group[0];
      d_lft_child->d_owner = d_owner;
      d_rht_child->d_owner = d_owner;
      return;
   }

   d_lft_child->d_group.resizeArray( d_group.size() );
   d_rht_child->d_group.resizeArray( d_group.size() );

#ifdef DEBUG_CHECK_ASSERTIONS
   /*
    * Only owner process should be here.
    */
   if ( d_common->rank != d_owner ) {
      TBOX_ERROR("Library error!");
   }
   TBOX_ASSERT( d_recv_msg.size() == 4*d_group.size() );
#endif


   int *lft_overlap = &d_recv_msg[0];
   int *rht_overlap = &d_recv_msg[1];


   const int imyself = findOwnerInGroup( d_common->rank, d_group );
   lft_overlap[imyself*4] = d_lft_child->d_overlap;
   rht_overlap[imyself*4] = d_rht_child->d_overlap;


   int *lft_criteria = NULL;
   int *rht_criteria = NULL;
   switch ( d_common->owner_mode ) {
   case SINGLE_OWNER:
      lft_criteria = &d_recv_msg[0];
      rht_criteria = &d_recv_msg[1];
      lft_criteria[imyself*4] = tbox::MathUtilities<int>::getMax();
      rht_criteria[imyself*4] = tbox::MathUtilities<int>::getMax();
      break;
   case MOST_OVERLAP:
      lft_criteria = &d_recv_msg[0];
      rht_criteria = &d_recv_msg[1];
      lft_criteria[imyself*4] = d_lft_child->d_overlap;
      rht_criteria[imyself*4] = d_rht_child->d_overlap;
      break;
   case FEWEST_OWNED:
      lft_criteria = &d_recv_msg[2];
      rht_criteria = &d_recv_msg[2];
      lft_criteria[imyself*4] = -d_common->num_nodes_owned;
      rht_criteria[imyself*4] = -d_common->num_nodes_owned;
      break;
   case LEAST_ACTIVE:
      lft_criteria = &d_recv_msg[3];
      rht_criteria = &d_recv_msg[3];
      lft_criteria[imyself*4] = -d_common->num_nodes_active;
      rht_criteria[imyself*4] = -d_common->num_nodes_active;
      break;
   default:
      TBOX_ERROR("LIBRARY error");
      break;
   }

   int n_lft=0;
   int n_rht=0;

   int lft_owner_score = tbox::MathUtilities<int>::getMin();
   int rht_owner_score = tbox::MathUtilities<int>::getMin();


   /*
    * Loop through the group to see which process should participate
    * on the left/right sides.  Also see which process should be the
    * owner of the left/right sides.  For efficiency in some searches
    * through d_groups, make sure that d_group is ordered.
    */
   for ( int i=0; i<d_group.size(); ++i ) {
      int i4 = i*4;
      if ( lft_overlap[i4] != 0 ) {
         d_lft_child->d_group[n_lft++] = d_group[i];
         if ( lft_criteria[i4] > lft_owner_score ) {
            d_lft_child->d_owner = d_group[i];
            lft_owner_score = lft_criteria[i4];
         }
      }
      if ( rht_overlap[i4] != 0 ) {
         d_rht_child->d_group[n_rht++] = d_group[i];
         if ( rht_criteria[i4] > rht_owner_score ) {
            d_rht_child->d_owner = d_group[i];
            rht_owner_score = rht_criteria[i4];
         }
      }
   }

   d_lft_child->d_group.resizeArray(n_lft);
   d_rht_child->d_group.resizeArray(n_rht);

#ifdef DEBUG_CHECK_ASSERTIONS
   // Recall that only the owner should execute this code.
   TBOX_ASSERT( d_lft_child->d_owner >= 0 );
   TBOX_ASSERT( d_lft_child->d_group.size() > 0 );
   TBOX_ASSERT( d_lft_child->d_group.size() <= d_group.size() );
   TBOX_ASSERT( d_common->owner_mode == SINGLE_OWNER ||
                ( ( d_lft_child->d_overlap == 0 ) !=
                  inGroup(d_lft_child->d_group) ) );
   TBOX_ASSERT( d_rht_child->d_owner >= 0 );
   TBOX_ASSERT( d_rht_child->d_group.size() > 0 );
   TBOX_ASSERT( d_rht_child->d_group.size() <= d_group.size() );
   TBOX_ASSERT( d_common->owner_mode == SINGLE_OWNER ||
                ( ( d_rht_child->d_overlap == 0 ) !=
                  inGroup(d_rht_child->d_group) ) );
   if ( d_common->owner_mode == SINGLE_OWNER ) {
      TBOX_ASSERT( inGroup(d_lft_child->d_group, d_owner) );
      TBOX_ASSERT( inGroup(d_rht_child->d_group, d_owner) );
   }
   for ( int i=0; i<d_group.size(); ++i ) {
      TBOX_ASSERT( i==0 || d_group[i] > d_group[i-1] );
      TBOX_ASSERT( ( lft_overlap[i*4] > 0 ||
                     ( d_group[i] == d_lft_child->d_owner ) )
                   == inGroup(d_lft_child->d_group, d_group[i]) );
      TBOX_ASSERT( ( rht_overlap[i*4] > 0 ||
                     ( d_group[i] == d_rht_child->d_owner ) )
                   == inGroup(d_rht_child->d_group, d_group[i]) );
   }
#endif

   return;
}




/*
***********************************************************************
* Form child groups and determine owners, taking advantage            *
* of level box data.  If the DNBG approach for box overlap            *
* computation is used, the levels may not have all its boxes          *
* stored locally, in which case, this method cannot be used.          *
*                                                                     *
* Each child group contains the processes that own patches            *
* overlaps with the box of that child.  The process having            *
* the greatest overlap becomes the owner of the child (except         *
* in single-owner mode).                                              *
***********************************************************************
*/
template<int DIM>
void AsyncBergerRigoutsosNode<DIM>::formChildGroupsUsingLevelBoxes()
{
   const hier::BoxArray<DIM> &level_boxes = d_common->level->getBoxes();
   const hier::ProcessorMapping &proc_map = d_common->level->getProcessorMapping();
   tbox::Pointer<hier::BoxTree<DIM> > box_tree = d_common->level->getBoxTree();

   int i;
   tbox::Array<int> overlap_box_indices;
   tbox::Array<int> overlap_counts( d_common->nproc );
   int owner_overlap;
   int group_size;
   hier::Box<DIM> grown_box;


   // Left side.

   for ( i=0; i<d_common->nproc; ++i ) overlap_counts[i] = 0;
   grown_box = d_lft_child->d_box;
   grown_box.grow(d_common->max_gcw);
   box_tree->findOverlapIndices( overlap_box_indices, grown_box );
   d_lft_child->d_overlap = 0;
   for ( i=0; i<overlap_box_indices.size(); ++i ) {
      int process = proc_map.getProcessorAssignment( overlap_box_indices[i] );
      const hier::Box<DIM> box = level_boxes( overlap_box_indices[i] );
      const hier::Box<DIM> intersection = box * grown_box;
      overlap_counts[process] += intersection.size();
   }
   d_lft_child->d_overlap = overlap_counts[d_common->rank];
   if ( d_common->owner_mode == SINGLE_OWNER ) {
      overlap_counts[d_owner] = tbox::MathUtilities<int>::getMax();
   }
   d_lft_child->d_group.setNull();
   d_lft_child->d_group.resizeArray(d_group.size());
   group_size = 0;
   owner_overlap = 0;
   for ( i=0; i<d_common->nproc; ++i ) {
      if ( overlap_counts[i] > 0 ) {
         d_lft_child->d_group[group_size++] = i;
         if ( overlap_counts[i] > owner_overlap ) {
            d_lft_child->d_owner = i;
            owner_overlap = overlap_counts[i];
         }
      }
   }
   d_lft_child->d_group.resizeArray( group_size );
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( d_lft_child->d_owner >= 0 );
   TBOX_ASSERT( d_lft_child->d_group.size() > 0 );
   TBOX_ASSERT( d_lft_child->d_group.size() <= d_group.size() );
   TBOX_ASSERT( d_common->owner_mode == SINGLE_OWNER ||
                ( ( d_lft_child->d_overlap == 0 ) !=
                  inGroup(d_lft_child->d_group) ) );
   if ( d_common->owner_mode == SINGLE_OWNER ) {
      TBOX_ASSERT( inGroup(d_lft_child->d_group, d_owner) );
   }
   for ( i=0; i<d_group.size(); ++i ) {
      TBOX_ASSERT( i==0 || d_group[i] > d_group[i-1] );
      TBOX_ASSERT( ( overlap_counts[d_group[i]] > 0 ||
                     ( d_group[i] == d_lft_child->d_owner ) )
                   == inGroup(d_lft_child->d_group, d_group[i]) );
   }
#endif


   // Right side.

   for ( i=0; i<d_common->nproc; ++i ) overlap_counts[i] = 0;
   grown_box = d_rht_child->d_box;
   grown_box.grow(d_common->max_gcw);
   box_tree->findOverlapIndices( overlap_box_indices, grown_box );
   d_rht_child->d_overlap = 0;
   for ( i=0; i<overlap_box_indices.size(); ++i ) {
      int process = proc_map.getProcessorAssignment( overlap_box_indices[i] );
      const hier::Box<DIM> box = level_boxes( overlap_box_indices[i] );
      const hier::Box<DIM> intersection = box * grown_box;
      overlap_counts[process] += intersection.size();
   }
   d_rht_child->d_overlap = overlap_counts[d_common->rank];
   if ( d_common->owner_mode == SINGLE_OWNER ) {
      overlap_counts[d_owner] = tbox::MathUtilities<int>::getMax();
   }
   d_rht_child->d_group.setNull();
   d_rht_child->d_group.resizeArray(d_group.size());
   group_size = 0;
   owner_overlap = 0;
   for ( i=0; i<d_common->nproc; ++i ) {
      if ( overlap_counts[i] > 0 ) {
         d_rht_child->d_group[group_size++] = i;
         if ( overlap_counts[i] > owner_overlap ) {
            d_rht_child->d_owner = i;
            owner_overlap = overlap_counts[i];
         }
      }
   }
   d_rht_child->d_group.resizeArray( group_size );
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( d_rht_child->d_owner >= 0 );
   TBOX_ASSERT( d_rht_child->d_group.size() > 0 );
   TBOX_ASSERT( d_rht_child->d_group.size() <= d_group.size() );
   TBOX_ASSERT( d_common->owner_mode == SINGLE_OWNER ||
                ( ( d_rht_child->d_overlap == 0 ) !=
                  inGroup(d_rht_child->d_group) ) );
   if ( d_common->owner_mode == SINGLE_OWNER ) {
      TBOX_ASSERT( inGroup(d_rht_child->d_group, d_owner) );
   }
   for ( i=0; i<d_group.size(); ++i ) {
      TBOX_ASSERT( i==0 || d_group[i] > d_group[i-1] );
      TBOX_ASSERT( ( overlap_counts[d_group[i]] > 0 ||
                     ( d_group[i] == d_rht_child->d_owner ) )
                   == inGroup(d_rht_child->d_group, d_group[i]) );
   }
#endif

   return;
}




/*
*************************************************************************
*                                                                       *
* Compute overlaps between the new graph node and nodes on              *
* the tagged level, saving that data in the form of edges.              *
*                                                                       *
* Note that the edge data may be duplicated in two objects.             *
* - tag_cnect_new stores the edges organized around each node           *
*   in the tagged level.  For each node on the tagged level,            *
*   we store a container of neighbors on the new layer.                 *
* - new_cnect_tags stores the edges organized around each NEW node.     *
*   For each new node we store a container of neighbors on the          *
*   tagged level.                                                       *
*                                                                       *
* If compute_edges > 0, we store tag_cnect_new.                         *
*                                                                       *
* If compute_edges > 1, we also compute new_cnect_tags.                 *
* The data in new_cnect_tags are                                        *
* computed by the particant processes but eventually stored on the      *
* owners of the new nodes, so their computation requires caching        *
* the edge data in edge_messages for sending to the appropriate         *
* processes later.                                                      *
*                                                                       *
*************************************************************************
*/
template<int DIM>
void AsyncBergerRigoutsosNode<DIM>::computeNewGraphEdges()
{
   d_common->t_compute_new_graph_edges->start();
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( d_common->compute_edges > 0 );
   TBOX_ASSERT( d_node.getLocalIndex() >= 0 );
   TBOX_ASSERT( boxAccepted() );
   TBOX_ASSERT( d_box_acceptance != accepted_by_dropout_bcast );
   TBOX_ASSERT( d_box.numberCells() >= d_common->min_box );
   /*
    * We should not compute nabrs if we got the node
    * by a dropout broadcast because we already know
    * there is no overlap!
    */
   TBOX_ASSERT( d_box_acceptance != accepted_by_dropout_bcast );
#endif

   GraphNabrContainer nabr_container;

   // Create an expanded box for intersection check.
   hier::Box<DIM> grown_box = d_box;
   grown_box.grow( d_common->max_gcw );


   /*
    * On the owner process, we store the neighbors of the new node.
    * This data is NOT required on other processes.
    */
   GraphNabrContainer *nabrs_of_new_node = NULL;
   if ( d_common->rank == d_owner ) {
      nabrs_of_new_node = &(d_common->new_cnect_tag[d_node.getLocalIndex()]);
   }

   // Data to send to d_owner regarding new edges found by local process.
  std::vector<int> *edge_message = NULL;
   if ( d_common->compute_edges > 1 && d_common->rank != d_owner ) {
      /*
       * Will have to send to d_owner the edges found locally for
       * graph node d_node.
       * Label the id of the new node and the (yet unknown) number
       * of edge found for it.
       *
       * The message to be sent to owner is appended the following
       * data:
       * - index of new node
       * - number of edges found for the new node
       * - index of nodes on the tagged level overlapping new node.
       */
      edge_message = &d_common->edge_messages[d_owner];
      edge_message->insert( edge_message->end(), d_node.getLocalIndex() );
      edge_message->insert( edge_message->end(), 0 );
   }


   const int index_of_counter =
      ( edge_message != NULL ? edge_message->size() : 0 ) - 1;
   const int ints_per_node = GraphNode::commBufferSize();


   const GraphNodeContainer &nodes_on_tagged_level =
      ((const hier::LayerNodeSet<DIM> &)d_common->tag_node_set).getNodeContainer();

   for ( typename GraphNodeContainer::const_iterator ni=nodes_on_tagged_level.begin();
         ni!=nodes_on_tagged_level.end(); ++ni ) {

      const GraphNode &tag_node = *ni;

      hier::Box<DIM> intersection = tag_node.getBox() * grown_box;

      if ( ! intersection.empty() ) {

         // Add d_node as a neighbor of tag_node.
         GraphNabrContainer &new_nabrs_of_tag_node =
            d_common->tag_cnect_new[tag_node.getLocalIndex()];
         new_nabrs_of_tag_node.insert( d_node );

         if ( nabrs_of_new_node != NULL ) {
            // Owner adds tag_node as a neighbor of d_node.
            (*nabrs_of_new_node).insert( tag_node );
         }

         if ( edge_message != NULL ) {
            /* Non-owners put found edge in the message
             * to (eventually) send to d_owner.
             */
            edge_message->insert( edge_message->end(), ints_per_node, 0 );
            int *ptr = &(*edge_message)[edge_message->size()-ints_per_node];
            tag_node.putToIntBuffer(ptr);
            ++(*edge_message)[index_of_counter];
         }

      }
   }

   if ( d_common->compute_edges > 1 &&
        d_common->rank == d_owner ) {
      /*
       * If box was accepted, the owner should remember
       * which process will be sending edge data.
       * Update the list of edge senders to make sure
       * it includes all processes in the group.
       * We use this list in shareNewEdgesWithOwners to
       * tell us which processors are sending us new edges.
       * The edge senders are the participants of the group.
       */
      d_common->edge_senders.insert( d_group.getPointer(),
                                     d_group.getPointer()+d_group.size() );
   }

   d_common->t_compute_new_graph_edges->stop();
   return;
}




/*
**********************************************************************
*                                                                    *
* Send new edges found by local process to owners of the new nodes   *
* associated with those edges.  Receive similar data from other      *
* processes.                                                         *
*                                                                    *
* Messages to be sent out were placed in d_common->edge_messages by  *
* computeNewGraphEdges().  This method sends out these messages      *
* and receives anticipated messages from processes listed in         *
* d_common->edge_senders.  Received messages are unpacked to get     *
* data on new edges.                                                 *
*                                                                    *
**********************************************************************
*/
template<int DIM>
void AsyncBergerRigoutsosNode<DIM>::shareNewEdgesWithOwners()
{
#if defined(HAVE_MPI)
   d_common->t_share_new_edges->start();
   IntSet edge_senders = d_common->edge_senders;
   std::map<int,std::vector<int> > &edge_messages = d_common->edge_messages;
   Connectivity &new_cnect_tag = d_common->new_cnect_tag;

   const int ints_per_node = GraphNode::commBufferSize();

   int ierr;
   tbox::SAMRAI_MPI::status mpi_status;




   // Nonblocking send of edge data.
   d_common->t_share_new_edges_send->start();
   tbox::Array<tbox::SAMRAI_MPI::request> mpi_request(edge_messages.size());
   std::map<int,std::vector<int> >::iterator send_i;
   int nsend = 0;
   for ( send_i = edge_messages.begin(), nsend = 0;
         send_i != edge_messages.end();
         ++send_i, ++nsend ) {
      const int &owner = (*send_i).first;
      std::vector<int> &msg = (*send_i).second;
      ierr = MPI_Isend( &msg[0],
                        msg.size(),
                        MPI_INT,
                        owner,
                        d_common->tag_upper_bound,
                        d_common->mpi_communicator,
                        &mpi_request[nsend] );
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT( ierr == MPI_SUCCESS );
#else
      NULL_USE(ierr);
#endif
   }
   d_common->t_share_new_edges_send->stop();

   {
      /*
       * The rest of this method assumes current process is NOT
       * in edge_senders, so remove it.  For efficiency, method
       * computeNewGraphEdges() (which created the edge senders)
       * did not remove it.
       */
      IntSet::iterator local = edge_senders.find(d_common->rank);
      if ( local != edge_senders.end() ) {
         edge_senders.erase(local);
      }
   }

   /*
    * Create set recved_from which is to contain ranks of
    * processes from which we've received the expected edge data.
    * The while loop goes until all expected messages have
    * been received from edge_senders.
    *
    * In the while loop:
    *    - Probe for an incomming message.
    *    - Determine its size allocate memory for receiving the message.
    *    - Receive the message.
    *    - Get edge data from the message.
    */
   IntSet recved_from;
   while ( recved_from.size() < edge_senders.size() ) {

      d_common->t_share_new_edges_recv->start();
      ierr = MPI_Probe( MPI_ANY_SOURCE,
                        d_common->tag_upper_bound,
                        d_common->mpi_communicator,
                        &mpi_status );
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT( ierr == MPI_SUCCESS );
#endif

      const int sender = mpi_status.MPI_SOURCE;
      int mesg_size = -1;
      MPI_Get_count( &mpi_status, MPI_INT, &mesg_size );
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT( edge_senders.find(sender) != edge_senders.end() );
      TBOX_ASSERT( recved_from.find(sender) == recved_from.end() );
      TBOX_ASSERT( mesg_size >= 0 );
#endif

      tbox::Array<int> buf(mesg_size);
      int *ptr = buf.getPointer();
      ierr = MPI_Recv( ptr,
                       mesg_size,
                       MPI_INT,
                       sender,
                       d_common->tag_upper_bound,
                       d_common->mpi_communicator,
                       &mpi_status );
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT( ierr == MPI_SUCCESS );
#endif
      d_common->t_share_new_edges_recv->stop();

      d_common->t_share_new_edges_unpack->start();
      int consumed = 0;
      while ( ptr < buf.getPointer() + buf.size() ) {
         int new_node_id = *(ptr++);
         int n_new_edges = *(ptr++);
         typename Connectivity::iterator nn = new_cnect_tag.find(new_node_id);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT( nn != new_cnect_tag.end() );
#endif
         GraphNabrContainer &nabr_container = (*nn).second;
         for ( int n=0; n<n_new_edges; ++n ) {
            GraphNode node;
            node.getFromIntBuffer(ptr);
            ptr += ints_per_node;
            nabr_container.insert(node);
         }
         consumed += 2 + n_new_edges*ints_per_node;
      }
      recved_from.insert(sender);
      d_common->t_share_new_edges_unpack->stop();
   }

   if ( nsend > 0 ) {
      // Make sure all nonblocking sends completed.
      d_common->t_share_new_edges_send->start();
      tbox::Array<tbox::SAMRAI_MPI::status> mpi_statuses(edge_messages.size());
      ierr = MPI_Waitall( edge_messages.size(),
                          mpi_request.getPointer(),
                          mpi_statuses.getPointer() );
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT( ierr == MPI_SUCCESS );
#endif
      d_common->t_share_new_edges_send->stop();
   }


#if 0
   // This is an expensive check, for serious debugging only!
   plog << "Edge sharing out/in: "
        << edge_messages.size() << ' ' << edge_senders.size() << "\n";

   LayerEdgeSet<DIM> tag_to_new;
   LayerEdgeSet<DIM> new_to_tag;
   tag_to_new.initialize( LayerEdgeSet<DIM>::DISTRIBUTED,
                          d_common->tag_node_set,
                          d_common->tag_node_set.getRefinementRatio(),
                          d_common->max_gcw,
                          &d_common->tag_cnect_new );
   new_to_tag.initialize( LayerEdgeSet<DIM>::DISTRIBUTED,
                          d_common->new_node_set,
                          d_common->new_node_set.getRefinementRatio(),
                          d_common->max_gcw,
                          &d_common->new_cnect_tag );

   tbox::plog << "======= Tagged layer edge set =======\n";
   tag_to_new.printClassData(tbox::plog, 2);
   tbox::plog << "======= New layer edge set =======\n";
   new_to_tag.printClassData(tbox::plog, 2);
   tbox::plog << "======= Checking graph accuracy =======\n";

   tag_to_new.checkConnectivity( d_common->new_node_set );
   new_to_tag.checkConnectivity( d_common->tag_node_set );
#endif

   d_common->t_share_new_edges->stop();

#endif
   return;
}



/*
********************************************************************
* Utility methods.                                                 *
********************************************************************
*/


template<int DIM>
tbox::List<tbox::RelaunchableJob*>::Iterator
AsyncBergerRigoutsosNode<DIM>::inRelaunchQueue(
   tbox::RelaunchableJob *node_ptr ) const
{
   tbox::List<tbox::RelaunchableJob*>::Iterator
      li(d_common->job_relauncher.getRelaunchQueue());
   if ( li ) {
      for ( ; li; li++ ) {
         if ( node_ptr == *li ) {
            break;
         }
      }
   }
   return li;
}




template<int DIM>
int AsyncBergerRigoutsosNode<DIM>::getHistogramBufferSize(
   const hier::Box<DIM> &box ) const
{
   int d, size=box.numberCells(0);
   for ( d=1; d<DIM; ++d ) {
      size += box.numberCells(d);
   }
   return size;
}



template<int DIM>
int *AsyncBergerRigoutsosNode<DIM>::putHistogramToBuffer(
   int *buffer )
{
   for ( int d=0; d<DIM; ++d ) {
      d_histogram[d].resizeArray( d_box.numberCells(d) );
      memcpy( buffer,
              d_histogram[d].getPointer(),
              d_box.numberCells(d)*sizeof(int) );
      buffer += d_box.numberCells(d);
   }
   return buffer;
}



template<int DIM>
int *AsyncBergerRigoutsosNode<DIM>::getHistogramFromBuffer(
   int *buffer )
{
   for ( int d=0; d<DIM; ++d ) {
      d_histogram[d].resizeArray( d_box.numberCells(d) );
      memcpy( d_histogram[d].getPointer(),
              buffer,
              d_box.numberCells(d)*sizeof(int) );
      buffer += d_box.numberCells(d);
   }
   return buffer;
}




template<int DIM>
int *AsyncBergerRigoutsosNode<DIM>::putBoxToBuffer(
   const hier::Box<DIM> &box, int *buffer ) const
{
   const hier::IntVector<DIM> &l = box.lower();
   const hier::IntVector<DIM> &u = box.upper();
   int d;
   for ( d=0; d<DIM; ++d ) {
      *(buffer++) = l(d);
      *(buffer++) = u(d);
   }
   return buffer;
}




template<int DIM>
int *AsyncBergerRigoutsosNode<DIM>::getBoxFromBuffer(
   hier::Box<DIM> &box, int *buffer ) const
{
   hier::IntVector<DIM> &l = box.lower();
   hier::IntVector<DIM> &u = box.upper();
   int d;
   for ( d=0; d<DIM; ++d ) {
      l(d) = *(buffer++);
      u(d) = *(buffer++);
   }
   return buffer;
}




/*
***********************************************************************
* Put in dropouts things that are in main_group but                   *
* not in sub_group.                                                   *
*                                                                     *
* Assume that sub_group is a subset of elements in main_group.        *
* Assume that sub_group and main_group are sorted in ascending order. *
*                                                                     *
* Assume add_root is NOT in the dropout and add it anyway.            *
***********************************************************************
*/
template<int DIM>
void AsyncBergerRigoutsosNode<DIM>::computeDropoutGroup(
   const tbox::Array<int> &main_group,
   const tbox::Array<int> &sub_group,
   tbox::Array<int> &dropout_group,
   int add_root ) const
{

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( main_group.size() >= sub_group.size() );
#endif

   dropout_group.setNull();
   dropout_group.resizeArray (main_group.size());

   int i, j, k=0;
   dropout_group[k++] = add_root;
   for ( i=0, j=0; i<main_group.size(); ++i ) {
      if ( main_group[i] != sub_group[j] ) {
         dropout_group[k++] = main_group[i];
      }
      else {
         ++j;
         if ( j == sub_group.size() ) {
            // No more in the sub_group so the rest of main_group
            // goes in dropout_group.
            for ( i=i+1; i<main_group.size(); ++i, ++k ) {
               dropout_group[k] = main_group[i];
            }
         }
      }
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( j = sub_group.size() );
#endif
   dropout_group.resizeArray(k);
   return;
}




/*
**********************************************************************
* Determine if the given rank is in the given group.                 *
**********************************************************************
*/
template<int DIM>
bool AsyncBergerRigoutsosNode<DIM>::inGroup( tbox::Array<int> &group,
                                              int rank ) const
{
   if ( rank < 0 ) rank = d_common->rank;
   int i;
   for ( i=0; i<group.size(); ++i ) {
      if ( rank == group[i] ) {
         return true;
      }
   }
   return false;
}




/*
**********************************************************************
* Heuristically determine the "best" tree degree for                 *
* the communication group.                                           *
* Use binary tree for less than 2^4 processes.                       *
* After that, for each 2^3 multipling of group size,                 *
* increase tree degree by 1.                                         *
**********************************************************************
*/
template<int DIM>
int AsyncBergerRigoutsosNode<DIM>::computeCommunicationTreeDegree(
   int group_size ) const
{
   int tree_deg = 2;
   int shifted_size = group_size >> 3;
   while ( shifted_size > 0 ) {
      shifted_size >>= 3;
      ++tree_deg;
   }
   return tree_deg;
}




template<int DIM>
int AsyncBergerRigoutsosNode<DIM>::findOwnerInGroup(
   int owner, const tbox::Array<int> &group ) const
{
   for ( int i=0; i<group.size(); ++i ) {
      if ( group[i] == owner ) return i;
   }
   return -1;
}
/*
**********************************************************************
* Claim a unique tag from the processor's available tag pool.        *
* Check that the pool is not overused.                               *
**********************************************************************
*/
template<int DIM>
void AsyncBergerRigoutsosNode<DIM>::claimMPITag()
{
#ifdef DEBUG_CHECK_ASSERTIONS
   /*
    * Each dendogram node should claim no more than one MPI tag
    * so make sure it does not already have one.
    */
   TBOX_ASSERT( d_mpi_tag < 0 );
#endif
   d_mpi_tag = d_common->available_mpi_tag;
   d_common->available_mpi_tag = d_mpi_tag + total_phase_tags;
   if ( d_mpi_tag+total_phase_tags-1 >
        d_common->tag_upper_bound/(d_common->nproc)*(d_common->rank+1) ) {
      /*
       * Each process is alloted tag_upper_bound/(d_common->nproc)
       * tag values.  If it needs more than this, it will encroach
       * on the tag pool of the next process and may lead to using
       * non-unique tags.
       */
      TBOX_ERROR("Out of MPI tag values need to ensure that\n"
                 <<"messages are properly differentiated."
                 <<"\nd_mpi_tag = " << d_mpi_tag
                 <<"\ntag_upper_bound = " << d_common->tag_upper_bound
                 <<"\nmber of nodes = " << d_common->nproc
                 <<"\nmax tag required = " << d_mpi_tag+total_phase_tags-1
                 <<"\nmax tag available = "
                 << d_common->tag_upper_bound/(d_common->nproc)*(d_common->rank+1));
      /*
       * It is probably safe to recycle tags if we run out of MPI tags.
       * This is not implemented because thus far, there is no need for it.
       * Recycling is starting over from the initial tag set aside for the
       * local process.  To make sure that recycled tags are not still
       * in use, we should claim a new (or recycled) tag for the dropout
       * broadcast phase.  This is because descendant nodes may recycle
       * the current claimed tag before this phase starts.  All other
       * phases are not interupted by descendant communications, so we
       * are assured that their tag is not doubly claimed.
       */
   }
   return;
}




/*
**********************************************************************
* Convert an integer value to BoxAcceptance.                         *
* This is needed because the compiler cannot                         *
* cast an integer to an enum type.                                   *
**********************************************************************
*/
template<int DIM>
typename AsyncBergerRigoutsosNode<DIM>::BoxAcceptance
AsyncBergerRigoutsosNode<DIM>::intToBoxAcceptance( int i ) const
{
   switch (i) {
   case undetermined              : return undetermined;
   case rejected_by_calculation   : return rejected_by_calculation;
   case accepted_by_calculation   : return accepted_by_calculation;
   case rejected_by_owner         : return rejected_by_owner;
   case accepted_by_owner         : return accepted_by_owner;
   case rejected_by_recombination : return rejected_by_recombination;
   case accepted_by_recombination : return accepted_by_recombination;
   case rejected_by_dropout_bcast : return rejected_by_dropout_bcast;
   case accepted_by_dropout_bcast : return accepted_by_dropout_bcast;
   default:
      TBOX_ERROR("Library error: bad BoxAcceptance data of " << i << ".\n");
   }
   return undetermined;
}




template<int DIM>
void AsyncBergerRigoutsosNode<DIM>::printClassData( std::ostream &os,
                                                     int detail_level ) const
{
   os << "ID              " << d_pos << " owner=" << d_owner << " box=" << d_box
      ;
   if ( detail_level > 0 ) {
      os << "\nfamily          " << (d_parent==NULL?0:d_parent->d_pos)
                                 << ' ' << (d_lft_child?(d_lft_child->d_pos):-1)
                                 << ' ' << (d_rht_child?(d_rht_child->d_pos):-1)
         ;
   }
   if ( detail_level > 1 ) {
      os << "\nthis           " << this
         << "\ngeneration     " << d_generation << " place=" << ((d_pos%2)?'r':'l')
         << "\nnode           " << d_node
         << "\nbox_acceptance " << ( d_box_acceptance ? d_box_acceptance : '-' )
         << "\noverlap        " << d_overlap
         << "\ngroup          " << d_group.size() << ':'
         ;
      int i;
      for ( i=0; i<d_group.size(); ++i ) { os << ' ' << d_group[i]; }
      os << std::endl;
   }
   return;
}






/*
**********************************************************************
*                                                                    *
* Methods for setting algorithm parameters before running.           *
*                                                                    *
**********************************************************************
*/

template<int DIM>
void AsyncBergerRigoutsosNode<DIM>::setMaxGhostCellWidth(
   const hier::IntVector<DIM> &max_gcw )
{
   d_common->max_gcw = max_gcw;
   return;
}


template<int DIM>
void AsyncBergerRigoutsosNode<DIM>::setOwnerMode( const std::string &mode )
{
   if ( mode == "SINGLE_OWNER" ) {
     d_common->owner_mode = SINGLE_OWNER;
   }
   else if ( mode == "MOST_OVERLAP" ) {
     d_common->owner_mode = MOST_OVERLAP;
   }
   else if ( mode == "FEWEST_OWNED" ) {
     d_common->owner_mode = FEWEST_OWNED;
   }
   else if ( mode == "LEAST_ACTIVE" ) {
     d_common->owner_mode = LEAST_ACTIVE;
   }
   else {
     TBOX_ERROR("AsyncBergerRigoutsosNode: Unrecognized owner mode request: "
                << mode);
   }
   return;
}




template<int DIM>
void AsyncBergerRigoutsosNode<DIM>::setUseLevelBoxes( bool flag )
{
   d_common->use_level_boxes = flag;
   return;
}


template<int DIM>
void AsyncBergerRigoutsosNode<DIM>::setComputeEdges( int flag )
{
   d_common->compute_edges = flag;
   return;
}


template<int DIM>
void AsyncBergerRigoutsosNode<DIM>::setLogNodeHistory( bool flag )
{
   d_common->log_node_history = flag;
   return;
}








/*
**********************************************************************
*                                                                    *
* Methods for collecting data on dendogram after it completes.       *
*                                                                    *
**********************************************************************
*/

template<int DIM>
int AsyncBergerRigoutsosNode<DIM>::getMaxNodes() const
{
   return d_common->max_nodes_allocated;
}

template<int DIM>
int AsyncBergerRigoutsosNode<DIM>::getMaxGeneration() const
{
   return d_common->max_generation;
}

template<int DIM>
int AsyncBergerRigoutsosNode<DIM>::getMaxOwnership() const
{
   return d_common->max_nodes_owned;
}

template<int DIM>
double AsyncBergerRigoutsosNode<DIM>::getAvgNumberOfCont() const
{
   if ( d_common->num_nodes_completed > 0 ) {
      return (double)d_common->num_conts_to_complete/d_common->num_nodes_completed;
   }
   return 0;
}

template<int DIM>
int AsyncBergerRigoutsosNode<DIM>::getMaxNumberOfCont() const
{
   return d_common->max_conts_to_complete;
}

template<int DIM>
int AsyncBergerRigoutsosNode<DIM>::getNumBoxesGenerated() const
{
   return d_common->num_boxes_generated;
}




/*
**********************************************************************
*                                                                    *
* Methods for accessing output from BR algorithm.                    *
*                                                                    *
**********************************************************************
*/

template<int DIM>
const hier::LayerNodeSet<DIM> &AsyncBergerRigoutsosNode<DIM>::getNewNodes() const
{
#ifdef DEBUG_CHECK_ASSERTION
   if ( d_wait_phase != completed ) {
      TBOX_ERROR("Cannot get results until algorithm completes running.");
   }
#endif
   return d_common->new_node_set;
}



template<int DIM>
const typename hier::LayerEdgeSet<DIM>::Connectivity
&AsyncBergerRigoutsosNode<DIM>::getNewCnect() const
{
#ifdef DEBUG_CHECK_ASSERTION
   if ( d_wait_phase != completed ) {
      TBOX_ERROR("Cannot get results until algorithm completes running.");
   }
#endif
  return d_common->tag_cnect_new;
}



/*
**********************************************************************
*                                                                    *
* Implementations for pure virtual methods defined by                *
* tbox::RelaunchableJob.                                             *
*                                                                    *
**********************************************************************
*/

template<int DIM>
void AsyncBergerRigoutsosNode<DIM>::continueJob() {
   continueAlgorithm();
   return;
}

template<int DIM>
tbox::RelaunchableJob::JobState AsyncBergerRigoutsosNode<DIM>::getJobState() {
   tbox::RelaunchableJob::JobState job_state = JOB_IS_COMPLETED;
   switch (d_wait_phase) {
   case reduce_histogram:
   case bcast_acceptability:
   case gather_grouping_criteria:
   case bcast_child_groups:
   case bcast_to_dropouts:
      job_state = COMMUNICATION_WAIT;
      break;
   case to_be_launched:
   case run_children:
      job_state = NONCOMMUNICATION_WAIT;
      break;
   case completed:
      job_state = JOB_IS_COMPLETED;
      break;
   default:
      TBOX_ERROR("Library error");
   }
   return job_state;
}

template<int DIM>
tbox::AsyncCommGroup *AsyncBergerRigoutsosNode<DIM>::getCommunicationGroup() {
   return d_comm_group;
}


}
}
#endif
