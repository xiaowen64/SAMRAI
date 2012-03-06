/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Node in asynchronous Berger-Rigoutsos dendogram
 *
 ************************************************************************/

namespace SAMRAI {
namespace mesh {

SAMRAI_INLINE_KEYWORD
void
BergerRigoutsosNode::setClusteringParameters(
   const int tag_data_index,
   const int tag_val,
   const hier::IntVector min_box,
   const double efficiency_tol,
   const double combine_tol,
   const hier::IntVector& max_box_size,
   const double max_lap_cut_from_center,
   const bool laplace_cut_long_dir_only)
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS2(d_dim, min_box, max_box_size);

   d_common->tag_data_index = tag_data_index;
   d_common->tag_val = tag_val;
   d_common->min_box = min_box;
   d_common->efficiency_tol = efficiency_tol;
   d_common->combine_tol = combine_tol;
   d_common->max_box_size = max_box_size;
   d_common->max_lap_cut_from_center = max_lap_cut_from_center;
   d_common->laplace_cut_long_dir_only = laplace_cut_long_dir_only;
}

SAMRAI_INLINE_KEYWORD
void
BergerRigoutsosNode::runChildren_start()
{
   /*
    * Children were created to store temporary data
    * and determine participation. Now, run them.
    */

#ifdef DEBUG_CHECK_ASSERTIONS
   /*
    * Should only be here if box is rejected based on calculation.
    */
   TBOX_ASSERT(d_box_acceptance == rejected_by_calculation ||
      d_box_acceptance == rejected_by_owner);
#endif

   d_lft_child->d_wait_phase = to_be_launched;
   d_rht_child->d_wait_phase = to_be_launched;

   /*
    * Queue the children so they get executed.
    * Put them at the front so that in synchronous
    * mode, they can complete first before moving
    * to another task (important in synchronous mode).
    * It also does not hurt to put children at the
    * front of the queue because they have
    * immediate computation (compute histogram)
    * to perform.  Put the left child in front
    * of the right to more closely match the
    * progression of the recursive BR (not essential).
    */
   d_common->relaunch_queue.push_front(d_rht_child);
   d_common->relaunch_queue.push_front(d_lft_child);
}

SAMRAI_INLINE_KEYWORD
bool
BergerRigoutsosNode::gatherGroupingCriteria_check()
{
   if (d_group.size() == 1) {
      return true;
   }
   d_comm_group->checkGather();
   /*
    * Do nothing yet with the overlap data d_recv_msg.
    * We extract it in formChildGroups().
    */
   return d_comm_group->isDone();
}

/*
 ********************************************************************
 * Utility methods.
 ********************************************************************
 */

SAMRAI_INLINE_KEYWORD
// tbox::List<BergerRigoutsosNode*>::Iterator
std::list<BergerRigoutsosNode *>::const_iterator
BergerRigoutsosNode::inRelaunchQueue(
   BergerRigoutsosNode* node_ptr) const
{
   std::list<BergerRigoutsosNode *>::const_iterator li =
      std::find(d_common->relaunch_queue.begin(),
         d_common->relaunch_queue.end(),
         node_ptr);
   return li;
}

SAMRAI_INLINE_KEYWORD
int
BergerRigoutsosNode::getHistogramBufferSize(
   const hier::Box& box) const
{
   int size = box.numberCells(0);
   for (int d = 1; d < d_dim.getValue(); ++d) {
      size += box.numberCells(d);
   }
   return size;
}

SAMRAI_INLINE_KEYWORD
int*
BergerRigoutsosNode::putHistogramToBuffer(
   int* buffer)
{
   for (int d = 0; d < d_dim.getValue(); ++d) {
      d_histogram[d].resize(d_box.numberCells(d), BAD_INTEGER);
      memcpy(buffer,
         &d_histogram[d][0],
         d_box.numberCells(d) * sizeof(int));
      buffer += d_box.numberCells(d);
   }
   return buffer;
}

SAMRAI_INLINE_KEYWORD
int*
BergerRigoutsosNode::getHistogramFromBuffer(
   int* buffer)
{
   for (unsigned int d = 0; d < d_dim.getValue(); ++d) {
      TBOX_ASSERT((int)d_histogram[d].size() == d_box.numberCells(d));
      // d_histogram[d].resizeArray( d_box.numberCells(d) );
      memcpy(&d_histogram[d][0],
         buffer,
         d_box.numberCells(d) * sizeof(int));
      buffer += d_box.numberCells(d);
   }
   return buffer;
}

SAMRAI_INLINE_KEYWORD
int*
BergerRigoutsosNode::putBoxToBuffer(
   const hier::Box& box,
   int* buffer) const
{
   const hier::IntVector& l = box.lower();
   const hier::IntVector& u = box.upper();
   for (int d = 0; d < d_dim.getValue(); ++d) {
      *(buffer++) = l(d);
      *(buffer++) = u(d);
   }
   return buffer;
}

SAMRAI_INLINE_KEYWORD
int*
BergerRigoutsosNode::getBoxFromBuffer(
   hier::Box& box,
   int* buffer) const
{
   hier::IntVector& l = box.lower();
   hier::IntVector& u = box.upper();
   for (int d = 0; d < d_dim.getValue(); ++d) {
      l(d) = *(buffer++);
      u(d) = *(buffer++);
   }
   return buffer;
}

SAMRAI_INLINE_KEYWORD
bool
BergerRigoutsosNode::boxAccepted() const
{
   return bool(d_box_acceptance >= 0 && d_box_acceptance % 2);
}

SAMRAI_INLINE_KEYWORD
bool
BergerRigoutsosNode::boxRejected() const
{
   return bool(d_box_acceptance >= 0 && d_box_acceptance % 2 == 0);
}

SAMRAI_INLINE_KEYWORD
bool
BergerRigoutsosNode::boxHasNoTag() const
{
   return bool(d_box_acceptance == -1);
}

/*
 **********************************************************************
 * Determine if the given rank is in the given group.
 **********************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool
BergerRigoutsosNode::inGroup(
   VectorOfInts& group,
   int rank) const
{
   if (rank < 0) {
      rank = d_common->rank;
   }
   for (size_t i = 0; i < group.size(); ++i) {
      if (rank == group[i]) {
         return true;
      }
   }
   return false;
}

/*
 **********************************************************************
 * Heuristically determine the "best" tree degree for
 * the communication group.
 * Use binary tree for less than 2^4 processes.
 * After that, for each 2^3 multiplying of group size,
 * increase tree degree by 1.
 **********************************************************************
 */
SAMRAI_INLINE_KEYWORD
int
BergerRigoutsosNode::computeCommunicationTreeDegree(
   int group_size) const
{
   int tree_deg = 2;
   int shifted_size = group_size >> 3;
   while (shifted_size > 0) {
      shifted_size >>= 3;
      ++tree_deg;
   }
   return tree_deg;
}

SAMRAI_INLINE_KEYWORD
int
BergerRigoutsosNode::findOwnerInGroup(
   int owner,
   const VectorOfInts& group) const
{
   for (unsigned int i = 0; i < group.size(); ++i) {
      if (group[i] == owner) {
         return i;
      }
   }
   return -1;
}

SAMRAI_INLINE_KEYWORD
void
BergerRigoutsosNode::setLogNodeHistory(
   bool flag)
{
   d_common->log_node_history = flag;
}

/*
 **********************************************************************
 *
 * Methods for collecting data on dendogram after it completes.
 *
 **********************************************************************
 */

SAMRAI_INLINE_KEYWORD
int
BergerRigoutsosNode::getNumTags() const
{
   return d_common->num_tags_in_all_nodes;
}

SAMRAI_INLINE_KEYWORD
int
BergerRigoutsosNode::getMaxTagsOwned() const
{
   return d_common->max_tags_owned;
}

SAMRAI_INLINE_KEYWORD
int
BergerRigoutsosNode::getMaxNodes() const
{
   return d_common->max_nodes_allocated;
}

SAMRAI_INLINE_KEYWORD
int
BergerRigoutsosNode::getMaxGeneration() const
{
   return d_common->max_generation;
}

SAMRAI_INLINE_KEYWORD
int
BergerRigoutsosNode::getMaxOwnership() const
{
   return d_common->max_nodes_owned;
}

SAMRAI_INLINE_KEYWORD
double
BergerRigoutsosNode::getAvgNumberOfCont() const
{
   if (d_common->num_nodes_completed > 0) {
      return (double)d_common->num_conts_to_complete
             / d_common->num_nodes_completed;
   }
   return 0;
}

SAMRAI_INLINE_KEYWORD
int
BergerRigoutsosNode::getMaxNumberOfCont() const
{
   return d_common->max_conts_to_complete;
}

SAMRAI_INLINE_KEYWORD
int
BergerRigoutsosNode::getNumBoxesGenerated() const
{
   return d_common->num_boxes_generated;
}

}
}
