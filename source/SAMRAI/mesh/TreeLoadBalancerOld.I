/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Scalable load balancer using tree algorithm.
 *
 ************************************************************************/

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <cassert>
#endif
#endif

namespace SAMRAI {
namespace mesh {



/*
**************************************************************************
**************************************************************************
*/
SAMRAI_INLINE_KEYWORD
void
TreeLoadBalancerOld::setUniformWorkload(
   int level_number)
{
   d_workload_data_id[level_number] = -1;
}



/*
 *************************************************************************
 * Compute the total combined penalty associated with box cutting.
 * Any kind of binary function can be used.  I am just experimenting
 * with various types and weights.
 *
 * All input penalties should all be non-dimensional.
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
double
TreeLoadBalancerOld::combinedBreakingPenalty(
   double balance_penalty,
   double surface_penalty,
   double slender_penalty) const
{
   double combined_penalty =
      d_balance_penalty_wt * balance_penalty * balance_penalty
      + d_surface_penalty_wt * surface_penalty * surface_penalty
      + d_slender_penalty_wt * slender_penalty * slender_penalty;
   return combined_penalty;
}



/*
 *************************************************************************
 * Return non-dimensional volume-weighted balance penalty for
 * two containers of boxes and how imbalanced they are.
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
double
TreeLoadBalancerOld::computeBalancePenalty(
   const std::vector<hier::Box>& a,
   const std::vector<hier::Box>& b,
   double imbalance) const
{
   NULL_USE(a);
   NULL_USE(b);
   return tbox::MathUtilities<double>::Abs(imbalance);
}



/*
 *************************************************************************
 * Return non-dimensional volume-weighted balance penalty for
 * two containers of boxes and how imbalanced they are.
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
double
TreeLoadBalancerOld::computeBalancePenalty(
   const TransitSet& a,
   const TransitSet& b,
   double imbalance) const
{
   NULL_USE(a);
   NULL_USE(b);
   return tbox::MathUtilities<double>::Abs(imbalance);
}



/*
 *************************************************************************
 * Return non-dimensional volume-weighted balance penalty for
 * a box and how imbalanced it is.
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
double
TreeLoadBalancerOld::computeBalancePenalty(
   const hier::Box& a,
   double imbalance) const
{
   NULL_USE(a);
   return tbox::MathUtilities<double>::Abs(imbalance);
}



/*
 *************************************************************************
 *
 * Private utility functions to determine parameter values for level.
 *
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
int
TreeLoadBalancerOld::getWorkloadDataId(
   int level_number) const
{

   TBOX_ASSERT(level_number >= 0);

   int wrk_id = (level_number < d_workload_data_id.getSize() ?
                 d_workload_data_id[level_number] :
                 d_master_workload_data_id);

   return wrk_id;
}



/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
TreeLoadBalancerOld::printStatistics(
   std::ostream& s) const
{
   BalanceUtilities::gatherAndReportLoadBalance(d_load_stat,
      tbox::SAMRAI_MPI::getSAMRAIWorld(),
      s);
}



/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
double
TreeLoadBalancerOld::computeLoad(
   const hier::Box& box,
   const hier::Box& restriction) const
{
   /*
    * Currently only for uniform loads, where the load is equal
    * to the number of cells.  For non-uniform loads, this method
    * needs the patch data index for the load.  It would summ up
    * the individual cell loads in the overlap region.
    */
   return double((box * restriction).size());
}



/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
double
TreeLoadBalancerOld::computeLoad(
   const hier::Box& mapped_box) const
{
   /*
    * Currently only for uniform loads, where the load is equal
    * to the number of cells.  For non-uniform loads, this method
    * needs the patch data index for the load.  It would summ up
    * the individual cell loads in the cell.
    */
   return double(mapped_box.size());
}



/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
int
TreeLoadBalancerOld::BoxInTransit::getOwnerRank() const
{
   return d_box.getOwnerRank();
}



/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
hier::LocalId
TreeLoadBalancerOld::BoxInTransit::getLocalId() const
{
   return d_box.getLocalId();
}



/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
hier::Box&
TreeLoadBalancerOld::BoxInTransit::getBox()
{
   return d_box;
}



/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
const hier::Box&
TreeLoadBalancerOld::BoxInTransit::getBox() const
{
   return d_box;
}



/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
int
TreeLoadBalancerOld::BoxInTransit::commBufferSize() const
{
   const tbox::Dimension& dim(d_box.getDim());
   return 2 * hier::Box::commBufferSize(dim) + 2 +
      static_cast<int>(d_proc_hist.size());
}



/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
const TreeLoadBalancerOld::BoxInTransit&
TreeLoadBalancerOld::BoxInTransit::operator = (
   const BoxInTransit& other)
{
   d_box = other.d_box;
   d_orig_box = other.d_orig_box;
   d_load = other.d_load;
   d_proc_hist = other.d_proc_hist;
   return *this;
}



/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
int*
TreeLoadBalancerOld::BoxInTransit::putToIntBuffer(
   int* buffer,
   bool skip_last_owner) const
{
   const int box_comm_buf_size = hier::Box::commBufferSize(d_box.getDim());
   d_box.putToIntBuffer(buffer);
   buffer += box_comm_buf_size;
   d_orig_box.putToIntBuffer(buffer);
   buffer += box_comm_buf_size;
   *(buffer++) = d_load;
   if (skip_last_owner) {
      TBOX_ASSERT( !d_proc_hist.empty() );
   }
   *(buffer++) = static_cast<int>(d_proc_hist.size()-skip_last_owner);
   for (unsigned int i = 0; i < d_proc_hist.size()-skip_last_owner; ++i) {
      *(buffer++) = d_proc_hist[i];
   }
   return buffer;
}



/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
const int*
TreeLoadBalancerOld::BoxInTransit::getFromIntBuffer(
   const int* buffer)
{
   const int box_comm_buf_size = hier::Box::commBufferSize(d_box.getDim());
   d_box.getFromIntBuffer(buffer);
   buffer += box_comm_buf_size;
   d_orig_box.getFromIntBuffer(buffer);
   buffer += box_comm_buf_size;
   d_load = *(buffer++);
   d_proc_hist.clear();
   d_proc_hist.insert(d_proc_hist.end(), *(buffer++), 0);
   for (unsigned int i = 0; i < d_proc_hist.size(); ++i) {
      d_proc_hist[i] = *(buffer++);
   }
   return buffer;
}



/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
TreeLoadBalancerOld::BoxInTransit::packForPreviousOwner(
   std::map<int,std::vector<int> > &outgoing_messages ) const
{
   const int prev_owner = d_proc_hist.back();
   std::vector<int> &msg(outgoing_messages[prev_owner]);
   const int cbs1 = commBufferSize()-1;
   msg.insert(msg.end(), cbs1, 0);
   putToIntBuffer(&msg[msg.size() - cbs1], true);
}



/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool
TreeLoadBalancerOld::BoxInTransitMoreLoad::operator () (
   const TreeLoadBalancerOld::BoxInTransit& a,
   const TreeLoadBalancerOld::BoxInTransit& b) const
{
   if (a.getBox().size() != b.getBox().size()) {
      return a.d_load > b.d_load;
   }
   return a.d_box.getId() < b.d_box.getId();
}

}
}
