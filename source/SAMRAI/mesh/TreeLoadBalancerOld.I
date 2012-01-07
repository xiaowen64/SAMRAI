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
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
TreeLoadBalancerOld::BoxInTransit::BoxInTransit(
   const tbox::Dimension &dim) :
   d_box(dim),
   d_orig_box(dim)
{
}



/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
TreeLoadBalancerOld::BoxInTransit::BoxInTransit(
   const hier::Box& origin):
   d_box(origin),
   d_orig_box(origin),
   d_load(origin.size()),
   d_proc_hist()
{
}



/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
TreeLoadBalancerOld::BoxInTransit::BoxInTransit(
   const int *&ptr,
   const tbox::Dimension &dim ):
   d_box(dim),
   d_orig_box(dim),
   d_load(0),
   d_proc_hist(0)
{
   getFromIntBuffer(ptr);
   ptr += commBufferSize();
}



/*
 *************************************************************************
 * Construct a new BoxInTransit with the history of an existing box.
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
TreeLoadBalancerOld::BoxInTransit::BoxInTransit(
   const BoxInTransit& other,
   const hier::Box& box,
   int rank,
   hier::LocalId local_id):
   d_box(box, local_id, rank),
   d_orig_box(other.d_orig_box),
   d_load(d_box.size()),
   d_proc_hist(other.d_proc_hist)
{
   if (rank != other.getOwnerRank()) {
      d_proc_hist.push_back(other.getOwnerRank());
   }
}



/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
int TreeLoadBalancerOld::BoxInTransit::getOwnerRank() const {
   return d_box.getOwnerRank();
}



/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
hier::LocalId TreeLoadBalancerOld::BoxInTransit::getLocalId() const {
   return d_box.getLocalId();
}



/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
hier::Box& TreeLoadBalancerOld::BoxInTransit::getBox() {
   return d_box;
}



/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
const hier::Box& TreeLoadBalancerOld::BoxInTransit::getBox() const {
   return d_box;
}



/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
int TreeLoadBalancerOld::BoxInTransit::commBufferSize() const
{
   const tbox::Dimension& dim(d_box.getDim());
   return 2 * hier::Box::commBufferSize(dim) + 2 + static_cast<int>(d_proc_hist.size());
}



/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
const TreeLoadBalancerOld::BoxInTransit& TreeLoadBalancerOld::BoxInTransit::operator = (
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
int * TreeLoadBalancerOld::BoxInTransit::putToIntBuffer(
   int* buffer,
   bool skip_last_owner) const
{
   const tbox::Dimension& dim(d_box.getDim());
   d_box.putToIntBuffer(buffer);
   buffer += hier::Box::commBufferSize(dim);
   d_orig_box.putToIntBuffer(buffer);
   buffer += hier::Box::commBufferSize(dim);
   *(buffer++) = d_load;
   if (skip_last_owner) { TBOX_ASSERT( !d_proc_hist.empty() ); }
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
const int * TreeLoadBalancerOld::BoxInTransit::getFromIntBuffer(
   const int* buffer)
{
   const tbox::Dimension& dim(d_box.getDim());
   d_box.getFromIntBuffer(buffer);
   buffer += hier::Box::commBufferSize(dim);
   d_orig_box.getFromIntBuffer(buffer);
   buffer += hier::Box::commBufferSize(dim);
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
void TreeLoadBalancerOld::BoxInTransit::packForPreviousOwner(
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
bool TreeLoadBalancerOld::BoxInTransitMoreLoad::operator () (
   const TreeLoadBalancerOld::BoxInTransit& a,
   const TreeLoadBalancerOld::BoxInTransit& b) const
{
   if (a.getBox().size() != b.getBox().size()) {
      return a.d_load > b.d_load;
   }
   return a.d_box.getId() < b.d_box.getId();
}



/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
double TreeLoadBalancerOld::computeLoad(
   const hier::Box& box,
   const hier::Box& restriction) const
{
   /*
    * Currently only for uniform loads, where the load is equal
    * to the number of cells.  For non-uniform loads, this method
    * needs the patch data index for the load.  It would summ up
    * the individual cell loads in the overlap region.
    */
   hier::Box overlap = box * restriction;
   return double(overlap.size());
}



/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
double TreeLoadBalancerOld::computeLoad(
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

}
}
