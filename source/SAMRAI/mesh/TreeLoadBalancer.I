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
TreeLoadBalancer::BoxInTransit::BoxInTransit(
   const tbox::Dimension &dim) :
   box(dim),
   orig_box(dim)
{
}



/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
TreeLoadBalancer::BoxInTransit::BoxInTransit(
   const hier::Box& origin):
   box(origin),
   orig_box(origin),
   load(origin.size()),
   proc_hist()
{
}



/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
TreeLoadBalancer::BoxInTransit::BoxInTransit(
   const int *&ptr,
   const tbox::Dimension &dim ):
   box(dim),
   orig_box(dim),
   load(0),
   proc_hist(0)
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
TreeLoadBalancer::BoxInTransit::BoxInTransit(
   const BoxInTransit& other,
   const hier::Box& box,
   int rank,
   LocalId local_id):
   box(box, local_id, rank, other.orig_box.getBlockId()),
   orig_box(other.orig_box),
   load(box.size()),
   proc_hist(other.proc_hist)
{
   if (rank != other.getOwnerRank()) {
      proc_hist.push_back(other.getOwnerRank());
   }
}



/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
int TreeLoadBalancer::BoxInTransit::getOwnerRank() const {
   return box.getOwnerRank();
}



/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
hier::LocalId TreeLoadBalancer::BoxInTransit::getLocalId() const {
   return box.getLocalId();
}



/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
hier::Box& TreeLoadBalancer::BoxInTransit::getBox() {
   return box;
}



/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
const hier::Box& TreeLoadBalancer::BoxInTransit::getBox() const {
   return box;
}



/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
int TreeLoadBalancer::BoxInTransit::commBufferSize() const
{
   const tbox::Dimension& dim(box.getDim());
   return 2 * hier::Box::commBufferSize(dim) + 2 + static_cast<int>(proc_hist.size());
}



/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
const TreeLoadBalancer::BoxInTransit& TreeLoadBalancer::BoxInTransit::operator = (
   const BoxInTransit& other)
{
   box = other.box;
   orig_box = other.orig_box;
   load = other.load;
   proc_hist = other.proc_hist;
   return *this;
}



/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
int * TreeLoadBalancer::BoxInTransit::putToIntBuffer(
   int* buffer,
   bool skip_last_owner) const
{
   const tbox::Dimension& dim(box.getDim());
   box.putToIntBuffer(buffer);
   buffer += hier::Box::commBufferSize(dim);
   orig_box.putToIntBuffer(buffer);
   buffer += hier::Box::commBufferSize(dim);
   *(buffer++) = load;
   if (skip_last_owner) { TBOX_ASSERT( !proc_hist.empty() ); }
   *(buffer++) = static_cast<int>(proc_hist.size()-skip_last_owner);
   for (unsigned int i = 0; i < proc_hist.size()-skip_last_owner; ++i) {
      *(buffer++) = proc_hist[i];
   }
   return buffer;
}



/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
const int * TreeLoadBalancer::BoxInTransit::getFromIntBuffer(
   const int* buffer)
{
   const tbox::Dimension& dim(box.getDim());
   box.getFromIntBuffer(buffer);
   buffer += hier::Box::commBufferSize(dim);
   orig_box.getFromIntBuffer(buffer);
   buffer += hier::Box::commBufferSize(dim);
   load = *(buffer++);
   proc_hist.clear();
   proc_hist.insert(proc_hist.end(), *(buffer++), 0);
   for (unsigned int i = 0; i < proc_hist.size(); ++i) {
      proc_hist[i] = *(buffer++);
   }
   return buffer;
}



/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void TreeLoadBalancer::BoxInTransit::packForPreviousOwner(
   std::map<int,std::vector<int> > &outgoing_messages ) const
{
   const int prev_owner = proc_hist.back();
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
bool TreeLoadBalancer::BoxInTransitMoreLoad::operator () (
   const TreeLoadBalancer::BoxInTransit& a,
   const TreeLoadBalancer::BoxInTransit& b) const
{
   if (a.getBox().size() != b.getBox().size()) {
      return a.load > b.load;
   }
   return a.box.getId() < b.box.getId();
}



/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
double TreeLoadBalancer::computeLoad(
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
double TreeLoadBalancer::computeLoad(
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
