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
   d_box(dim),
   d_orig_box(dim)
{
}



/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
TreeLoadBalancer::BoxInTransit::BoxInTransit(
   const hier::Box& origin):
   d_box(origin),
   d_orig_box(origin),
   d_load(origin.size())
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
   d_box(dim),
   d_orig_box(dim),
   d_load(0)
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
   hier::LocalId local_id):
   d_box(box, local_id, rank),
   d_orig_box(other.d_orig_box),
   d_load(d_box.size())
{
}



/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
int TreeLoadBalancer::BoxInTransit::getOwnerRank() const {
   return d_box.getOwnerRank();
}



/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
hier::LocalId TreeLoadBalancer::BoxInTransit::getLocalId() const {
   return d_box.getLocalId();
}



/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
hier::Box& TreeLoadBalancer::BoxInTransit::getBox() {
   return d_box;
}



/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
const hier::Box& TreeLoadBalancer::BoxInTransit::getBox() const {
   return d_box;
}



/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
int TreeLoadBalancer::BoxInTransit::commBufferSize() const
{
   const tbox::Dimension& dim(d_box.getDim());
   return 2 * hier::Box::commBufferSize(dim) + 1;
}



/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
const TreeLoadBalancer::BoxInTransit& TreeLoadBalancer::BoxInTransit::operator = (
   const BoxInTransit& other)
{
   d_box = other.d_box;
   d_orig_box = other.d_orig_box;
   d_load = other.d_load;
   return *this;
}



/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
int * TreeLoadBalancer::BoxInTransit::putToIntBuffer(
   int* buffer) const
{
   const tbox::Dimension& dim(d_box.getDim());
   d_box.putToIntBuffer(buffer);
   buffer += hier::Box::commBufferSize(dim);
   d_orig_box.putToIntBuffer(buffer);
   buffer += hier::Box::commBufferSize(dim);
   *(buffer++) = d_load;
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
   const tbox::Dimension& dim(d_box.getDim());
   d_box.getFromIntBuffer(buffer);
   buffer += hier::Box::commBufferSize(dim);
   d_orig_box.getFromIntBuffer(buffer);
   buffer += hier::Box::commBufferSize(dim);
   d_load = *(buffer++);
   return buffer;
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
      return a.d_load > b.d_load;
   }
   return a.d_box.getId() < b.d_box.getId();
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
