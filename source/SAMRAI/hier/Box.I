/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Box representing a portion of the AMR index space 
 *
 ************************************************************************/

#include "SAMRAI/tbox/MathUtilities.h"

namespace SAMRAI {
namespace hier {

SAMRAI_INLINE_KEYWORD
const tbox::Dimension& Box::getDim() const
{
   return d_lo.getDim();
}

SAMRAI_INLINE_KEYWORD
Box::Box():
   d_lo(tbox::Dimension::getInvalidDimension(), tbox::MathUtilities<int>::getMax()),
   d_hi(tbox::Dimension::getInvalidDimension(), tbox::MathUtilities<int>::getMin())
{
}

SAMRAI_INLINE_KEYWORD
Box::Box(
   const tbox::Dimension& dim):
   d_lo(dim, tbox::MathUtilities<int>::getMax()),
   d_hi(dim, tbox::MathUtilities<int>::getMin())
{
   TBOX_DIM_ASSERT_CHECK_DIM(dim);
}

SAMRAI_INLINE_KEYWORD
Box::Box(
   const Index& lower,
   const Index& upper):
   d_lo(lower),
   d_hi(upper)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(lower, upper);
}

SAMRAI_INLINE_KEYWORD
Box::Box(
   const Box& box):d_lo(box.d_lo),
   d_hi(box.d_hi)
{
}

SAMRAI_INLINE_KEYWORD
Box::Box(
   const tbox::DatabaseBox& box):
   d_lo(box.getDim(), tbox::MathUtilities<int>::getMax()),
   d_hi(box.getDim(), tbox::MathUtilities<int>::getMin())
{
   set_Box_from_DatabaseBox(box);
}

SAMRAI_INLINE_KEYWORD
Box::~Box()
{
}


SAMRAI_INLINE_KEYWORD
Box& Box::Box_from_DatabaseBox(
   const tbox::DatabaseBox& box)
{
   set_Box_from_DatabaseBox(box);
   return *this;
}

SAMRAI_INLINE_KEYWORD
Box& Box::operator = (
   const Box& box)
{

   /*
    * Allow assignment of to an uninitialized box but
    * not from an uninitialized box.
    * This is needed for tbox::Array and vector.
    */
   if (getDim().isValid()) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);
   } else {
      TBOX_DIM_ASSERT_CHECK_DIM(box.getDim());
   }

   d_lo = box.d_lo;
   d_hi = box.d_hi;
   return *this;
}

SAMRAI_INLINE_KEYWORD
Box& Box::operator = (
   const tbox::DatabaseBox& box)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   return Box_from_DatabaseBox(box);
}

SAMRAI_INLINE_KEYWORD
Index& Box::lower()
{
   return d_lo;
}

SAMRAI_INLINE_KEYWORD
Index& Box::upper()
{
   return d_hi;
}

SAMRAI_INLINE_KEYWORD
const Index& Box::lower() const
{
   return d_lo;
}

SAMRAI_INLINE_KEYWORD
const Index& Box::upper() const
{
   return d_hi;
}

SAMRAI_INLINE_KEYWORD
int& Box::lower(
   const int i)
{
   return d_lo(i);
}

SAMRAI_INLINE_KEYWORD
int& Box::upper(
   const int i)
{
   return d_hi(i);
}

SAMRAI_INLINE_KEYWORD
const int& Box::lower(
   const int i) const
{
   return d_lo(i);
}

SAMRAI_INLINE_KEYWORD
const int& Box::upper(
   const int i) const
{
   return d_hi(i);
}

SAMRAI_INLINE_KEYWORD
bool Box::empty() const
{
   for (int i = 0; i < getDim().getValue(); i++) {
      if (d_hi(i) < d_lo(i)) return true;
   }

   return false;
}

SAMRAI_INLINE_KEYWORD
bool Box::isEmpty() const
{
   for (int i = 0; i < getDim().getValue(); i++) {
      if (d_hi(i) < d_lo(i)) return true;
   }

   return false;
}

SAMRAI_INLINE_KEYWORD
void Box::setEmpty()
{
   const tbox::Dimension& dim(getDim());

   d_lo = Index(dim, tbox::MathUtilities<int>::getMax());
   d_hi = Index(dim, tbox::MathUtilities<int>::getMin());
}

SAMRAI_INLINE_KEYWORD
IntVector Box::numberCells() const
{
   if (empty()) {
      return IntVector::getZero(getDim());
   } else {
      return d_hi - d_lo + 1;
   }
}

SAMRAI_INLINE_KEYWORD
int Box::numberCells(
   const int i) const
{
   if (empty()) {
      return 0;
   } else {
      return d_hi(i) - d_lo(i) + 1;
   }
}

SAMRAI_INLINE_KEYWORD
int Box::size() const
{
   int mysize = 0;
   if (!empty()) {
      mysize = 1;
      for (int i = 0; i < getDim().getValue(); i++)
         mysize *= (d_hi(i) - d_lo(i) + 1);
   }

   return mysize;
}

SAMRAI_INLINE_KEYWORD
int Box::offset(
   const Index& p) const
{
   int myoffset = 0;

   for (int i = getDim().getValue() - 1; i > 0; i--) {
      myoffset = (d_hi(i - 1) - d_lo(i - 1) + 1) * (p(i) - d_lo(i) + myoffset);
   }
   myoffset += p(0) - d_lo(0);

   return myoffset;
}

SAMRAI_INLINE_KEYWORD
Index Box::index(
   const int offset) const
{
   TBOX_ASSERT(offset >= 0);
   TBOX_ASSERT(offset <= size());

   hier::IntVector n(getDim());
   hier::IntVector index(getDim());

   n = numberCells();

   int remainder = offset;

   for (int d = getDim().getValue() - 1; d > -1; d--) {
      /* Compute the stride for indexing */
      int stride = 1;
      for (int stride_dim = 0; stride_dim < d; stride_dim++) {
         stride *= n[stride_dim];
      }

      /* Compute the local index */
      index[d] = remainder / stride;
      remainder -= index[d] * stride;

      /* Compute the global index */
      index[getDim().getValue()] += lower(d);
   }

   hier::Index idx(index);

   return idx;
}

SAMRAI_INLINE_KEYWORD
bool Box::contains(
   const Index& p) const
{
   for (int i = 0; i < getDim().getValue(); i++) {
      if ((p(i) < d_lo(i)) || (p(i) > d_hi(i))) return false;
   }
   return true;
}

SAMRAI_INLINE_KEYWORD
bool Box::contains(
   const Box& box) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   if (box.empty()) return true;

   if (!contains(box.lower())) return false;

   if (!contains(box.upper())) return false;

   return true;
}

SAMRAI_INLINE_KEYWORD
Box::operator tbox::DatabaseBox ()
{
   return DatabaseBox_from_Box();
}

SAMRAI_INLINE_KEYWORD
Box::operator tbox::DatabaseBox () const
{
   return DatabaseBox_from_Box();
}

SAMRAI_INLINE_KEYWORD
int Box::operator == (
   const Box& box) const
{
   return ((d_lo == box.d_lo) && (d_hi == box.d_hi)) || (empty() && box.empty());
}

SAMRAI_INLINE_KEYWORD
int Box::operator != (
   const Box& box) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   return !(*this == box);
}

SAMRAI_INLINE_KEYWORD
Box Box::operator * (
   const Box& box) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   Box both(*this);
   both.d_lo.max(box.d_lo);
   both.d_hi.min(box.d_hi);
   return both;
}

SAMRAI_INLINE_KEYWORD
bool Box::intersects(
   const Box& box) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   for (int i = 0; i < getDim().getValue(); i++) {
      if (tbox::MathUtilities<int>::Max(d_lo(i), box.d_lo(i)) >
          tbox::MathUtilities<int>::Min(d_hi(i), box.d_hi(i))) {
         return false;
      }
   }
   return true;
}

SAMRAI_INLINE_KEYWORD
Box Box::operator + (
   const Box& box) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   Box bbox(*this);
   bbox += box;
   return bbox;
}

SAMRAI_INLINE_KEYWORD
void Box::grow(
   const IntVector& ghosts)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, ghosts);

   if (!empty()) {
      d_lo -= ghosts;
      d_hi += ghosts;
   }
}

SAMRAI_INLINE_KEYWORD
void Box::grow(
   const int direction,
   const int ghosts)
{
   TBOX_ASSERT((direction >= 0) && (direction < getDim().getValue()));

   if (!empty()) {
      d_lo(direction) -= ghosts;
      d_hi(direction) += ghosts;
   }
}

SAMRAI_INLINE_KEYWORD
void Box::growLower(
   const IntVector& ghosts)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, ghosts);

   if (!empty()) {
      d_lo -= ghosts;
   }
}

SAMRAI_INLINE_KEYWORD
void Box::growLower(
   const int direction,
   const int ghosts)
{
   TBOX_ASSERT((direction >= 0) && (direction < getDim().getValue()));

   if (!empty()) {
      d_lo(direction) -= ghosts;
   }
}

SAMRAI_INLINE_KEYWORD
void Box::growUpper(
   const IntVector& ghosts)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, ghosts);

   if (!empty()) {
      d_hi += ghosts;
   }
}

SAMRAI_INLINE_KEYWORD
void Box::growUpper(
   const int direction,
   const int ghosts)
{
   TBOX_ASSERT((direction >= 0) && (direction < getDim().getValue()));

   if (!empty()) {
      d_hi(direction) += ghosts;
   }
}

SAMRAI_INLINE_KEYWORD
void Box::lengthen(
   const int direction,
   const int ghosts)
{
   TBOX_ASSERT((direction >= 0) && (direction < getDim().getValue()));

   if (!empty()) {
      if (ghosts > 0) {
         d_hi(direction) += ghosts;
      } else {
         d_lo(direction) += ghosts;
      }
   }
}

SAMRAI_INLINE_KEYWORD
void Box::shorten(
   const int direction,
   const int ghosts)
{
   TBOX_ASSERT((direction >= 0) && (direction < getDim().getValue()));

   if (!empty()) {
      if (ghosts > 0) {
         d_hi(direction) -= ghosts;
      } else {
         d_lo(direction) -= ghosts;
      }
   }
}

SAMRAI_INLINE_KEYWORD
void Box::shift(
   const IntVector& offset)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, offset);

   d_lo += offset;
   d_hi += offset;
}

SAMRAI_INLINE_KEYWORD
void Box::shift(
   const int direction,
   const int offset)
{
   TBOX_ASSERT((direction >= 0) && (direction < getDim().getValue()));

   d_lo(direction) += offset;
   d_hi(direction) += offset;
}

SAMRAI_INLINE_KEYWORD
int Box::coarsen(
   const int index,
   const int ratio)
{
   return index < 0 ? (index + 1) / ratio - 1 : index / ratio;
}

SAMRAI_INLINE_KEYWORD
void Box::coarsen(
   const IntVector& ratio)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, ratio);

   for (int i = 0; i < getDim().getValue(); i++) {
      d_lo(i) = coarsen(d_lo(i), ratio(i));
      d_hi(i) = coarsen(d_hi(i), ratio(i));
   }
}

SAMRAI_INLINE_KEYWORD
void Box::refine(
   const IntVector& ratio)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, ratio);

   bool negative_ratio = false;
   for (int d = 0; d < getDim().getValue(); d++) {
      if (ratio(d) < 0) {
         negative_ratio = true;
         break;
      }
   }

   if (!negative_ratio) {
      d_lo *= ratio;
      d_hi = d_hi * ratio + (ratio - 1);
   } else {
      for (int i = 0; i < getDim().getValue(); i++) {
         if (ratio(i) > 0) {
            d_lo(i) *= ratio(i);
            d_hi(i) = d_hi(i) * ratio(i) + (ratio(i) - 1);
         } else {
            d_lo(i) = coarsen(d_lo(i), -ratio(i));
            d_hi(i) = coarsen(d_hi(i), -ratio(i));
         }
      }
   }
}

SAMRAI_INLINE_KEYWORD
Box Box::grow(
   const Box& box,
   const IntVector& ghosts)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(box, ghosts);

   Box tmp = box;
   tmp.grow(ghosts);
   return tmp;
}

SAMRAI_INLINE_KEYWORD
Box Box::shift(
   const Box& box,
   const IntVector& ghosts)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(box, ghosts);

   return Box(box.lower() + ghosts, box.upper() + ghosts);
}

SAMRAI_INLINE_KEYWORD
Box Box::refine(
   const Box& box,
   const IntVector& ratio)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(box, ratio);

   Box tmp = box;
   tmp.refine(ratio);
   return tmp;
}

SAMRAI_INLINE_KEYWORD
Box Box::coarsen(
   const Box& box,
   const IntVector& ratio)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(box, ratio);

   Box tmp = box;
   tmp.coarsen(ratio);
   return tmp;
}

SAMRAI_INLINE_KEYWORD
BoxIterator::BoxIterator():
   d_index(tbox::Dimension::getInvalidDimension())
{
}

SAMRAI_INLINE_KEYWORD
BoxIterator::BoxIterator(
   const Box& box):
   d_index(box.lower()),
   d_box(box)
{
   TBOX_DIM_ASSERT_CHECK_DIM(box.getDim());
}

SAMRAI_INLINE_KEYWORD
BoxIterator::BoxIterator(
   const BoxIterator& iter):
   d_index(iter.d_index),
   d_box(iter.d_box)
{
}

SAMRAI_INLINE_KEYWORD
void BoxIterator::initialize(
   const Box& box)
{
   TBOX_DIM_ASSERT_CHECK_DIM(box.getDim());

   d_index = box.lower();
   d_box = box;
}

SAMRAI_INLINE_KEYWORD
BoxIterator&
BoxIterator::operator = (
   const BoxIterator& iter)
{
   d_index = iter.d_index;
   d_box = iter.d_box;
   return *this;
}

SAMRAI_INLINE_KEYWORD
BoxIterator::~BoxIterator()
{
}

SAMRAI_INLINE_KEYWORD
const Index& BoxIterator::operator * () const
{
   return d_index;
}

SAMRAI_INLINE_KEYWORD
const Index& BoxIterator::operator () () const
{
   return d_index;
}

SAMRAI_INLINE_KEYWORD
BoxIterator::operator bool () const
{
   for (int i = 0; i < d_index.getDim().getValue(); i++) {
      if (d_index(i) > d_box.upper(i)) return false;
   }
   return true;
}

#ifndef LACKS_BOOL_VOID_RESOLUTION
SAMRAI_INLINE_KEYWORD
BoxIterator::operator const void* () const
{
   return BoxIterator::operator bool () ? this : NULL;
}
#endif

SAMRAI_INLINE_KEYWORD
bool BoxIterator::operator ! () const
{
   return !BoxIterator::operator bool ();
}

SAMRAI_INLINE_KEYWORD
void BoxIterator::operator ++ (
   int)
{
   d_index(0)++;
   for (int i = 0; i < (d_index.getDim().getValue() - 1); i++) {
      if (d_index(i) > d_box.upper(i)) {
         d_index(i) = d_box.lower(i);
         d_index(i + 1)++;
      } else
         break;
   }
}

SAMRAI_INLINE_KEYWORD
bool BoxIterator::operator == (
   const BoxIterator& iter) const
{
   return d_index == iter.d_index;
}

SAMRAI_INLINE_KEYWORD
bool BoxIterator::operator != (
   const BoxIterator& iter) const
{
   return d_index != iter.d_index;
}

SAMRAI_INLINE_KEYWORD
const Box& Box::getEmptyBox(
   const tbox::Dimension& dim) {
   TBOX_DIM_ASSERT_CHECK_DIM(dim);
   return *(s_emptys[dim.getValue() - 1]);
}

SAMRAI_INLINE_KEYWORD
const Box& Box::getUniverse(
   const tbox::Dimension& dim) {
   TBOX_DIM_ASSERT_CHECK_DIM(dim);
   return *(s_universes[dim.getValue() - 1]);
}

}
}
