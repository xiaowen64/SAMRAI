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
#include "SAMRAI/hier/BoxId.h"
#include "SAMRAI/hier/PeriodicShiftCatalog.h"

namespace SAMRAI {
namespace hier {

SAMRAI_INLINE_KEYWORD
const tbox::Dimension&
Box::getDim() const
{
   return d_lo.getDim();
}

SAMRAI_INLINE_KEYWORD
void
Box::initialize(
   const Box& box,
   const LocalId& local_id,
   const int owner,
   const PeriodicId& periodic_id)
{
   d_lo = box.d_lo;
   d_hi = box.d_hi;
   d_block_id = box.d_block_id;
   if (!d_id_locked) {
      d_id.initialize(local_id, owner, periodic_id);
   } else {
      TBOX_ERROR("Attempted to change BoxId that is locked in an ordered BoxContainer.");
   }
}

SAMRAI_INLINE_KEYWORD
void
Box::setId(const BoxId& box_id) {
   d_id = box_id;
}

SAMRAI_INLINE_KEYWORD
const BoxId&
Box::getId() const {
   return d_id;
}

SAMRAI_INLINE_KEYWORD
int
Box::getOwnerRank() const
{
   return d_id.getOwnerRank();
}

SAMRAI_INLINE_KEYWORD
const LocalId&
Box::getLocalId() const
{
   return d_id.getLocalId();
}
SAMRAI_INLINE_KEYWORD
const GlobalId&
Box::getGlobalId() const
{
   return d_id.getGlobalId();
}

SAMRAI_INLINE_KEYWORD
void
Box::setBlockId(const BlockId& block_id)
{
   d_block_id = block_id;
}

SAMRAI_INLINE_KEYWORD
const BlockId&
Box::getBlockId() const
{
   return d_block_id;
}

SAMRAI_INLINE_KEYWORD
const PeriodicId&
Box::getPeriodicId() const
{
   return d_id.getPeriodicId();
}

SAMRAI_INLINE_KEYWORD
bool
Box::isPeriodicImage() const
{
   return d_id.isPeriodicImage();
}

SAMRAI_INLINE_KEYWORD
bool
Box::isIdEqual(const Box& other) const
{
   return d_id == other.d_id;
}

SAMRAI_INLINE_KEYWORD
int
Box::commBufferSize(
   const tbox::Dimension& dim)
{
   return BoxId::commBufferSize() + (2 * dim.getValue()) + 1;
}

SAMRAI_INLINE_KEYWORD
Box&
Box::Box_from_DatabaseBox(
   const tbox::DatabaseBox& box)
{
   set_Box_from_DatabaseBox(box);
   return *this;
}

SAMRAI_INLINE_KEYWORD
Box&
Box::operator = (
   const tbox::DatabaseBox& box)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);
#ifdef BOX_TELEMETRY
   // Increment the cumulative assigned count only.
   ++s_cumulative_assigned_ct;
#endif

   return Box_from_DatabaseBox(box);
}

SAMRAI_INLINE_KEYWORD
Index&
Box::lower()
{
   return d_lo;
}

SAMRAI_INLINE_KEYWORD
Index&
Box::upper()
{
   return d_hi;
}

SAMRAI_INLINE_KEYWORD
const Index&
Box::lower() const
{
   return d_lo;
}

SAMRAI_INLINE_KEYWORD
const Index&
Box::upper() const
{
   return d_hi;
}

SAMRAI_INLINE_KEYWORD
int&
Box::lower(
   const int i)
{
   return d_lo(i);
}

SAMRAI_INLINE_KEYWORD
int&
Box::upper(
   const int i)
{
   return d_hi(i);
}

SAMRAI_INLINE_KEYWORD
const int&
Box::lower(
   const int i) const
{
   return d_lo(i);
}

SAMRAI_INLINE_KEYWORD
const int&
Box::upper(
   const int i) const
{
   return d_hi(i);
}

SAMRAI_INLINE_KEYWORD
bool
Box::empty() const
{
   for (int i = 0; i < getDim().getValue(); i++) {
      if (d_hi(i) < d_lo(i)) {
         return true;
      }
   }

   return false;
}

SAMRAI_INLINE_KEYWORD
bool
Box::isEmpty() const
{
   for (int i = 0; i < getDim().getValue(); i++) {
      if (d_hi(i) < d_lo(i)) {
         return true;
      }
   }

   return false;
}

SAMRAI_INLINE_KEYWORD
void
Box::setEmpty()
{
   const tbox::Dimension& dim(getDim());

   d_lo = Index(dim, tbox::MathUtilities<int>::getMax());
   d_hi = Index(dim, tbox::MathUtilities<int>::getMin());
}

SAMRAI_INLINE_KEYWORD
IntVector
Box::numberCells() const
{
   if (empty()) {
      return IntVector::getZero(getDim());
   } else {
      return d_hi - d_lo + 1;
   }
}

SAMRAI_INLINE_KEYWORD
int
Box::numberCells(
   const int i) const
{
   if (empty()) {
      return 0;
   } else {
      return d_hi(i) - d_lo(i) + 1;
   }
}

SAMRAI_INLINE_KEYWORD
int
Box::size() const
{
   int mysize = 0;
   if (!empty()) {
      mysize = 1;
      for (int i = 0; i < getDim().getValue(); i++) {
         mysize *= (d_hi(i) - d_lo(i) + 1);
      }
   }

   return mysize;
}

SAMRAI_INLINE_KEYWORD
int
Box::offset(
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
bool
Box::contains(
   const Index& p) const
{
   for (int i = 0; i < getDim().getValue(); i++) {
      if ((p(i) < d_lo(i)) || (p(i) > d_hi(i))) {
         return false;
      }
   }
   return true;
}

SAMRAI_INLINE_KEYWORD
Box::operator tbox::DatabaseBox()
{
   return DatabaseBox_from_Box();
}

SAMRAI_INLINE_KEYWORD
Box::operator tbox::DatabaseBox() const
{
   return DatabaseBox_from_Box();
}

SAMRAI_INLINE_KEYWORD
bool
Box::isSpatiallyEqual(
   const Box& box) const
{
   return ((d_lo == box.d_lo) && (d_hi == box.d_hi) && (d_block_id == box.d_block_id)) || (empty() && box.empty());
}

SAMRAI_INLINE_KEYWORD
void
Box::grow(
   const IntVector& ghosts)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, ghosts);

   if (!empty()) {
      d_lo -= ghosts;
      d_hi += ghosts;
   }
}

SAMRAI_INLINE_KEYWORD
void
Box::grow(
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
void
Box::growLower(
   const IntVector& ghosts)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, ghosts);

   if (!empty()) {
      d_lo -= ghosts;
   }
}

SAMRAI_INLINE_KEYWORD
void
Box::growLower(
   const int direction,
   const int ghosts)
{
   TBOX_ASSERT((direction >= 0) && (direction < getDim().getValue()));

   if (!empty()) {
      d_lo(direction) -= ghosts;
   }
}

SAMRAI_INLINE_KEYWORD
void
Box::growUpper(
   const IntVector& ghosts)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, ghosts);

   if (!empty()) {
      d_hi += ghosts;
   }
}

SAMRAI_INLINE_KEYWORD
void
Box::growUpper(
   const int direction,
   const int ghosts)
{
   TBOX_ASSERT((direction >= 0) && (direction < getDim().getValue()));

   if (!empty()) {
      d_hi(direction) += ghosts;
   }
}

SAMRAI_INLINE_KEYWORD
void
Box::shift(
   const IntVector& offset)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, offset);

   d_lo += offset;
   d_hi += offset;
}

SAMRAI_INLINE_KEYWORD
void
Box::shift(
   const int direction,
   const int offset)
{
   TBOX_ASSERT((direction >= 0) && (direction < getDim().getValue()));

   d_lo(direction) += offset;
   d_hi(direction) += offset;
}

SAMRAI_INLINE_KEYWORD
int
Box::coarsen(
   const int index,
   const int ratio)
{
   return index < 0 ? (index + 1) / ratio - 1 : index / ratio;
}

SAMRAI_INLINE_KEYWORD
void
Box::coarsen(
   const IntVector& ratio)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, ratio);

   for (int i = 0; i < getDim().getValue(); i++) {
      d_lo(i) = coarsen(d_lo(i), ratio(i));
      d_hi(i) = coarsen(d_hi(i), ratio(i));
   }
}

SAMRAI_INLINE_KEYWORD
Box
Box::grow(
   const Box& box,
   const IntVector& ghosts)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(box, ghosts);

   Box tmp = box;
   tmp.grow(ghosts);
   return tmp;
}

SAMRAI_INLINE_KEYWORD
Box
Box::shift(
   const Box& box,
   const IntVector& ghosts)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(box, ghosts);

   return Box(box.lower() + ghosts, box.upper() + ghosts, box.d_block_id);
}

SAMRAI_INLINE_KEYWORD
Box
Box::refine(
   const Box& box,
   const IntVector& ratio)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(box, ratio);

   Box tmp = box;
   tmp.refine(ratio);
   return tmp;
}

SAMRAI_INLINE_KEYWORD
Box
Box::coarsen(
   const Box& box,
   const IntVector& ratio)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(box, ratio);

   Box tmp = box;
   tmp.coarsen(ratio);
   return tmp;
}

SAMRAI_INLINE_KEYWORD
void
BoxIterator::initialize(
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
const Index&
BoxIterator::operator * () const
{
   return d_index;
}

SAMRAI_INLINE_KEYWORD
const Index&
BoxIterator::operator () () const
{
   return d_index;
}

SAMRAI_INLINE_KEYWORD
BoxIterator::operator bool () const
{
   for (int i = 0; i < d_index.getDim().getValue(); i++) {
      if (d_index(i) > d_box.upper(i)) {
         return false;
      }
   }
   return true;
}

#ifndef LACKS_BOOL_VOID_RESOLUTION
SAMRAI_INLINE_KEYWORD
BoxIterator::operator const void * () const
{
   return BoxIterator::operator bool () ? this : NULL;
}
#endif

SAMRAI_INLINE_KEYWORD
bool
BoxIterator::operator ! () const
{
   return !BoxIterator::operator bool ();
}

SAMRAI_INLINE_KEYWORD
void
BoxIterator::operator ++ (
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
bool
BoxIterator::operator == (
   const BoxIterator& iter) const
{
   return d_index == iter.d_index;
}

SAMRAI_INLINE_KEYWORD
bool
BoxIterator::operator != (
   const BoxIterator& iter) const
{
   return d_index != iter.d_index;
}

SAMRAI_INLINE_KEYWORD
const Box&
Box::getEmptyBox(
   const tbox::Dimension& dim) {
   TBOX_DIM_ASSERT_CHECK_DIM(dim);
   return *(s_emptys[dim.getValue() - 1]);
}

SAMRAI_INLINE_KEYWORD
const Box&
Box::getUniverse(
   const tbox::Dimension& dim) {
   TBOX_DIM_ASSERT_CHECK_DIM(dim);
   return *(s_universes[dim.getValue() - 1]);
}

SAMRAI_INLINE_KEYWORD
void
Box::lockId() {
   d_id_locked = true;
}

}
}
