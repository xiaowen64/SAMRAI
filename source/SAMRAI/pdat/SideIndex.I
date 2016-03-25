/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace pdat {

SAMRAI_INLINE_KEYWORD
SideIndex::SideIndex(
   const tbox::Dimension& dim):
   hier::Index(dim)
{
}

SAMRAI_INLINE_KEYWORD
SideIndex::SideIndex(
   const hier::Index& rhs,
   const int axis,
   const int side):
   hier::Index(rhs),
   d_axis(axis)
{
   (*this)(d_axis) += side;
}

SAMRAI_INLINE_KEYWORD
SideIndex::SideIndex(
   const SideIndex& rhs):
   hier::Index(rhs),
   d_axis(rhs.d_axis)
{
}

SAMRAI_INLINE_KEYWORD
SideIndex& SideIndex::operator = (
   const SideIndex& rhs)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   hier::Index::operator = (rhs);
   d_axis = rhs.d_axis;
   return *this;
}

SAMRAI_INLINE_KEYWORD
SideIndex::~SideIndex()
{
}

SAMRAI_INLINE_KEYWORD
int SideIndex::getAxis() const
{
   return d_axis;
}

SAMRAI_INLINE_KEYWORD
void SideIndex::setAxis(
   const int axis)
{
   d_axis = axis;
}

SAMRAI_INLINE_KEYWORD
hier::Index SideIndex::toCell(
   const int side) const
{
   const tbox::Dimension& dim(getDim());

   hier::Index index(dim);

   for (int i = 0; i < dim.getValue(); i++) {
      index(i) = (*this)(i);
   }

   index(d_axis) += (side - 1);

   return index;
}

SAMRAI_INLINE_KEYWORD
SideIndex& SideIndex::operator += (
   const hier::IntVector& rhs)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   hier::Index::operator += (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
SideIndex SideIndex::operator + (
   const hier::IntVector& rhs) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   SideIndex tmp = *this;
   tmp += rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
SideIndex& SideIndex::operator += (
   const int rhs)
{
   hier::Index::operator += (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
SideIndex SideIndex::operator + (
   const int rhs) const
{
   SideIndex tmp = *this;
   tmp += rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
SideIndex& SideIndex::operator -= (
   const hier::IntVector& rhs)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   hier::Index::operator -= (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
SideIndex SideIndex::operator - (
   const hier::IntVector& rhs) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   SideIndex tmp = *this;
   tmp -= rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
SideIndex& SideIndex::operator -= (
   const int rhs)
{
   hier::Index::operator -= (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
SideIndex SideIndex::operator - (
   const int rhs) const
{
   SideIndex tmp = *this;
   tmp -= rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
SideIndex& SideIndex::operator *= (
   const hier::IntVector& rhs)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   hier::Index::operator *= (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
SideIndex SideIndex::operator * (
   const hier::IntVector& rhs) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   SideIndex tmp = *this;
   tmp *= rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
SideIndex& SideIndex::operator *= (
   const int rhs)
{
   hier::Index::operator *= (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
SideIndex SideIndex::operator * (
   const int rhs) const
{
   SideIndex tmp = *this;
   tmp *= rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
bool SideIndex::operator == (
   const SideIndex& rhs) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   return ((hier::Index *)this)->operator == (rhs) && (d_axis == rhs.d_axis);
}

SAMRAI_INLINE_KEYWORD
bool SideIndex::operator != (
   const SideIndex& rhs) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   return ((hier::Index *)this)->operator != (rhs) || (d_axis != rhs.d_axis);
}

}
}
