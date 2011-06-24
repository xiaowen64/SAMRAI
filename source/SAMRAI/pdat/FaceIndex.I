/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   hier 
 *
 ************************************************************************/

#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace pdat {

SAMRAI_INLINE_KEYWORD
FaceIndex::FaceIndex(
   const tbox::Dimension& dim):
   hier::Index(dim)
{
}

SAMRAI_INLINE_KEYWORD
FaceIndex::FaceIndex(
   const hier::Index& rhs,
   const int axis,
   const int face):
   hier::Index(rhs),
   d_axis(axis)
{
   (*this)(0) = rhs(d_axis) + face;
   for (int i = 1; i < getDim().getValue(); i++) {
      (*this)(i) = rhs((d_axis + i) % getDim().getValue());
   }
}

SAMRAI_INLINE_KEYWORD
FaceIndex::FaceIndex(
   const FaceIndex& rhs):
   hier::Index(rhs),
   d_axis(rhs.d_axis)
{
}

SAMRAI_INLINE_KEYWORD
FaceIndex& FaceIndex::operator = (
   const FaceIndex& rhs)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   hier::Index::operator = (rhs);
   d_axis = rhs.d_axis;
   return *this;
}

SAMRAI_INLINE_KEYWORD
FaceIndex::~FaceIndex()
{
}

SAMRAI_INLINE_KEYWORD
int FaceIndex::getAxis() const
{
   return d_axis;
}

SAMRAI_INLINE_KEYWORD
void FaceIndex::setAxis(
   const int axis)
{
   d_axis = axis;
}

SAMRAI_INLINE_KEYWORD
hier::Index FaceIndex::toCell(
   const int face) const
{
   hier::Index index(getDim());
   index(d_axis) = (*this)(0) + face - 1;
   for (int i = 1; i < getDim().getValue(); i++) {
      index((d_axis + i) % getDim().getValue()) = (*this)(i);
   }
   return index;
}

SAMRAI_INLINE_KEYWORD
FaceIndex& FaceIndex::operator += (
   const hier::IntVector& rhs)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   hier::Index::operator += (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
FaceIndex FaceIndex::operator + (
   const hier::IntVector& rhs) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   FaceIndex tmp = *this;
   tmp += rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
FaceIndex& FaceIndex::operator += (
   const int rhs)
{
   hier::Index::operator += (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
FaceIndex FaceIndex::operator + (
   const int rhs) const
{
   FaceIndex tmp = *this;
   tmp += rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
FaceIndex& FaceIndex::operator -= (
   const hier::IntVector& rhs)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   hier::Index::operator -= (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
FaceIndex FaceIndex::operator - (
   const hier::IntVector& rhs) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   FaceIndex tmp = *this;
   tmp -= rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
FaceIndex& FaceIndex::operator -= (
   const int rhs)
{
   hier::Index::operator -= (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
FaceIndex FaceIndex::operator - (
   const int rhs) const
{
   FaceIndex tmp = *this;
   tmp -= rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
FaceIndex& FaceIndex::operator *= (
   const hier::IntVector& rhs)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   hier::Index::operator *= (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
FaceIndex FaceIndex::operator * (
   const hier::IntVector& rhs) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   FaceIndex tmp = *this;
   tmp *= rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
FaceIndex& FaceIndex::operator *= (
   const int rhs)
{
   hier::Index::operator *= (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
FaceIndex FaceIndex::operator * (
   const int rhs) const
{
   FaceIndex tmp = *this;
   tmp *= rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
bool FaceIndex::operator == (
   const FaceIndex& rhs) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   return ((hier::Index *)this)->operator == (rhs) && (d_axis == rhs.d_axis);
}

SAMRAI_INLINE_KEYWORD
bool FaceIndex::operator != (
   const FaceIndex& rhs) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   return ((hier::Index *)this)->operator != (rhs) || (d_axis != rhs.d_axis);
}

}
}
