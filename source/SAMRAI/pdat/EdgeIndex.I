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
EdgeIndex&
EdgeIndex::operator = (
   const EdgeIndex& rhs)
{
   hier::Index::operator = (rhs);
   d_axis = rhs.d_axis;
   return *this;
}

SAMRAI_INLINE_KEYWORD
int
EdgeIndex::getAxis() const
{
   return d_axis;
}

SAMRAI_INLINE_KEYWORD
void
EdgeIndex::setAxis(
   const int axis)
{
   d_axis = axis;
}

SAMRAI_INLINE_KEYWORD
EdgeIndex&
EdgeIndex::operator += (
   const hier::IntVector& rhs)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   hier::Index::operator += (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
EdgeIndex
EdgeIndex::operator + (
   const hier::IntVector& rhs) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   EdgeIndex tmp = *this;
   tmp += rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
EdgeIndex&
EdgeIndex::operator += (
   const int rhs)
{
   hier::Index::operator += (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
EdgeIndex
EdgeIndex::operator + (
   const int rhs) const
{
   EdgeIndex tmp = *this;
   tmp += rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
EdgeIndex&
EdgeIndex::operator -= (
   const hier::IntVector& rhs)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   hier::Index::operator -= (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
EdgeIndex
EdgeIndex::operator - (
   const hier::IntVector& rhs) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   EdgeIndex tmp = *this;
   tmp -= rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
EdgeIndex&
EdgeIndex::operator -= (
   const int rhs)
{
   hier::Index::operator -= (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
EdgeIndex
EdgeIndex::operator - (
   const int rhs) const
{
   EdgeIndex tmp = *this;
   tmp -= rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
EdgeIndex&
EdgeIndex::operator *= (
   const hier::IntVector& rhs)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   hier::Index::operator *= (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
EdgeIndex
EdgeIndex::operator * (
   const hier::IntVector& rhs) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   EdgeIndex tmp = *this;
   tmp *= rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
EdgeIndex&
EdgeIndex::operator *= (
   const int rhs)
{
   hier::Index::operator *= (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
EdgeIndex
EdgeIndex::operator * (
   const int rhs) const
{
   EdgeIndex tmp = *this;
   tmp *= rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
bool
EdgeIndex::operator == (
   const EdgeIndex& rhs) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   return ((hier::Index *)this)->operator == (rhs) && (d_axis == rhs.d_axis);
}

SAMRAI_INLINE_KEYWORD
bool
EdgeIndex::operator != (
   const EdgeIndex& rhs) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   return ((hier::Index *)this)->operator != (rhs) || (d_axis != rhs.d_axis);
}

}
}
