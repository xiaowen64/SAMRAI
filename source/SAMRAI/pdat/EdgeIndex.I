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
EdgeIndex::EdgeIndex(
   const tbox::Dimension& dim):
   hier::Index(dim)
{
}

SAMRAI_INLINE_KEYWORD
EdgeIndex::EdgeIndex(
   const hier::Index& rhs,
   const int axis,
   const int edge):
   hier::Index(rhs),
   d_axis(axis)
{
   if (getDim() > tbox::Dimension(1)) {
      (*this)((d_axis + 1) % getDim().getValue()) += edge % 2;
   }
   for (int j = 2; j < getDim().getValue(); j++) {
      (*this)((d_axis + j) % getDim().getValue()) += (edge / (1 << (j - 1))) % 2;
   }
}

SAMRAI_INLINE_KEYWORD
EdgeIndex::EdgeIndex(
   const EdgeIndex& rhs):
   hier::Index(rhs),
   d_axis(rhs.d_axis)
{
}

SAMRAI_INLINE_KEYWORD
EdgeIndex& EdgeIndex::operator = (
   const EdgeIndex& rhs)
{
   hier::Index::operator = (rhs);
   d_axis = rhs.d_axis;
   return *this;
}

SAMRAI_INLINE_KEYWORD
EdgeIndex::~EdgeIndex()
{
}

SAMRAI_INLINE_KEYWORD
int EdgeIndex::getAxis() const
{
   return d_axis;
}

SAMRAI_INLINE_KEYWORD
void EdgeIndex::setAxis(
   const int axis)
{
   d_axis = axis;
}

SAMRAI_INLINE_KEYWORD
hier::Index EdgeIndex::toCell(
   const int edge) const
{
   const tbox::Dimension& dim(getDim());
   hier::Index index(dim);

   for (int i = 0; i < dim.getValue(); i++) {
      index(i) = (*this)(i);
   }

   if (dim > tbox::Dimension(1)) {
      index((d_axis + 1) % dim.getValue()) += ((edge % 2) - 1);
   }
   for (int j = 2; j < dim.getValue(); j++) {
      index((d_axis + j) % dim.getValue()) += (((edge / (2 << (j - 1))) % 2) - 1);
   }
   return index;
}

SAMRAI_INLINE_KEYWORD
EdgeIndex& EdgeIndex::operator += (
   const hier::IntVector& rhs)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   hier::Index::operator += (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
EdgeIndex EdgeIndex::operator + (
   const hier::IntVector& rhs) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   EdgeIndex tmp = *this;
   tmp += rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
EdgeIndex& EdgeIndex::operator += (
   const int rhs)
{
   hier::Index::operator += (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
EdgeIndex EdgeIndex::operator + (
   const int rhs) const
{
   EdgeIndex tmp = *this;
   tmp += rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
EdgeIndex& EdgeIndex::operator -= (
   const hier::IntVector& rhs)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   hier::Index::operator -= (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
EdgeIndex EdgeIndex::operator - (
   const hier::IntVector& rhs) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   EdgeIndex tmp = *this;
   tmp -= rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
EdgeIndex& EdgeIndex::operator -= (
   const int rhs)
{
   hier::Index::operator -= (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
EdgeIndex EdgeIndex::operator - (
   const int rhs) const
{
   EdgeIndex tmp = *this;
   tmp -= rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
EdgeIndex& EdgeIndex::operator *= (
   const hier::IntVector& rhs)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   hier::Index::operator *= (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
EdgeIndex EdgeIndex::operator * (
   const hier::IntVector& rhs) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   EdgeIndex tmp = *this;
   tmp *= rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
EdgeIndex& EdgeIndex::operator *= (
   const int rhs)
{
   hier::Index::operator *= (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
EdgeIndex EdgeIndex::operator * (
   const int rhs) const
{
   EdgeIndex tmp = *this;
   tmp *= rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
bool EdgeIndex::operator == (
   const EdgeIndex& rhs) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   return ((hier::Index *)this)->operator == (rhs) && (d_axis == rhs.d_axis);
}

SAMRAI_INLINE_KEYWORD
bool EdgeIndex::operator != (
   const EdgeIndex& rhs) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   return ((hier::Index *)this)->operator != (rhs) || (d_axis != rhs.d_axis);
}

}
}
