/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Interface for the AMR Index object
 *
 ************************************************************************/

#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"

namespace SAMRAI {
namespace hier {

SAMRAI_INLINE_KEYWORD
Index&
Index::operator = (
   const Index& rhs)
{
   /*
    * Allow assignment of to an uninitialized
    * but do not allow assignment from an
    * uninitialized.
    */
   if (getDim().isValid()) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);
   } else {
      TBOX_DIM_ASSERT_CHECK_DIM(rhs.getDim());
   }

   IntVector::operator = (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
Index&
Index::operator = (
   const IntVector& rhs)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   IntVector::operator = (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
Index&
Index::operator += (
   const IntVector& rhs)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   IntVector::operator += (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
Index
Index::operator + (
   const IntVector& rhs) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   Index tmp = *this;
   tmp += rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
Index&
Index::operator += (
   const int rhs)
{
   TBOX_DIM_ASSERT_CHECK_DIM(getDim());

   IntVector::operator += (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
Index
Index::operator + (
   const int rhs) const
{
   TBOX_DIM_ASSERT_CHECK_DIM(getDim());

   Index tmp = *this;
   tmp += rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
Index&
Index::operator -= (
   const IntVector& rhs)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   IntVector::operator -= (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
Index
Index::operator - (
   const IntVector& rhs) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   Index tmp = *this;
   tmp -= rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
Index&
Index::operator -= (
   const int rhs)
{
   TBOX_DIM_ASSERT_CHECK_DIM(getDim());

   IntVector::operator -= (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
Index
Index::operator - (
   const int rhs) const
{
   TBOX_DIM_ASSERT_CHECK_DIM(getDim());

   Index tmp = *this;
   tmp -= rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
Index&
Index::operator *= (
   const IntVector& rhs)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   IntVector::operator *= (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
Index
Index::operator * (
   const IntVector& rhs) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   Index tmp = *this;
   tmp *= rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
Index&
Index::operator *= (
   const int rhs)
{
   TBOX_DIM_ASSERT_CHECK_DIM(getDim());

   IntVector::operator *= (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
Index
Index::operator * (
   const int rhs) const
{
   TBOX_DIM_ASSERT_CHECK_DIM(getDim());

   Index tmp = *this;
   tmp *= rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
Index&
Index::operator /= (
   const IntVector& rhs)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   IntVector::operator /= (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
Index
Index::operator / (
   const IntVector& rhs) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   Index tmp = *this;
   tmp /= rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
Index&
Index::operator /= (
   const int rhs)
{
   TBOX_DIM_ASSERT_CHECK_DIM(getDim());

   IntVector::operator /= (rhs);
   return *this;
}

SAMRAI_INLINE_KEYWORD
Index
Index::operator / (
   const int rhs) const
{
   TBOX_DIM_ASSERT_CHECK_DIM(getDim());

   Index tmp = *this;
   tmp /= rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
int
Index::coarsen(
   const int index,
   const int ratio)
{
   return index < 0 ? (index + 1) / ratio - 1 : index / ratio;
}

SAMRAI_INLINE_KEYWORD
Index&
Index::coarsen(
   const IntVector& ratio)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, ratio);

   for (int d = 0; d < getDim().getValue(); ++d) {
      (*this)(d) = coarsen((*this)(d), ratio(d));
   }
   return *this;
}

SAMRAI_INLINE_KEYWORD
const Index&
Index::getZeroIndex(
   const tbox::Dimension& dim) {
   return *(s_zeros[dim.getValue() - 1]);
}

SAMRAI_INLINE_KEYWORD
const Index&
Index::getOneIndex(
   const tbox::Dimension& dim) {
   return *(s_ones[dim.getValue() - 1]);
}

SAMRAI_INLINE_KEYWORD
const Index&
Index::getMinIndex(
   const tbox::Dimension& dim) {
   return *(s_mins[dim.getValue() - 1]);
}

SAMRAI_INLINE_KEYWORD
const Index&
Index::getMaxIndex(
   const tbox::Dimension& dim) {
   return *(s_maxs[dim.getValue() - 1]);
}

SAMRAI_INLINE_KEYWORD
Index
Index::coarsen(
   const Index& index,
   const IntVector& ratio)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(index, ratio);

   tbox::Dimension dim(index.getDim());

   Index tmp(dim);
   for (int d = 0; d < dim.getValue(); ++d) {
      tmp(d) = coarsen(index(d), ratio(d));
   }
   return tmp;
}

}
}
