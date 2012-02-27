/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   A n-dimensional integer vector
 *
 ************************************************************************/

#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"

namespace SAMRAI {
namespace hier {

SAMRAI_INLINE_KEYWORD
const tbox::Dimension&
IntVector::getDim() const
{
   return d_dim;
}

SAMRAI_INLINE_KEYWORD
IntVector& IntVector::operator = (
   const IntVector& rhs)
{
   /*
    * Allow assignment of to an uninitialized but do not allow
    * assignment from an uninitialized.
    */
   if (d_dim.isValid()) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);
   } else {
      TBOX_DIM_ASSERT_CHECK_DIM(rhs.getDim());
      d_dim = rhs.getDim();
   }

   for (int i = 0; i < d_dim.getValue(); i++)
      d_vector[i] = rhs.d_vector[i];

#ifdef DEBUG_INITIALIZE_UNDEFINED
   for (int i = d_dim.getValue(); i < tbox::Dimension::MAXIMUM_DIMENSION_VALUE; i++) {
      d_vector[i] = tbox::MathUtilities<int>::getMin();
   }
#endif

   return *this;
}

SAMRAI_INLINE_KEYWORD
int&
IntVector::operator [] (
   const int i)
{
   TBOX_DIM_ASSERT_CHECK_DIM(d_dim);
   TBOX_ASSERT(i >= 0 && i < d_dim.getValue());
   return d_vector[i];
}

SAMRAI_INLINE_KEYWORD
const int&
IntVector::operator [] (
   const int i) const
{
   TBOX_DIM_ASSERT_CHECK_DIM(d_dim);
   TBOX_ASSERT(i >= 0 && i < d_dim.getValue());
   return d_vector[i];
}

SAMRAI_INLINE_KEYWORD
int&
IntVector::operator () (
   const int i)
{
   TBOX_DIM_ASSERT_CHECK_DIM(d_dim);
   TBOX_ASSERT(i >= 0 && i < d_dim.getValue());
   return d_vector[i];
}

SAMRAI_INLINE_KEYWORD
const int&
IntVector::operator () (
   const int i) const
{
   TBOX_DIM_ASSERT_CHECK_DIM(d_dim);
   TBOX_ASSERT(i >= 0 && i < d_dim.getValue());
   return d_vector[i];
}

SAMRAI_INLINE_KEYWORD
IntVector&
IntVector::operator += (
   const IntVector& rhs)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   for (int i = 0; i < d_dim.getValue(); i++)
      d_vector[i] += rhs.d_vector[i];

   return *this;
}

SAMRAI_INLINE_KEYWORD
IntVector
IntVector::operator + (
   const IntVector& rhs) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   IntVector tmp(*this);
   tmp += rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
IntVector&
IntVector::operator += (
   const int rhs)
{
   TBOX_DIM_ASSERT_CHECK_DIM(d_dim);

   for (int i = 0; i < d_dim.getValue(); i++)
      d_vector[i] += rhs;

   return *this;
}

SAMRAI_INLINE_KEYWORD
IntVector
IntVector::operator + (
   const int rhs) const
{
   TBOX_DIM_ASSERT_CHECK_DIM(d_dim);

   IntVector tmp = *this;
   tmp += rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
IntVector&
IntVector::operator -= (
   const IntVector& rhs)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   for (int i = 0; i < d_dim.getValue(); i++)
      d_vector[i] -= rhs.d_vector[i];

   return *this;
}

SAMRAI_INLINE_KEYWORD
IntVector
IntVector::operator - (
   const IntVector& rhs) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   IntVector tmp = *this;
   tmp -= rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
IntVector&
IntVector::operator -= (
   const int rhs)
{
   TBOX_DIM_ASSERT_CHECK_DIM(d_dim);

   for (int i = 0; i < d_dim.getValue(); i++)
      d_vector[i] -= rhs;

   return *this;
}

SAMRAI_INLINE_KEYWORD
IntVector
IntVector::operator - (
   const int rhs) const
{
   TBOX_DIM_ASSERT_CHECK_DIM(d_dim);

   IntVector tmp = *this;
   tmp -= rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
IntVector&
IntVector::operator *= (
   const IntVector& rhs)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   for (int i = 0; i < d_dim.getValue(); i++)
      d_vector[i] *= rhs.d_vector[i];

   return *this;
}

SAMRAI_INLINE_KEYWORD
IntVector
IntVector::operator * (
   const IntVector& rhs) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   IntVector tmp = *this;
   tmp *= rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
IntVector&
IntVector::operator *= (
   const int rhs)
{
   TBOX_DIM_ASSERT_CHECK_DIM(d_dim);

   for (int i = 0; i < d_dim.getValue(); i++)
      d_vector[i] *= rhs;

   return *this;
}

SAMRAI_INLINE_KEYWORD
IntVector
IntVector::operator * (
   const int rhs) const
{
   TBOX_DIM_ASSERT_CHECK_DIM(d_dim);

   IntVector tmp = *this;
   tmp *= rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
IntVector&
IntVector::operator /= (
   const IntVector& rhs)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   for (int i = 0; i < d_dim.getValue(); i++)
      d_vector[i] /= rhs.d_vector[i];

   return *this;
}

SAMRAI_INLINE_KEYWORD
void
IntVector::ceilingDivide(
   const IntVector& denominator)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, denominator);

   for (int i = 0; i < getDim().getValue(); i++) {
      /*
       * This is the formula for integer divide, rounding away from
       * zero.  It is meant as an extension of the ceilingDivide quotient of
       * 2 positive integers.
       *
       * The ceilingDivide is the integer divide plus 0, -1 or 1 representing
       * the results of rounding.
       * - Add zero if there's no remainder to round.
       * - Round remainder to 1 if numerator and denominator has same sign.
       * - Round remainder to -1 if numerator and denominator has opposite sign.
       */
      d_vector[i] = (d_vector[i] / denominator[i])
         + ((d_vector[i] % denominator[i]) ?
            ((d_vector[i] > 0) == (denominator[i] > 0) ? 1 : -1) : 0);
   }
}

SAMRAI_INLINE_KEYWORD
IntVector
IntVector::ceilingDivide(
   const IntVector& numerator,
   const IntVector& denominator)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(numerator, denominator);

   IntVector rval(numerator.getDim());
   for (int i = 0; i < numerator.getDim().getValue(); i++) {
      /*
       * This is the formula for integer divide, rounding away from
       * zero.  It is meant as an extension of the ceilingDivide quotient of
       * 2 positive integers.
       *
       * The ceilingDivide is the integer divide plus 0, -1 or 1 representing
       * the results of rounding.
       * - Add zero if there's no remainder to round.
       * - Round remainder to 1 if numerator and denominator has same sign.
       * - Round remainder to -1 if numerator and denominator has opposite sign.
       */
      rval[i] = (numerator[i] / denominator[i])
         + ((numerator[i] % denominator[i]) ?
            ((numerator[i] > 0) == (denominator[i] > 0) ? 1 : -1) : 0);
   }
   return rval;
}

SAMRAI_INLINE_KEYWORD
IntVector
IntVector::operator / (
   const IntVector& rhs) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   IntVector tmp = *this;
   tmp /= rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
IntVector&
IntVector::operator /= (
   const int rhs)
{
   TBOX_DIM_ASSERT_CHECK_DIM(d_dim);

   for (int i = 0; i < d_dim.getValue(); i++)
      d_vector[i] /= rhs;

   return *this;
}

SAMRAI_INLINE_KEYWORD
IntVector
IntVector::operator / (
   const int rhs) const
{
   TBOX_DIM_ASSERT_CHECK_DIM(d_dim);

   IntVector tmp = *this;
   tmp /= rhs;
   return tmp;
}

SAMRAI_INLINE_KEYWORD
IntVector
IntVector::operator - () const
{
   TBOX_DIM_ASSERT_CHECK_DIM(d_dim);

   IntVector tmp(d_dim);
   for (int i = 0; i < d_dim.getValue(); i++)
      tmp.d_vector[i] = -d_vector[i];

   return tmp;
}

SAMRAI_INLINE_KEYWORD
bool
IntVector::operator == (
   int rhs) const
{
   TBOX_DIM_ASSERT_CHECK_DIM(d_dim);

   bool result = true;

   for (int i = 0; result && (i < d_dim.getValue()); i++)
      result = d_vector[i] == rhs;

   return result;
}

SAMRAI_INLINE_KEYWORD
bool
IntVector::operator != (
   int rhs) const
{
   TBOX_DIM_ASSERT_CHECK_DIM(d_dim);

   bool result = true;

   for (int i = 0; result && (i < d_dim.getValue()); i++)
      result = d_vector[i] != rhs;

   return result;
}

SAMRAI_INLINE_KEYWORD
bool
IntVector::operator == (
   const IntVector& rhs) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   bool result = true;

   for (int i = 0; result && (i < d_dim.getValue()); i++)
      result = result && (d_vector[i] == rhs.d_vector[i]);

   return result;
}

SAMRAI_INLINE_KEYWORD
bool
IntVector::operator != (
   const IntVector& rhs) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   return !(*this == rhs);
}

SAMRAI_INLINE_KEYWORD
void
IntVector::min(
   IntVector const& rhs)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   for (int i = 0; i < d_dim.getValue(); i++)
      if (rhs.d_vector[i] < d_vector[i])
         d_vector[i] = rhs.d_vector[i];

}

SAMRAI_INLINE_KEYWORD
bool
IntVector::operator < (
   const IntVector& rhs) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   bool result = true;

   for (int i = 0; result && (i < d_dim.getValue()); i++)
      result = result && (d_vector[i] < rhs.d_vector[i]);

   return result;
}

SAMRAI_INLINE_KEYWORD
bool
IntVector::operator <= (
   const IntVector& rhs) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   bool result = true;

   for (int i = 0; result && (i < d_dim.getValue()); i++)
      result = result && (d_vector[i] <= rhs.d_vector[i]);

   return result;

}

SAMRAI_INLINE_KEYWORD
bool
IntVector::operator > (
   const IntVector& rhs) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   bool result = true;

   for (int i = 0; result && (i < d_dim.getValue()); i++)
      result = result && (d_vector[i] > rhs.d_vector[i]);

   return result;
}

SAMRAI_INLINE_KEYWORD
bool
IntVector::operator >= (
   const IntVector& rhs) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   bool result = true;

   for (int i = 0; result && (i < d_dim.getValue()); i++)
      result = result && (d_vector[i] >= rhs.d_vector[i]);

   return result;
}

SAMRAI_INLINE_KEYWORD
int
IntVector::min() const
{
   TBOX_DIM_ASSERT_CHECK_DIM(d_dim);

   int min = d_vector[0];

   for (int i = 1; i < d_dim.getValue(); i++)
      if (d_vector[i] < min)
         min = d_vector[i];

   return min;
}

SAMRAI_INLINE_KEYWORD
void
IntVector::max(
   const IntVector& rhs)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   for (int i = 0; i < d_dim.getValue(); i++)
      if (rhs.d_vector[i] > d_vector[i])
         d_vector[i] = rhs.d_vector[i];
}

SAMRAI_INLINE_KEYWORD
int
IntVector::max() const
{
   TBOX_DIM_ASSERT_CHECK_DIM(d_dim);

   int max = d_vector[0];

   for (int i = 1; i < d_dim.getValue(); i++)
      if (d_vector[i] > max)
         max = d_vector[i];

   return max;
}

SAMRAI_INLINE_KEYWORD
IntVector
IntVector::min(
   const IntVector& a,
   const IntVector& b)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(a, b);

   IntVector tmp = a;
   tmp.min(b);
   return tmp;
}

SAMRAI_INLINE_KEYWORD
IntVector
IntVector::max(
   const IntVector& a,
   const IntVector& b)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(a, b);

   IntVector tmp = a;
   tmp.max(b);
   return tmp;
}

SAMRAI_INLINE_KEYWORD
int
IntVector::getProduct() const
{
   int prod = 1;
   for (int i = 0; i < d_dim.getValue(); i++) {
      prod *= d_vector[i];
   }
   return prod;
}

SAMRAI_INLINE_KEYWORD
const IntVector&
IntVector::getZero(
   const tbox::Dimension& dim) {
   return *(s_zeros[dim.getValue() - 1]);
}

SAMRAI_INLINE_KEYWORD
const IntVector&
IntVector::getOne(
   const tbox::Dimension& dim) {
   return *(s_ones[dim.getValue() - 1]);
}

}
}
