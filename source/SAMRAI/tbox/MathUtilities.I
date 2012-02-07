/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Utilities class to access common POSIX constants and math ops
 *
 ************************************************************************/

namespace SAMRAI {
namespace tbox {

template<class TYPE>
SAMRAI_INLINE_KEYWORD
TYPE MathUtilities<TYPE>::getZero()
{
   return s_zero;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
TYPE MathUtilities<TYPE>::getOne()
{
   return s_one;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
TYPE MathUtilities<TYPE>::getSignalingNaN()
{
   return s_signaling_nan;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
TYPE MathUtilities<TYPE>::getMax()
{
   return s_max;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
TYPE MathUtilities<TYPE>::getMin()
{
   return s_min;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
TYPE MathUtilities<TYPE>::getEpsilon()
{
   return s_epsilon;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
bool MathUtilities<TYPE>::isNaN(
   const TYPE& value)
{
   NULL_USE(value);
   return false;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
bool MathUtilities<TYPE>::equalEps(
   const TYPE& a,
   const TYPE& b)
{
   return a == b;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
TYPE MathUtilities<TYPE>::Min(
   TYPE a,
   TYPE b)
{
   return a < b ? a : b;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
TYPE MathUtilities<TYPE>::Max(
   TYPE a,
   TYPE b)
{
   return a > b ? a : b;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
TYPE MathUtilities<TYPE>::Abs(
   TYPE value)
{
   return value;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
TYPE MathUtilities<TYPE>::round(
   TYPE x)
{
   return x;
}

}
}
