/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   IEEE routines to set up handlers and get signaling NaNs
 *
 ************************************************************************/

#include "SAMRAI/tbox/MathUtilities.h"

namespace SAMRAI {
namespace tbox {

/*
 * Routines to get signaling NaN value.
 */

SAMRAI_INLINE_KEYWORD
float
IEEE::getSignalingFloatNaN()
{
   return MathUtilities<float>::getSignalingNaN();
}

SAMRAI_INLINE_KEYWORD
double
IEEE::getSignalingNaN()
{
   return MathUtilities<double>::getSignalingNaN();
}

SAMRAI_INLINE_KEYWORD
dcomplex
IEEE::getSignalingComplexNaN()
{
   return MathUtilities<dcomplex>::getSignalingNaN();
}

/*
 * Routines to set signaling NaN value.
 */

SAMRAI_INLINE_KEYWORD
void
IEEE::setNaN(
   float& f)
{
   f = MathUtilities<float>::getSignalingNaN();
}

SAMRAI_INLINE_KEYWORD
void
IEEE::setNaN(
   double& d)
{
   d = MathUtilities<double>::getSignalingNaN();
}

SAMRAI_INLINE_KEYWORD
void
IEEE::setNaN(
   dcomplex& dc)
{
   dc = MathUtilities<dcomplex>::getSignalingNaN();
}

/*
 * Routines to check against signaling NaN value.
 */

SAMRAI_INLINE_KEYWORD
bool
IEEE::isNaN(
   const float& f)
{
   return MathUtilities<float>::isNaN(f);
}

SAMRAI_INLINE_KEYWORD
bool
IEEE::isNaN(
   const double& d)
{
   return MathUtilities<double>::isNaN(d);
}

SAMRAI_INLINE_KEYWORD
bool
IEEE::isNaN(
   const dcomplex& dc)
{
   return MathUtilities<dcomplex>::isNaN(dc);
}

/*
 *************************************************************************
 *
 * Routines to initialize arrays to signaling NaNs.
 *
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
void
IEEE::initializeArrayToSignalingNaN(
   Array<float>& array)
{
   MathUtilities<float>::setArrayToSignalingNaN(array);
}

SAMRAI_INLINE_KEYWORD
void
IEEE::initializeArrayToSignalingNaN(
   Array<double>& array)
{
   MathUtilities<double>::setArrayToSignalingNaN(array);
}

SAMRAI_INLINE_KEYWORD
void
IEEE::initializeArrayToSignalingNaN(
   Array<dcomplex>& array)
{
   MathUtilities<dcomplex>::setArrayToSignalingNaN(array);
}

SAMRAI_INLINE_KEYWORD
void
IEEE::initializeArrayToSignalingNaN(
   float* array,
   int n)
{
   MathUtilities<float>::setArrayToSignalingNaN(array, n);
}

SAMRAI_INLINE_KEYWORD
void
IEEE::initializeArrayToSignalingNaN(
   double* array,
   int n)
{
   MathUtilities<double>::setArrayToSignalingNaN(array, n);
}

SAMRAI_INLINE_KEYWORD
void
IEEE::initializeArrayToSignalingNaN(
   dcomplex* array,
   int n)
{
   MathUtilities<dcomplex>::setArrayToSignalingNaN(array, n);
}

}
}
