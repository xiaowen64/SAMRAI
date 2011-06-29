/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Utilities class to access common POSIX constants and math ops 
 *
 ************************************************************************/

#ifndef included_tbox_MathUtilities_C
#define included_tbox_MathUtilities_C

#include "SAMRAI/tbox/MathUtilities.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/tbox/MathUtilities.I"
#endif

namespace SAMRAI {
namespace tbox {

/*
 *************************************************************************
 *                                                                       *
 * Routines to initialize arrays to signaling NaNs.                      *
 *                                                                       *
 *************************************************************************
 */

template<class TYPE>
void MathUtilities<TYPE>::setArrayToSignalingNaN(
   Array<TYPE>& array)
{
   for (int i = 0; i < array.getSize(); i++) {
      array[i] = getSignalingNaN();
   }
}

template<class TYPE>
void MathUtilities<TYPE>::setArrayToSignalingNaN(
   TYPE* array,
   int n)
{
   for (int i = 0; i < n; i++) {
      array[i] = getSignalingNaN();
   }
}

/*
 *************************************************************************
 *                                                                       *
 * Routines to initialize arrays to max value for type.                  *
 *                                                                       *
 *************************************************************************
 */

template<class TYPE>
void MathUtilities<TYPE>::setArrayToMax(
   Array<TYPE>& array)
{
   for (int i = 0; i < array.getSize(); i++) {
      array[i] = getMax();
   }
}

template<class TYPE>
void MathUtilities<TYPE>::setArrayToMax(
   TYPE* array,
   int n)
{
   for (int i = 0; i < n; i++) {
      array[i] = getMax();
   }
}

/*
 *************************************************************************
 *                                                                       *
 * Routines to initialize arrays to min value for type.                  *
 *                                                                       *
 *************************************************************************
 */

template<class TYPE>
void MathUtilities<TYPE>::setArrayToMin(
   Array<TYPE>& array)
{
   for (int i = 0; i < array.getSize(); i++) {
      array[i] = getMin();
   }
}

template<class TYPE>
void MathUtilities<TYPE>::setArrayToMin(
   TYPE* array,
   int n)
{
   for (int i = 0; i < n; i++) {
      array[i] = getMin();
   }
}

/*
 *************************************************************************
 *                                                                       *
 * Routines to initialize arrays to epsilon value for type.              *
 *                                                                       *
 *************************************************************************
 */

template<class TYPE>
void MathUtilities<TYPE>::setArrayToEpsilon(
   Array<TYPE>& array)
{
   for (int i = 0; i < array.getSize(); i++) {
      array[i] = getEpsilon();
   }
}

template<class TYPE>
void MathUtilities<TYPE>::setArrayToEpsilon(
   TYPE* array,
   int n)
{
   for (int i = 0; i < n; i++) {
      array[i] = getEpsilon();
   }
}

}
}

#endif
