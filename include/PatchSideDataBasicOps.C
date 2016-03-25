//
// File:	PatchSideDataBasicOps.C
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Basic templated side-centered patch data operations.
//

#ifndef included_math_PatchSideDataBasicOps_C
#define included_math_PatchSideDataBasicOps_C

#include "tbox/MathUtilities.h"
#include "PatchSideDataBasicOps.h"
#include "SideGeometry.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

namespace SAMRAI {
    namespace math {

template<int DIM, class TYPE>
PatchSideDataBasicOps<DIM,TYPE>::PatchSideDataBasicOps()
{
}

template<int DIM, class TYPE>
PatchSideDataBasicOps<DIM,TYPE>::~PatchSideDataBasicOps()
{
}

/*
*************************************************************************
*                                                                       *
* The const constructor and assignment operator are not actually used   *
* but are defined here for compilers that require an implementation for *
* every declaration.                                                    *
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
PatchSideDataBasicOps<DIM,TYPE>::PatchSideDataBasicOps(
   const PatchSideDataBasicOps<DIM,TYPE>& foo)
{
   (void) foo;  // not implemented (but needed by some compilers)
}

template<int DIM, class TYPE>
void PatchSideDataBasicOps<DIM,TYPE>::operator=(
   const PatchSideDataBasicOps<DIM,TYPE>& foo)
{
   (void) foo;  // not implemented (but needed by some compilers)
}

/*
*************************************************************************
*                                                                       *
* General basic templated opertions for side data.                      *
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
void PatchSideDataBasicOps<DIM,TYPE>::scale(
   tbox::Pointer< pdat::SideData<DIM,TYPE> >& dst,
   const TYPE& alpha,
   const tbox::Pointer< pdat::SideData<DIM,TYPE> >& src,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!dst.isNull() && !src.isNull());
   assert(dst->getDirectionVector() == src->getDirectionVector());
#endif
   const hier::IntVector<DIM>& directions = dst->getDirectionVector();
   for (int d = 0; d < DIM; d++) {
      if (directions(d)) {
         const hier::Box<DIM> side_box = pdat::SideGeometry<DIM>::toSideBox(box, d);
         d_array_ops.scale(dst->getArrayData(d),
                           alpha, src->getArrayData(d),
                           side_box);
      }
   }
}

template<int DIM, class TYPE>
void PatchSideDataBasicOps<DIM,TYPE>::addScalar(
   tbox::Pointer< pdat::SideData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::SideData<DIM,TYPE> >& src,
   const TYPE& alpha,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!dst.isNull() && !src.isNull());
   assert(dst->getDirectionVector() == src->getDirectionVector());
#endif
   const hier::IntVector<DIM>& directions = dst->getDirectionVector();
   for (int d = 0; d < DIM; d++) {
      if (directions(d)) {
         const hier::Box<DIM> side_box = pdat::SideGeometry<DIM>::toSideBox(box, d);
         d_array_ops.addScalar(dst->getArrayData(d),
                               src->getArrayData(d), alpha,
                               side_box);
      }
   }
}

template<int DIM, class TYPE>
void PatchSideDataBasicOps<DIM,TYPE>::add(
   tbox::Pointer< pdat::SideData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::SideData<DIM,TYPE> >& src1,
   const tbox::Pointer< pdat::SideData<DIM,TYPE> >& src2,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!dst.isNull() && !src1.isNull() && !src2.isNull());
   assert(dst->getDirectionVector() == src1->getDirectionVector());
   assert(dst->getDirectionVector() == src2->getDirectionVector());
#endif
   const hier::IntVector<DIM>& directions = dst->getDirectionVector(); 
   for (int d = 0; d < DIM; d++) {
      if (directions(d)) {      
         const hier::Box<DIM> side_box = pdat::SideGeometry<DIM>::toSideBox(box, d);
         d_array_ops.add(dst->getArrayData(d),
                         src1->getArrayData(d), src2->getArrayData(d),
                         side_box);
      }
   }
}

template<int DIM, class TYPE>
void PatchSideDataBasicOps<DIM,TYPE>::subtract(
   tbox::Pointer< pdat::SideData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::SideData<DIM,TYPE> >& src1,
   const tbox::Pointer< pdat::SideData<DIM,TYPE> >& src2,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!dst.isNull() && !src1.isNull() && !src2.isNull());
   assert(dst->getDirectionVector() == src1->getDirectionVector());
   assert(dst->getDirectionVector() == src2->getDirectionVector());
#endif
   const hier::IntVector<DIM>& directions = dst->getDirectionVector();
   for (int d = 0; d < DIM; d++) {
      if (directions(d)) {
         const hier::Box<DIM> side_box = pdat::SideGeometry<DIM>::toSideBox(box, d);
         d_array_ops.subtract(dst->getArrayData(d),
                              src1->getArrayData(d), src2->getArrayData(d),
                              side_box);
      }
   }
}

template<int DIM, class TYPE>
void PatchSideDataBasicOps<DIM,TYPE>::multiply(
   tbox::Pointer< pdat::SideData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::SideData<DIM,TYPE> >& src1,
   const tbox::Pointer< pdat::SideData<DIM,TYPE> >& src2,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!dst.isNull() && !src1.isNull() && !src2.isNull());
   assert(dst->getDirectionVector() == src1->getDirectionVector());
   assert(dst->getDirectionVector() == src2->getDirectionVector());
#endif
   const hier::IntVector<DIM>& directions = dst->getDirectionVector();
   for (int d = 0; d < DIM; d++) {
      if (directions(d)) {
         const hier::Box<DIM> side_box = pdat::SideGeometry<DIM>::toSideBox(box, d);
         d_array_ops.multiply(dst->getArrayData(d),
                              src1->getArrayData(d), src2->getArrayData(d),
                              side_box);
      }
   }
}

template<int DIM, class TYPE>
void PatchSideDataBasicOps<DIM,TYPE>::divide(
   tbox::Pointer< pdat::SideData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::SideData<DIM,TYPE> >& src1,
   const tbox::Pointer< pdat::SideData<DIM,TYPE> >& src2,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!dst.isNull() && !src1.isNull() && !src2.isNull());
   assert(dst->getDirectionVector() == src1->getDirectionVector());
   assert(dst->getDirectionVector() == src2->getDirectionVector());
#endif
   const hier::IntVector<DIM>& directions = dst->getDirectionVector();
   for (int d = 0; d < DIM; d++) {
      if (directions(d)) {
         const hier::Box<DIM> side_box = pdat::SideGeometry<DIM>::toSideBox(box, d);
         d_array_ops.divide(dst->getArrayData(d),
                            src1->getArrayData(d), src2->getArrayData(d),
                            side_box);
      }
   }
}

template<int DIM, class TYPE>
void PatchSideDataBasicOps<DIM,TYPE>::reciprocal(
   tbox::Pointer< pdat::SideData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::SideData<DIM,TYPE> >& src,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!dst.isNull() && !src.isNull());
   assert(dst->getDirectionVector() == src->getDirectionVector());
#endif
   const hier::IntVector<DIM>& directions = dst->getDirectionVector();
   for (int d = 0; d < DIM; d++) {
      if (directions(d)) {
         const hier::Box<DIM> side_box = pdat::SideGeometry<DIM>::toSideBox(box, d);
         d_array_ops.reciprocal(dst->getArrayData(d),
                                src->getArrayData(d),
                                side_box);
      }
   }
}

template<int DIM, class TYPE>
void PatchSideDataBasicOps<DIM,TYPE>::linearSum(
   tbox::Pointer< pdat::SideData<DIM,TYPE> >& dst,
   const TYPE& alpha,
   const tbox::Pointer< pdat::SideData<DIM,TYPE> >& src1,
   const TYPE& beta,
   const tbox::Pointer< pdat::SideData<DIM,TYPE> >& src2,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!dst.isNull() && !src1.isNull() && !src2.isNull());
   assert(dst->getDirectionVector() == src1->getDirectionVector());
   assert(dst->getDirectionVector() == src2->getDirectionVector());
#endif
   const hier::IntVector<DIM>& directions = dst->getDirectionVector();
   for (int d = 0; d < DIM; d++) {
      if (directions(d)) {
         const hier::Box<DIM> side_box = pdat::SideGeometry<DIM>::toSideBox(box, d);
         d_array_ops.linearSum(dst->getArrayData(d),
                               alpha, src1->getArrayData(d),
                               beta, src2->getArrayData(d),
                               side_box);
      }
   }
}

template<int DIM, class TYPE>
void PatchSideDataBasicOps<DIM,TYPE>::axpy(
   tbox::Pointer< pdat::SideData<DIM,TYPE> >& dst,
   const TYPE& alpha,
   const tbox::Pointer< pdat::SideData<DIM,TYPE> >& src1,
   const tbox::Pointer< pdat::SideData<DIM,TYPE> >& src2,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!dst.isNull() && !src1.isNull() && !src2.isNull());
   assert(dst->getDirectionVector() == src1->getDirectionVector());
   assert(dst->getDirectionVector() == src2->getDirectionVector());
#endif
   const hier::IntVector<DIM>& directions = dst->getDirectionVector();
   for (int d = 0; d < DIM; d++) {
      if (directions(d)) {
         const hier::Box<DIM> side_box = pdat::SideGeometry<DIM>::toSideBox(box, d);
         d_array_ops.axpy(dst->getArrayData(d),
                          alpha, src1->getArrayData(d),
                          src2->getArrayData(d),
                          side_box);
      }
   }
}

template<int DIM, class TYPE>
void PatchSideDataBasicOps<DIM,TYPE>::axmy(
   tbox::Pointer< pdat::SideData<DIM,TYPE> >& dst,
   const TYPE& alpha,
   const tbox::Pointer< pdat::SideData<DIM,TYPE> >& src1,
   const tbox::Pointer< pdat::SideData<DIM,TYPE> >& src2,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!dst.isNull() && !src1.isNull() && !src2.isNull());
   assert(dst->getDirectionVector() == src1->getDirectionVector());
   assert(dst->getDirectionVector() == src2->getDirectionVector());
#endif
   const hier::IntVector<DIM>& directions = dst->getDirectionVector();
   for (int d = 0; d < DIM; d++) {
      if (directions(d)) {
         const hier::Box<DIM> side_box = pdat::SideGeometry<DIM>::toSideBox(box, d);
         d_array_ops.axmy(dst->getArrayData(d),
                          alpha, src1->getArrayData(d),
                          src2->getArrayData(d),
                          side_box);
      }
   }
}

template<int DIM, class TYPE>
void PatchSideDataBasicOps<DIM,TYPE>::setRandomValues(
   tbox::Pointer< pdat::SideData<DIM,TYPE> >& dst,
   const TYPE& width,
   const TYPE& low,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!dst.isNull());
#endif
   const hier::IntVector<DIM>& directions = dst->getDirectionVector();
   for (int d = 0; d < DIM; d++) {
      if (directions(d)) {
         const hier::Box<DIM> side_box = pdat::SideGeometry<DIM>::toSideBox(box, d);
         d_array_ops.setRandomValues(dst->getArrayData(d),
                                     width, low, side_box);
      }
   }
}

template<int DIM, class TYPE>
TYPE PatchSideDataBasicOps<DIM,TYPE>::min(
   const tbox::Pointer< pdat::SideData<DIM,TYPE> >& data,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!data.isNull());
#endif
    TYPE minval = tbox::MathUtilities<TYPE>::getMax();
   const hier::IntVector<DIM>& directions = data->getDirectionVector();
   for (int d = 0; d < DIM; d++) {
      if (directions(d)) {
         const hier::Box<DIM> side_box = pdat::SideGeometry<DIM>::toSideBox(box, d);
         minval = tbox::MathUtilities<TYPE>::Min(
                  minval, d_array_ops.min(data->getArrayData(d), side_box) );
      }
   }
   return(minval);
}

template<int DIM, class TYPE>
TYPE PatchSideDataBasicOps<DIM,TYPE>::max(
   const tbox::Pointer< pdat::SideData<DIM,TYPE> >& data,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!data.isNull());
#endif
   TYPE maxval = -tbox::MathUtilities<TYPE>::getMax();
   const hier::IntVector<DIM>& directions = data->getDirectionVector();
   for (int d = 0; d < DIM; d++) {
      if (directions(d)) {
         const hier::Box<DIM> side_box = pdat::SideGeometry<DIM>::toSideBox(box, d);
         maxval = tbox::MathUtilities<TYPE>::Max(
                  maxval, d_array_ops.max(data->getArrayData(d), side_box) );
      }
   }
   return(maxval);
}

}
}
#endif
