//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-0/source/mathops/face/PatchFaceDataBasicOps.C $
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2249 $
// Modified:	$LastChangedDate: 2008-07-03 08:17:20 -0700 (Thu, 03 Jul 2008) $
// Description:	Basic templated face-centered patch data operations.
//

#ifndef included_math_PatchFaceDataBasicOps_C
#define included_math_PatchFaceDataBasicOps_C

#include "tbox/MathUtilities.h"
#include "tbox/Utilities.h"
#include "PatchFaceDataBasicOps.h"
#include "FaceGeometry.h"

namespace SAMRAI {
    namespace math {

template<int DIM, class TYPE>
PatchFaceDataBasicOps<DIM,TYPE>::PatchFaceDataBasicOps()
{
}

template<int DIM, class TYPE>
PatchFaceDataBasicOps<DIM,TYPE>::~PatchFaceDataBasicOps()
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
PatchFaceDataBasicOps<DIM,TYPE>::PatchFaceDataBasicOps(
   const PatchFaceDataBasicOps<DIM,TYPE>& foo)
{
   NULL_USE(foo);
}

template<int DIM, class TYPE>
void PatchFaceDataBasicOps<DIM,TYPE>::operator=(
   const PatchFaceDataBasicOps<DIM,TYPE>& foo)
{
   NULL_USE(foo);
}

/*
*************************************************************************
*                                                                       *
* General basic templated opertions for face data.                      *
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
void PatchFaceDataBasicOps<DIM,TYPE>::scale(
   tbox::Pointer< pdat::FaceData<DIM,TYPE> >& dst,
   const TYPE& alpha,
   const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& src,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
#endif
   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> face_box = pdat::FaceGeometry<DIM>::toFaceBox(box, d);
      d_array_ops.scale(dst->getArrayData(d),
                        alpha, src->getArrayData(d),
                        face_box);
   }
}

template<int DIM, class TYPE>
void PatchFaceDataBasicOps<DIM,TYPE>::addScalar(
   tbox::Pointer< pdat::FaceData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& src,
   const TYPE& alpha,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
#endif
   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> face_box = pdat::FaceGeometry<DIM>::toFaceBox(box, d);
      d_array_ops.addScalar(dst->getArrayData(d),
                            src->getArrayData(d), alpha,
                            face_box);
   }
}

template<int DIM, class TYPE>
void PatchFaceDataBasicOps<DIM,TYPE>::add(
   tbox::Pointer< pdat::FaceData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& src1,
   const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& src2,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
#endif
   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> face_box = pdat::FaceGeometry<DIM>::toFaceBox(box, d);
      d_array_ops.add(dst->getArrayData(d),
                      src1->getArrayData(d), src2->getArrayData(d),
                      face_box);
   }
}

template<int DIM, class TYPE>
void PatchFaceDataBasicOps<DIM,TYPE>::subtract(
   tbox::Pointer< pdat::FaceData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& src1,
   const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& src2,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
#endif
   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> face_box = pdat::FaceGeometry<DIM>::toFaceBox(box, d);
      d_array_ops.subtract(dst->getArrayData(d),
                           src1->getArrayData(d), src2->getArrayData(d),
                           face_box);
   }
}

template<int DIM, class TYPE>
void PatchFaceDataBasicOps<DIM,TYPE>::multiply(
   tbox::Pointer< pdat::FaceData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& src1,
   const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& src2,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
#endif
   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> face_box = pdat::FaceGeometry<DIM>::toFaceBox(box, d);
      d_array_ops.multiply(dst->getArrayData(d),
                           src1->getArrayData(d), src2->getArrayData(d),
                           face_box);
   }
}

template<int DIM, class TYPE>
void PatchFaceDataBasicOps<DIM,TYPE>::divide(
   tbox::Pointer< pdat::FaceData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& src1,
   const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& src2,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
#endif
   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> face_box = pdat::FaceGeometry<DIM>::toFaceBox(box, d);
      d_array_ops.divide(dst->getArrayData(d),
                         src1->getArrayData(d), src2->getArrayData(d),
                         face_box);
   }
}

template<int DIM, class TYPE>
void PatchFaceDataBasicOps<DIM,TYPE>::reciprocal(
   tbox::Pointer< pdat::FaceData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& src,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
#endif
   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> face_box = pdat::FaceGeometry<DIM>::toFaceBox(box, d);
      d_array_ops.reciprocal(dst->getArrayData(d),
                             src->getArrayData(d),
                             face_box);
   }
}

template<int DIM, class TYPE>
void PatchFaceDataBasicOps<DIM,TYPE>::linearSum(
   tbox::Pointer< pdat::FaceData<DIM,TYPE> >& dst,
   const TYPE& alpha,
   const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& src1,
   const TYPE& beta,
   const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& src2,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
#endif
   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> face_box = pdat::FaceGeometry<DIM>::toFaceBox(box, d);
      d_array_ops.linearSum(dst->getArrayData(d),
                            alpha, src1->getArrayData(d),
                            beta, src2->getArrayData(d),
                            face_box);
   }
}

template<int DIM, class TYPE>
void PatchFaceDataBasicOps<DIM,TYPE>::axpy(
   tbox::Pointer< pdat::FaceData<DIM,TYPE> >& dst,
   const TYPE& alpha,
   const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& src1,
   const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& src2,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
#endif
   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> face_box = pdat::FaceGeometry<DIM>::toFaceBox(box, d);
      d_array_ops.axpy(dst->getArrayData(d),
                       alpha, src1->getArrayData(d),
                       src2->getArrayData(d),
                       face_box);
   }
}

template<int DIM, class TYPE>
void PatchFaceDataBasicOps<DIM,TYPE>::axmy(
   tbox::Pointer< pdat::FaceData<DIM,TYPE> >& dst,
   const TYPE& alpha,
   const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& src1,
   const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& src2,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
#endif
   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> face_box = pdat::FaceGeometry<DIM>::toFaceBox(box, d);
      d_array_ops.axmy(dst->getArrayData(d),
                       alpha, src1->getArrayData(d),
                       src2->getArrayData(d),
                       face_box);
   }
}

template<int DIM, class TYPE>
void PatchFaceDataBasicOps<DIM,TYPE>::setRandomValues(
   tbox::Pointer< pdat::FaceData<DIM,TYPE> >& dst,
   const TYPE& width,
   const TYPE& low,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull());
#endif
   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> face_box = pdat::FaceGeometry<DIM>::toFaceBox(box, d);
      d_array_ops.setRandomValues(dst->getArrayData(d),
                                  width, low, face_box);
   }
}

template<int DIM, class TYPE>
TYPE PatchFaceDataBasicOps<DIM,TYPE>::min(
   const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& data,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
#endif
   TYPE minval = tbox::MathUtilities<TYPE>::getMax();
   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> face_box = pdat::FaceGeometry<DIM>::toFaceBox(box, d);
      minval = tbox::MathUtilities<TYPE>::Min(
                  minval, d_array_ops.min(data->getArrayData(d), face_box) );
   }
   return(minval);
}

template<int DIM, class TYPE>
TYPE PatchFaceDataBasicOps<DIM,TYPE>::max(
   const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& data,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
#endif
   TYPE maxval = -tbox::MathUtilities<TYPE>::getMax();
   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> face_box = pdat::FaceGeometry<DIM>::toFaceBox(box, d);
      maxval = tbox::MathUtilities<TYPE>::Max(
                  maxval, d_array_ops.max(data->getArrayData(d), face_box) );
   }
   return(maxval);
}

}
}
#endif
