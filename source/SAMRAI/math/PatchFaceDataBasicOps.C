/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Basic templated face-centered patch data operations. 
 *
 ************************************************************************/

#ifndef included_math_PatchFaceDataBasicOps_C
#define included_math_PatchFaceDataBasicOps_C

#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/math/PatchFaceDataBasicOps.h"
#include "SAMRAI/pdat/FaceGeometry.h"

namespace SAMRAI {
namespace math {

template<class TYPE>
PatchFaceDataBasicOps<TYPE>::PatchFaceDataBasicOps()
{
}

template<class TYPE>
PatchFaceDataBasicOps<TYPE>::~PatchFaceDataBasicOps()
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

template<class TYPE>
PatchFaceDataBasicOps<TYPE>::PatchFaceDataBasicOps(
   const PatchFaceDataBasicOps<TYPE>& foo)
{
   NULL_USE(foo);
}

template<class TYPE>
void PatchFaceDataBasicOps<TYPE>::operator = (
   const PatchFaceDataBasicOps<TYPE>& foo)
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

template<class TYPE>
void PatchFaceDataBasicOps<TYPE>::scale(
   tbox::Pointer<pdat::FaceData<TYPE> >& dst,
   const TYPE& alpha,
   const tbox::Pointer<pdat::FaceData<TYPE> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst, *src, box);

   const tbox::Dimension& dim(dst->getDim());

   for (int d = 0; d < dim.getValue(); d++) {
      const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
      d_array_ops.scale(dst->getArrayData(d),
         alpha, src->getArrayData(d),
         face_box);
   }
}

template<class TYPE>
void PatchFaceDataBasicOps<TYPE>::addScalar(
   tbox::Pointer<pdat::FaceData<TYPE> >& dst,
   const tbox::Pointer<pdat::FaceData<TYPE> >& src,
   const TYPE& alpha,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst, *src, box);

   const tbox::Dimension& dim(dst->getDim());

   for (int d = 0; d < dim.getValue(); d++) {
      const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
      d_array_ops.addScalar(dst->getArrayData(d),
         src->getArrayData(d), alpha,
         face_box);
   }
}

template<class TYPE>
void PatchFaceDataBasicOps<TYPE>::add(
   tbox::Pointer<pdat::FaceData<TYPE> >& dst,
   const tbox::Pointer<pdat::FaceData<TYPE> >& src1,
   const tbox::Pointer<pdat::FaceData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst, *src1, *src2, box);

   const tbox::Dimension& dim(dst->getDim());

   for (int d = 0; d < dim.getValue(); d++) {
      const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
      d_array_ops.add(dst->getArrayData(d),
         src1->getArrayData(d), src2->getArrayData(d),
         face_box);
   }
}

template<class TYPE>
void PatchFaceDataBasicOps<TYPE>::subtract(
   tbox::Pointer<pdat::FaceData<TYPE> >& dst,
   const tbox::Pointer<pdat::FaceData<TYPE> >& src1,
   const tbox::Pointer<pdat::FaceData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst, *src1, *src2, box);

   const tbox::Dimension& dim(dst->getDim());

   for (int d = 0; d < dim.getValue(); d++) {
      const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
      d_array_ops.subtract(dst->getArrayData(d),
         src1->getArrayData(d), src2->getArrayData(d),
         face_box);
   }
}

template<class TYPE>
void PatchFaceDataBasicOps<TYPE>::multiply(
   tbox::Pointer<pdat::FaceData<TYPE> >& dst,
   const tbox::Pointer<pdat::FaceData<TYPE> >& src1,
   const tbox::Pointer<pdat::FaceData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst, *src1, *src2, box);

   const tbox::Dimension& dim(dst->getDim());

   for (int d = 0; d < dim.getValue(); d++) {
      const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
      d_array_ops.multiply(dst->getArrayData(d),
         src1->getArrayData(d), src2->getArrayData(d),
         face_box);
   }
}

template<class TYPE>
void PatchFaceDataBasicOps<TYPE>::divide(
   tbox::Pointer<pdat::FaceData<TYPE> >& dst,
   const tbox::Pointer<pdat::FaceData<TYPE> >& src1,
   const tbox::Pointer<pdat::FaceData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst, *src1, *src2, box);

   const tbox::Dimension& dim(dst->getDim());

   for (int d = 0; d < dim.getValue(); d++) {
      const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
      d_array_ops.divide(dst->getArrayData(d),
         src1->getArrayData(d), src2->getArrayData(d),
         face_box);
   }
}

template<class TYPE>
void PatchFaceDataBasicOps<TYPE>::reciprocal(
   tbox::Pointer<pdat::FaceData<TYPE> >& dst,
   const tbox::Pointer<pdat::FaceData<TYPE> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst, *src, box);

   const tbox::Dimension& dim(dst->getDim());

   for (int d = 0; d < dim.getValue(); d++) {
      const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
      d_array_ops.reciprocal(dst->getArrayData(d),
         src->getArrayData(d),
         face_box);
   }
}

template<class TYPE>
void PatchFaceDataBasicOps<TYPE>::linearSum(
   tbox::Pointer<pdat::FaceData<TYPE> >& dst,
   const TYPE& alpha,
   const tbox::Pointer<pdat::FaceData<TYPE> >& src1,
   const TYPE& beta,
   const tbox::Pointer<pdat::FaceData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst, *src1, *src2, box);

   const tbox::Dimension& dim(dst->getDim());

   for (int d = 0; d < dim.getValue(); d++) {
      const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
      d_array_ops.linearSum(dst->getArrayData(d),
         alpha, src1->getArrayData(d),
         beta, src2->getArrayData(d),
         face_box);
   }
}

template<class TYPE>
void PatchFaceDataBasicOps<TYPE>::axpy(
   tbox::Pointer<pdat::FaceData<TYPE> >& dst,
   const TYPE& alpha,
   const tbox::Pointer<pdat::FaceData<TYPE> >& src1,
   const tbox::Pointer<pdat::FaceData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst, *src1, *src2, box);

   const tbox::Dimension& dim(dst->getDim());

   for (int d = 0; d < dim.getValue(); d++) {
      const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
      d_array_ops.axpy(dst->getArrayData(d),
         alpha, src1->getArrayData(d),
         src2->getArrayData(d),
         face_box);
   }
}

template<class TYPE>
void PatchFaceDataBasicOps<TYPE>::axmy(
   tbox::Pointer<pdat::FaceData<TYPE> >& dst,
   const TYPE& alpha,
   const tbox::Pointer<pdat::FaceData<TYPE> >& src1,
   const tbox::Pointer<pdat::FaceData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst, *src1, *src2, box);

   const tbox::Dimension& dim(dst->getDim());

   for (int d = 0; d < dim.getValue(); d++) {
      const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
      d_array_ops.axmy(dst->getArrayData(d),
         alpha, src1->getArrayData(d),
         src2->getArrayData(d),
         face_box);
   }
}

template<class TYPE>
void PatchFaceDataBasicOps<TYPE>::setRandomValues(
   tbox::Pointer<pdat::FaceData<TYPE> >& dst,
   const TYPE& width,
   const TYPE& low,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS2(*dst, box);

   const tbox::Dimension& dim(dst->getDim());

   for (int d = 0; d < dim.getValue(); d++) {
      const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
      d_array_ops.setRandomValues(dst->getArrayData(d),
         width, low, face_box);
   }
}

template<class TYPE>
TYPE PatchFaceDataBasicOps<TYPE>::min(
   const tbox::Pointer<pdat::FaceData<TYPE> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(!data.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   const tbox::Dimension& dim(data->getDim());

   TYPE minval = tbox::MathUtilities<TYPE>::getMax();
   for (int d = 0; d < dim.getValue(); d++) {
      const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
      minval = tbox::MathUtilities<TYPE>::Min(
            minval, d_array_ops.min(data->getArrayData(d), face_box));
   }
   return minval;
}

template<class TYPE>
TYPE PatchFaceDataBasicOps<TYPE>::max(
   const tbox::Pointer<pdat::FaceData<TYPE> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(!data.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   const tbox::Dimension& dim(data->getDim());

   TYPE maxval = -tbox::MathUtilities<TYPE>::getMax();
   for (int d = 0; d < dim.getValue(); d++) {
      const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
      maxval = tbox::MathUtilities<TYPE>::Max(
            maxval, d_array_ops.max(data->getArrayData(d), face_box));
   }
   return maxval;
}

}
}
#endif
