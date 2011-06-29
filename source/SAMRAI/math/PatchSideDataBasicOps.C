/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Basic templated side-centered patch data operations. 
 *
 ************************************************************************/

#ifndef included_math_PatchSideDataBasicOps_C
#define included_math_PatchSideDataBasicOps_C

#include "SAMRAI/math/PatchSideDataBasicOps.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/pdat/SideGeometry.h"

namespace SAMRAI {
namespace math {

template<class TYPE>
PatchSideDataBasicOps<TYPE>::PatchSideDataBasicOps()
{
}

template<class TYPE>
PatchSideDataBasicOps<TYPE>::~PatchSideDataBasicOps()
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
PatchSideDataBasicOps<TYPE>::PatchSideDataBasicOps(
   const PatchSideDataBasicOps<TYPE>& foo)
{
   NULL_USE(foo);
}

template<class TYPE>
void PatchSideDataBasicOps<TYPE>::operator = (
   const PatchSideDataBasicOps<TYPE>& foo)
{
   NULL_USE(foo);
}

/*
 *************************************************************************
 *                                                                       *
 * General basic templated opertions for side data.                      *
 *                                                                       *
 *************************************************************************
 */

template<class TYPE>
void PatchSideDataBasicOps<TYPE>::scale(
   tbox::Pointer<pdat::SideData<TYPE> >& dst,
   const TYPE& alpha,
   const tbox::Pointer<pdat::SideData<TYPE> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
   TBOX_ASSERT(dst->getDirectionVector() == src->getDirectionVector());
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst, *src, box);

   const tbox::Dimension& dim(dst->getDim());

   const hier::IntVector& directions = dst->getDirectionVector();
   for (int d = 0; d < dim.getValue(); d++) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         d_array_ops.scale(dst->getArrayData(d),
            alpha, src->getArrayData(d),
            side_box);
      }
   }
}

template<class TYPE>
void PatchSideDataBasicOps<TYPE>::addScalar(
   tbox::Pointer<pdat::SideData<TYPE> >& dst,
   const tbox::Pointer<pdat::SideData<TYPE> >& src,
   const TYPE& alpha,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
   TBOX_ASSERT(dst->getDirectionVector() == src->getDirectionVector());
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst, *src, box);

   const tbox::Dimension& dim(dst->getDim());

   const hier::IntVector& directions = dst->getDirectionVector();
   for (int d = 0; d < dim.getValue(); d++) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         d_array_ops.addScalar(dst->getArrayData(d),
            src->getArrayData(d), alpha,
            side_box);
      }
   }
}

template<class TYPE>
void PatchSideDataBasicOps<TYPE>::add(
   tbox::Pointer<pdat::SideData<TYPE> >& dst,
   const tbox::Pointer<pdat::SideData<TYPE> >& src1,
   const tbox::Pointer<pdat::SideData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
   TBOX_ASSERT(dst->getDirectionVector() == src1->getDirectionVector());
   TBOX_ASSERT(dst->getDirectionVector() == src2->getDirectionVector());
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst, *src1, *src2, box);

   const tbox::Dimension& dim(dst->getDim());

   const hier::IntVector& directions = dst->getDirectionVector();
   for (int d = 0; d < dim.getValue(); d++) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         d_array_ops.add(dst->getArrayData(d),
            src1->getArrayData(d), src2->getArrayData(d),
            side_box);
      }
   }
}

template<class TYPE>
void PatchSideDataBasicOps<TYPE>::subtract(
   tbox::Pointer<pdat::SideData<TYPE> >& dst,
   const tbox::Pointer<pdat::SideData<TYPE> >& src1,
   const tbox::Pointer<pdat::SideData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
   TBOX_ASSERT(dst->getDirectionVector() == src1->getDirectionVector());
   TBOX_ASSERT(dst->getDirectionVector() == src2->getDirectionVector());
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst, *src1, *src2, box);

   const tbox::Dimension& dim(dst->getDim());

   const hier::IntVector& directions = dst->getDirectionVector();
   for (int d = 0; d < dim.getValue(); d++) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         d_array_ops.subtract(dst->getArrayData(d),
            src1->getArrayData(d), src2->getArrayData(d),
            side_box);
      }
   }
}

template<class TYPE>
void PatchSideDataBasicOps<TYPE>::multiply(
   tbox::Pointer<pdat::SideData<TYPE> >& dst,
   const tbox::Pointer<pdat::SideData<TYPE> >& src1,
   const tbox::Pointer<pdat::SideData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
   TBOX_ASSERT(dst->getDirectionVector() == src1->getDirectionVector());
   TBOX_ASSERT(dst->getDirectionVector() == src2->getDirectionVector());
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst, *src1, *src2, box);

   const tbox::Dimension& dim(dst->getDim());

   const hier::IntVector& directions = dst->getDirectionVector();
   for (int d = 0; d < dim.getValue(); d++) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         d_array_ops.multiply(dst->getArrayData(d),
            src1->getArrayData(d), src2->getArrayData(d),
            side_box);
      }
   }
}

template<class TYPE>
void PatchSideDataBasicOps<TYPE>::divide(
   tbox::Pointer<pdat::SideData<TYPE> >& dst,
   const tbox::Pointer<pdat::SideData<TYPE> >& src1,
   const tbox::Pointer<pdat::SideData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
   TBOX_ASSERT(dst->getDirectionVector() == src1->getDirectionVector());
   TBOX_ASSERT(dst->getDirectionVector() == src2->getDirectionVector());
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst, *src1, *src2, box);

   const tbox::Dimension& dim(dst->getDim());

   const hier::IntVector& directions = dst->getDirectionVector();
   for (int d = 0; d < dim.getValue(); d++) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         d_array_ops.divide(dst->getArrayData(d),
            src1->getArrayData(d), src2->getArrayData(d),
            side_box);
      }
   }
}

template<class TYPE>
void PatchSideDataBasicOps<TYPE>::reciprocal(
   tbox::Pointer<pdat::SideData<TYPE> >& dst,
   const tbox::Pointer<pdat::SideData<TYPE> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
   TBOX_ASSERT(dst->getDirectionVector() == src->getDirectionVector());
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst, *src, box);

   const tbox::Dimension& dim(dst->getDim());

   const hier::IntVector& directions = dst->getDirectionVector();
   for (int d = 0; d < dim.getValue(); d++) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         d_array_ops.reciprocal(dst->getArrayData(d),
            src->getArrayData(d),
            side_box);
      }
   }
}

template<class TYPE>
void PatchSideDataBasicOps<TYPE>::linearSum(
   tbox::Pointer<pdat::SideData<TYPE> >& dst,
   const TYPE& alpha,
   const tbox::Pointer<pdat::SideData<TYPE> >& src1,
   const TYPE& beta,
   const tbox::Pointer<pdat::SideData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
   TBOX_ASSERT(dst->getDirectionVector() == src1->getDirectionVector());
   TBOX_ASSERT(dst->getDirectionVector() == src2->getDirectionVector());
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst, *src1, *src2, box);

   const tbox::Dimension& dim(dst->getDim());

   const hier::IntVector& directions = dst->getDirectionVector();
   for (int d = 0; d < dim.getValue(); d++) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         d_array_ops.linearSum(dst->getArrayData(d),
            alpha, src1->getArrayData(d),
            beta, src2->getArrayData(d),
            side_box);
      }
   }
}

template<class TYPE>
void PatchSideDataBasicOps<TYPE>::axpy(
   tbox::Pointer<pdat::SideData<TYPE> >& dst,
   const TYPE& alpha,
   const tbox::Pointer<pdat::SideData<TYPE> >& src1,
   const tbox::Pointer<pdat::SideData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
   TBOX_ASSERT(dst->getDirectionVector() == src1->getDirectionVector());
   TBOX_ASSERT(dst->getDirectionVector() == src2->getDirectionVector());
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst, *src1, *src2, box);

   const tbox::Dimension& dim(dst->getDim());

   const hier::IntVector& directions = dst->getDirectionVector();
   for (int d = 0; d < dim.getValue(); d++) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         d_array_ops.axpy(dst->getArrayData(d),
            alpha, src1->getArrayData(d),
            src2->getArrayData(d),
            side_box);
      }
   }
}

template<class TYPE>
void PatchSideDataBasicOps<TYPE>::axmy(
   tbox::Pointer<pdat::SideData<TYPE> >& dst,
   const TYPE& alpha,
   const tbox::Pointer<pdat::SideData<TYPE> >& src1,
   const tbox::Pointer<pdat::SideData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
   TBOX_ASSERT(dst->getDirectionVector() == src1->getDirectionVector());
   TBOX_ASSERT(dst->getDirectionVector() == src2->getDirectionVector());
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst, *src1, *src2, box);

   const tbox::Dimension& dim(dst->getDim());

   const hier::IntVector& directions = dst->getDirectionVector();
   for (int d = 0; d < dim.getValue(); d++) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         d_array_ops.axmy(dst->getArrayData(d),
            alpha, src1->getArrayData(d),
            src2->getArrayData(d),
            side_box);
      }
   }
}

template<class TYPE>
void PatchSideDataBasicOps<TYPE>::setRandomValues(
   tbox::Pointer<pdat::SideData<TYPE> >& dst,
   const TYPE& width,
   const TYPE& low,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS2(*dst, box);

   const tbox::Dimension& dim(dst->getDim());

   const hier::IntVector& directions = dst->getDirectionVector();
   for (int d = 0; d < dim.getValue(); d++) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         d_array_ops.setRandomValues(dst->getArrayData(d),
            width, low, side_box);
      }
   }
}

template<class TYPE>
TYPE PatchSideDataBasicOps<TYPE>::min(
   const tbox::Pointer<pdat::SideData<TYPE> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(!data.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   const tbox::Dimension& dim(data->getDim());

   TYPE minval = tbox::MathUtilities<TYPE>::getMax();
   const hier::IntVector& directions = data->getDirectionVector();
   for (int d = 0; d < dim.getValue(); d++) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         minval = tbox::MathUtilities<TYPE>::Min(
               minval, d_array_ops.min(data->getArrayData(d), side_box));
      }
   }
   return minval;
}

template<class TYPE>
TYPE PatchSideDataBasicOps<TYPE>::max(
   const tbox::Pointer<pdat::SideData<TYPE> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(!data.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   const tbox::Dimension& dim(data->getDim());

   TYPE maxval = -tbox::MathUtilities<TYPE>::getMax();
   const hier::IntVector& directions = data->getDirectionVector();
   for (int d = 0; d < dim.getValue(); d++) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         maxval = tbox::MathUtilities<TYPE>::Max(
               maxval, d_array_ops.max(data->getArrayData(d), side_box));
      }
   }
   return maxval;
}

}
}
#endif
