/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Basic templated edge-centered patch data operations.
 *
 ************************************************************************/

#ifndef included_math_PatchEdgeDataBasicOps_C
#define included_math_PatchEdgeDataBasicOps_C

#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/math/PatchEdgeDataBasicOps.h"
#include "SAMRAI/pdat/EdgeGeometry.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace math {

template<class TYPE>
PatchEdgeDataBasicOps<TYPE>::PatchEdgeDataBasicOps()
{
}

template<class TYPE>
PatchEdgeDataBasicOps<TYPE>::~PatchEdgeDataBasicOps()
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
PatchEdgeDataBasicOps<TYPE>::PatchEdgeDataBasicOps(
   const PatchEdgeDataBasicOps<TYPE>& foo)
{
   NULL_USE(foo);
}

template<class TYPE>
void PatchEdgeDataBasicOps<TYPE>::operator = (
   const PatchEdgeDataBasicOps<TYPE>& foo)
{
   NULL_USE(foo);
}

/*
 *************************************************************************
 *                                                                       *
 * General basic templated operations for edge data.                     *
 *                                                                       *
 *************************************************************************
 */

template<class TYPE>
void PatchEdgeDataBasicOps<TYPE>::scale(
   tbox::Pointer<pdat::EdgeData<TYPE> >& dst,
   const TYPE& alpha,
   const tbox::Pointer<pdat::EdgeData<TYPE> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst, *src, box);

   const tbox::Dimension& dim(box.getDim());

   for (int d = 0; d < dim.getValue(); d++) {
      const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
      d_array_ops.scale(dst->getArrayData(d),
         alpha, src->getArrayData(d),
         edge_box);
   }
}

template<class TYPE>
void PatchEdgeDataBasicOps<TYPE>::addScalar(
   tbox::Pointer<pdat::EdgeData<TYPE> >& dst,
   const tbox::Pointer<pdat::EdgeData<TYPE> >& src,
   const TYPE& alpha,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst, *src, box);

   const tbox::Dimension& dim(box.getDim());

   for (int d = 0; d < dim.getValue(); d++) {
      const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
      d_array_ops.addScalar(dst->getArrayData(d),
         src->getArrayData(d), alpha,
         edge_box);
   }
}

template<class TYPE>
void PatchEdgeDataBasicOps<TYPE>::add(
   tbox::Pointer<pdat::EdgeData<TYPE> >& dst,
   const tbox::Pointer<pdat::EdgeData<TYPE> >& src1,
   const tbox::Pointer<pdat::EdgeData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst, *src1, *src2, box);

   const tbox::Dimension& dim(box.getDim());

   for (int d = 0; d < dim.getValue(); d++) {
      const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
      d_array_ops.add(dst->getArrayData(d),
         src1->getArrayData(d), src2->getArrayData(d),
         edge_box);
   }
}

template<class TYPE>
void PatchEdgeDataBasicOps<TYPE>::subtract(
   tbox::Pointer<pdat::EdgeData<TYPE> >& dst,
   const tbox::Pointer<pdat::EdgeData<TYPE> >& src1,
   const tbox::Pointer<pdat::EdgeData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst, *src1, *src2, box);

   const tbox::Dimension& dim(box.getDim());

   for (int d = 0; d < dim.getValue(); d++) {
      const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
      d_array_ops.subtract(dst->getArrayData(d),
         src1->getArrayData(d), src2->getArrayData(d),
         edge_box);
   }
}

template<class TYPE>
void PatchEdgeDataBasicOps<TYPE>::multiply(
   tbox::Pointer<pdat::EdgeData<TYPE> >& dst,
   const tbox::Pointer<pdat::EdgeData<TYPE> >& src1,
   const tbox::Pointer<pdat::EdgeData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst, *src1, *src2, box);

   const tbox::Dimension& dim(box.getDim());

   for (int d = 0; d < dim.getValue(); d++) {
      const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
      d_array_ops.multiply(dst->getArrayData(d),
         src1->getArrayData(d), src2->getArrayData(d),
         edge_box);
   }
}

template<class TYPE>
void PatchEdgeDataBasicOps<TYPE>::divide(
   tbox::Pointer<pdat::EdgeData<TYPE> >& dst,
   const tbox::Pointer<pdat::EdgeData<TYPE> >& src1,
   const tbox::Pointer<pdat::EdgeData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst, *src1, *src2, box);

   const tbox::Dimension& dim(box.getDim());

   for (int d = 0; d < dim.getValue(); d++) {
      const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
      d_array_ops.divide(dst->getArrayData(d),
         src1->getArrayData(d), src2->getArrayData(d),
         edge_box);
   }
}

template<class TYPE>
void PatchEdgeDataBasicOps<TYPE>::reciprocal(
   tbox::Pointer<pdat::EdgeData<TYPE> >& dst,
   const tbox::Pointer<pdat::EdgeData<TYPE> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst, *src, box);

   const tbox::Dimension& dim(box.getDim());

   for (int d = 0; d < dim.getValue(); d++) {
      const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
      d_array_ops.reciprocal(dst->getArrayData(d),
         src->getArrayData(d),
         edge_box);
   }
}

template<class TYPE>
void PatchEdgeDataBasicOps<TYPE>::linearSum(
   tbox::Pointer<pdat::EdgeData<TYPE> >& dst,
   const TYPE& alpha,
   const tbox::Pointer<pdat::EdgeData<TYPE> >& src1,
   const TYPE& beta,
   const tbox::Pointer<pdat::EdgeData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst, *src1, *src2, box);

   const tbox::Dimension& dim(box.getDim());

   for (int d = 0; d < dim.getValue(); d++) {
      const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
      d_array_ops.linearSum(dst->getArrayData(d),
         alpha, src1->getArrayData(d),
         beta, src2->getArrayData(d),
         edge_box);
   }
}

template<class TYPE>
void PatchEdgeDataBasicOps<TYPE>::axpy(
   tbox::Pointer<pdat::EdgeData<TYPE> >& dst,
   const TYPE& alpha,
   const tbox::Pointer<pdat::EdgeData<TYPE> >& src1,
   const tbox::Pointer<pdat::EdgeData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst, *src1, *src2, box);

   const tbox::Dimension& dim(box.getDim());

   for (int d = 0; d < dim.getValue(); d++) {
      const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
      d_array_ops.axpy(dst->getArrayData(d),
         alpha, src1->getArrayData(d),
         src2->getArrayData(d),
         edge_box);
   }
}

template<class TYPE>
void PatchEdgeDataBasicOps<TYPE>::axmy(
   tbox::Pointer<pdat::EdgeData<TYPE> >& dst,
   const TYPE& alpha,
   const tbox::Pointer<pdat::EdgeData<TYPE> >& src1,
   const tbox::Pointer<pdat::EdgeData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst, *src1, *src2, box);

   const tbox::Dimension& dim(box.getDim());

   for (int d = 0; d < dim.getValue(); d++) {
      const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
      d_array_ops.axmy(dst->getArrayData(d),
         alpha, src1->getArrayData(d),
         src2->getArrayData(d),
         edge_box);
   }
}

template<class TYPE>
void PatchEdgeDataBasicOps<TYPE>::setRandomValues(
   tbox::Pointer<pdat::EdgeData<TYPE> >& dst,
   const TYPE& width,
   const TYPE& low,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS2(*dst, box);

   const tbox::Dimension& dim(box.getDim());

   for (int d = 0; d < dim.getValue(); d++) {
      const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
      d_array_ops.setRandomValues(dst->getArrayData(d),
         width, low, edge_box);
   }
}

template<class TYPE>
TYPE PatchEdgeDataBasicOps<TYPE>::min(
   const tbox::Pointer<pdat::EdgeData<TYPE> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(!data.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   const tbox::Dimension& dim(box.getDim());
   TYPE minval = tbox::MathUtilities<TYPE>::getMax();

   for (int d = 0; d < dim.getValue(); d++) {
      const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
      minval = tbox::MathUtilities<TYPE>::Min(
            minval, d_array_ops.min(data->getArrayData(d), edge_box));
   }
   return minval;
}

template<class TYPE>
TYPE PatchEdgeDataBasicOps<TYPE>::max(
   const tbox::Pointer<pdat::EdgeData<TYPE> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(!data.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   const tbox::Dimension& dim(box.getDim());

   TYPE maxval = -tbox::MathUtilities<TYPE>::getMax();
   for (int d = 0; d < dim.getValue(); d++) {
      const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
      maxval = tbox::MathUtilities<TYPE>::Max(
            maxval, d_array_ops.max(data->getArrayData(d), edge_box));
   }
   return maxval;
}

}
}
#endif
