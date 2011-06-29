/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Basic templated node-centered patch data operations. 
 *
 ************************************************************************/

#ifndef included_math_PatchNodeDataBasicOps_C
#define included_math_PatchNodeDataBasicOps_C

#include "SAMRAI/math/PatchNodeDataBasicOps.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/pdat/NodeGeometry.h"

namespace SAMRAI {
namespace math {

template<class TYPE>
PatchNodeDataBasicOps<TYPE>::PatchNodeDataBasicOps()
{
}

template<class TYPE>
PatchNodeDataBasicOps<TYPE>::~PatchNodeDataBasicOps()
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
PatchNodeDataBasicOps<TYPE>::PatchNodeDataBasicOps(
   const PatchNodeDataBasicOps<TYPE>& foo)
{
   NULL_USE(foo);
}

template<class TYPE>
void PatchNodeDataBasicOps<TYPE>::operator = (
   const PatchNodeDataBasicOps<TYPE>& foo)
{
   NULL_USE(foo);
}

/*
 *************************************************************************
 *                                                                       *
 * Generic basic templated operations for node-centered patch data.      *
 *                                                                       *
 *************************************************************************
 */

template<class TYPE>
void PatchNodeDataBasicOps<TYPE>::scale(
   tbox::Pointer<pdat::NodeData<TYPE> >& dst,
   const TYPE& alpha,
   const tbox::Pointer<pdat::NodeData<TYPE> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst, *src, box);

   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   d_array_ops.scale(dst->getArrayData(),
      alpha, src->getArrayData(),
      node_box);
}

template<class TYPE>
void PatchNodeDataBasicOps<TYPE>::addScalar(
   tbox::Pointer<pdat::NodeData<TYPE> >& dst,
   const tbox::Pointer<pdat::NodeData<TYPE> >& src,
   const TYPE& alpha,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst, *src, box);

   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   d_array_ops.addScalar(dst->getArrayData(),
      src->getArrayData(), alpha,
      node_box);
}

template<class TYPE>
void PatchNodeDataBasicOps<TYPE>::add(
   tbox::Pointer<pdat::NodeData<TYPE> >& dst,
   const tbox::Pointer<pdat::NodeData<TYPE> >& src1,
   const tbox::Pointer<pdat::NodeData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst, *src1, *src2, box);

   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   d_array_ops.add(dst->getArrayData(),
      src1->getArrayData(), src2->getArrayData(),
      node_box);
}

template<class TYPE>
void PatchNodeDataBasicOps<TYPE>::subtract(
   tbox::Pointer<pdat::NodeData<TYPE> >& dst,
   const tbox::Pointer<pdat::NodeData<TYPE> >& src1,
   const tbox::Pointer<pdat::NodeData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst, *src1, *src2, box);

   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   d_array_ops.subtract(dst->getArrayData(),
      src1->getArrayData(), src2->getArrayData(),
      node_box);
}

template<class TYPE>
void PatchNodeDataBasicOps<TYPE>::multiply(
   tbox::Pointer<pdat::NodeData<TYPE> >& dst,
   const tbox::Pointer<pdat::NodeData<TYPE> >& src1,
   const tbox::Pointer<pdat::NodeData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst, *src1, *src2, box);

   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   d_array_ops.multiply(dst->getArrayData(),
      src1->getArrayData(), src2->getArrayData(),
      node_box);
}

template<class TYPE>
void PatchNodeDataBasicOps<TYPE>::divide(
   tbox::Pointer<pdat::NodeData<TYPE> >& dst,
   const tbox::Pointer<pdat::NodeData<TYPE> >& src1,
   const tbox::Pointer<pdat::NodeData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst, *src1, *src2, box);

   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   d_array_ops.divide(dst->getArrayData(),
      src1->getArrayData(), src2->getArrayData(),
      node_box);
}

template<class TYPE>
void PatchNodeDataBasicOps<TYPE>::reciprocal(
   tbox::Pointer<pdat::NodeData<TYPE> >& dst,
   const tbox::Pointer<pdat::NodeData<TYPE> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst, *src, box);

   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   d_array_ops.reciprocal(dst->getArrayData(),
      src->getArrayData(),
      node_box);
}

template<class TYPE>
void PatchNodeDataBasicOps<TYPE>::linearSum(
   tbox::Pointer<pdat::NodeData<TYPE> >& dst,
   const TYPE& alpha,
   const tbox::Pointer<pdat::NodeData<TYPE> >& src1,
   const TYPE& beta,
   const tbox::Pointer<pdat::NodeData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst, *src1, *src2, box);

   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   d_array_ops.linearSum(dst->getArrayData(),
      alpha, src1->getArrayData(),
      beta, src2->getArrayData(),
      node_box);
}

template<class TYPE>
void PatchNodeDataBasicOps<TYPE>::axpy(
   tbox::Pointer<pdat::NodeData<TYPE> >& dst,
   const TYPE& alpha,
   const tbox::Pointer<pdat::NodeData<TYPE> >& src1,
   const tbox::Pointer<pdat::NodeData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst, *src1, *src2, box);

   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   d_array_ops.axpy(dst->getArrayData(),
      alpha, src1->getArrayData(),
      src2->getArrayData(),
      node_box);
}

template<class TYPE>
void PatchNodeDataBasicOps<TYPE>::axmy(
   tbox::Pointer<pdat::NodeData<TYPE> >& dst,
   const TYPE& alpha,
   const tbox::Pointer<pdat::NodeData<TYPE> >& src1,
   const tbox::Pointer<pdat::NodeData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst, *src1, *src2, box);

   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   d_array_ops.axmy(dst->getArrayData(),
      alpha, src1->getArrayData(),
      src2->getArrayData(),
      node_box);
}

template<class TYPE>
TYPE PatchNodeDataBasicOps<TYPE>::min(
   const tbox::Pointer<pdat::NodeData<TYPE> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(!data.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   return d_array_ops.min(data->getArrayData(), node_box);
}

template<class TYPE>
TYPE PatchNodeDataBasicOps<TYPE>::max(
   const tbox::Pointer<pdat::NodeData<TYPE> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(!data.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   return d_array_ops.max(data->getArrayData(), node_box);
}

template<class TYPE>
void PatchNodeDataBasicOps<TYPE>::setRandomValues(
   tbox::Pointer<pdat::NodeData<TYPE> >& dst,
   const TYPE& width,
   const TYPE& low,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS2(*dst, box);

   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   d_array_ops.setRandomValues(dst->getArrayData(),
      width, low, node_box);
}

}
}
#endif
