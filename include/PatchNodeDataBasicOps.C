//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-1/source/mathops/node/PatchNodeDataBasicOps.C $
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1846 $
// Modified:	$LastChangedDate: 2008-01-11 09:51:05 -0800 (Fri, 11 Jan 2008) $
// Description:	Basic templated node-centered patch data operations.
//

#ifndef included_math_PatchNodeDataBasicOps_C
#define included_math_PatchNodeDataBasicOps_C

#include "PatchNodeDataBasicOps.h"
#include "NodeGeometry.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#include "tbox/Utilities.h"
#endif

namespace SAMRAI {
    namespace math {

template<int DIM, class TYPE>
PatchNodeDataBasicOps<DIM,TYPE>::PatchNodeDataBasicOps()
{
}

template<int DIM, class TYPE>
PatchNodeDataBasicOps<DIM,TYPE>::~PatchNodeDataBasicOps()
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
PatchNodeDataBasicOps<DIM,TYPE>::PatchNodeDataBasicOps(
   const PatchNodeDataBasicOps<DIM,TYPE>& foo)
{
   NULL_USE(foo);
}

template<int DIM, class TYPE>
void PatchNodeDataBasicOps<DIM,TYPE>::operator=(
   const PatchNodeDataBasicOps<DIM,TYPE>& foo)
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

template<int DIM, class TYPE>
void PatchNodeDataBasicOps<DIM,TYPE>::scale(
   tbox::Pointer< pdat::NodeData<DIM,TYPE> >& dst,
   const TYPE& alpha,
   const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& src,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
#endif
   const hier::Box<DIM> node_box = pdat::NodeGeometry<DIM>::toNodeBox(box);
   d_array_ops.scale(dst->getArrayData(),
                     alpha, src->getArrayData(),
                     node_box);
}

template<int DIM, class TYPE>
void PatchNodeDataBasicOps<DIM,TYPE>::addScalar(
   tbox::Pointer< pdat::NodeData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& src,
   const TYPE& alpha,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
#endif
   const hier::Box<DIM> node_box = pdat::NodeGeometry<DIM>::toNodeBox(box);
   d_array_ops.addScalar(dst->getArrayData(),
                         src->getArrayData(), alpha,
                         node_box);
}

template<int DIM, class TYPE>
void PatchNodeDataBasicOps<DIM,TYPE>::add(
   tbox::Pointer< pdat::NodeData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& src1,
   const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& src2,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
#endif
   const hier::Box<DIM> node_box = pdat::NodeGeometry<DIM>::toNodeBox(box);
   d_array_ops.add(dst->getArrayData(),
                   src1->getArrayData(), src2->getArrayData(),
                   node_box);
}

template<int DIM, class TYPE>
void PatchNodeDataBasicOps<DIM,TYPE>::subtract(
   tbox::Pointer< pdat::NodeData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& src1,
   const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& src2,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
#endif
   const hier::Box<DIM> node_box = pdat::NodeGeometry<DIM>::toNodeBox(box);
   d_array_ops.subtract(dst->getArrayData(),
                        src1->getArrayData(), src2->getArrayData(),
                        node_box);
}

template<int DIM, class TYPE>
void PatchNodeDataBasicOps<DIM,TYPE>::multiply(
   tbox::Pointer< pdat::NodeData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& src1,
   const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& src2,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
#endif
   const hier::Box<DIM> node_box = pdat::NodeGeometry<DIM>::toNodeBox(box);
   d_array_ops.multiply(dst->getArrayData(),
                        src1->getArrayData(), src2->getArrayData(),
                        node_box);
}

template<int DIM, class TYPE>
void PatchNodeDataBasicOps<DIM,TYPE>::divide(
   tbox::Pointer< pdat::NodeData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& src1,
   const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& src2,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
#endif
   const hier::Box<DIM> node_box = pdat::NodeGeometry<DIM>::toNodeBox(box);
   d_array_ops.divide(dst->getArrayData(),
                      src1->getArrayData(), src2->getArrayData(),
                      node_box);
}

template<int DIM, class TYPE>
void PatchNodeDataBasicOps<DIM,TYPE>::reciprocal(
   tbox::Pointer< pdat::NodeData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& src,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
#endif
   const hier::Box<DIM> node_box = pdat::NodeGeometry<DIM>::toNodeBox(box);
   d_array_ops.reciprocal(dst->getArrayData(),
                          src->getArrayData(),
                          node_box);
}

template<int DIM, class TYPE>
void PatchNodeDataBasicOps<DIM,TYPE>::linearSum(
   tbox::Pointer< pdat::NodeData<DIM,TYPE> >& dst,
   const TYPE& alpha,
   const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& src1,
   const TYPE& beta,
   const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& src2,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
#endif
   const hier::Box<DIM> node_box = pdat::NodeGeometry<DIM>::toNodeBox(box);
   d_array_ops.linearSum(dst->getArrayData(),
                         alpha, src1->getArrayData(),
                         beta, src2->getArrayData(),
                         node_box);
}

template<int DIM, class TYPE>
void PatchNodeDataBasicOps<DIM,TYPE>::axpy(
   tbox::Pointer< pdat::NodeData<DIM,TYPE> >& dst,
   const TYPE& alpha,
   const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& src1,
   const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& src2,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
#endif
   const hier::Box<DIM> node_box = pdat::NodeGeometry<DIM>::toNodeBox(box);
   d_array_ops.axpy(dst->getArrayData(),
                    alpha, src1->getArrayData(),
                    src2->getArrayData(),
                    node_box);
}

template<int DIM, class TYPE>
void PatchNodeDataBasicOps<DIM,TYPE>::axmy(
   tbox::Pointer< pdat::NodeData<DIM,TYPE> >& dst,
   const TYPE& alpha,
   const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& src1,
   const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& src2,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
#endif
   const hier::Box<DIM> node_box = pdat::NodeGeometry<DIM>::toNodeBox(box);
   d_array_ops.axmy(dst->getArrayData(),
                    alpha, src1->getArrayData(),
                    src2->getArrayData(),
                    node_box);
}

template<int DIM, class TYPE>
TYPE PatchNodeDataBasicOps<DIM,TYPE>::min(
   const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& data,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
#endif
   const hier::Box<DIM> node_box = pdat::NodeGeometry<DIM>::toNodeBox(box);
   return( d_array_ops.min(data->getArrayData(), node_box) );
}

template<int DIM, class TYPE>
TYPE PatchNodeDataBasicOps<DIM,TYPE>::max(
   const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& data,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
#endif
   const hier::Box<DIM> node_box = pdat::NodeGeometry<DIM>::toNodeBox(box);
   return( d_array_ops.max(data->getArrayData(), node_box) );
}

template<int DIM, class TYPE>
void PatchNodeDataBasicOps<DIM,TYPE>::setRandomValues(
   tbox::Pointer< pdat::NodeData<DIM,TYPE> >& dst,
   const TYPE& width,
   const TYPE& low,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull());
#endif
   const hier::Box<DIM> node_box = pdat::NodeGeometry<DIM>::toNodeBox(box);
   d_array_ops.setRandomValues(dst->getArrayData(),
                               width, low, node_box);
}

}
}
#endif
