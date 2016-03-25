//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-1/source/mathops/edge/PatchEdgeDataBasicOps.C $
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1846 $
// Modified:	$LastChangedDate: 2008-01-11 09:51:05 -0800 (Fri, 11 Jan 2008) $
// Description:	Basic templated edge-centered patch data operations.
//

#ifndef included_math_PatchEdgeDataBasicOps_C
#define included_math_PatchEdgeDataBasicOps_C

#include "tbox/MathUtilities.h"
#include "PatchEdgeDataBasicOps.h"
#include "EdgeGeometry.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#include "tbox/Utilities.h"
#endif

namespace SAMRAI {
    namespace math {

template<int DIM, class TYPE>
PatchEdgeDataBasicOps<DIM,TYPE>::PatchEdgeDataBasicOps()
{
}

template<int DIM, class TYPE>
PatchEdgeDataBasicOps<DIM,TYPE>::~PatchEdgeDataBasicOps()
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
PatchEdgeDataBasicOps<DIM,TYPE>::PatchEdgeDataBasicOps(
   const PatchEdgeDataBasicOps<DIM,TYPE>& foo)
{
   NULL_USE(foo); 
}

template<int DIM, class TYPE>
void PatchEdgeDataBasicOps<DIM,TYPE>::operator=(
   const PatchEdgeDataBasicOps<DIM,TYPE>& foo)
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

template<int DIM, class TYPE>
void PatchEdgeDataBasicOps<DIM,TYPE>::scale(
   tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& dst,
   const TYPE& alpha,
   const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& src,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
#endif
   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> edge_box = pdat::EdgeGeometry<DIM>::toEdgeBox(box, d);
      d_array_ops.scale(dst->getArrayData(d),
                        alpha, src->getArrayData(d),
                        edge_box);
   }
}

template<int DIM, class TYPE>
void PatchEdgeDataBasicOps<DIM,TYPE>::addScalar(
   tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& src,
   const TYPE& alpha,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
#endif
   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> edge_box = pdat::EdgeGeometry<DIM>::toEdgeBox(box, d);
      d_array_ops.addScalar(dst->getArrayData(d),
                            src->getArrayData(d), alpha,
                            edge_box);
   }
}

template<int DIM, class TYPE>
void PatchEdgeDataBasicOps<DIM,TYPE>::add(
   tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& src1,
   const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& src2,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
#endif
   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> edge_box = pdat::EdgeGeometry<DIM>::toEdgeBox(box, d);
      d_array_ops.add(dst->getArrayData(d),
                      src1->getArrayData(d), src2->getArrayData(d),
                      edge_box);
   }
}

template<int DIM, class TYPE>
void PatchEdgeDataBasicOps<DIM,TYPE>::subtract(
   tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& src1,
   const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& src2,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
#endif
   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> edge_box = pdat::EdgeGeometry<DIM>::toEdgeBox(box, d);
      d_array_ops.subtract(dst->getArrayData(d),
                           src1->getArrayData(d), src2->getArrayData(d),
                           edge_box);
   }
}

template<int DIM, class TYPE>
void PatchEdgeDataBasicOps<DIM,TYPE>::multiply(
   tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& src1,
   const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& src2,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
#endif
   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> edge_box = pdat::EdgeGeometry<DIM>::toEdgeBox(box, d);
      d_array_ops.multiply(dst->getArrayData(d),
                           src1->getArrayData(d), src2->getArrayData(d),
                           edge_box);
   }
}

template<int DIM, class TYPE>
void PatchEdgeDataBasicOps<DIM,TYPE>::divide(
   tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& src1,
   const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& src2,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
#endif
   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> edge_box = pdat::EdgeGeometry<DIM>::toEdgeBox(box, d);
      d_array_ops.divide(dst->getArrayData(d),
                         src1->getArrayData(d), src2->getArrayData(d),
                         edge_box);
   }
}

template<int DIM, class TYPE>
void PatchEdgeDataBasicOps<DIM,TYPE>::reciprocal(
   tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& src,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
#endif
   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> edge_box = pdat::EdgeGeometry<DIM>::toEdgeBox(box, d);
      d_array_ops.reciprocal(dst->getArrayData(d),
                             src->getArrayData(d),
                             edge_box);
   }
}

template<int DIM, class TYPE>
void PatchEdgeDataBasicOps<DIM,TYPE>::linearSum(
   tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& dst,
   const TYPE& alpha,
   const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& src1,
   const TYPE& beta,
   const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& src2,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
#endif
   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> edge_box = pdat::EdgeGeometry<DIM>::toEdgeBox(box, d);
      d_array_ops.linearSum(dst->getArrayData(d),
                            alpha, src1->getArrayData(d),
                            beta, src2->getArrayData(d),
                            edge_box);
   }
}

template<int DIM, class TYPE>
void PatchEdgeDataBasicOps<DIM,TYPE>::axpy(
   tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& dst,
   const TYPE& alpha,
   const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& src1,
   const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& src2,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
#endif
   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> edge_box = pdat::EdgeGeometry<DIM>::toEdgeBox(box, d);
      d_array_ops.axpy(dst->getArrayData(d),
                       alpha, src1->getArrayData(d),
                       src2->getArrayData(d),
                       edge_box);
   }
}

template<int DIM, class TYPE>
void PatchEdgeDataBasicOps<DIM,TYPE>::axmy(
   tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& dst,
   const TYPE& alpha,
   const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& src1,
   const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& src2,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
#endif
   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> edge_box = pdat::EdgeGeometry<DIM>::toEdgeBox(box, d);
      d_array_ops.axmy(dst->getArrayData(d),
                       alpha, src1->getArrayData(d),
                       src2->getArrayData(d),
                       edge_box);
   }
}

template<int DIM, class TYPE>
void PatchEdgeDataBasicOps<DIM,TYPE>::setRandomValues(
   tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& dst,
   const TYPE& width,
   const TYPE& low,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull());
#endif
   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> edge_box = pdat::EdgeGeometry<DIM>::toEdgeBox(box, d);
      d_array_ops.setRandomValues(dst->getArrayData(d),
                                  width, low, edge_box);
   }
}

template<int DIM, class TYPE>
TYPE PatchEdgeDataBasicOps<DIM,TYPE>::min(
   const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& data,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
#endif
   TYPE minval = tbox::MathUtilities<TYPE>::getMax();
   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> edge_box = pdat::EdgeGeometry<DIM>::toEdgeBox(box, d);
      minval = tbox::MathUtilities<TYPE>::Min(
                  minval, d_array_ops.min(data->getArrayData(d), edge_box) );
   }
   return(minval);
}

template<int DIM, class TYPE>
TYPE PatchEdgeDataBasicOps<DIM,TYPE>::max(
   const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& data,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
#endif
   TYPE maxval = -tbox::MathUtilities<TYPE>::getMax();
   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> edge_box = pdat::EdgeGeometry<DIM>::toEdgeBox(box, d);
      maxval = tbox::MathUtilities<TYPE>::Max(
                  maxval, d_array_ops.max(data->getArrayData(d), edge_box) );
   }
   return(maxval);
}

}
}
#endif
