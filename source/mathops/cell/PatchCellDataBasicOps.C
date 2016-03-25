//
// File:	PatchCellDataBasicOps.C
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Basic templated cell-centered patch data operations.
//

#ifndef included_math_PatchCellDataBasicOps_C
#define included_math_PatchCellDataBasicOps_C

#include "PatchCellDataBasicOps.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

namespace SAMRAI {
    namespace math {

template<int DIM, class TYPE>
PatchCellDataBasicOps<DIM,TYPE>::PatchCellDataBasicOps()
{
}

template<int DIM, class TYPE>
PatchCellDataBasicOps<DIM,TYPE>::~PatchCellDataBasicOps()
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
PatchCellDataBasicOps<DIM,TYPE>::PatchCellDataBasicOps(
   const PatchCellDataBasicOps<DIM,TYPE>& foo)
{
   (void) foo;  // not implemented (but needed by some compilers)
}

template<int DIM, class TYPE>
void PatchCellDataBasicOps<DIM,TYPE>::operator=(
   const PatchCellDataBasicOps<DIM,TYPE>& foo)
{
   (void) foo;  // not implemented (but needed by some compilers)
}

/*
*************************************************************************
*                                                                       *
* Generic templated basic operations for cell-centered patch data.      *
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
void PatchCellDataBasicOps<DIM,TYPE>::scale(
   tbox::Pointer< pdat::CellData<DIM,TYPE> >& dst,
   const TYPE& alpha,
   const tbox::Pointer< pdat::CellData<DIM,TYPE> >& src,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!dst.isNull() && !src.isNull());
#endif
   d_array_ops.scale(dst->getArrayData(),
                     alpha, src->getArrayData(),
                     box);
}

template<int DIM, class TYPE>
void PatchCellDataBasicOps<DIM,TYPE>::addScalar(
   tbox::Pointer< pdat::CellData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::CellData<DIM,TYPE> >& src,
   const TYPE& alpha,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!dst.isNull() && !src.isNull());
#endif
   d_array_ops.addScalar(dst->getArrayData(),
                         src->getArrayData(), alpha,
                         box);
}

template<int DIM, class TYPE>
void PatchCellDataBasicOps<DIM,TYPE>::add(
   tbox::Pointer< pdat::CellData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::CellData<DIM,TYPE> >& src1,
   const tbox::Pointer< pdat::CellData<DIM,TYPE> >& src2,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!dst.isNull() && !src1.isNull() && !src2.isNull());
#endif
   d_array_ops.add(dst->getArrayData(),
                   src1->getArrayData(), src2->getArrayData(),
                   box);
}

template<int DIM, class TYPE>
void PatchCellDataBasicOps<DIM,TYPE>::subtract(
   tbox::Pointer< pdat::CellData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::CellData<DIM,TYPE> >& src1,
   const tbox::Pointer< pdat::CellData<DIM,TYPE> >& src2,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!dst.isNull() && !src1.isNull() && !src2.isNull());
#endif
   d_array_ops.subtract(dst->getArrayData(),
                        src1->getArrayData(), src2->getArrayData(),
                        box);
}

template<int DIM, class TYPE>
void PatchCellDataBasicOps<DIM,TYPE>::multiply(
   tbox::Pointer< pdat::CellData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::CellData<DIM,TYPE> >& src1,
   const tbox::Pointer< pdat::CellData<DIM,TYPE> >& src2,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!dst.isNull() && !src1.isNull() && !src2.isNull());
#endif
   d_array_ops.multiply(dst->getArrayData(),
                        src1->getArrayData(), src2->getArrayData(),
                        box);
}

template<int DIM, class TYPE>
void PatchCellDataBasicOps<DIM,TYPE>::divide(
   tbox::Pointer< pdat::CellData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::CellData<DIM,TYPE> >& src1,
   const tbox::Pointer< pdat::CellData<DIM,TYPE> >& src2,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!dst.isNull() && !src1.isNull() && !src2.isNull());
#endif
   d_array_ops.divide(dst->getArrayData(),
                      src1->getArrayData(), src2->getArrayData(),
                      box);
}

template<int DIM, class TYPE>
void PatchCellDataBasicOps<DIM,TYPE>::reciprocal(
   tbox::Pointer< pdat::CellData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::CellData<DIM,TYPE> >& src,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!dst.isNull() && !src.isNull());
#endif
   d_array_ops.reciprocal(dst->getArrayData(),
                          src->getArrayData(),
                          box);
}

template<int DIM, class TYPE>
void PatchCellDataBasicOps<DIM,TYPE>::linearSum(
   tbox::Pointer< pdat::CellData<DIM,TYPE> >& dst,
   const TYPE& alpha,
   const tbox::Pointer< pdat::CellData<DIM,TYPE> >& src1,
   const TYPE& beta,
   const tbox::Pointer< pdat::CellData<DIM,TYPE> >& src2,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!dst.isNull() && !src1.isNull() && !src2.isNull());
#endif
   d_array_ops.linearSum(dst->getArrayData(),
                         alpha, src1->getArrayData(),
                         beta, src2->getArrayData(),
                         box);
}

template<int DIM, class TYPE>
void PatchCellDataBasicOps<DIM,TYPE>::axpy(
   tbox::Pointer< pdat::CellData<DIM,TYPE> >& dst,
   const TYPE& alpha,
   const tbox::Pointer< pdat::CellData<DIM,TYPE> >& src1,
   const tbox::Pointer< pdat::CellData<DIM,TYPE> >& src2,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!dst.isNull() && !src1.isNull() && !src2.isNull());
#endif
   d_array_ops.axpy(dst->getArrayData(),
                    alpha, src1->getArrayData(),
                    src2->getArrayData(),
                    box);
}

template<int DIM, class TYPE>
void PatchCellDataBasicOps<DIM,TYPE>::axmy(
   tbox::Pointer< pdat::CellData<DIM,TYPE> >& dst,
   const TYPE& alpha,
   const tbox::Pointer< pdat::CellData<DIM,TYPE> >& src1,
   const tbox::Pointer< pdat::CellData<DIM,TYPE> >& src2,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!dst.isNull() && !src1.isNull() && !src2.isNull());
#endif
   d_array_ops.axmy(dst->getArrayData(),
                    alpha, src1->getArrayData(),
                    src2->getArrayData(),
                    box);
}

template<int DIM, class TYPE>
TYPE PatchCellDataBasicOps<DIM,TYPE>::min(
   const tbox::Pointer< pdat::CellData<DIM,TYPE> >& data,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!data.isNull());
#endif
   return( d_array_ops.min(data->getArrayData(), box) );
}

template<int DIM, class TYPE>
TYPE PatchCellDataBasicOps<DIM,TYPE>::max(
   const tbox::Pointer< pdat::CellData<DIM,TYPE> >& data,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!data.isNull());
#endif
   return( d_array_ops.max(data->getArrayData(), box) );
}

template<int DIM, class TYPE>
void PatchCellDataBasicOps<DIM,TYPE>::setRandomValues(
   tbox::Pointer< pdat::CellData<DIM,TYPE> >& dst,
   const TYPE& width,
   const TYPE& low,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!dst.isNull());
#endif
   d_array_ops.setRandomValues(dst->getArrayData(),
                               width, low, box);
}

}
}
#endif
