/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Basic templated cell-centered patch data operations. 
 *
 ************************************************************************/

#ifndef included_math_PatchCellDataBasicOps_C
#define included_math_PatchCellDataBasicOps_C

#include "SAMRAI/math/PatchCellDataBasicOps.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace math {

template<class TYPE>
PatchCellDataBasicOps<TYPE>::PatchCellDataBasicOps()
{
}

template<class TYPE>
PatchCellDataBasicOps<TYPE>::~PatchCellDataBasicOps()
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
PatchCellDataBasicOps<TYPE>::PatchCellDataBasicOps(
   const PatchCellDataBasicOps<TYPE>& foo)
{
   NULL_USE(foo);
}

template<class TYPE>
void PatchCellDataBasicOps<TYPE>::operator = (
   const PatchCellDataBasicOps<TYPE>& foo)
{
   NULL_USE(foo);
}

/*
 *************************************************************************
 *                                                                       *
 * Generic templated basic operations for cell-centered patch data.      *
 *                                                                       *
 *************************************************************************
 */

template<class TYPE>
void PatchCellDataBasicOps<TYPE>::scale(
   tbox::Pointer<pdat::CellData<TYPE> >& dst,
   const TYPE& alpha,
   const tbox::Pointer<pdat::CellData<TYPE> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst, *src, box);

   d_array_ops.scale(dst->getArrayData(),
      alpha, src->getArrayData(),
      box);
}

template<class TYPE>
void PatchCellDataBasicOps<TYPE>::addScalar(
   tbox::Pointer<pdat::CellData<TYPE> >& dst,
   const tbox::Pointer<pdat::CellData<TYPE> >& src,
   const TYPE& alpha,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst, *src, box);

   d_array_ops.addScalar(dst->getArrayData(),
      src->getArrayData(), alpha,
      box);
}

template<class TYPE>
void PatchCellDataBasicOps<TYPE>::add(
   tbox::Pointer<pdat::CellData<TYPE> >& dst,
   const tbox::Pointer<pdat::CellData<TYPE> >& src1,
   const tbox::Pointer<pdat::CellData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst, *src1, *src2, box);

   d_array_ops.add(dst->getArrayData(),
      src1->getArrayData(), src2->getArrayData(),
      box);
}

template<class TYPE>
void PatchCellDataBasicOps<TYPE>::subtract(
   tbox::Pointer<pdat::CellData<TYPE> >& dst,
   const tbox::Pointer<pdat::CellData<TYPE> >& src1,
   const tbox::Pointer<pdat::CellData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst, *src1, *src2, box);

   d_array_ops.subtract(dst->getArrayData(),
      src1->getArrayData(), src2->getArrayData(),
      box);
}

template<class TYPE>
void PatchCellDataBasicOps<TYPE>::multiply(
   tbox::Pointer<pdat::CellData<TYPE> >& dst,
   const tbox::Pointer<pdat::CellData<TYPE> >& src1,
   const tbox::Pointer<pdat::CellData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst, *src1, *src2, box);

   d_array_ops.multiply(dst->getArrayData(),
      src1->getArrayData(), src2->getArrayData(),
      box);
}

template<class TYPE>
void PatchCellDataBasicOps<TYPE>::divide(
   tbox::Pointer<pdat::CellData<TYPE> >& dst,
   const tbox::Pointer<pdat::CellData<TYPE> >& src1,
   const tbox::Pointer<pdat::CellData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst, *src1, *src2, box);

   d_array_ops.divide(dst->getArrayData(),
      src1->getArrayData(), src2->getArrayData(),
      box);
}

template<class TYPE>
void PatchCellDataBasicOps<TYPE>::reciprocal(
   tbox::Pointer<pdat::CellData<TYPE> >& dst,
   const tbox::Pointer<pdat::CellData<TYPE> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst, *src, box);

   d_array_ops.reciprocal(dst->getArrayData(),
      src->getArrayData(),
      box);
}

template<class TYPE>
void PatchCellDataBasicOps<TYPE>::linearSum(
   tbox::Pointer<pdat::CellData<TYPE> >& dst,
   const TYPE& alpha,
   const tbox::Pointer<pdat::CellData<TYPE> >& src1,
   const TYPE& beta,
   const tbox::Pointer<pdat::CellData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst, *src1, *src2, box);

   d_array_ops.linearSum(dst->getArrayData(),
      alpha, src1->getArrayData(),
      beta, src2->getArrayData(),
      box);
}

template<class TYPE>
void PatchCellDataBasicOps<TYPE>::axpy(
   tbox::Pointer<pdat::CellData<TYPE> >& dst,
   const TYPE& alpha,
   const tbox::Pointer<pdat::CellData<TYPE> >& src1,
   const tbox::Pointer<pdat::CellData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst, *src1, *src2, box);

   d_array_ops.axpy(dst->getArrayData(),
      alpha, src1->getArrayData(),
      src2->getArrayData(),
      box);
}

template<class TYPE>
void PatchCellDataBasicOps<TYPE>::axmy(
   tbox::Pointer<pdat::CellData<TYPE> >& dst,
   const TYPE& alpha,
   const tbox::Pointer<pdat::CellData<TYPE> >& src1,
   const tbox::Pointer<pdat::CellData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src1.isNull() && !src2.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst, *src1, *src2, box);

   d_array_ops.axmy(dst->getArrayData(),
      alpha, src1->getArrayData(),
      src2->getArrayData(),
      box);
}

template<class TYPE>
TYPE PatchCellDataBasicOps<TYPE>::min(
   const tbox::Pointer<pdat::CellData<TYPE> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(!data.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   return d_array_ops.min(data->getArrayData(), box);
}

template<class TYPE>
TYPE PatchCellDataBasicOps<TYPE>::max(
   const tbox::Pointer<pdat::CellData<TYPE> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(!data.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   return d_array_ops.max(data->getArrayData(), box);
}

template<class TYPE>
void PatchCellDataBasicOps<TYPE>::setRandomValues(
   tbox::Pointer<pdat::CellData<TYPE> >& dst,
   const TYPE& width,
   const TYPE& low,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS2(*dst, box);

   d_array_ops.setRandomValues(dst->getArrayData(),
      width, low, box);
}

}
}
#endif
