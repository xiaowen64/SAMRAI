/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Templated miscellaneous operations for real node-centered data. 
 *
 ************************************************************************/

#ifndef included_math_PatchNodeDataMiscellaneousOpsReal_C
#define included_math_PatchNodeDataMiscellaneousOpsReal_C

#include "SAMRAI/math/PatchNodeDataMiscellaneousOpsReal.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/pdat/NodeGeometry.h"

namespace SAMRAI {
namespace math {

template<class TYPE>
PatchNodeDataMiscellaneousOpsReal<TYPE>::PatchNodeDataMiscellaneousOpsReal()
{
}

template<class TYPE>
PatchNodeDataMiscellaneousOpsReal<TYPE>::~PatchNodeDataMiscellaneousOpsReal()
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
PatchNodeDataMiscellaneousOpsReal<TYPE>::PatchNodeDataMiscellaneousOpsReal(
   const PatchNodeDataMiscellaneousOpsReal<TYPE>& foo)
{
   NULL_USE(foo);
}

template<class TYPE>
void PatchNodeDataMiscellaneousOpsReal<TYPE>::operator = (
   const PatchNodeDataMiscellaneousOpsReal<TYPE>& foo)
{
   NULL_USE(foo);
}

/*
 *************************************************************************
 *                                                                       *
 * Templated miscellaneous opertions for real node-centered data.        *
 *                                                                       *
 *************************************************************************
 */

template<class TYPE>
int PatchNodeDataMiscellaneousOpsReal<TYPE>::computeConstrProdPos(
   const tbox::Pointer<pdat::NodeData<TYPE> >& data1,
   const tbox::Pointer<pdat::NodeData<TYPE> >& data2,
   const hier::Box& box,
   const tbox::Pointer<pdat::NodeData<double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data1.isNull() && !data2.isNull());
#endif
   int retval;
   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   if (cvol.isNull()) {
      retval = d_array_ops.computeConstrProdPos(data1->getArrayData(),
            data2->getArrayData(),
            node_box);
   } else {
      retval = d_array_ops.computeConstrProdPosWithControlVolume(
            data1->getArrayData(),
            data2->getArrayData(),
            cvol->getArrayData(),
            node_box);
   }
   return retval;
}

template<class TYPE>
void PatchNodeDataMiscellaneousOpsReal<TYPE>::compareToScalar(
   tbox::Pointer<pdat::NodeData<TYPE> >& dst,
   const tbox::Pointer<pdat::NodeData<TYPE> >& src,
   const TYPE& alpha,
   const hier::Box& box,
   const tbox::Pointer<pdat::NodeData<double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
#endif
   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   if (cvol.isNull()) {
      d_array_ops.compareToScalar(dst->getArrayData(),
         src->getArrayData(),
         alpha,
         node_box);
   } else {
      d_array_ops.compareToScalarWithControlVolume(dst->getArrayData(),
         src->getArrayData(),
         alpha,
         cvol->getArrayData(),
         node_box);
   }
}

template<class TYPE>
int PatchNodeDataMiscellaneousOpsReal<TYPE>::testReciprocal(
   tbox::Pointer<pdat::NodeData<TYPE> >& dst,
   const tbox::Pointer<pdat::NodeData<TYPE> >& src,
   const hier::Box& box,
   const tbox::Pointer<pdat::NodeData<double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
#endif
   int retval;
   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   if (cvol.isNull()) {
      retval = d_array_ops.testReciprocal(dst->getArrayData(),
            src->getArrayData(),
            node_box);
   } else {
      retval = d_array_ops.testReciprocalWithControlVolume(
            dst->getArrayData(),
            src->getArrayData(),
            cvol->getArrayData(),
            node_box);
   }
   return retval;
}

template<class TYPE>
TYPE PatchNodeDataMiscellaneousOpsReal<TYPE>::maxPointwiseDivide(
   const tbox::Pointer<pdat::NodeData<TYPE> >& numer,
   const tbox::Pointer<pdat::NodeData<TYPE> >& denom,
   const hier::Box& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!numer.isNull() && !denom.isNull());
#endif
   TYPE retval;
   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   retval = d_array_ops.maxPointwiseDivide(numer->getArrayData(),
         denom->getArrayData(),
         node_box);
   return retval;
}

template<class TYPE>
TYPE PatchNodeDataMiscellaneousOpsReal<TYPE>::minPointwiseDivide(
   const tbox::Pointer<pdat::NodeData<TYPE> >& numer,
   const tbox::Pointer<pdat::NodeData<TYPE> >& denom,
   const hier::Box& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!numer.isNull() && !denom.isNull());
#endif
   TYPE retval;
   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   retval = d_array_ops.minPointwiseDivide(numer->getArrayData(),
         denom->getArrayData(),
         node_box);
   return retval;
}

}
}
#endif
