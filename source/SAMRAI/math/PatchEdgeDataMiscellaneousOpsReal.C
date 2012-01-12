/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Templated miscellaneous operations for real edge-centered data.
 *
 ************************************************************************/

#ifndef included_math_PatchEdgeDataMiscellaneousOpsReal_C
#define included_math_PatchEdgeDataMiscellaneousOpsReal_C

#include "SAMRAI/math/PatchEdgeDataMiscellaneousOpsReal.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/pdat/EdgeGeometry.h"

namespace SAMRAI {
namespace math {

template<class TYPE>
PatchEdgeDataMiscellaneousOpsReal<TYPE>::PatchEdgeDataMiscellaneousOpsReal()
{
}

template<class TYPE>
PatchEdgeDataMiscellaneousOpsReal<TYPE>::~PatchEdgeDataMiscellaneousOpsReal()
{
}

/*
 *************************************************************************
 *
 * The const constructor and assignment operator are not actually used
 * but are defined here for compilers that require an implementation for
 * every declaration.
 *
 *************************************************************************
 */

template<class TYPE>
PatchEdgeDataMiscellaneousOpsReal<TYPE>::PatchEdgeDataMiscellaneousOpsReal(
   const PatchEdgeDataMiscellaneousOpsReal<TYPE>& foo)
{
   NULL_USE(foo);
}

template<class TYPE>
void PatchEdgeDataMiscellaneousOpsReal<TYPE>::operator = (
   const PatchEdgeDataMiscellaneousOpsReal<TYPE>& foo)
{
   NULL_USE(foo);
}

/*
 *************************************************************************
 *
 * Templated miscellaneous opertions for real edge-centered data.
 *
 *************************************************************************
 */

template<class TYPE>
int PatchEdgeDataMiscellaneousOpsReal<TYPE>::computeConstrProdPos(
   const tbox::Pointer<pdat::EdgeData<TYPE> >& data1,
   const tbox::Pointer<pdat::EdgeData<TYPE> >& data2,
   const hier::Box& box,
   const tbox::Pointer<pdat::EdgeData<double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(data1 && data2);
#endif
   const tbox::Dimension& dim(data1->getDim());

   int retval = 1;
   if (!cvol) {
      for (int d = 0; d < dim.getValue(); d++) {
         const hier::Box edge_box =
            pdat::EdgeGeometry::toEdgeBox(box, d);
         retval = tbox::MathUtilities<int>::Min(retval,
               d_array_ops.computeConstrProdPos(
                  data1->getArrayData(d),
                  data2->getArrayData(d),
                  edge_box));
      }
   } else {
      for (int d = 0; d < dim.getValue(); d++) {
         const hier::Box edge_box =
            pdat::EdgeGeometry::toEdgeBox(box, d);
         retval = tbox::MathUtilities<int>::Min(retval,
               d_array_ops.computeConstrProdPosWithControlVolume(
                  data1->getArrayData(d),
                  data2->getArrayData(d),
                  cvol->getArrayData(d),
                  edge_box));
      }
   }
   return retval;
}

template<class TYPE>
void PatchEdgeDataMiscellaneousOpsReal<TYPE>::compareToScalar(
   tbox::Pointer<pdat::EdgeData<TYPE> >& dst,
   const tbox::Pointer<pdat::EdgeData<TYPE> >& src,
   const TYPE& alpha,
   const hier::Box& box,
   const tbox::Pointer<pdat::EdgeData<double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(dst && src);
#endif
   const tbox::Dimension& dim(dst->getDim());

   if (!cvol) {
      for (int d = 0; d < dim.getValue(); d++) {
         const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
         d_array_ops.compareToScalar(dst->getArrayData(d),
            src->getArrayData(d),
            alpha,
            edge_box);
      }
   } else {
      for (int d = 0; d < dim.getValue(); d++) {
         const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
         d_array_ops.compareToScalarWithControlVolume(dst->getArrayData(d),
            src->getArrayData(d),
            alpha,
            cvol->getArrayData(d),
            edge_box);
      }
   }
}

template<class TYPE>
int PatchEdgeDataMiscellaneousOpsReal<TYPE>::testReciprocal(
   tbox::Pointer<pdat::EdgeData<TYPE> >& dst,
   const tbox::Pointer<pdat::EdgeData<TYPE> >& src,
   const hier::Box& box,
   const tbox::Pointer<pdat::EdgeData<double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(dst && src);
#endif
   const tbox::Dimension& dim(dst->getDim());

   int retval = 1;
   if (!cvol) {
      for (int d = 0; d < dim.getValue(); d++) {
         const hier::Box edge_box =
            pdat::EdgeGeometry::toEdgeBox(box, d);
         retval = tbox::MathUtilities<int>::Min(retval,
               d_array_ops.testReciprocal(
                  dst->getArrayData(d),
                  src->getArrayData(d),
                  edge_box));
      }
   } else {
      for (int d = 0; d < dim.getValue(); d++) {
         const hier::Box edge_box =
            pdat::EdgeGeometry::toEdgeBox(box, d);
         retval = tbox::MathUtilities<int>::Min(retval,
               d_array_ops.testReciprocalWithControlVolume(
                  dst->getArrayData(d),
                  src->getArrayData(d),
                  cvol->getArrayData(d),
                  edge_box));
      }
   }
   return retval;
}

template<class TYPE>
TYPE PatchEdgeDataMiscellaneousOpsReal<TYPE>::maxPointwiseDivide(
   const tbox::Pointer<pdat::EdgeData<TYPE> >& numer,
   const tbox::Pointer<pdat::EdgeData<TYPE> >& denom,
   const hier::Box& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(numer && denom);
#endif
   const tbox::Dimension& dim(numer->getDim());

   TYPE retval = 0.0;
   for (int d = 0; d < dim.getValue(); d++) {
      const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
      TYPE dirval = d_array_ops.maxPointwiseDivide(numer->getArrayData(d),
            denom->getArrayData(d),
            edge_box);
      retval = tbox::MathUtilities<TYPE>::Max(retval, dirval);
   }
   return retval;
}

template<class TYPE>
TYPE PatchEdgeDataMiscellaneousOpsReal<TYPE>::minPointwiseDivide(
   const tbox::Pointer<pdat::EdgeData<TYPE> >& numer,
   const tbox::Pointer<pdat::EdgeData<TYPE> >& denom,
   const hier::Box& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(numer && denom);
#endif
   const tbox::Dimension& dim(numer->getDim());

   TYPE retval = 0.0;
   for (int d = 0; d < dim.getValue(); d++) {
      const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
      TYPE dirval = d_array_ops.minPointwiseDivide(numer->getArrayData(d),
            denom->getArrayData(d),
            edge_box);
      retval = tbox::MathUtilities<TYPE>::Min(retval, dirval);
   }
   return retval;
}

}
}
#endif
