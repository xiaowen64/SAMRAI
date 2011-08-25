/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Templated miscellaneous operations for real face-centered data.
 *
 ************************************************************************/

#ifndef included_math_PatchFaceDataMiscellaneousOpsReal_C
#define included_math_PatchFaceDataMiscellaneousOpsReal_C

#include "SAMRAI/math/PatchFaceDataMiscellaneousOpsReal.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/pdat/FaceGeometry.h"

namespace SAMRAI {
namespace math {

template<class TYPE>
PatchFaceDataMiscellaneousOpsReal<TYPE>::PatchFaceDataMiscellaneousOpsReal()
{
}

template<class TYPE>
PatchFaceDataMiscellaneousOpsReal<TYPE>::~PatchFaceDataMiscellaneousOpsReal()
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
PatchFaceDataMiscellaneousOpsReal<TYPE>::PatchFaceDataMiscellaneousOpsReal(
   const PatchFaceDataMiscellaneousOpsReal<TYPE>& foo)
{
   NULL_USE(foo);
}

template<class TYPE>
void PatchFaceDataMiscellaneousOpsReal<TYPE>::operator = (
   const PatchFaceDataMiscellaneousOpsReal<TYPE>& foo)
{
   NULL_USE(foo);
}

/*
 *************************************************************************
 *                                                                       *
 * Templated miscellaneous opertions for real face-centered data.        *
 *                                                                       *
 *************************************************************************
 */

template<class TYPE>
int PatchFaceDataMiscellaneousOpsReal<TYPE>::computeConstrProdPos(
   const tbox::Pointer<pdat::FaceData<TYPE> >& data1,
   const tbox::Pointer<pdat::FaceData<TYPE> >& data2,
   const hier::Box& box,
   const tbox::Pointer<pdat::FaceData<double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data1.isNull() && !data2.isNull());
#endif
   const tbox::Dimension& dim(data1->getDim());

   int retval = 1;
   if (cvol.isNull()) {
      for (int d = 0; d < dim.getValue(); d++) {
         const hier::Box face_box =
            pdat::FaceGeometry::toFaceBox(box, d);
         retval = tbox::MathUtilities<int>::Min(retval,
               d_array_ops.computeConstrProdPos(
                  data1->getArrayData(d),
                  data2->getArrayData(d),
                  face_box));
      }
   } else {
      for (int d = 0; d < dim.getValue(); d++) {
         const hier::Box face_box =
            pdat::FaceGeometry::toFaceBox(box, d);
         retval = tbox::MathUtilities<int>::Min(retval,
               d_array_ops.computeConstrProdPosWithControlVolume(
                  data1->getArrayData(d),
                  data2->getArrayData(d),
                  cvol->getArrayData(d),
                  face_box));
      }
   }
   return retval;
}

template<class TYPE>
void PatchFaceDataMiscellaneousOpsReal<TYPE>::compareToScalar(
   tbox::Pointer<pdat::FaceData<TYPE> >& dst,
   const tbox::Pointer<pdat::FaceData<TYPE> >& src,
   const TYPE& alpha,
   const hier::Box& box,
   const tbox::Pointer<pdat::FaceData<double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
#endif
   const tbox::Dimension& dim(dst->getDim());

   if (cvol.isNull()) {
      for (int d = 0; d < dim.getValue(); d++) {
         const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
         d_array_ops.compareToScalar(dst->getArrayData(d),
            src->getArrayData(d),
            alpha,
            face_box);
      }
   } else {
      for (int d = 0; d < dim.getValue(); d++) {
         const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
         d_array_ops.compareToScalarWithControlVolume(dst->getArrayData(d),
            src->getArrayData(d),
            alpha,
            cvol->getArrayData(d),
            face_box);
      }
   }
}

template<class TYPE>
int PatchFaceDataMiscellaneousOpsReal<TYPE>::testReciprocal(
   tbox::Pointer<pdat::FaceData<TYPE> >& dst,
   const tbox::Pointer<pdat::FaceData<TYPE> >& src,
   const hier::Box& box,
   const tbox::Pointer<pdat::FaceData<double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
#endif
   const tbox::Dimension& dim(dst->getDim());

   int retval = 1;
   if (cvol.isNull()) {
      for (int d = 0; d < dim.getValue(); d++) {
         const hier::Box face_box =
            pdat::FaceGeometry::toFaceBox(box, d);
         retval = tbox::MathUtilities<int>::Min(retval,
               d_array_ops.testReciprocal(
                  dst->getArrayData(d),
                  src->getArrayData(d),
                  face_box));
      }
   } else {
      for (int d = 0; d < dim.getValue(); d++) {
         const hier::Box face_box =
            pdat::FaceGeometry::toFaceBox(box, d);
         retval = tbox::MathUtilities<int>::Min(retval,
               d_array_ops.testReciprocalWithControlVolume(
                  dst->getArrayData(d),
                  src->getArrayData(d),
                  cvol->getArrayData(d),
                  face_box));
      }
   }
   return retval;
}

template<class TYPE>
TYPE PatchFaceDataMiscellaneousOpsReal<TYPE>::maxPointwiseDivide(
   const tbox::Pointer<pdat::FaceData<TYPE> >& numer,
   const tbox::Pointer<pdat::FaceData<TYPE> >& denom,
   const hier::Box& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!numer.isNull() && !denom.isNull());
#endif
   const tbox::Dimension& dim(numer->getDim());

   TYPE retval = 0.0;
   for (int d = 0; d < dim.getValue(); d++) {
      const hier::Box face_box =
         pdat::FaceGeometry::toFaceBox(box, d);
      TYPE dirval = d_array_ops.maxPointwiseDivide(numer->getArrayData(d),
            denom->getArrayData(d),
            face_box);
      retval = tbox::MathUtilities<TYPE>::Max(retval, dirval);
   }
   return retval;
}

template<class TYPE>
TYPE PatchFaceDataMiscellaneousOpsReal<TYPE>::minPointwiseDivide(
   const tbox::Pointer<pdat::FaceData<TYPE> >& numer,
   const tbox::Pointer<pdat::FaceData<TYPE> >& denom,
   const hier::Box& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!numer.isNull() && !denom.isNull());
#endif
   const tbox::Dimension& dim(numer->getDim());

   TYPE retval = 0.0;
   for (int d = 0; d < dim.getValue(); d++) {
      const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
      TYPE dirval = d_array_ops.minPointwiseDivide(numer->getArrayData(d),
            denom->getArrayData(d),
            face_box);
      retval = tbox::MathUtilities<TYPE>::Min(retval, dirval);
   }
   return retval;
}

}
}
#endif
