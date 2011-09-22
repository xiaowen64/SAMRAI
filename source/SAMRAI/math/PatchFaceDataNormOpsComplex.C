/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Norm operations for complex face-centered patch data.
 *
 ************************************************************************/

#ifndef included_math_PatchFaceDataNormOpsComplex_C
#define included_math_PatchFaceDataNormOpsComplex_C

#include "SAMRAI/math/PatchFaceDataNormOpsComplex.h"
#include "SAMRAI/pdat/FaceGeometry.h"
#include "SAMRAI/tbox/MathUtilities.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#include "SAMRAI/tbox/Utilities.h"
#endif

namespace SAMRAI {
namespace math {

PatchFaceDataNormOpsComplex::PatchFaceDataNormOpsComplex()
{
}

PatchFaceDataNormOpsComplex::~PatchFaceDataNormOpsComplex()
{
}

/*
 *************************************************************************
 *
 * Compute the number of data entries on a patch in the given box.
 *
 *************************************************************************
 */

int PatchFaceDataNormOpsComplex::numberOfEntries(
   const tbox::Pointer<pdat::FaceData<dcomplex> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(!data.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   const tbox::Dimension& dim(box.getDim());

   int retval = 0;
   const hier::Box ibox = box * data->getGhostBox();
   for (int d = 0; d < dim.getValue(); d++) {
      const hier::Box dbox = pdat::FaceGeometry::toFaceBox(ibox, d);
      retval += (dbox.size() * data->getDepth());
   }
   return retval;
}

/*
 *************************************************************************
 *
 * Norm operations for complex face-centered data.
 *
 *************************************************************************
 */

double PatchFaceDataNormOpsComplex::sumControlVolumes(
   const tbox::Pointer<pdat::FaceData<dcomplex> >& data,
   const tbox::Pointer<pdat::FaceData<double> >& cvol,
   const hier::Box& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull() && !cvol.isNull());
#endif
   const tbox::Dimension& dim(box.getDim());

   double retval = 0.0;
   for (int d = 0; d < dim.getValue(); d++) {
      const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
      retval += d_array_ops.sumControlVolumes(data->getArrayData(d),
            cvol->getArrayData(d),
            face_box);
   }
   return retval;
}

void PatchFaceDataNormOpsComplex::abs(
   tbox::Pointer<pdat::FaceData<double> >& dst,
   const tbox::Pointer<pdat::FaceData<dcomplex> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst, *src, box);

   const tbox::Dimension& dim(box.getDim());

   for (int d = 0; d < dim.getValue(); d++) {
      const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
      d_array_ops.abs(dst->getArrayData(d),
         src->getArrayData(d),
         face_box);
   }
}

double PatchFaceDataNormOpsComplex::L1Norm(
   const tbox::Pointer<pdat::FaceData<dcomplex> >& data,
   const hier::Box& box,
   const tbox::Pointer<pdat::FaceData<double> > cvol) const
{
   TBOX_ASSERT(!data.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   const tbox::Dimension& dim(box.getDim());

   double retval = 0.0;
   if (cvol.isNull()) {
      for (int d = 0; d < dim.getValue(); d++) {
         const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
         retval += d_array_ops.L1Norm(data->getArrayData(d), face_box);
      }
   } else {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*data, *cvol);

      for (int d = 0; d < dim.getValue(); d++) {
         const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
         retval += d_array_ops.L1NormWithControlVolume(data->getArrayData(d),
               cvol->getArrayData(d),
               face_box);
      }
   }
   return retval;
}

double PatchFaceDataNormOpsComplex::L2Norm(
   const tbox::Pointer<pdat::FaceData<dcomplex> >& data,
   const hier::Box& box,
   const tbox::Pointer<pdat::FaceData<double> > cvol) const
{
   TBOX_ASSERT(!data.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   const tbox::Dimension& dim(box.getDim());

   double retval = 0.0;
   if (cvol.isNull()) {
      for (int d = 0; d < dim.getValue(); d++) {
         const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
         double aval = d_array_ops.L2Norm(data->getArrayData(d), face_box);
         retval += aval * aval;
      }
   } else {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*data, *cvol);

      for (int d = 0; d < dim.getValue(); d++) {
         const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
         double aval = d_array_ops.L2NormWithControlVolume(
               data->getArrayData(d),
               cvol->getArrayData(d),
               face_box);
         retval += aval * aval;
      }
   }
   return sqrt(retval);
}

double PatchFaceDataNormOpsComplex::weightedL2Norm(
   const tbox::Pointer<pdat::FaceData<dcomplex> >& data,
   const tbox::Pointer<pdat::FaceData<dcomplex> >& weight,
   const hier::Box& box,
   const tbox::Pointer<pdat::FaceData<double> > cvol) const
{
   TBOX_ASSERT(!data.isNull() && !weight.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS3(*data, *weight, box);

   const tbox::Dimension& dim(box.getDim());

   double retval = 0.0;
   if (cvol.isNull()) {
      for (int d = 0; d < dim.getValue(); d++) {
         const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
         double aval = d_array_ops.weightedL2Norm(data->getArrayData(d),
               weight->getArrayData(d),
               face_box);
         retval += aval * aval;
      }
   } else {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*data, *cvol);

      for (int d = 0; d < dim.getValue(); d++) {
         const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
         double aval = d_array_ops.weightedL2NormWithControlVolume(
               data->getArrayData(d),
               weight->getArrayData(d),
               cvol->getArrayData(d),
               face_box);
         retval += aval * aval;
      }
   }
   return sqrt(retval);
}

double PatchFaceDataNormOpsComplex::RMSNorm(
   const tbox::Pointer<pdat::FaceData<dcomplex> >& data,
   const hier::Box& box,
   const tbox::Pointer<pdat::FaceData<double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
#endif
   double retval = L2Norm(data, box, cvol);
   if (cvol.isNull()) {
      retval /= sqrt((double)numberOfEntries(data, box));
   } else {
      retval /= sqrt(sumControlVolumes(data, cvol, box));
   }
   return retval;
}

double PatchFaceDataNormOpsComplex::weightedRMSNorm(
   const tbox::Pointer<pdat::FaceData<dcomplex> >& data,
   const tbox::Pointer<pdat::FaceData<dcomplex> >& weight,
   const hier::Box& box,
   const tbox::Pointer<pdat::FaceData<double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull() && !weight.isNull());
#endif
   double retval = weightedL2Norm(data, weight, box, cvol);
   if (cvol.isNull()) {
      retval /= sqrt((double)numberOfEntries(data, box));
   } else {
      retval /= sqrt(sumControlVolumes(data, cvol, box));
   }
   return retval;
}

double PatchFaceDataNormOpsComplex::maxNorm(
   const tbox::Pointer<pdat::FaceData<dcomplex> >& data,
   const hier::Box& box,
   const tbox::Pointer<pdat::FaceData<double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
#endif
   const tbox::Dimension& dim(box.getDim());

   double retval = 0.0;
   if (cvol.isNull()) {
      for (int d = 0; d < dim.getValue(); d++) {
         const hier::Box face_box =
            pdat::FaceGeometry::toFaceBox(box, d);
         retval = tbox::MathUtilities<double>::Max(retval,
               d_array_ops.maxNorm(data->getArrayData(d), face_box));
      }
   } else {
      for (int d = 0; d < dim.getValue(); d++) {
         const hier::Box face_box =
            pdat::FaceGeometry::toFaceBox(box, d);
         retval = tbox::MathUtilities<double>::Max(retval,
               d_array_ops.maxNormWithControlVolume(
                  data->getArrayData(d), cvol->getArrayData(d), face_box));
      }
   }
   return retval;
}

dcomplex PatchFaceDataNormOpsComplex::dot(
   const tbox::Pointer<pdat::FaceData<dcomplex> >& data1,
   const tbox::Pointer<pdat::FaceData<dcomplex> >& data2,
   const hier::Box& box,
   const tbox::Pointer<pdat::FaceData<double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data1.isNull() && !data2.isNull());
#endif
   const tbox::Dimension& dim(box.getDim());

   dcomplex retval = dcomplex(0.0, 0.0);
   if (cvol.isNull()) {
      for (int d = 0; d < dim.getValue(); d++) {
         const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
         retval += d_array_ops.dot(data1->getArrayData(d),
               data2->getArrayData(d),
               face_box);
      }
   } else {
      for (int d = 0; d < dim.getValue(); d++) {
         const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
         retval += d_array_ops.dotWithControlVolume(
               data1->getArrayData(d),
               data2->getArrayData(d),
               cvol->getArrayData(d),
               face_box);
      }
   }
   return retval;
}

dcomplex PatchFaceDataNormOpsComplex::integral(
   const tbox::Pointer<pdat::FaceData<dcomplex> >& data,
   const hier::Box& box,
   const tbox::Pointer<pdat::FaceData<double> > vol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
#endif
   const tbox::Dimension& dim(box.getDim());

   dcomplex retval = dcomplex(0.0, 0.0);

   for (int d = 0; d < dim.getValue(); d++) {
      const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
      retval += d_array_ops.integral(data->getArrayData(d),
            vol->getArrayData(d),
            face_box);
   }

   return retval;
}

}
}
#endif
