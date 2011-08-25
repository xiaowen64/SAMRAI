/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Norm operations for complex side-centered patch data.
 *
 ************************************************************************/

#ifndef included_math_PatchSideDataNormOpsComplex_C
#define included_math_PatchSideDataNormOpsComplex_C

#include "SAMRAI/math/PatchSideDataNormOpsComplex.h"
#include "SAMRAI/pdat/SideGeometry.h"
#include "SAMRAI/tbox/MathUtilities.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#include "SAMRAI/tbox/Utilities.h"
#endif

namespace SAMRAI {
namespace math {

PatchSideDataNormOpsComplex::PatchSideDataNormOpsComplex()
{
}

PatchSideDataNormOpsComplex::~PatchSideDataNormOpsComplex()
{
}

/*
 *************************************************************************
 *                                                                       *
 * Compute the number of data entries on a patch in the given box.       *
 *                                                                       *
 *************************************************************************
 */

int PatchSideDataNormOpsComplex::numberOfEntries(
   const tbox::Pointer<pdat::SideData<dcomplex> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(!data.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   const tbox::Dimension& dim(box.getDim());

   int retval = 0;
   const hier::Box ibox = box * data->getGhostBox();
   const hier::IntVector& directions = data->getDirectionVector();
   for (int d = 0; d < dim.getValue(); d++) {
      if (directions(d)) {
         const hier::Box dbox = pdat::SideGeometry::toSideBox(ibox, d);
         retval += (dbox.size() * data->getDepth());
      }
   }
   return retval;
}

/*
 *************************************************************************
 *                                                                       *
 * Norm operations for complex side-centered data.                       *
 *                                                                       *
 *************************************************************************
 */

double PatchSideDataNormOpsComplex::sumControlVolumes(
   const tbox::Pointer<pdat::SideData<dcomplex> >& data,
   const tbox::Pointer<pdat::SideData<double> >& cvol,
   const hier::Box& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull() && !cvol.isNull());
#endif
   double retval = 0.0;
   const hier::IntVector& directions = data->getDirectionVector();
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(directions ==
      hier::IntVector::min(directions, cvol->getDirectionVector()));
#endif
   const tbox::Dimension& dim(box.getDim());

   for (int d = 0; d < dim.getValue(); d++) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         retval += d_array_ops.sumControlVolumes(data->getArrayData(d),
               cvol->getArrayData(d),
               side_box);
      }
   }
   return retval;
}

void PatchSideDataNormOpsComplex::abs(
   tbox::Pointer<pdat::SideData<double> >& dst,
   const tbox::Pointer<pdat::SideData<dcomplex> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
   TBOX_ASSERT(dst->getDirectionVector() == src->getDirectionVector());
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst, *src, box);

   const tbox::Dimension& dim(box.getDim());

   const hier::IntVector& directions = dst->getDirectionVector();
   for (int d = 0; d < dim.getValue(); d++) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         d_array_ops.abs(dst->getArrayData(d),
            src->getArrayData(d),
            side_box);
      }
   }
}

double PatchSideDataNormOpsComplex::L1Norm(
   const tbox::Pointer<pdat::SideData<dcomplex> >& data,
   const hier::Box& box,
   const tbox::Pointer<pdat::SideData<double> > cvol) const
{
   TBOX_ASSERT(!data.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   const tbox::Dimension& dim(box.getDim());

   double retval = 0.0;
   const hier::IntVector& directions = data->getDirectionVector();
   if (cvol.isNull()) {
      for (int d = 0; d < dim.getValue(); d++) {
         if (directions(d)) {
            const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
            retval += d_array_ops.L1Norm(data->getArrayData(d), side_box);
         }
      }
   } else {
      TBOX_ASSERT(directions ==
         hier::IntVector::min(directions, cvol->getDirectionVector()));
      TBOX_DIM_ASSERT_CHECK_ARGS2(*data, *cvol);

      for (int d = 0; d < dim.getValue(); d++) {
         if (directions(d)) {
            const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
            retval += d_array_ops.L1NormWithControlVolume(data->getArrayData(d),
                  cvol->getArrayData(d),
                  side_box);
         }
      }
   }
   return retval;
}

double PatchSideDataNormOpsComplex::L2Norm(
   const tbox::Pointer<pdat::SideData<dcomplex> >& data,
   const hier::Box& box,
   const tbox::Pointer<pdat::SideData<double> > cvol) const
{
   TBOX_ASSERT(!data.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   const tbox::Dimension& dim(box.getDim());

   double retval = 0.0;
   const hier::IntVector& directions = data->getDirectionVector();
   if (cvol.isNull()) {
      for (int d = 0; d < dim.getValue(); d++) {
         if (directions(d)) {
            const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
            double aval = d_array_ops.L2Norm(data->getArrayData(d), side_box);
            retval += aval * aval;
         }
      }
   } else {
      TBOX_ASSERT(directions ==
         hier::IntVector::min(directions, cvol->getDirectionVector()));
      TBOX_DIM_ASSERT_CHECK_ARGS2(*data, *cvol);

      for (int d = 0; d < dim.getValue(); d++) {
         if (directions(d)) {
            const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
            double aval = d_array_ops.L2NormWithControlVolume(
                  data->getArrayData(d),
                  cvol->getArrayData(d),
                  side_box);
            retval += aval * aval;
         }
      }
   }
   return sqrt(retval);
}

double PatchSideDataNormOpsComplex::weightedL2Norm(
   const tbox::Pointer<pdat::SideData<dcomplex> >& data,
   const tbox::Pointer<pdat::SideData<dcomplex> >& weight,
   const hier::Box& box,
   const tbox::Pointer<pdat::SideData<double> > cvol) const
{
   TBOX_ASSERT(!data.isNull() && !weight.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS3(*data, *weight, box);

   const tbox::Dimension& dim(box.getDim());

   double retval = 0.0;
   const hier::IntVector& directions = data->getDirectionVector();
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(directions ==
      hier::IntVector::min(directions, weight->getDirectionVector()));
#endif
   if (cvol.isNull()) {
      for (int d = 0; d < dim.getValue(); d++) {
         if (directions(d)) {
            const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
            double aval = d_array_ops.weightedL2Norm(data->getArrayData(d),
                  weight->getArrayData(d),
                  side_box);
            retval += aval * aval;
         }
      }
   } else {
      TBOX_ASSERT(directions ==
         hier::IntVector::min(directions, cvol->getDirectionVector()));
      TBOX_DIM_ASSERT_CHECK_ARGS2(*data, *cvol);

      for (int d = 0; d < dim.getValue(); d++) {
         if (directions(d)) {
            const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
            double aval = d_array_ops.weightedL2NormWithControlVolume(
                  data->getArrayData(d),
                  weight->getArrayData(d),
                  cvol->getArrayData(d),
                  side_box);
            retval += aval * aval;
         }
      }
   }
   return sqrt(retval);
}

double PatchSideDataNormOpsComplex::RMSNorm(
   const tbox::Pointer<pdat::SideData<dcomplex> >& data,
   const hier::Box& box,
   const tbox::Pointer<pdat::SideData<double> > cvol) const
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

double PatchSideDataNormOpsComplex::weightedRMSNorm(
   const tbox::Pointer<pdat::SideData<dcomplex> >& data,
   const tbox::Pointer<pdat::SideData<dcomplex> >& weight,
   const hier::Box& box,
   const tbox::Pointer<pdat::SideData<double> > cvol) const
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

double PatchSideDataNormOpsComplex::maxNorm(
   const tbox::Pointer<pdat::SideData<dcomplex> >& data,
   const hier::Box& box,
   const tbox::Pointer<pdat::SideData<double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
#endif
   const tbox::Dimension& dim(box.getDim());

   double retval = 0.0;
   const hier::IntVector& directions = data->getDirectionVector();
   if (cvol.isNull()) {
      for (int d = 0; d < dim.getValue(); d++) {
         if (directions(d)) {
            const hier::Box side_box =
               pdat::SideGeometry::toSideBox(box, d);
            retval = tbox::MathUtilities<double>::Max(retval,
                  d_array_ops.maxNorm(data->getArrayData(d), side_box));
         }
      }
   } else {
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(directions ==
         hier::IntVector::min(directions, cvol->getDirectionVector()));
#endif
      for (int d = 0; d < dim.getValue(); d++) {
         if (directions(d)) {
            const hier::Box side_box =
               pdat::SideGeometry::toSideBox(box, d);
            retval = tbox::MathUtilities<double>::Max(retval,
                  d_array_ops.maxNormWithControlVolume(
                     data->getArrayData(d),
                     cvol->getArrayData(d), side_box));
         }
      }
   }
   return retval;
}

dcomplex PatchSideDataNormOpsComplex::dot(
   const tbox::Pointer<pdat::SideData<dcomplex> >& data1,
   const tbox::Pointer<pdat::SideData<dcomplex> >& data2,
   const hier::Box& box,
   const tbox::Pointer<pdat::SideData<double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data1.isNull() && !data2.isNull());
   TBOX_ASSERT(data1->getDirectionVector() == data2->getDirectionVector());
#endif
   const tbox::Dimension& dim(box.getDim());

   dcomplex retval = dcomplex(0.0, 0.0);
   const hier::IntVector& directions = data1->getDirectionVector();
   if (cvol.isNull()) {
      for (int d = 0; d < dim.getValue(); d++) {
         if (directions(d)) {
            const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
            retval += d_array_ops.dot(data1->getArrayData(d),
                  data2->getArrayData(d),
                  side_box);
         }
      }
   } else {
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(directions ==
         hier::IntVector::min(directions, cvol->getDirectionVector()));
#endif
      for (int d = 0; d < dim.getValue(); d++) {
         if (directions(d)) {
            const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
            retval += d_array_ops.dotWithControlVolume(
                  data1->getArrayData(d),
                  data2->getArrayData(d),
                  cvol->getArrayData(d),
                  side_box);
         }
      }
   }
   return retval;
}

dcomplex PatchSideDataNormOpsComplex::integral(
   const tbox::Pointer<pdat::SideData<dcomplex> >& data,
   const hier::Box& box,
   const tbox::Pointer<pdat::SideData<double> > vol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
#endif
   const tbox::Dimension& dim(box.getDim());

   dcomplex retval = dcomplex(0.0, 0.0);
   const hier::IntVector& directions = data->getDirectionVector();

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(directions ==
      hier::IntVector::min(directions, vol->getDirectionVector()));
#endif
   for (int d = 0; d < dim.getValue(); d++) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         retval += d_array_ops.integral(
               data->getArrayData(d),
               vol->getArrayData(d),
               side_box);
      }
   }

   return retval;
}

}
}
#endif
