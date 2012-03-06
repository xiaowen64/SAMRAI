/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Norm operations for complex edge-centered patch data.
 *
 ************************************************************************/

#ifndef included_math_PatchEdgeDataNormOpsComplex_C
#define included_math_PatchEdgeDataNormOpsComplex_C

#include "SAMRAI/math/PatchEdgeDataNormOpsComplex.h"
#include "SAMRAI/pdat/EdgeGeometry.h"
#include "SAMRAI/tbox/MathUtilities.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#include "SAMRAI/tbox/Utilities.h"
#endif

#ifndef SAMRAI_INLINE
#include "SAMRAI/math/PatchEdgeDataNormOpsComplex.I"
#endif

namespace SAMRAI {
namespace math {

PatchEdgeDataNormOpsComplex::PatchEdgeDataNormOpsComplex()
{
}

PatchEdgeDataNormOpsComplex::~PatchEdgeDataNormOpsComplex()
{
}

/*
 *************************************************************************
 *
 * Norm operations for complex edge-centered data.
 *
 *************************************************************************
 */

double
PatchEdgeDataNormOpsComplex::L1Norm(
   const boost::shared_ptr<pdat::EdgeData<dcomplex> >& data,
   const hier::Box& box,
   const boost::shared_ptr<pdat::EdgeData<double> >& cvol) const
{
   TBOX_ASSERT(data);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   int dimVal = box.getDim().getValue();

   double retval = 0.0;
   if (!cvol) {
      for (int d = 0; d < dimVal; d++) {
         const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
         retval += d_array_ops.L1Norm(data->getArrayData(d), edge_box);
      }
   } else {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*data, *cvol);

      for (int d = 0; d < dimVal; d++) {
         const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
         retval += d_array_ops.L1NormWithControlVolume(data->getArrayData(d),
               cvol->getArrayData(d),
               edge_box);
      }
   }
   return retval;
}

double
PatchEdgeDataNormOpsComplex::L2Norm(
   const boost::shared_ptr<pdat::EdgeData<dcomplex> >& data,
   const hier::Box& box,
   const boost::shared_ptr<pdat::EdgeData<double> >& cvol) const
{
   TBOX_ASSERT(data);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   int dimVal = box.getDim().getValue();

   double retval = 0.0;
   if (!cvol) {
      for (int d = 0; d < dimVal; d++) {
         const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
         double aval = d_array_ops.L2Norm(data->getArrayData(d), edge_box);
         retval += aval * aval;
      }
   } else {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*data, *cvol);

      for (int d = 0; d < dimVal; d++) {
         const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
         double aval = d_array_ops.L2NormWithControlVolume(
               data->getArrayData(d),
               cvol->getArrayData(d),
               edge_box);
         retval += aval * aval;
      }
   }
   return sqrt(retval);
}

double
PatchEdgeDataNormOpsComplex::weightedL2Norm(
   const boost::shared_ptr<pdat::EdgeData<dcomplex> >& data,
   const boost::shared_ptr<pdat::EdgeData<dcomplex> >& weight,
   const hier::Box& box,
   const boost::shared_ptr<pdat::EdgeData<double> >& cvol) const
{
   TBOX_ASSERT(data && weight);
   TBOX_DIM_ASSERT_CHECK_ARGS3(*data, *weight, box);

   int dimVal = box.getDim().getValue();

   double retval = 0.0;
   if (!cvol) {
      for (int d = 0; d < dimVal; d++) {
         const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
         double aval = d_array_ops.weightedL2Norm(data->getArrayData(d),
               weight->getArrayData(d),
               edge_box);
         retval += aval * aval;
      }
   } else {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*data, *cvol);

      for (int d = 0; d < dimVal; d++) {
         const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
         double aval = d_array_ops.weightedL2NormWithControlVolume(
               data->getArrayData(d),
               weight->getArrayData(d),
               cvol->getArrayData(d),
               edge_box);
         retval += aval * aval;
      }
   }
   return sqrt(retval);
}

double
PatchEdgeDataNormOpsComplex::RMSNorm(
   const boost::shared_ptr<pdat::EdgeData<dcomplex> >& data,
   const hier::Box& box,
   const boost::shared_ptr<pdat::EdgeData<double> >& cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(data);
#endif
   double retval = L2Norm(data, box, cvol);
   if (!cvol) {
      retval /= sqrt((double)numberOfEntries(data, box));
   } else {
      retval /= sqrt(sumControlVolumes(data, cvol, box));
   }
   return retval;
}

double
PatchEdgeDataNormOpsComplex::weightedRMSNorm(
   const boost::shared_ptr<pdat::EdgeData<dcomplex> >& data,
   const boost::shared_ptr<pdat::EdgeData<dcomplex> >& weight,
   const hier::Box& box,
   const boost::shared_ptr<pdat::EdgeData<double> >& cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(data && weight);
#endif
   double retval = weightedL2Norm(data, weight, box, cvol);
   if (!cvol) {
      retval /= sqrt((double)numberOfEntries(data, box));
   } else {
      retval /= sqrt(sumControlVolumes(data, cvol, box));
   }
   return retval;
}

double
PatchEdgeDataNormOpsComplex::maxNorm(
   const boost::shared_ptr<pdat::EdgeData<dcomplex> >& data,
   const hier::Box& box,
   const boost::shared_ptr<pdat::EdgeData<double> >& cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(data);
#endif
   int dimVal = box.getDim().getValue();

   double retval = 0.0;
   if (!cvol) {
      for (int d = 0; d < dimVal; d++) {
         const hier::Box edge_box =
            pdat::EdgeGeometry::toEdgeBox(box, d);
         retval = tbox::MathUtilities<double>::Max(retval,
               d_array_ops.maxNorm(data->getArrayData(d), edge_box));
      }
   } else {
      for (int d = 0; d < dimVal; d++) {
         const hier::Box edge_box =
            pdat::EdgeGeometry::toEdgeBox(box, d);
         retval = tbox::MathUtilities<double>::Max(retval,
               d_array_ops.maxNormWithControlVolume(
                  data->getArrayData(d), cvol->getArrayData(d), edge_box));
      }
   }
   return retval;
}

dcomplex
PatchEdgeDataNormOpsComplex::dot(
   const boost::shared_ptr<pdat::EdgeData<dcomplex> >& data1,
   const boost::shared_ptr<pdat::EdgeData<dcomplex> >& data2,
   const hier::Box& box,
   const boost::shared_ptr<pdat::EdgeData<double> >& cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(data1 && data2);
#endif
   int dimVal = box.getDim().getValue();

   dcomplex retval = dcomplex(0.0, 0.0);
   if (!cvol) {
      for (int d = 0; d < dimVal; d++) {
         const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
         retval += d_array_ops.dot(data1->getArrayData(d),
               data2->getArrayData(d),
               edge_box);
      }
   } else {
      for (int d = 0; d < dimVal; d++) {
         const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
         retval += d_array_ops.dotWithControlVolume(
               data1->getArrayData(d),
               data2->getArrayData(d),
               cvol->getArrayData(d),
               edge_box);
      }
   }
   return retval;
}

}
}
#endif
