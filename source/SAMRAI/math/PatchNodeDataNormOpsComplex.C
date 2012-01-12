/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Norm operations for complex node-centered patch data.
 *
 ************************************************************************/

#ifndef included_math_PatchNodeDataNormOpsComplex_C
#define included_math_PatchNodeDataNormOpsComplex_C

#include "SAMRAI/math/PatchNodeDataNormOpsComplex.h"
#include "SAMRAI/pdat/NodeGeometry.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#include "SAMRAI/tbox/Utilities.h"
#endif

namespace SAMRAI {
namespace math {

PatchNodeDataNormOpsComplex::PatchNodeDataNormOpsComplex()
{
}

PatchNodeDataNormOpsComplex::~PatchNodeDataNormOpsComplex()
{
}

/*
 *************************************************************************
 *
 * Compute the number of data entries on a patch in the given box.
 *
 *************************************************************************
 */

int PatchNodeDataNormOpsComplex::numberOfEntries(
   const tbox::Pointer<pdat::NodeData<dcomplex> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(data);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   const hier::Box ibox =
      pdat::NodeGeometry::toNodeBox(box * data->getGhostBox());
   int retval = ibox.size() * data->getDepth();
   return retval;
}

/*
 *************************************************************************
 *
 * Norm operations for complex node-centered data.
 *
 *************************************************************************
 */

double PatchNodeDataNormOpsComplex::sumControlVolumes(
   const tbox::Pointer<pdat::NodeData<dcomplex> >& data,
   const tbox::Pointer<pdat::NodeData<double> >& cvol,
   const hier::Box& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(data && cvol);
#endif
   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   return d_array_ops.sumControlVolumes(data->getArrayData(),
      cvol->getArrayData(),
      node_box);
}

void PatchNodeDataNormOpsComplex::abs(
   tbox::Pointer<pdat::NodeData<double> >& dst,
   const tbox::Pointer<pdat::NodeData<dcomplex> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst, *src, box);

   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   d_array_ops.abs(dst->getArrayData(),
      src->getArrayData(),
      node_box);
}

double PatchNodeDataNormOpsComplex::L1Norm(
   const tbox::Pointer<pdat::NodeData<dcomplex> >& data,
   const hier::Box& box,
   const tbox::Pointer<pdat::NodeData<double> > cvol) const
{
   TBOX_ASSERT(data);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   double retval;
   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   if (!cvol) {
      retval = d_array_ops.L1Norm(data->getArrayData(), node_box);
   } else {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*data, *cvol);

      retval = d_array_ops.L1NormWithControlVolume(data->getArrayData(),
            cvol->getArrayData(),
            node_box);
   }
   return retval;
}

double PatchNodeDataNormOpsComplex::L2Norm(
   const tbox::Pointer<pdat::NodeData<dcomplex> >& data,
   const hier::Box& box,
   const tbox::Pointer<pdat::NodeData<double> > cvol) const
{
   TBOX_ASSERT(data);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   double retval;
   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   if (!cvol) {
      retval = d_array_ops.L2Norm(data->getArrayData(), node_box);
   } else {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*data, *cvol);

      retval = d_array_ops.L2NormWithControlVolume(data->getArrayData(),
            cvol->getArrayData(),
            node_box);
   }
   return retval;
}

double PatchNodeDataNormOpsComplex::weightedL2Norm(
   const tbox::Pointer<pdat::NodeData<dcomplex> >& data,
   const tbox::Pointer<pdat::NodeData<dcomplex> >& weight,
   const hier::Box& box,
   const tbox::Pointer<pdat::NodeData<double> > cvol) const
{
   TBOX_ASSERT(data && weight);
   TBOX_DIM_ASSERT_CHECK_ARGS3(*data, *weight, box);

   double retval;
   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   if (!cvol) {
      retval = d_array_ops.weightedL2Norm(data->getArrayData(),
            weight->getArrayData(),
            node_box);
   } else {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*data, *cvol);

      retval = d_array_ops.weightedL2NormWithControlVolume(
            data->getArrayData(),
            weight->getArrayData(),
            cvol->getArrayData(),
            node_box);
   }
   return retval;
}

double PatchNodeDataNormOpsComplex::RMSNorm(
   const tbox::Pointer<pdat::NodeData<dcomplex> >& data,
   const hier::Box& box,
   const tbox::Pointer<pdat::NodeData<double> > cvol) const
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

double PatchNodeDataNormOpsComplex::weightedRMSNorm(
   const tbox::Pointer<pdat::NodeData<dcomplex> >& data,
   const tbox::Pointer<pdat::NodeData<dcomplex> >& weight,
   const hier::Box& box,
   const tbox::Pointer<pdat::NodeData<double> > cvol) const
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

double PatchNodeDataNormOpsComplex::maxNorm(
   const tbox::Pointer<pdat::NodeData<dcomplex> >& data,
   const hier::Box& box,
   const tbox::Pointer<pdat::NodeData<double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(data);
#endif
   double retval;
   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   if (!cvol) {
      retval = d_array_ops.maxNorm(data->getArrayData(), node_box);
   } else {
      retval = d_array_ops.maxNormWithControlVolume(data->getArrayData(),
            cvol->getArrayData(),
            node_box);
   }
   return retval;
}

dcomplex PatchNodeDataNormOpsComplex::dot(
   const tbox::Pointer<pdat::NodeData<dcomplex> >& data1,
   const tbox::Pointer<pdat::NodeData<dcomplex> >& data2,
   const hier::Box& box,
   const tbox::Pointer<pdat::NodeData<double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(data1 && data2);
#endif
   dcomplex retval;
   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   if (!cvol) {
      retval = d_array_ops.dot(data1->getArrayData(),
            data2->getArrayData(),
            node_box);
   } else {
      retval = d_array_ops.dotWithControlVolume(
            data1->getArrayData(),
            data2->getArrayData(),
            cvol->getArrayData(),
            node_box);
   }
   return retval;
}

dcomplex PatchNodeDataNormOpsComplex::integral(
   const tbox::Pointer<pdat::NodeData<dcomplex> >& data,
   const hier::Box& box,
   const tbox::Pointer<pdat::NodeData<double> > vol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(data);
#endif
   dcomplex retval;
   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);

   retval = d_array_ops.integral(
         data->getArrayData(),
         vol->getArrayData(),
         node_box);

   return retval;
}

}
}
#endif
