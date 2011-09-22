/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Norm operations for complex cell-centered patch data.
 *
 ************************************************************************/

#ifndef included_math_PatchCellDataNormOpsComplex_C
#define included_math_PatchCellDataNormOpsComplex_C

#include "SAMRAI/math/PatchCellDataNormOpsComplex.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#include "SAMRAI/tbox/Utilities.h"
#endif

namespace SAMRAI {
namespace math {

PatchCellDataNormOpsComplex::PatchCellDataNormOpsComplex()
{
}

PatchCellDataNormOpsComplex::~PatchCellDataNormOpsComplex()
{
}

/*
 *************************************************************************
 *
 * Compute the number of data entries on a patch in the given box.
 *
 *************************************************************************
 */

int PatchCellDataNormOpsComplex::numberOfEntries(
   const tbox::Pointer<pdat::CellData<dcomplex> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(!data.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   const hier::Box ibox = box * data->getGhostBox();
   int retval = ibox.size() * data->getDepth();
   return retval;
}

/*
 *************************************************************************
 *
 * Norm operations for complex cell-centered data.
 *
 *************************************************************************
 */

double PatchCellDataNormOpsComplex::sumControlVolumes(
   const tbox::Pointer<pdat::CellData<dcomplex> >& data,
   const tbox::Pointer<pdat::CellData<double> >& cvol,
   const hier::Box& box) const
{
   TBOX_ASSERT(!data.isNull() && !cvol.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS3(*data, *cvol, box);

   return d_array_ops.sumControlVolumes(data->getArrayData(),
      cvol->getArrayData(),
      box);
}

void PatchCellDataNormOpsComplex::abs(
   tbox::Pointer<pdat::CellData<double> >& dst,
   const tbox::Pointer<pdat::CellData<dcomplex> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst, *src, box);

   d_array_ops.abs(dst->getArrayData(),
      src->getArrayData(),
      box);
}

double PatchCellDataNormOpsComplex::L1Norm(
   const tbox::Pointer<pdat::CellData<dcomplex> >& data,
   const hier::Box& box,
   const tbox::Pointer<pdat::CellData<double> > cvol) const
{
   TBOX_ASSERT(!data.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   double retval;
   if (cvol.isNull()) {
      retval = d_array_ops.L1Norm(data->getArrayData(), box);
   } else {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*data, *cvol);

      retval = d_array_ops.L1NormWithControlVolume(data->getArrayData(),
            cvol->getArrayData(),
            box);
   }
   return retval;
}

double PatchCellDataNormOpsComplex::L2Norm(
   const tbox::Pointer<pdat::CellData<dcomplex> >& data,
   const hier::Box& box,
   const tbox::Pointer<pdat::CellData<double> > cvol) const
{
   TBOX_ASSERT(!data.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   double retval;
   if (cvol.isNull()) {
      retval = d_array_ops.L2Norm(data->getArrayData(), box);
   } else {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*data, *cvol);

      retval = d_array_ops.L2NormWithControlVolume(data->getArrayData(),
            cvol->getArrayData(),
            box);
   }
   return retval;
}

double PatchCellDataNormOpsComplex::weightedL2Norm(
   const tbox::Pointer<pdat::CellData<dcomplex> >& data,
   const tbox::Pointer<pdat::CellData<dcomplex> >& weight,
   const hier::Box& box,
   const tbox::Pointer<pdat::CellData<double> > cvol) const
{
   TBOX_ASSERT(!data.isNull() && !weight.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS3(*data, *weight, box);

   double retval;
   if (cvol.isNull()) {
      retval = d_array_ops.weightedL2Norm(data->getArrayData(),
            weight->getArrayData(),
            box);
   } else {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*data, *cvol);

      retval = d_array_ops.weightedL2NormWithControlVolume(
            data->getArrayData(),
            weight->getArrayData(),
            cvol->getArrayData(),
            box);
   }
   return retval;
}

double PatchCellDataNormOpsComplex::RMSNorm(
   const tbox::Pointer<pdat::CellData<dcomplex> >& data,
   const hier::Box& box,
   const tbox::Pointer<pdat::CellData<double> > cvol) const
{
   TBOX_ASSERT(!data.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   double retval = L2Norm(data, box, cvol);
   if (cvol.isNull()) {
      retval /= sqrt((double)numberOfEntries(data, box));
   } else {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*data, *cvol);

      retval /= sqrt(sumControlVolumes(data, cvol, box));
   }
   return retval;
}

double PatchCellDataNormOpsComplex::weightedRMSNorm(
   const tbox::Pointer<pdat::CellData<dcomplex> >& data,
   const tbox::Pointer<pdat::CellData<dcomplex> >& weight,
   const hier::Box& box,
   const tbox::Pointer<pdat::CellData<double> > cvol) const
{
   TBOX_ASSERT(!data.isNull() && !weight.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS3(*data, *weight, box);

   double retval = weightedL2Norm(data, weight, box, cvol);
   if (cvol.isNull()) {
      retval /= sqrt((double)numberOfEntries(data, box));
   } else {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*data, *cvol);

      retval /= sqrt(sumControlVolumes(data, cvol, box));
   }
   return retval;
}

double PatchCellDataNormOpsComplex::maxNorm(
   const tbox::Pointer<pdat::CellData<dcomplex> >& data,
   const hier::Box& box,
   const tbox::Pointer<pdat::CellData<double> > cvol) const
{
   TBOX_ASSERT(!data.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   double retval;
   if (cvol.isNull()) {
      retval = d_array_ops.maxNorm(data->getArrayData(), box);
   } else {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*data, *cvol);

      retval = d_array_ops.maxNormWithControlVolume(data->getArrayData(),
            cvol->getArrayData(),
            box);
   }
   return retval;
}

dcomplex PatchCellDataNormOpsComplex::dot(
   const tbox::Pointer<pdat::CellData<dcomplex> >& data1,
   const tbox::Pointer<pdat::CellData<dcomplex> >& data2,
   const hier::Box& box,
   const tbox::Pointer<pdat::CellData<double> > cvol) const
{
   TBOX_ASSERT(!data1.isNull() && !data2.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS3(*data1, *data2, box);

   dcomplex retval;
   if (cvol.isNull()) {
      retval = d_array_ops.dot(data1->getArrayData(),
            data2->getArrayData(),
            box);
   } else {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*data1, *cvol);

      retval = d_array_ops.dotWithControlVolume(
            data1->getArrayData(),
            data2->getArrayData(),
            cvol->getArrayData(),
            box);
   }
   return retval;
}

dcomplex PatchCellDataNormOpsComplex::integral(
   const tbox::Pointer<pdat::CellData<dcomplex> >& data,
   const hier::Box& box,
   const tbox::Pointer<pdat::CellData<double> > vol) const
{
   TBOX_ASSERT(!data.isNull());
   TBOX_ASSERT(!vol.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS3(*data, box, *vol);

   dcomplex retval = dcomplex(0.0, 0.0);

   retval = d_array_ops.integral(data->getArrayData(),
         vol->getArrayData(),
         box);

   return retval;
}

}
}
#endif
