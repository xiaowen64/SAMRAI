/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Norm operations for complex cell-centered patch data.
 *
 ************************************************************************/

namespace SAMRAI {
namespace math {

/*
 *************************************************************************
 *
 * Compute the number of data entries on a patch in the given box.
 *
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
int
PatchCellDataNormOpsComplex::numberOfEntries(
   const boost::shared_ptr<pdat::CellData<dcomplex> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(data);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   return ((box * data->getGhostBox()).size()) * data->getDepth();
}

/*
 *************************************************************************
 *
 * Norm operations for complex cell-centered data.
 *
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
double
PatchCellDataNormOpsComplex::sumControlVolumes(
   const boost::shared_ptr<pdat::CellData<dcomplex> >& data,
   const boost::shared_ptr<pdat::CellData<double> >& cvol,
   const hier::Box& box) const
{
   TBOX_ASSERT(data && cvol);
   TBOX_DIM_ASSERT_CHECK_ARGS3(*data, *cvol, box);

   return d_array_ops.sumControlVolumes(data->getArrayData(),
      cvol->getArrayData(),
      box);
}

SAMRAI_INLINE_KEYWORD
void
PatchCellDataNormOpsComplex::abs(
   const boost::shared_ptr<pdat::CellData<double> >& dst,
   const boost::shared_ptr<pdat::CellData<dcomplex> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst, *src, box);

   d_array_ops.abs(dst->getArrayData(),
      src->getArrayData(),
      box);
}

SAMRAI_INLINE_KEYWORD
dcomplex
PatchCellDataNormOpsComplex::integral(
   const boost::shared_ptr<pdat::CellData<dcomplex> >& data,
   const hier::Box& box,
   const boost::shared_ptr<pdat::CellData<double> >& vol) const
{
   TBOX_ASSERT(data);
   TBOX_ASSERT(vol);
   TBOX_DIM_ASSERT_CHECK_ARGS3(*data, box, *vol);

   return d_array_ops.integral(data->getArrayData(), vol->getArrayData(), box);
}

}
}
