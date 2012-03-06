/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Operations for integer cell-centered patch data.
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
PatchCellDataOpsInteger::numberOfEntries(
   const boost::shared_ptr<pdat::CellData<int> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(data);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   return ((box * data->getGhostBox()).size()) * data->getDepth();
}

SAMRAI_INLINE_KEYWORD
void
PatchCellDataOpsInteger::printData(
   const boost::shared_ptr<pdat::CellData<int> >& data,
   const hier::Box& box,
   std::ostream& s) const
{
   TBOX_ASSERT(data);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   s << "Data box = " << box << std::endl;
   data->print(box, s);
   s << "\n";
}

SAMRAI_INLINE_KEYWORD
void
PatchCellDataOpsInteger::copyData(
   const boost::shared_ptr<pdat::CellData<int> >& dst,
   const boost::shared_ptr<pdat::CellData<int> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst, *src, box);

   dst->getArrayData().copy(src->getArrayData(), box);
}

SAMRAI_INLINE_KEYWORD
void
PatchCellDataOpsInteger::setToScalar(
   const boost::shared_ptr<pdat::CellData<int> >& dst,
   const int& alpha,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*dst, box);

   dst->fillAll(alpha, box);
}

SAMRAI_INLINE_KEYWORD
void
PatchCellDataOpsInteger::abs(
   const boost::shared_ptr<pdat::CellData<int> >& dst,
   const boost::shared_ptr<pdat::CellData<int> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst, *src, box);

   d_array_ops.abs(dst->getArrayData(),
      src->getArrayData(),
      box);
}

}
}
