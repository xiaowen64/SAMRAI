/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Operations for integer cell-centered patch data.
 *
 ************************************************************************/

#ifndef included_math_PatchCellDataOpsInteger_C
#define included_math_PatchCellDataOpsInteger_C

#include "SAMRAI/math/PatchCellDataOpsInteger.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#include "SAMRAI/tbox/Utilities.h"
#endif

namespace SAMRAI {
namespace math {

PatchCellDataOpsInteger::PatchCellDataOpsInteger()
{
}

PatchCellDataOpsInteger::~PatchCellDataOpsInteger()
{
}

/*
 *************************************************************************
 *
 * Compute the number of data entries on a patch in the given box.
 *
 *************************************************************************
 */

int PatchCellDataOpsInteger::numberOfEntries(
   const tbox::Pointer<pdat::CellData<int> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(data);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   const hier::Box ibox = box * data->getGhostBox();
   int retval = ibox.size() * data->getDepth();
   return retval;
}

/*
 *************************************************************************
 *
 * General operations for integer cell-centered patch data.
 *
 *************************************************************************
 */

void PatchCellDataOpsInteger::swapData(
   tbox::Pointer<hier::Patch> patch,
   const int data1_id,
   const int data2_id) const
{
   TBOX_ASSERT(patch);

   tbox::Pointer<pdat::CellData<int> > d1(
      patch->getPatchData(data1_id),
      tbox::__dynamic_cast_tag());
   tbox::Pointer<pdat::CellData<int> > d2(
      patch->getPatchData(data2_id),
      tbox::__dynamic_cast_tag());
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(d1 && d2);
   TBOX_ASSERT(d1->getDepth() && d2->getDepth());
   TBOX_ASSERT(d1->getBox().isSpatiallyEqual(d2->getBox()));
   TBOX_ASSERT(d1->getGhostBox().isSpatiallyEqual(d2->getGhostBox()));
#endif
   patch->setPatchData(data1_id, d2);
   patch->setPatchData(data2_id, d1);
}

void PatchCellDataOpsInteger::printData(
   const tbox::Pointer<pdat::CellData<int> >& data,
   const hier::Box& box,
   std::ostream& s) const
{
   TBOX_ASSERT(data);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   s << "Data box = " << box << std::endl;
   data->print(box, s);
   s << "\n";
}

void PatchCellDataOpsInteger::copyData(
   tbox::Pointer<pdat::CellData<int> >& dst,
   const tbox::Pointer<pdat::CellData<int> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst, *src, box);

   (dst->getArrayData()).copy(src->getArrayData(), box);
}

void PatchCellDataOpsInteger::setToScalar(
   tbox::Pointer<pdat::CellData<int> >& dst,
   const int& alpha,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*dst, box);

   dst->fillAll(alpha, box);
}

void PatchCellDataOpsInteger::abs(
   tbox::Pointer<pdat::CellData<int> >& dst,
   const tbox::Pointer<pdat::CellData<int> >& src,
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
#endif
