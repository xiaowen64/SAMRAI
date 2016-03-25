/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Operations for integer side-centered patch data.
 *
 ************************************************************************/

#ifndef included_math_PatchSideDataOpsInteger_C
#define included_math_PatchSideDataOpsInteger_C

#include "SAMRAI/math/PatchSideDataOpsInteger.h"
#include "SAMRAI/pdat/SideGeometry.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#include "SAMRAI/tbox/Utilities.h"
#endif

namespace SAMRAI {
namespace math {

PatchSideDataOpsInteger::PatchSideDataOpsInteger()
{
}

PatchSideDataOpsInteger::~PatchSideDataOpsInteger()
{
}

/*
 *************************************************************************
 *
 * Compute the number of data entries on a patch in the given box.
 *
 *************************************************************************
 */

int PatchSideDataOpsInteger::numberOfEntries(
   const tbox::Pointer<pdat::SideData<int> >& data,
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
 *
 * General operations for integer side-centered patch data.
 *
 *************************************************************************
 */

void PatchSideDataOpsInteger::swapData(
   tbox::Pointer<hier::Patch> patch,
   const int data1_id,
   const int data2_id) const
{
   TBOX_ASSERT(!patch.isNull());

   tbox::Pointer<pdat::SideData<int> > d1 = patch->getPatchData(data1_id);
   tbox::Pointer<pdat::SideData<int> > d2 = patch->getPatchData(data2_id);

   TBOX_ASSERT(!d1.isNull() && !d2.isNull());
   TBOX_ASSERT(d1->getDepth() && d2->getDepth());
   TBOX_ASSERT(d1->getBox().isSpatiallyEqual(d2->getBox()));
   TBOX_ASSERT(d1->getDirectionVector() == d2->getDirectionVector());
   TBOX_ASSERT(d1->getGhostBox().isSpatiallyEqual(d2->getGhostBox()));

   patch->setPatchData(data1_id, d2);
   patch->setPatchData(data2_id, d1);
}

void PatchSideDataOpsInteger::printData(
   const tbox::Pointer<pdat::SideData<int> >& data,
   const hier::Box& box,
   std::ostream& s) const
{
   TBOX_ASSERT(!data.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   s << "Data box = " << box << std::endl;
   data->print(box, s);
   s << "\n";
}

void PatchSideDataOpsInteger::copyData(
   tbox::Pointer<pdat::SideData<int> >& dst,
   const tbox::Pointer<pdat::SideData<int> >& src,
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
         (dst->getArrayData(d)).copy(src->getArrayData(d), side_box);
      }
   }
}

void PatchSideDataOpsInteger::setToScalar(
   tbox::Pointer<pdat::SideData<int> >& dst,
   const int& alpha,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS2(*dst, box);

   dst->fillAll(alpha, box);
}

void PatchSideDataOpsInteger::abs(
   tbox::Pointer<pdat::SideData<int> >& dst,
   const tbox::Pointer<pdat::SideData<int> >& src,
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

}
}
#endif
