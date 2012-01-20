/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Operations for integer face-centered patch data.
 *
 ************************************************************************/

#ifndef included_math_PatchFaceDataOpsInteger_C
#define included_math_PatchFaceDataOpsInteger_C

#include "SAMRAI/math/PatchFaceDataOpsInteger.h"
#include "SAMRAI/pdat/FaceGeometry.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#include "SAMRAI/tbox/Utilities.h"
#endif

namespace SAMRAI {
namespace math {

PatchFaceDataOpsInteger::PatchFaceDataOpsInteger()
{
}

PatchFaceDataOpsInteger::~PatchFaceDataOpsInteger()
{
}

/*
 *************************************************************************
 *
 * Compute the number of data entries on a patch in the given box.
 *
 *************************************************************************
 */

int PatchFaceDataOpsInteger::numberOfEntries(
   const boost::shared_ptr<pdat::FaceData<int> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(data);
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
 * General operations for integer face-centered patch data.
 *
 *************************************************************************
 */

void PatchFaceDataOpsInteger::swapData(
   boost::shared_ptr<hier::Patch> patch,
   const int data1_id,
   const int data2_id) const
{
   TBOX_ASSERT(patch);

   boost::shared_ptr<pdat::FaceData<int> > d1(
      patch->getPatchData(data1_id),
      boost::detail::dynamic_cast_tag());
   boost::shared_ptr<pdat::FaceData<int> > d2(
      patch->getPatchData(data2_id),
      boost::detail::dynamic_cast_tag());
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(d1 && d2);
   TBOX_ASSERT(d1->getDepth() && d2->getDepth());
   TBOX_ASSERT(d1->getBox().isSpatiallyEqual(d2->getBox()));
   TBOX_ASSERT(d1->getGhostBox().isSpatiallyEqual(d2->getGhostBox()));
#endif
   patch->setPatchData(data1_id, d2);
   patch->setPatchData(data2_id, d1);
}

void PatchFaceDataOpsInteger::printData(
   const boost::shared_ptr<pdat::FaceData<int> >& data,
   const hier::Box& box,
   std::ostream& s) const
{
   TBOX_ASSERT(data);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   s << "Data box = " << box << std::endl;
   data->print(box, s);
   s << "\n";
}

void PatchFaceDataOpsInteger::copyData(
   boost::shared_ptr<pdat::FaceData<int> >& dst,
   const boost::shared_ptr<pdat::FaceData<int> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst, *src, box);

   const tbox::Dimension& dim(box.getDim());

   for (int d = 0; d < dim.getValue(); d++) {
      const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
      (dst->getArrayData(d)).copy(src->getArrayData(d), face_box);
   }
}

void PatchFaceDataOpsInteger::setToScalar(
   boost::shared_ptr<pdat::FaceData<int> >& dst,
   const int& alpha,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*dst, box);

   dst->fillAll(alpha, box);
}

void PatchFaceDataOpsInteger::abs(
   boost::shared_ptr<pdat::FaceData<int> >& dst,
   const boost::shared_ptr<pdat::FaceData<int> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst, *src, box);

   const tbox::Dimension& dim(box.getDim());

   for (int d = 0; d < dim.getValue(); d++) {
      const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
      d_array_ops.abs(dst->getArrayData(d),
         src->getArrayData(d),
         face_box);
   }
}

}
}
#endif
