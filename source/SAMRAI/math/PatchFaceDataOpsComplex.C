/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Operations for complex face-centered patch data.
 *
 ************************************************************************/

#ifndef included_math_PatchFaceDataOpsComplex_C
#define included_math_PatchFaceDataOpsComplex_C

#include "SAMRAI/math/PatchFaceDataOpsComplex.h"
#include "SAMRAI/pdat/FaceGeometry.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#include "SAMRAI/tbox/Utilities.h"
#endif

namespace SAMRAI {
namespace math {

PatchFaceDataOpsComplex::PatchFaceDataOpsComplex()
{
}

PatchFaceDataOpsComplex::~PatchFaceDataOpsComplex()
{
}

/*
 *************************************************************************
 *                                                                       *
 * General operations for complex face-centered patch data.              *
 *                                                                       *
 *************************************************************************
 */

void PatchFaceDataOpsComplex::swapData(
   tbox::Pointer<hier::Patch> patch,
   const int data1_id,
   const int data2_id) const
{
   TBOX_ASSERT(!patch.isNull());

   tbox::Pointer<pdat::FaceData<dcomplex> > d1 = patch->getPatchData(data1_id);
   tbox::Pointer<pdat::FaceData<dcomplex> > d2 = patch->getPatchData(data2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d1.isNull() && !d2.isNull());
   TBOX_ASSERT(d1->getDepth() && d2->getDepth());
   TBOX_ASSERT(d1->getBox().isSpatiallyEqual(d2->getBox()));
   TBOX_ASSERT(d1->getGhostBox().isSpatiallyEqual(d2->getGhostBox()));
#endif
   patch->setPatchData(data1_id, d2);
   patch->setPatchData(data2_id, d1);
}

void PatchFaceDataOpsComplex::printData(
   const tbox::Pointer<pdat::FaceData<dcomplex> >& data,
   const hier::Box& box,
   std::ostream& s) const
{
   TBOX_ASSERT(!data.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   s << "Data box = " << box << std::endl;
   data->print(box, s);
   s << "\n";
}

void PatchFaceDataOpsComplex::copyData(
   tbox::Pointer<pdat::FaceData<dcomplex> >& dst,
   const tbox::Pointer<pdat::FaceData<dcomplex> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst, *src, box);

   const tbox::Dimension& dim(box.getDim());

   for (int d = 0; d < dim.getValue(); d++) {
      const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
      (dst->getArrayData(d)).copy(src->getArrayData(d), face_box);
   }
}

void PatchFaceDataOpsComplex::setToScalar(
   tbox::Pointer<pdat::FaceData<dcomplex> >& dst,
   const dcomplex& alpha,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS2(*dst, box);

   dst->fillAll(alpha, box);
}

}
}
#endif
