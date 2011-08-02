/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Templated operations for real side-centered patch data. 
 *
 ************************************************************************/

#ifndef included_math_PatchSideDataOpsReal_C
#define included_math_PatchSideDataOpsReal_C

#include "SAMRAI/math/PatchSideDataOpsReal.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/pdat/SideGeometry.h"

namespace SAMRAI {
namespace math {

template<class TYPE>
PatchSideDataOpsReal<TYPE>::PatchSideDataOpsReal()
{
}

#if 0
/*
 * This was moved into the header due to what looks like bug in the
 * XLC compiler.
 */
template<class TYPE>
PatchSideDataOpsReal<TYPE>::~PatchSideDataOpsReal()
{
}
#endif

/*
 *************************************************************************
 *                                                                       *
 * The const constructor and assignment operator are not actually used   *
 * but are defined here for compilers that require an implementation for *
 * every declaration.                                                    *
 *                                                                       *
 *************************************************************************
 */

template<class TYPE>
PatchSideDataOpsReal<TYPE>::PatchSideDataOpsReal(
   const PatchSideDataOpsReal<TYPE>& foo)
{
   NULL_USE(foo);
}

template<class TYPE>
void PatchSideDataOpsReal<TYPE>::operator = (
   const PatchSideDataOpsReal<TYPE>& foo)
{
   NULL_USE(foo);
}

/*
 *************************************************************************
 *                                                                       *
 * General templated operations for real side-centered patch data.       *
 *                                                                       *
 *************************************************************************
 */

template<class TYPE>
void PatchSideDataOpsReal<TYPE>::swapData(
   tbox::Pointer<hier::Patch> patch,
   const int data1_id,
   const int data2_id) const
{
   TBOX_ASSERT(!patch.isNull());

   tbox::Pointer<pdat::SideData<TYPE> > d1 = patch->getPatchData(data1_id);
   tbox::Pointer<pdat::SideData<TYPE> > d2 = patch->getPatchData(data2_id);

   TBOX_ASSERT(!d1.isNull() && !d2.isNull());
   TBOX_ASSERT(d1->getDepth() && d2->getDepth());
   TBOX_ASSERT(d1->getBox().isSpatiallyEqual(d2->getBox()));
   TBOX_ASSERT(d1->getDirectionVector() == d2->getDirectionVector());
   TBOX_ASSERT(d1->getGhostBox().isSpatiallyEqual(d2->getGhostBox()));

   patch->setPatchData(data1_id, d2);
   patch->setPatchData(data2_id, d1);
}

template<class TYPE>
void PatchSideDataOpsReal<TYPE>::printData(
   const tbox::Pointer<pdat::SideData<TYPE> >& data,
   const hier::Box& box,
   std::ostream& s) const
{
   TBOX_ASSERT(!data.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   s << "Data box = " << box << std::endl;
   data->print(box, s);
   s << "\n";
}

template<class TYPE>
void PatchSideDataOpsReal<TYPE>::copyData(
   tbox::Pointer<pdat::SideData<TYPE> >& dst,
   const tbox::Pointer<pdat::SideData<TYPE> >& src,
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

template<class TYPE>
void PatchSideDataOpsReal<TYPE>::setToScalar(
   tbox::Pointer<pdat::SideData<TYPE> >& dst,
   const TYPE& alpha,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS2(*dst, box);

   dst->fillAll(alpha, box);
}

}
}
#endif
