/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Operations for complex node-centered patch data.
 *
 ************************************************************************/

#ifndef included_math_PatchNodeDataOpsComplex_C
#define included_math_PatchNodeDataOpsComplex_C

#include "SAMRAI/math/PatchNodeDataOpsComplex.h"
#include "SAMRAI/pdat/NodeGeometry.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#include "SAMRAI/tbox/Utilities.h"
#endif

namespace SAMRAI {
namespace math {

PatchNodeDataOpsComplex::PatchNodeDataOpsComplex()
{
}

PatchNodeDataOpsComplex::~PatchNodeDataOpsComplex()
{
}

/*
 *************************************************************************
 *
 * General operations for complex node-centered patch data.
 *
 *************************************************************************
 */

void PatchNodeDataOpsComplex::swapData(
   tbox::Pointer<hier::Patch> patch,
   const int data1_id,
   const int data2_id) const
{
   TBOX_ASSERT(patch);

   tbox::Pointer<pdat::NodeData<dcomplex> > d1(
      patch->getPatchData(data1_id),
      tbox::__dynamic_cast_tag());
   tbox::Pointer<pdat::NodeData<dcomplex> > d2(
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

void PatchNodeDataOpsComplex::printData(
   const tbox::Pointer<pdat::NodeData<dcomplex> >& data,
   const hier::Box& box,
   std::ostream& s) const
{
   TBOX_ASSERT(data);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   s << "Data box = " << box << std::endl;
   data->print(box, s);
   s << "\n";
}

void PatchNodeDataOpsComplex::copyData(
   tbox::Pointer<pdat::NodeData<dcomplex> >& dst,
   const tbox::Pointer<pdat::NodeData<dcomplex> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst, *src, box);

   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   (dst->getArrayData()).copy(src->getArrayData(), node_box);
}

void PatchNodeDataOpsComplex::setToScalar(
   tbox::Pointer<pdat::NodeData<dcomplex> >& dst,
   const dcomplex& alpha,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*dst, box);

   dst->fillAll(alpha, box);
}

}
}
#endif
