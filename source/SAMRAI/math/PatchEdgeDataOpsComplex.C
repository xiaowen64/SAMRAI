/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Operations for complex edge-centered patch data. 
 *
 ************************************************************************/

#ifndef included_math_PatchEdgeDataOpsComplex_C
#define included_math_PatchEdgeDataOpsComplex_C

#include "SAMRAI/math/PatchEdgeDataOpsComplex.h"
#include "SAMRAI/pdat/EdgeGeometry.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#include "SAMRAI/tbox/Utilities.h"
#endif

namespace SAMRAI {
namespace math {

PatchEdgeDataOpsComplex::PatchEdgeDataOpsComplex()
{
}

PatchEdgeDataOpsComplex::~PatchEdgeDataOpsComplex()
{
}

/*
 *************************************************************************
 *                                                                       *
 * General operations for complex edge-centered patch data.              *
 *                                                                       *
 *************************************************************************
 */

void PatchEdgeDataOpsComplex::swapData(
   tbox::Pointer<hier::Patch> patch,
   const int data1_id,
   const int data2_id) const
{
   TBOX_ASSERT(!patch.isNull());

   tbox::Pointer<pdat::EdgeData<dcomplex> > d1 = patch->getPatchData(data1_id);
   tbox::Pointer<pdat::EdgeData<dcomplex> > d2 = patch->getPatchData(data2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d1.isNull() && !d2.isNull());
   TBOX_ASSERT(d1->getDepth() && d2->getDepth());
   TBOX_ASSERT(d1->getBox() == d2->getBox());
   TBOX_ASSERT(d1->getGhostBox() == d2->getGhostBox());
#endif
   patch->setPatchData(data1_id, d2);
   patch->setPatchData(data2_id, d1);
}

void PatchEdgeDataOpsComplex::printData(
   const tbox::Pointer<pdat::EdgeData<dcomplex> >& data,
   const hier::Box& box,
   std::ostream& s) const
{
   TBOX_ASSERT(!data.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   s << "Data box = " << box << std::endl;
   data->print(box, s);
   s << "\n";
}

void PatchEdgeDataOpsComplex::copyData(
   tbox::Pointer<pdat::EdgeData<dcomplex> >& dst,
   const tbox::Pointer<pdat::EdgeData<dcomplex> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst, *src, box);

   const tbox::Dimension& dim(box.getDim());

   for (int d = 0; d < dim.getValue(); d++) {
      const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
      (dst->getArrayData(d)).copy(src->getArrayData(d), edge_box);
   }
}

void PatchEdgeDataOpsComplex::setToScalar(
   tbox::Pointer<pdat::EdgeData<dcomplex> >& dst,
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
