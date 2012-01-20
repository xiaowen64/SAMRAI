/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Operations for complex cell-centered patch data.
 *
 ************************************************************************/

#ifndef included_math_PatchCellDataOpsComplex_C
#define included_math_PatchCellDataOpsComplex_C

#include "SAMRAI/math/PatchCellDataOpsComplex.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#include "SAMRAI/tbox/Utilities.h"
#endif

namespace SAMRAI {
namespace math {

PatchCellDataOpsComplex::PatchCellDataOpsComplex()
{
}

PatchCellDataOpsComplex::~PatchCellDataOpsComplex()
{
}

/*
 *************************************************************************
 *
 * General operations for complex cell-centered patch data.
 *
 *************************************************************************
 */

void PatchCellDataOpsComplex::swapData(
   boost::shared_ptr<hier::Patch> patch,
   const int data1_id,
   const int data2_id) const
{
   TBOX_ASSERT(patch);

   boost::shared_ptr<pdat::CellData<dcomplex> > d1(
      patch->getPatchData(data1_id),
      boost::detail::dynamic_cast_tag());
   boost::shared_ptr<pdat::CellData<dcomplex> > d2(
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

void PatchCellDataOpsComplex::printData(
   const boost::shared_ptr<pdat::CellData<dcomplex> >& data,
   const hier::Box& box,
   std::ostream& s) const
{
   TBOX_ASSERT(data);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   s << "Data box = " << box << std::endl;
   data->print(box, s);
   s << "\n";
}

void PatchCellDataOpsComplex::copyData(
   boost::shared_ptr<pdat::CellData<dcomplex> >& dst,
   const boost::shared_ptr<pdat::CellData<dcomplex> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst, *src, box);

   (dst->getArrayData()).copy(src->getArrayData(), box);
}

void PatchCellDataOpsComplex::setToScalar(
   boost::shared_ptr<pdat::CellData<dcomplex> >& dst,
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
