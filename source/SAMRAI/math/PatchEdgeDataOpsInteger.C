/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Operations for integer edge-centered patch data.
 *
 ************************************************************************/

#ifndef included_math_PatchEdgeDataOpsInteger_C
#define included_math_PatchEdgeDataOpsInteger_C

#include "SAMRAI/math/PatchEdgeDataOpsInteger.h"
#include "SAMRAI/pdat/EdgeGeometry.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#include "SAMRAI/tbox/Utilities.h"
#endif

#ifndef SAMRAI_INLINE
#include "SAMRAI/math/PatchEdgeDataOpsInteger.I"
#endif

namespace SAMRAI {
namespace math {

PatchEdgeDataOpsInteger::PatchEdgeDataOpsInteger()
{
}

PatchEdgeDataOpsInteger::~PatchEdgeDataOpsInteger()
{
}

/*
 *************************************************************************
 *
 * General operations for integer edge-centered patch data.
 *
 *************************************************************************
 */

void
PatchEdgeDataOpsInteger::swapData(
   const boost::shared_ptr<hier::Patch>& patch,
   const int data1_id,
   const int data2_id) const
{
   TBOX_ASSERT(patch);

   boost::shared_ptr<pdat::EdgeData<int> > d1(
      patch->getPatchData(data1_id),
      boost::detail::dynamic_cast_tag());
   boost::shared_ptr<pdat::EdgeData<int> > d2(
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

}
}
#endif
