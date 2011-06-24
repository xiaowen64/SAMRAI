/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Constant refine operator for cell-centered double data on
 *                a Moving mesh. 
 *
 ************************************************************************/

#ifndef included_geom_SkeletonRefine_C
#define included_geom_SkeletonRefine_C

#include "SAMRAI/geom/SkeletonRefine.h"

#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace geom {

SkeletonRefine::SkeletonRefine(
   const tbox::Dimension& dim):
   hier::RefineOperator(dim, "SKELETON_REFINE")
{
}

SkeletonRefine::~SkeletonRefine()
{
}

bool SkeletonRefine::findRefineOperator(
   const tbox::Pointer<hier::Variable>& var,
   const std::string& op_name) const
{
   NULL_USE(var);
   if (op_name == getOperatorName()) {
      return true;
   } else {
      return false;
   }
}

int SkeletonRefine::getOperatorPriority() const
{
   return 0;
}

hier::IntVector
SkeletonRefine::getStencilWidth() const {
   return hier::IntVector::getZero(getDim());
}

void SkeletonRefine::refine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const int dst_component,
   const int src_component,
   const hier::BoxOverlap& fine_overlap,
   const hier::IntVector& ratio) const
{
   //no operation for the empty refine operator
   NULL_USE(fine);
   NULL_USE(coarse);
   NULL_USE(dst_component);
   NULL_USE(src_component);
   NULL_USE(fine_overlap);
   NULL_USE(ratio);
}

}
}
#endif
