/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Weighted averaging operator for cell-centered double data on
 *                a Moving mesh.
 *
 ************************************************************************/

#ifndef included_geom_SkeletonCoarsen_C
#define included_geom_SkeletonCoarsen_C

#include "SAMRAI/geom/SkeletonCoarsen.h"

#include "SAMRAI/tbox/Utilities.h"

#include <float.h>
#include <math.h>

namespace SAMRAI {
namespace geom {

// using namespace std;

SkeletonCoarsen::SkeletonCoarsen(
   const tbox::Dimension& dim):
   hier::CoarsenOperator(dim, "SKELETON_COARSEN")
{
}

SkeletonCoarsen::~SkeletonCoarsen()
{
}

bool SkeletonCoarsen::findCoarsenOperator(
   const boost::shared_ptr<hier::Variable>& var,
   const std::string& op_name) const
{
   NULL_USE(var);
   if (op_name == getOperatorName()) {
      return true;
   } else {
      return false;
   }
}

int SkeletonCoarsen::getOperatorPriority() const
{
   return 0;
}

hier::IntVector SkeletonCoarsen::getStencilWidth() const {
   return hier::IntVector::getZero(getDim());
}

void SkeletonCoarsen::coarsen(
   hier::Patch& coarse,
   const hier::Patch& fine,
   const int dst_component,
   const int src_component,
   const hier::Box& coarse_box,
   const hier::IntVector& ratio) const
{
   //no operation for skeleton coarsen operator
   NULL_USE(coarse);
   NULL_USE(fine);
   NULL_USE(dst_component);
   NULL_USE(src_component);
   NULL_USE(coarse_box);
   NULL_USE(ratio);
}

}
}
#endif
