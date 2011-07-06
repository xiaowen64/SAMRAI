/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   AMR hierarchy generation and regridding routines. 
 *
 ************************************************************************/

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <cassert>
#endif
#endif

namespace SAMRAI {
namespace mesh {

SAMRAI_INLINE_KEYWORD
bool
GriddingAlgorithm::errorEstimationUsesTimeIntegration()
const
{
   return d_tag_init_strategy->usesTimeIntegration();
}

SAMRAI_INLINE_KEYWORD
tbox::Pointer<TagAndInitializeStrategy>
GriddingAlgorithm::getTagAndInitializeStrategy() const
{
   return d_tag_init_strategy;
}

SAMRAI_INLINE_KEYWORD
tbox::Pointer<LoadBalanceStrategy>
GriddingAlgorithm::getLoadBalanceStrategy() const
{
   return d_load_balancer;
}

SAMRAI_INLINE_KEYWORD
tbox::Pointer<LoadBalanceStrategy>
GriddingAlgorithm::getLoadBalanceStrategy0() const
{
   return d_load_balancer0;
}

SAMRAI_INLINE_KEYWORD
double
GriddingAlgorithm::getEfficiencyTolerance(
   const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((level_number >= 0) && (level_number < d_hierarchy->getMaxNumberOfLevels()));
#endif
   int size = d_efficiency_tolerance.getSize();
   return (level_number < size)
          ? d_efficiency_tolerance[level_number]
          : d_efficiency_tolerance[size - 1];
}

SAMRAI_INLINE_KEYWORD
double
GriddingAlgorithm::getCombineEfficiency(
   const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((level_number >= 0) && (level_number < d_hierarchy->getMaxNumberOfLevels()));
#endif
   int size = d_combine_efficiency.getSize();
   return (level_number < size)
          ? d_combine_efficiency[level_number]
          : d_combine_efficiency[size - 1];
}

}
}
