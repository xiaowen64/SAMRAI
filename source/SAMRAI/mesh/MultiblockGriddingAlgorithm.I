/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   AMR hierarchy generation and regridding routines. 
 *
 ************************************************************************/

namespace SAMRAI {
namespace mesh {

SAMRAI_INLINE_KEYWORD
bool
MultiblockGriddingAlgorithm::errorEstimationUsesTimeIntegration()
const
{
   return d_tag_init_strategy->usesTimeIntegration();
}

SAMRAI_INLINE_KEYWORD
tbox::Pointer<TagAndInitializeStrategy>
MultiblockGriddingAlgorithm::getTagAndInitializeStrategy() const
{
   return d_tag_init_strategy;
}

SAMRAI_INLINE_KEYWORD
tbox::Pointer<LoadBalanceStrategy>
MultiblockGriddingAlgorithm::getLoadBalanceStrategy() const
{
   return d_load_balancer;
}

SAMRAI_INLINE_KEYWORD
const hier::IntVector&
MultiblockGriddingAlgorithm::getRatioToCoarserLevel(
   const int level_number) const
{
   TBOX_ASSERT((level_number >= 0) && (level_number < d_mb_hierarchy->getMaxNumberOfLevels()));
   return d_mb_hierarchy->getRatioToCoarserLevel(level_number);
}

SAMRAI_INLINE_KEYWORD
double
MultiblockGriddingAlgorithm::getEfficiencyTolerance(
   const int level_number) const
{
   TBOX_ASSERT((level_number >= 0) && (level_number < d_mb_hierarchy->getMaxNumberOfLevels()));
   int size = d_efficiency_tolerance.getSize();
   return (level_number < size)
          ? d_efficiency_tolerance[level_number]
          : d_efficiency_tolerance[size - 1];
}

SAMRAI_INLINE_KEYWORD
double
MultiblockGriddingAlgorithm::getCombineEfficiency(
   const int level_number) const
{
   TBOX_ASSERT((level_number >= 0) && (level_number < d_mb_hierarchy->getMaxNumberOfLevels()));
   int size = d_combine_efficiency.getSize();
   return (level_number < size)
          ? d_combine_efficiency[level_number]
          : d_combine_efficiency[size - 1];
}

/*!
 * Return proper nesting buffer width for level.
 *
 * Level number level_number+1 must nest inside level level_number
 * by the width returned (except where the levels touch the domain
 * boundary).
 */
SAMRAI_INLINE_KEYWORD
int
MultiblockGriddingAlgorithm::getProperNestingBuffer(
   const int level_number) const
{
   TBOX_ASSERT((level_number >= 0) && (level_number < d_mb_hierarchy->getMaxNumberOfLevels()));
   return d_mb_hierarchy->getProperNestingBuffer(level_number);
}

SAMRAI_INLINE_KEYWORD
const hier::IntVector&
MultiblockGriddingAlgorithm::getSmallestPatchSize(
   const int level_number) const
{
   TBOX_ASSERT((level_number >= 0) && (level_number < d_mb_hierarchy->getMaxNumberOfLevels()));
   return d_mb_hierarchy->getSmallestPatchSize(level_number);
}

SAMRAI_INLINE_KEYWORD
const hier::IntVector&
MultiblockGriddingAlgorithm::getLargestPatchSize(
   const int level_number) const
{
   TBOX_ASSERT((level_number >= 0) && (level_number < d_mb_hierarchy->getMaxNumberOfLevels()));
   return d_mb_hierarchy->getLargestPatchSize(level_number);
}

SAMRAI_INLINE_KEYWORD
const std::string& MultiblockGriddingAlgorithm::getObjectName() const
{
   return d_object_name;
}

}
}
