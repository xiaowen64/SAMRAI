/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Time integration manager for AMR with local time stepping.
 *
 ************************************************************************/

namespace SAMRAI {
namespace algs {

SAMRAI_INLINE_KEYWORD
double TimeRefinementIntegrator::getIntegratorTime() const
{
   return d_integrator_time;
}

SAMRAI_INLINE_KEYWORD
double TimeRefinementIntegrator::getStartTime() const
{
   return d_start_time;
}

SAMRAI_INLINE_KEYWORD
double TimeRefinementIntegrator::getEndTime() const
{
   return d_end_time;
}

SAMRAI_INLINE_KEYWORD
int TimeRefinementIntegrator::getIntegratorStep() const
{
   return d_integrator_step;
}

SAMRAI_INLINE_KEYWORD
int TimeRefinementIntegrator::getMaxIntegratorSteps() const
{
   return d_max_integrator_steps;
}

SAMRAI_INLINE_KEYWORD
bool TimeRefinementIntegrator::stepsRemaining(
   const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((level_number >= 0) &&
      (level_number <= d_patch_hierarchy->getFinestLevelNumber()));
#endif
   return d_step_level[level_number] < d_max_steps_level[level_number];
}

SAMRAI_INLINE_KEYWORD
bool TimeRefinementIntegrator::stepsRemaining() const
{
   return d_integrator_step < d_max_integrator_steps;
}

SAMRAI_INLINE_KEYWORD
double TimeRefinementIntegrator::getLevelDtActual(
   const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((level_number >= 0) &&
      (level_number <= d_patch_hierarchy->getFinestLevelNumber()));
#endif
   return d_dt_actual_level[level_number];
}

SAMRAI_INLINE_KEYWORD
double TimeRefinementIntegrator::getLevelDtMax(
   const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((level_number >= 0) &&
      (level_number <= d_patch_hierarchy->getFinestLevelNumber()));
#endif
   return d_dt_max_level[level_number];
}

SAMRAI_INLINE_KEYWORD
double TimeRefinementIntegrator::getLevelSimTime(
   const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((level_number >= 0) &&
      (level_number <= d_patch_hierarchy->getFinestLevelNumber()));
#endif
   return d_level_sim_time[level_number];
}

SAMRAI_INLINE_KEYWORD
int TimeRefinementIntegrator::getLevelStep(
   const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((level_number >= 0) &&
      (level_number <= d_patch_hierarchy->getFinestLevelNumber()));
#endif
   return d_step_level[level_number];
}

SAMRAI_INLINE_KEYWORD
int TimeRefinementIntegrator::getLevelMaxSteps(
   const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((level_number >= 0) &&
      (level_number <= d_patch_hierarchy->getFinestLevelNumber()));
#endif
   return d_max_steps_level[level_number];
}

SAMRAI_INLINE_KEYWORD
const boost::shared_ptr<hier::PatchHierarchy>
TimeRefinementIntegrator::getPatchHierarchy() const
{
   return d_patch_hierarchy;
}

SAMRAI_INLINE_KEYWORD
boost::shared_ptr<TimeRefinementLevelStrategy>
TimeRefinementIntegrator::getLevelIntegrator() const
{
   return d_refine_level_integrator;
}

SAMRAI_INLINE_KEYWORD
boost::shared_ptr<mesh::GriddingAlgorithmStrategy>
TimeRefinementIntegrator::getGriddingAlgorithm() const
{
   return d_gridding_algorithm;
}

SAMRAI_INLINE_KEYWORD
bool TimeRefinementIntegrator::firstLevelStep(
   const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((level_number >= 0) &&
      (level_number <= d_patch_hierarchy->getFinestLevelNumber()));
#endif
   return d_step_level[level_number] <= 0;
}

SAMRAI_INLINE_KEYWORD
bool TimeRefinementIntegrator::lastLevelStep(
   const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((level_number >= 0) &&
      (level_number <= d_patch_hierarchy->getFinestLevelNumber()));
#endif
   return d_step_level[level_number] >= d_max_steps_level[level_number];
}

SAMRAI_INLINE_KEYWORD
void TimeRefinementIntegrator::setRegridInterval(
   const int regrid_interval)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_use_refined_timestepping);
#endif

   for (int i = 0; i < d_regrid_interval.getSize(); i++) {
      d_regrid_interval[i] = regrid_interval;
   }
}

SAMRAI_INLINE_KEYWORD
const std::string& TimeRefinementIntegrator::getObjectName() const
{
   return d_object_name;
}

}
}
