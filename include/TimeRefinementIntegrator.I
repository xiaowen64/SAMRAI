//
// File:        TimeRefinementIntegrator.I
// Package:     SAMRAI algorithms
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Time integration manager for AMR with local time stepping.
//

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace algs {

template<int DIM> inline 
double TimeRefinementIntegrator<DIM>::getIntegratorTime() const
{
   return (d_integrator_time);
}
 
template<int DIM> inline 
double TimeRefinementIntegrator<DIM>::getStartTime() const
{
   return (d_start_time);
}
 
template<int DIM> inline 
double TimeRefinementIntegrator<DIM>::getEndTime() const
{
   return (d_end_time);
}
 
template<int DIM> inline 
int TimeRefinementIntegrator<DIM>::getIntegratorStep() const
{
   return (d_integrator_step);
}
 
template<int DIM> inline 
int TimeRefinementIntegrator<DIM>::getMaxIntegratorSteps() const
{
   return (d_max_integrator_steps);
}
 
template<int DIM> inline 
bool TimeRefinementIntegrator<DIM>::stepsRemaining(
   const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( (level_number >= 0) &&
           (level_number <= d_patch_hierarchy->getFinestLevelNumber()) );
#endif
   return (d_step_level[level_number] < d_max_steps_level[level_number]);
}
 
template<int DIM> inline 
bool TimeRefinementIntegrator<DIM>::stepsRemaining() const
{
   return (d_integrator_step < d_max_integrator_steps);
}
 
template<int DIM> inline 
double TimeRefinementIntegrator<DIM>::getLevelDtActual(
   const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( (level_number >= 0) &&
           (level_number <= d_patch_hierarchy->getFinestLevelNumber()) );
#endif
   return (d_dt_actual_level[level_number]);
}
 
template<int DIM> inline 
double TimeRefinementIntegrator<DIM>::getLevelDtMax(
   const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( (level_number >= 0) &&
           (level_number <= d_patch_hierarchy->getFinestLevelNumber()) );
#endif
   return (d_dt_max_level[level_number]);
}
 
template<int DIM> inline 
double TimeRefinementIntegrator<DIM>::getLevelSimTime(
   const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( (level_number >= 0) &&
           (level_number <= d_patch_hierarchy->getFinestLevelNumber()) );
#endif
   return (d_level_sim_time[level_number]);
}
 
template<int DIM> inline 
int TimeRefinementIntegrator<DIM>::getLevelStep(
   const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( (level_number >= 0) &&
           (level_number <= d_patch_hierarchy->getFinestLevelNumber()) );
#endif
   return (d_step_level[level_number]);
}
 
template<int DIM> inline 
int TimeRefinementIntegrator<DIM>::getLevelMaxSteps(
   const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( (level_number >= 0) &&
           (level_number <= d_patch_hierarchy->getFinestLevelNumber()) );
#endif
   return (d_max_steps_level[level_number]);
}

template<int DIM> inline 
const tbox::Pointer< hier::BasePatchHierarchy<DIM> >
TimeRefinementIntegrator<DIM>::getPatchHierarchy() const
{
   return (d_patch_hierarchy);
}
 
template<int DIM> inline
tbox::Pointer<TimeRefinementLevelStrategy<DIM> >
TimeRefinementIntegrator<DIM>::getLevelIntegrator() const
{
   return (d_refine_level_integrator);
}

template<int DIM> inline 
tbox::Pointer< mesh::BaseGriddingAlgorithm<DIM> >
TimeRefinementIntegrator<DIM>::getGriddingAlgorithm() const
{
   return (d_gridding_algorithm);
}

template<int DIM> inline 
bool TimeRefinementIntegrator<DIM>::firstLevelStep(
   const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( (level_number >= 0) &&
           (level_number <= d_patch_hierarchy->getFinestLevelNumber()) );
#endif
   return ( d_step_level[level_number] <= 0 );
}
 
template<int DIM> inline 
bool TimeRefinementIntegrator<DIM>::lastLevelStep(
   const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( (level_number >= 0) &&
           (level_number <= d_patch_hierarchy->getFinestLevelNumber()) );
#endif
   return ( d_step_level[level_number] >= d_max_steps_level[level_number] );
}

template<int DIM> inline
void TimeRefinementIntegrator<DIM>::setRegridInterval(
   const int regrid_interval)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d_use_refined_timestepping);
#endif

   for (int i = 0; i < d_regrid_interval.getSize(); i++) {
      d_regrid_interval[i] = regrid_interval;
   }
}


#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
