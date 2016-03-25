//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/algorithm/time_refinement/TimeRefinementIntegrator.I $
// Package:     SAMRAI algorithms
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1917 $
// Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description: Time integration manager for AMR with local time stepping.
//


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
   TBOX_ASSERT( (level_number >= 0) &&
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
   TBOX_ASSERT( (level_number >= 0) &&
           (level_number <= d_patch_hierarchy->getFinestLevelNumber()) );
#endif
   return (d_dt_actual_level[level_number]);
}
 
template<int DIM> inline 
double TimeRefinementIntegrator<DIM>::getLevelDtMax(
   const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( (level_number >= 0) &&
           (level_number <= d_patch_hierarchy->getFinestLevelNumber()) );
#endif
   return (d_dt_max_level[level_number]);
}
 
template<int DIM> inline 
double TimeRefinementIntegrator<DIM>::getLevelSimTime(
   const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( (level_number >= 0) &&
           (level_number <= d_patch_hierarchy->getFinestLevelNumber()) );
#endif
   return (d_level_sim_time[level_number]);
}
 
template<int DIM> inline 
int TimeRefinementIntegrator<DIM>::getLevelStep(
   const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( (level_number >= 0) &&
           (level_number <= d_patch_hierarchy->getFinestLevelNumber()) );
#endif
   return (d_step_level[level_number]);
}
 
template<int DIM> inline 
int TimeRefinementIntegrator<DIM>::getLevelMaxSteps(
   const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( (level_number >= 0) &&
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
   TBOX_ASSERT( (level_number >= 0) &&
           (level_number <= d_patch_hierarchy->getFinestLevelNumber()) );
#endif
   return ( d_step_level[level_number] <= 0 );
}
 
template<int DIM> inline 
bool TimeRefinementIntegrator<DIM>::lastLevelStep(
   const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( (level_number >= 0) &&
           (level_number <= d_patch_hierarchy->getFinestLevelNumber()) );
#endif
   return ( d_step_level[level_number] >= d_max_steps_level[level_number] );
}

template<int DIM> inline
void TimeRefinementIntegrator<DIM>::setRegridInterval(
   const int regrid_interval)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_use_refined_timestepping);
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
