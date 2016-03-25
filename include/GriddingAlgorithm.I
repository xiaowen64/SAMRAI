//
// File:        GriddingAlgorithm.I
// Package:     SAMRAI mesh
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: AMR hierarchy generation and regridding routines.
//
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
    namespace mesh {

template<int DIM> inline
bool
GriddingAlgorithm<DIM>::errorEstimationUsesTimeIntegration() 
const
{
   return(d_tag_init_strategy->usesTimeIntegration());
}

template<int DIM> inline
int
GriddingAlgorithm<DIM>::getErrorCoarsenRatio() 
const
{
   return(d_tag_init_strategy->getErrorCoarsenRatio());
}

template<int DIM> inline 
bool 
GriddingAlgorithm<DIM>::levelCanBeRefined(const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(level_number >= 0);
#endif
   return(level_number < d_max_levels-1);
}

template<int DIM> inline
tbox::Pointer< TagAndInitializeStrategy<DIM> > 
GriddingAlgorithm<DIM>::getTagAndInitializeStrategy() const
{
   return(d_tag_init_strategy);
}

template<int DIM> inline
tbox::Pointer< LoadBalanceStrategy<DIM> >
GriddingAlgorithm<DIM>::getLoadBalanceStrategy() const
{
   return(d_load_balancer);
}

template<int DIM> inline 
int
GriddingAlgorithm<DIM>::getMaxLevels() const
{
   return(d_max_levels);
}

template<int DIM> inline 
const hier::IntVector<DIM>&
GriddingAlgorithm<DIM>::getRatioToCoarserLevel(
   const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( (level_number >= 0) && (level_number < d_max_levels) );
#endif
   return(d_ratio_to_coarser[level_number]);
}

template<int DIM> inline
double
GriddingAlgorithm<DIM>::getEfficiencyTolerance(
   const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( (level_number >= 0) && (level_number < d_max_levels - 1) );
#endif
   int size = d_efficiency_tolerance.getSize();
   return( (level_number < size)
           ? d_efficiency_tolerance[level_number]
           : d_efficiency_tolerance[size-1] );
}

template<int DIM> inline
double
GriddingAlgorithm<DIM>::getCombineEfficiency(
   const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( (level_number >= 0) && (level_number < d_max_levels - 1) );
#endif
   int size = d_combine_efficiency.getSize(); 
   return( (level_number < size)
           ? d_combine_efficiency[level_number]
           : d_combine_efficiency[size-1] );
}

template<int DIM> inline
int
GriddingAlgorithm<DIM>::getProperNestingBuffer(
   const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( (level_number >= 0) && (level_number < d_max_levels - 1) );
#endif
   int size = d_proper_nesting_buffer.getSize();
   return( (level_number < size)
           ? d_proper_nesting_buffer[level_number]
           : d_proper_nesting_buffer[size-1] );
}

template<int DIM> inline
const hier::IntVector<DIM>&
GriddingAlgorithm<DIM>::getSmallestPatchSize(
   const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( (level_number >= 0) && (level_number < d_max_levels) );
#endif
   int size = d_smallest_patch_size.getSize(); 
   return( (level_number < size)
           ? d_smallest_patch_size[level_number]
           : d_smallest_patch_size[size-1] );
}

template<int DIM> inline
const hier::IntVector<DIM>&
GriddingAlgorithm<DIM>::getLargestPatchSize(
   const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( (level_number >= 0) && (level_number < d_max_levels) );
#endif
   int size = d_largest_patch_size.getSize();
   return( (level_number < size)
           ? d_largest_patch_size[level_number]
           : d_largest_patch_size[size-1] );
}



#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
