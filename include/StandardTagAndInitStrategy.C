//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/mesh/gridding/StandardTagAndInitStrategy.C $
// Package:     SAMRAI mesh
// Copyright:   (c) 1997-2000 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Strategy interface for Richardson extrapolation error detection.
//

#ifndef included_mesh_StandardTagAndInitStrategy_C
#define included_mesh_StandardTagAndInitStrategy_C

#include "StandardTagAndInitStrategy.h"

#include "tbox/Utilities.h"

namespace SAMRAI {
    namespace mesh {

template<int DIM> 
StandardTagAndInitStrategy<DIM>::StandardTagAndInitStrategy() 
{
}

template<int DIM> 
StandardTagAndInitStrategy<DIM>::~StandardTagAndInitStrategy()
{
}

/*
*************************************************************************
*                                                                       *
* Default virtual function implementations.                             * 
*                                                                       *
*************************************************************************
*/

template<int DIM> void StandardTagAndInitStrategy<DIM>::applyGradientDetector(
   const tbox::Pointer< hier::BasePatchHierarchy<DIM> > hierarchy,
   const int level_number, 
   const double error_data_time, 
   const int tag_index, 
   const bool initial_time,
   const bool uses_richardson_extrapolation_too)
{
   NULL_USE(hierarchy);
   NULL_USE(level_number);
   NULL_USE(error_data_time);
   NULL_USE(tag_index);
   NULL_USE(initial_time);
   NULL_USE(uses_richardson_extrapolation_too);
   TBOX_WARNING("StandardTagAndInitStrategy<DIM>::applyGradientDetector()"
                << "\nNo class supplies a concrete implementation for "
                << "\nthis method.  The default abstract method (which "
                << "\ndoes no cell tagging) is executed" << std::endl);
}

template<int DIM> void StandardTagAndInitStrategy<DIM>::coarsenDataForRichardsonExtrapolation(
   const tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
   const int level_number,
   const tbox::Pointer< hier::PatchLevel<DIM> > coarser_level,
   const double coarsen_data_time,
   const bool before_advance)
{
   NULL_USE(hierarchy);
   NULL_USE(level_number);
   NULL_USE(coarser_level);
   NULL_USE(coarsen_data_time);
   NULL_USE(before_advance);
   TBOX_WARNING("StandardTagAndInitStrategy<DIM>::"
                << "coarsenDataForRichardsonExtrapolation()"
                << "\nNo class supplies a concrete implementation for "
                << "\nthis method.  The default abstract method (which "
                << "\ndoes nothing) is executed" << std::endl);
}

template<int DIM> void StandardTagAndInitStrategy<DIM>::applyRichardsonExtrapolation(
   const tbox::Pointer< hier::PatchLevel<DIM> > level, 
   const double error_data_time, 
   const int tag_index, 
   const double deltat, 
   const int error_coarsen_ratio, 
   const bool initial_time,
   const bool uses_gradient_detector_too)
{
   NULL_USE(level);
   NULL_USE(error_data_time);
   NULL_USE(tag_index);
   NULL_USE(deltat);
   NULL_USE(error_coarsen_ratio);
   NULL_USE(initial_time);
   NULL_USE(uses_gradient_detector_too);
   TBOX_WARNING("StandardTagAndInitStrategy<DIM>::"
                << "applyRichardsonExtrapolation()"
                << "\nNo class supplies a concrete implementation for "
                << "\nthis method.  The default abstract method (which "
                << "\ndoes nothing) is executed" << std::endl);
}

template<int DIM> double StandardTagAndInitStrategy<DIM>::getLevelDt(
   const tbox::Pointer< hier::BasePatchLevel<DIM> > level,
   const double dt_time,
   const bool initial_time)
{
   NULL_USE(level);
   NULL_USE(dt_time);
   NULL_USE(initial_time);
   TBOX_WARNING("StandardTagAndInitStrategy<DIM>::getLevelDt()"
                << "\nNo class supplies a concrete implementation for "
                << "\nthis method.  The default abstract method (which "
                << "\nsimply returns 0.) is executed" << std::endl);
   return(0.0);
}

template<int DIM> void StandardTagAndInitStrategy<DIM>::resetTimeDependentData(
   const tbox::Pointer< hier::BasePatchLevel<DIM> > level, 
   const double new_time, 
   const bool can_be_refined)
{
   NULL_USE(level);
   NULL_USE(new_time);
   NULL_USE(can_be_refined);
   TBOX_WARNING("StandardTagAndInitStrategy<DIM>::resetTimeDependentData()"
                << "\nNo class supplies a concrete implementation for "
                << "\nthis method.  The default abstract method (which "
                << "\ndoes nothing) is executed" << std::endl);
}

template<int DIM> double StandardTagAndInitStrategy<DIM>::advanceLevel(
   const tbox::Pointer< hier::BasePatchLevel<DIM> > level, 
   const tbox::Pointer< hier::BasePatchHierarchy<DIM> > hierarchy, 
   const double current_time, 
   const double new_time, 
   const bool first_step, 
   const bool last_step, 
   const bool regrid_advance)
{
   NULL_USE(level);
   NULL_USE(hierarchy);
   NULL_USE(current_time);
   NULL_USE(new_time);
   NULL_USE(first_step);
   NULL_USE(last_step);
   NULL_USE(regrid_advance);
   TBOX_WARNING("StandardTagAndInitStrategy<DIM>::advanceLevel()"
                << "\nNo class supplies a concrete implementation for "
                << "\nthis method.  The default abstract method (which "
                << "\ndoes nothing) is executed" << std::endl);
   return(0.0);
}

template<int DIM> void StandardTagAndInitStrategy<DIM>::resetDataToPreadvanceState(
   const tbox::Pointer< hier::BasePatchLevel<DIM> > level)
{
   NULL_USE(level);
   TBOX_WARNING("StandardTagAndInitStrategy<DIM>::"
                << "resetDataToPreadvanceState()"
                << "\nNo class supplies a concrete implementation for "
                << "\nthis method.  The default abstract method (which "
                << "\ndoes nothing) is executed" << std::endl);
}

}
}
#endif
