//
// File:        HyperbolicPatchStrategy.C
// Package:     SAMRAI algorithms
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Interface to patch routines for hyperbolic integration scheme.
//
 
#ifndef included_algs_HyperbolicPatchStrategy_C
#define included_algs_HyperbolicPatchStrategy_C

#include "HyperbolicPatchStrategy.h"

#include "tbox/Utilities.h"

namespace SAMRAI {
    namespace algs {

template<int DIM>  HyperbolicPatchStrategy<DIM>::HyperbolicPatchStrategy()
{
   d_data_context.setNull();
}
 
template<int DIM>  HyperbolicPatchStrategy<DIM>::~HyperbolicPatchStrategy()
{
}

/*
*************************************************************************
*                                                                       *
* Default virtual function implementations.                             * 
*                                                                       *
*************************************************************************
*/

template<int DIM> void HyperbolicPatchStrategy<DIM>::tagGradientDetectorCells(
   hier::Patch<DIM>& patch, 
   const double regrid_time, 
   const bool initial_error, 
   const int tag_index, 
   const bool uses_richardson_extrapolation_too)
{
   NULL_USE(patch);
   NULL_USE(regrid_time);
   NULL_USE(initial_error);
   NULL_USE(tag_index);
   NULL_USE(uses_richardson_extrapolation_too);
   TBOX_WARNING("HyperbolicPatchStrategy::tagGradientDetectorCells()"
                << "\nNo class supplies a concrete implementation for "
                << "\nthis method.  The default abstract method (which "
                << "\ndoes no cell tagging) is executed" << endl);
}


template<int DIM> void HyperbolicPatchStrategy<DIM>::tagRichardsonExtrapolationCells(
   hier::Patch<DIM>& patch,
   const int error_level_number,
   const tbox::Pointer<hier::VariableContext> coarsened_fine,
   const tbox::Pointer<hier::VariableContext> advanced_coarse,
   const double regrid_time,
   const double deltat, 
   const int error_coarsen_ratio,
   const bool initial_error, 
   const int tag_index,
   const bool uses_gradient_detector_too)
{
   NULL_USE(patch);
   NULL_USE(error_level_number);
   NULL_USE(coarsened_fine);
   NULL_USE(advanced_coarse);
   NULL_USE(regrid_time);
   NULL_USE(deltat);
   NULL_USE(error_coarsen_ratio);
   NULL_USE(initial_error);
   NULL_USE(tag_index);
   NULL_USE(uses_gradient_detector_too);
   TBOX_WARNING("HyperbolicPatchStrategy::tagRichardsonExtrapolationCells()"
                << "\nNo class supplies a concrete implementation for "
                << "\nthis method.  The default abstract method (which "
                << "\ndoes no cell tagging) is executed" << endl);
}
}
}
#endif
