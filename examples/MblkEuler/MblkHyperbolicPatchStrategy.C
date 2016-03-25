//
// File:        MblkHyperbolicPatchStrategy.C
// Package:     SAMRAI algorithms
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1917 $
// Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description: Interface to patch routines for hyperbolic integration scheme.
//

#include "MblkHyperbolicPatchStrategy.h"

#include "tbox/Utilities.h"

using namespace SAMRAI;

MblkHyperbolicPatchStrategy::MblkHyperbolicPatchStrategy()
{
   d_data_context.setNull();
}


MblkHyperbolicPatchStrategy::~MblkHyperbolicPatchStrategy()
{
}

/*
*************************************************************************
*                                                                       *
* Default virtual function implementations.                             * 
*                                                                       *
*************************************************************************
*/

void MblkHyperbolicPatchStrategy::tagGradientDetectorCells(
   hier::Patch<NDIM>& patch, 
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
   TBOX_WARNING("MblkHyperbolicPatchStrategy::tagGradientDetectorCells()"
                << "\nNo class supplies a concrete implementation for "
                << "\nthis method.  The default abstract method (which "
                << "\ndoes no cell tagging) is executed" << endl);
}


void MblkHyperbolicPatchStrategy::tagRichardsonExtrapolationCells(
   hier::Patch<NDIM>& patch,
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
   TBOX_WARNING("MblkHyperbolicPatchStrategy::tagRichardsonExtrapolationCells()"
                << "\nNo class supplies a concrete implementation for "
                << "\nthis method.  The default abstract method (which "
                << "\ndoes no cell tagging) is executed" << endl);
}

