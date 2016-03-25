//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-1/source/apputils/plotting/VisDerivedDataStrategy.C $
// Package:     SAMRAI application utilities
// Copyright:   (c) 1997-2003 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1848 $
// Modified:    $LastChangedDate: 2008-01-11 16:26:13 -0800 (Fri, 11 Jan 2008) $
// Description: Interface for writing user-defined data to either VisIt or
//              Vizamrai dump file
//

#ifndef included_appu_VisDerivedDataStrategy_C
#define included_appu_VisDerivedDataStrategy_C

#include "VisDerivedDataStrategy.h"

#include "tbox/Utilities.h"

namespace SAMRAI {
    namespace appu {

template<int DIM>  VisDerivedDataStrategy<DIM>::VisDerivedDataStrategy()
{
}

template<int DIM>  VisDerivedDataStrategy<DIM>::~VisDerivedDataStrategy()
{
}

template<int DIM>
bool VisDerivedDataStrategy<DIM>::packMixedDerivedDataIntoDoubleBuffer(
   double *buffer,
   std::vector<double>& mixbuffer,
   const hier::Patch<DIM>& patch,
   const hier::Box<DIM>& region,
   const std::string& variable_name,
   int   depth_index) const
{
   NULL_USE(buffer);
   NULL_USE(mixbuffer);
   NULL_USE(patch);
   NULL_USE(region);
   NULL_USE(variable_name);
   NULL_USE(depth_index);
   TBOX_ERROR ("VisDerivedDataStrategy<DIM>::"
               << "packMixedDerivedDataIntoDoubleBuffer()"
               << "\nNo class supplies a concrete implementation for "
               << "\nthis method.  The default abstract method (which "
               << "\ndoes nothing) is executed" << std::endl);
   return 0;
}


}
}

#endif
