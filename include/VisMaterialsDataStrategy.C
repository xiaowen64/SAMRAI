//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/apputils/plotting/VisMaterialsDataStrategy.C $
// Package:     SAMRAI application utilities
// Copyright:   (c) 1997-2003 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Interface for writing material related data to VisIt 
//              dump file
//

#ifndef included_appu_VisMaterialsDataStrategy_C
#define included_appu_VisMaterialsDataStrategy_C

#include "VisMaterialsDataStrategy.h"

#include "tbox/Utilities.h"

namespace SAMRAI {
    namespace appu {

template<int DIM>  VisMaterialsDataStrategy<DIM>::VisMaterialsDataStrategy()
{
}

template<int DIM>  VisMaterialsDataStrategy<DIM>::~VisMaterialsDataStrategy()
{
}

template<int DIM> int VisMaterialsDataStrategy<DIM>::packMaterialFractionsIntoDoubleBuffer(
      double *buffer,
      const hier::Patch<DIM>& patch,
      const hier::Box<DIM>& region,
      const std::string& material_name) const
{
   NULL_USE(buffer);
   NULL_USE(patch);
   NULL_USE(region);
   NULL_USE(material_name);
   TBOX_ERROR ("VisMaterialsDataStrategy<DIM>::"
               << "packMaterialFractionsIntoDoubleBuffer()"
               << "\nNo class supplies a concrete implementation for "
               << "\nthis method.  The default abstract method (which "
               << "\ndoes nothing) is executed" << std::endl);
   return 0;
}


template<int DIM> int VisMaterialsDataStrategy<DIM>::packSpeciesFractionsIntoDoubleBuffer(
      double *buffer,
      const hier::Patch<DIM>& patch,
      const hier::Box<DIM>& region,
      const std::string& material_name,
      const std::string& species_name) const
{
   NULL_USE(buffer);
   NULL_USE(patch);
   NULL_USE(region);
   NULL_USE(material_name);
   NULL_USE(species_name);
   TBOX_ERROR ("VisMaterialsDataStrategy<DIM>::"
               << "packSpeciesFractionsIntoDoubleBuffer()"
               << "\nNo class supplies a concrete implementation for "
               << "\nthis method.  The default abstract method (which "
               << "\ndoes nothing) is executed" << std::endl);
   return 0;
}



}
}

#endif
