//
// File:        VisMaterialsDataStrategy.C
// Package:     SAMRAI application utilities
// Copyright:   (c) 1997-2003 The Regents of the University of California
// Revision:    $Revision: 47 $
// Modified:    $Date: 2004-12-09 16:08:57 -0800 (Thu, 09 Dec 2004) $
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
      const string& material_name)
{
   NULL_USE(buffer);
   NULL_USE(patch);
   NULL_USE(region);
   NULL_USE(material_name);
   TBOX_ERROR ("VisMaterialsDataStrategy<DIM>::"
               << "packMaterialFractionsIntoDoubleBuffer()"
               << "\nNo class supplies a concrete implementation for "
               << "\nthis method.  The default abstract method (which "
               << "\ndoes nothing) is executed" << endl);
   return 0;
}


template<int DIM> int VisMaterialsDataStrategy<DIM>::packSpeciesFractionsIntoDoubleBuffer(
      double *buffer,
      const hier::Patch<DIM>& patch,
      const hier::Box<DIM>& region,
      const string& material_name,
      const string& species_name)
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
               << "\ndoes nothing) is executed" << endl);
   return 0;
}



}
}

#endif
