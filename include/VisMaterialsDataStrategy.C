//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-1/source/apputils/plotting/VisMaterialsDataStrategy.C $
// Package:     SAMRAI application utilities
// Copyright:   (c) 1997-2003 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1848 $
// Modified:    $LastChangedDate: 2008-01-11 16:26:13 -0800 (Fri, 11 Jan 2008) $
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

template<int DIM>
int VisMaterialsDataStrategy<DIM>::packMaterialFractionsIntoDoubleBuffer(
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

                                                                                
template<int DIM>
int VisMaterialsDataStrategy<DIM>::packMaterialFractionsIntoSparseBuffers(
       int* mat_list,
       std::vector<int>& mix_zones,
       std::vector<int>& mix_mat,
       std::vector<double>& vol_fracs,
       std::vector<int>& next_mat,
       const hier::Patch<NDIM>& patch,
       const hier::Box<NDIM>& region) const
{
   NULL_USE(mat_list);
   NULL_USE(mix_zones);
   NULL_USE(mix_mat);
   NULL_USE(vol_fracs);
   NULL_USE(next_mat);
   NULL_USE(patch);
   NULL_USE(region);
   TBOX_ERROR ("VisMaterialsDataStrategy<DIM>::"
               << "packSparseMaterialFractionsIntoDoubleBuffer()"
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
