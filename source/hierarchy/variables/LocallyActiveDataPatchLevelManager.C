//
// File:        LocallyActiveDataPatchLevelManager.C
// Package:     SAMRAI hierarchy 
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 481 $
// Modified:    $Date: 2005-07-21 13:50:43 -0700 (Thu, 21 Jul 2005) $
// Description: Class for managing locally-active data on a single patch level.
//

#ifndef included_hier_LocallyActiveDataPatchLevelManager_C
#define included_hier_LocallyActiveDataPatchLevelManager_C

#include "LocallyActiveDataPatchLevelManager.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

#include "LocallyActiveVariableDatabase.h"

#ifdef DEBUG_NO_INLINE
#include "LocallyActiveDataPatchLevelManager.I"
#endif

namespace SAMRAI {
    namespace hier {

/*
*************************************************************************
*                                                                       *
* Constructors.                                                         *
*                                                                       *
*************************************************************************
*/

template<int DIM>
LocallyActiveDataPatchLevelManager<DIM>::LocallyActiveDataPatchLevelManager()
{
   d_number_patches = -1;
}

template<int DIM>
LocallyActiveDataPatchLevelManager<DIM>::LocallyActiveDataPatchLevelManager(
   const hier::PatchLevel<DIM>& level)
{
   reset(level);
}

template<int DIM>
LocallyActiveDataPatchLevelManager<DIM>::LocallyActiveDataPatchLevelManager(
   const tbox::Pointer< hier::PatchLevel<DIM> > level)
{
   reset(level);
}

template<int DIM>
LocallyActiveDataPatchLevelManager<DIM>::~LocallyActiveDataPatchLevelManager()
{
   d_number_patches = 0;
   d_active_patch_data.resizeArray(0); 
}

/*
*************************************************************************
*                                                                       *
* Initialization functions.                                             *
*                                                                       *
*************************************************************************
*/

template<int DIM>
void LocallyActiveDataPatchLevelManager<DIM>::reset(
   const hier::PatchLevel<DIM>& level)
{
   d_patch_level = (hier::PatchLevel<DIM>*) &level; 
   d_number_patches = d_patch_level->getNumberOfPatches();
   d_active_patch_data.resizeArray(0);
   d_active_patch_data.resizeArray(d_number_patches);
   for (int ip = 0; ip < d_number_patches; ip++) {
      d_active_patch_data[ip] = new hier::ComponentSelector();
   }
}

template<int DIM>
void LocallyActiveDataPatchLevelManager<DIM>::reset(
   const tbox::Pointer< hier::PatchLevel<DIM> > level)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!level.isNull());
#endif
   d_patch_level = level;
   d_number_patches = d_patch_level->getNumberOfPatches();
   d_active_patch_data.resizeArray(0);
   d_active_patch_data.resizeArray(d_number_patches);
   for (int ip = 0; ip < d_number_patches; ip++) {
      d_active_patch_data[ip] = new hier::ComponentSelector();
   }
}

/*
*************************************************************************
*                                                                       *
* Accessory functions to check whether patch data is allocated on       *
* active patches, and to allocate/deallocate data on active patches.    *
*                                                                       *
*************************************************************************
*/

template<int DIM>
bool LocallyActiveDataPatchLevelManager<DIM>::checkAllocated(
   int patch_data_index) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d_patch_level.isNull());
#endif
   bool ret_val = true;
   Iterator ip = getIterator(patch_data_index);
   for ( ; ip; ip++) {
      ret_val &= 
         d_patch_level->getPatch(ip())->checkAllocated(patch_data_index);
   }
   return(ret_val);
}

template<int DIM>
void LocallyActiveDataPatchLevelManager<DIM>::allocatePatchData(
   int patch_data_index,
   double timestamp,
   tbox::Pointer<tbox::Arena> pool) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d_patch_level.isNull());
#endif 
   Iterator ip = getIterator(patch_data_index);
   for ( ; ip; ip++) {
      d_patch_level->getPatch(ip())->allocatePatchData(patch_data_index,
                                                       timestamp, pool);
   }
}

template<int DIM>
void LocallyActiveDataPatchLevelManager<DIM>::deallocatePatchData(
   int patch_data_index) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d_patch_level.isNull());
#endif
   Iterator ip = getIterator(patch_data_index);
   for ( ; ip; ip++) {
      d_patch_level->getPatch(ip())->deallocatePatchData(patch_data_index);
   }
}

/*
*************************************************************************
*                                                                       *
* Print all locally-active variable information to given output stream. *
*                                                                       *
*************************************************************************
*/

template<int DIM>
void LocallyActiveDataPatchLevelManager<DIM>::printClassData(ostream& os) const
{
   os << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
      << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
   os << "Printing LocallyActiveDataPatchLevelManager<DIM> information...";
   os << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
   os << "d_patch_level = " << (hier::PatchLevel<DIM>*)d_patch_level << endl;
   os << "level number = " << d_patch_level->getLevelNumber() << endl;
   os << "d_number_patches = " << d_number_patches << endl;
   for (int ip = 0; ip < d_number_patches; ip++) {
      const int bvsize = d_active_patch_data[ip]->getSize();
      os << "\nActive data ids on patch " << ip 
         << ":\n    ";
      int id = 0;
      bool found_first = false;
      while (!found_first && (id < bvsize)) {
         if (d_active_patch_data[ip]->isSet(id)) {
            found_first = true;
            os << id;
         }
         id++;
      }
      for ( ; id < bvsize; id++) {
         if (d_active_patch_data[ip]->isSet(id)) {
            os << " , " << id;
         }
      }
   }
}

/*
*************************************************************************
*                                                                       *
* Implementation of LocallyActiveDataPatchLevelIterator<DIM>.           *
*                                                                       *
*************************************************************************
*/

template<int DIM>
LocallyActiveVariableDatabase<DIM>*
   LocallyActiveDataPatchLevelIterator<DIM>::s_variable_database =
      LocallyActiveVariableDatabase<DIM>::getDatabase();

template<int DIM>
LocallyActiveDataPatchLevelIterator<DIM>::
LocallyActiveDataPatchLevelIterator(
   const tbox::Pointer< hier::Variable<DIM> > variable,
   const hier::PatchLevel<DIM>& pl)
{
   initialize(variable, pl);
   while ( notActivePatch(d_patch) ) {
      d_patch++;
   }
}

template<int DIM>
LocallyActiveDataPatchLevelIterator<DIM>::
LocallyActiveDataPatchLevelIterator(
   const tbox::Pointer< hier::Variable<DIM> > variable,
   const hier::PatchLevel<DIM>* pl)
{
   initialize(variable, pl);
   while ( notActivePatch(d_patch) ) {
      d_patch++;
   }
}

template<int DIM>
LocallyActiveDataPatchLevelIterator<DIM>::
LocallyActiveDataPatchLevelIterator(
   int patch_data_index,
   const hier::PatchLevel<DIM>* pl,
   const tbox::Array< tbox::Pointer<hier::ComponentSelector> >* active_data_indices)
{
   d_patch = 0;
   d_data_index = patch_data_index;
   d_active_patch_info = active_data_indices;
   d_number_patches = pl->getNumberOfPatches();
   d_mapping = &(pl->getProcessorMapping());
   while ( notActivePatch(d_patch) ) {
      d_patch++;
   }  
}

template<int DIM>
LocallyActiveDataPatchLevelIterator<DIM>::
~LocallyActiveDataPatchLevelIterator()
{
}

template<int DIM>
void LocallyActiveDataPatchLevelIterator<DIM>::initialize(
   const tbox::Pointer< hier::Variable<DIM> > variable,
   const hier::PatchLevel<DIM>& pl)
{
   d_patch = 0;
   d_data_index = s_variable_database->mapVariableToIndex(variable);
   d_active_patch_info = &(s_variable_database->
                           getLevelManager(pl)->d_active_patch_data);
   d_number_patches = pl.getNumberOfPatches();
   d_mapping = &(pl.getProcessorMapping());
}

template<int DIM>
void LocallyActiveDataPatchLevelIterator<DIM>::initialize(
   const tbox::Pointer< hier::Variable<DIM> > variable,
   const hier::PatchLevel<DIM>* pl)
{
   d_patch = 0;
   d_data_index = s_variable_database->mapVariableToIndex(variable);
   d_active_patch_info = &(s_variable_database->
                           getLevelManager(pl)->d_active_patch_data);
   d_number_patches = pl->getNumberOfPatches();
   d_mapping = &(pl->getProcessorMapping());
}

template<int DIM>
void LocallyActiveDataPatchLevelIterator<DIM>::operator++(int)
{
   do {
      d_patch++;
   } while ( notActivePatch(d_patch) );
}

}
}

#endif
