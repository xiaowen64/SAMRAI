//
// File:        LocallyActiveVariableDatabase.C
// Package:     SAMRAI hierarchy 
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 731 $
// Modified:    $Date: 2005-11-15 14:05:53 -0800 (Tue, 15 Nov 2005) $
// Description: Singleton database for variables defined on subset of hierarchy patches.
//

#ifndef included_hier_LocallyActiveVariableDatabase_C
#define included_hier_LocallyActiveVariableDatabase_C

#include "LocallyActiveVariableDatabase.h"

#include "tbox/MathUtilities.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif
#include "LocallyActiveDataPatchLevelManager.h"

#include <typeinfo>
using namespace std;

#ifdef DEBUG_NO_INLINE
#include "LocallyActiveVariableDatabase.I"
#endif

namespace SAMRAI {
    namespace hier {

#define PATCHLEVEL_ARRAY_SCRATCH_SPACE (10)
#define VARIABLE_ARRAY_SCRATCH_SPACE (100)
#define DESCRIPTOR_ARRAY_SCRATCH_SPACE (200)

#ifndef NULL
#defined NULL (0)
#endif

template<int DIM> LocallyActiveVariableDatabase<DIM>*
LocallyActiveVariableDatabase<DIM>::
   s_locally_active_variable_database_instance = NULL;

/*
*************************************************************************
*                                                                       *
* Static database member functions.                                     *
*                                                                       *
*************************************************************************
*/

template<int DIM>
LocallyActiveVariableDatabase<DIM>* 
LocallyActiveVariableDatabase<DIM>::getDatabase()
{
   if (!s_locally_active_variable_database_instance) {
      s_locally_active_variable_database_instance = 
         new LocallyActiveVariableDatabase<DIM>();
   }
   return(s_locally_active_variable_database_instance);
}

/*
*************************************************************************
*                                                                       *
* Database constructor and destructor.                                  *
*                                                                       *
*************************************************************************
*/

template<int DIM>
LocallyActiveVariableDatabase<DIM>::LocallyActiveVariableDatabase()
{
   d_locally_active_context = VariableDatabase<DIM>::getContext("LOCALLY_ACTIVE");
   d_num_variables_registered = 0;
   d_max_locally_active_variable_id = -1;

   registerSingletonSubclassInstance(this);
}

template<int DIM>
LocallyActiveVariableDatabase<DIM>::~LocallyActiveVariableDatabase()
{
   d_patch_level_active_data_manager.resizeArray(0);
}

/*
*************************************************************************
*                                                                       *
* Add variable to database if it doesn't already exist in the database. *
* If variable already exists in the database, do nothing.  Note that    *
* we check ensure that no two distinct variables can exist in the       *
* database with the same name.                                          *
*                                                                       *
* Note that each locally-active variable is also maintained in the      *
* standard variable database (superclass) instance.                     *
*                                                                       *
*************************************************************************
*/

template<int DIM>
int LocallyActiveVariableDatabase<DIM>::registerVariable(
   const tbox::Pointer< hier::Variable<DIM> > variable,
   const hier::IntVector<DIM>& ghosts)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!variable.isNull());
   assert(ghosts.min() >= 0);
#endif

   int desc_id = registerVariableAndContext(variable,
                                            d_locally_active_context,
                                            ghosts);

   int var_id = variable->getInstanceIdentifier();

   bool var_found = false;
   bool grow_array = false;

   if (var_id < d_locally_active_variables.getSize()) {
      var_found = !(d_locally_active_variables[var_id].isNull());
   } else {
      grow_array = true;
   }

   if (!var_found) {

      if (grow_array) {
         const int newsize = d_locally_active_variables.getSize() +
                             VARIABLE_ARRAY_SCRATCH_SPACE;
         d_locally_active_variables.resizeArray(newsize);
         d_locally_active_variable_descriptor_indices.resizeArray(newsize);
      }

      d_max_locally_active_variable_id = 
         tbox::MathUtilities<int>::Max(d_max_locally_active_variable_id, var_id);

      d_locally_active_variables[var_id] = variable;
      d_locally_active_variable_descriptor_indices[var_id] = desc_id;

      d_num_variables_registered++;
   }

   return(desc_id);

}

/*
*************************************************************************
*                                                                       *
* Return variable in database with given name.  If no such variable     *
* resides in database, return a null pointer.                           *
*                                                                       *
*************************************************************************
*/

template<int DIM>
tbox::Pointer< hier::Variable<DIM> >
LocallyActiveVariableDatabase<DIM>::getVariable(const string& name) const
{
   tbox::Pointer< hier::Variable<DIM> > variable(NULL);

   int var_id = VariableDatabase<DIM>::getVariableId(name);

   if ( (var_id != VariableDatabase<DIM>::idUndefined()) && 
        (d_locally_active_variables.getSize() > var_id) ) {
      variable = d_locally_active_variables[var_id];
   }

   return(variable);
}

/*
*************************************************************************
*                                                                       *
* Return true if variable with given name exists in database.           *
* Otherwise, return false.                                              *
*                                                                       *
*************************************************************************
*/

template<int DIM>
bool LocallyActiveVariableDatabase<DIM>::checkVariableExists(
   const string& name) const
{
   int var_id = VariableDatabase<DIM>::getVariableId(name);

   return( (var_id != VariableDatabase<DIM>::idUndefined()) &&
           (d_max_locally_active_variable_id > var_id) && 
           (!d_locally_active_variables[var_id].isNull()) );
}

/*
*************************************************************************
*                                                                       *
* Accessory functions for managing array of active data managers for    *
* patch levels in AMR hierarchy.                                        *
*                                                                       *
*************************************************************************
*/

template<int DIM>
tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> >
LocallyActiveVariableDatabase<DIM>::getLocallyActiveDataPatchLevelManager(
   const tbox::Pointer< hier::PatchLevel<DIM> > level)
{
   tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > ret_plman;

   if ( validLevel(level) ) {

      const int level_number = level->getLevelNumber();
      if (d_patch_level_active_data_manager.size() <= level_number) {
         const int newsize = d_patch_level_active_data_manager.getSize() +
                             PATCHLEVEL_ARRAY_SCRATCH_SPACE;
         d_patch_level_active_data_manager.resizeArray(newsize);
      }

      ret_plman = d_patch_level_active_data_manager[level_number]; 

      if (ret_plman.isNull()) {
         d_patch_level_active_data_manager[level_number] =
            new hier::LocallyActiveDataPatchLevelManager<DIM>(level);
         ret_plman = d_patch_level_active_data_manager[level_number];
      } else {
         if ( !(ret_plman->checkLevel(level)) ) {
            ret_plman.setNull(); 
            TBOX_ERROR(
               "LocallyActiveVariableDatabase<DIM>::getLocallyActiveDataPatchLevelManager"
               << " error..."
               << "\n argument level with level number " << level->getLevelNumber()
               << " is inconsistent with the current manager for that level number."
               << endl);
         }
      }

   }

   return(ret_plman);
}

template<int DIM>
void
LocallyActiveVariableDatabase<DIM>::resetLocallyActiveDataPatchLevelManager(
   const tbox::Pointer< hier::PatchLevel<DIM> > level)
{
   if ( validLevel(level) ) {

      const int level_number = level->getLevelNumber();
      if (d_patch_level_active_data_manager.size() <= level_number) {
         const int newsize = d_patch_level_active_data_manager.getSize() +
                             PATCHLEVEL_ARRAY_SCRATCH_SPACE;
         d_patch_level_active_data_manager.resizeArray(newsize);
      }

      if ( !d_patch_level_active_data_manager[level_number].isNull() ) {
         d_patch_level_active_data_manager[level_number].setNull();
      }

      d_patch_level_active_data_manager[level_number] =
         new hier::LocallyActiveDataPatchLevelManager<DIM>(level);

   }
}

/*
*************************************************************************
*                                                                       *
* Print all context, variable, and descriptor index data                *
* contained in database to given output stream.                         *
*                                                                       *
*************************************************************************
*/

template<int DIM>
void LocallyActiveVariableDatabase<DIM>::printClassData(ostream& os) const
{
   int i;
   os << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
      << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
   os << "Printing LocallyActiveVariableDatabase<DIM> information...";
   os << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
   os << "Variable Context shared by all variables in locally-active database:" << endl;
   if (!d_locally_active_context.isNull()) {
      os << "   Context name = " << d_locally_active_context->getName() << endl;
      os << "   Context identifier = " << d_locally_active_context->getIndex() << endl;
   } else {
         os << " : NO SHARED VARIABLE CONTEXT IN DATABASE";
   }
   os << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
      << endl << flush;
   os << d_num_variables_registered 
      << " variables registered with locally-active database..." << endl;
   os << "These variables and their corresponding patch descriptor indices are:";
   for (i = 0; i <= d_max_locally_active_variable_id; i++) {
      if (!d_locally_active_variables[i].isNull()) {
         os << "\nVariable instance = " << i;
         os << "\n";
         os << "   Variable name = " 
            << d_locally_active_variables[i]->getName();
         os << "\n   Variable type = " 
            << typeid(*(d_locally_active_variables[i])).name();
         os << "\n   Patch descriptor id = " 
            << d_locally_active_variable_descriptor_indices[i];
      }
   }
   os << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
      << endl << flush;
   os << "Active patch data information by variable:";
   for (i = 0; i <= d_max_locally_active_variable_id; i++) {
      if (!d_locally_active_variables[i].isNull()) {
         os << "\nVariable name = " << d_locally_active_variables[i]->getName();
         const int var_data_index = d_locally_active_variable_descriptor_indices[i];
         for (int ln = 0; ln < d_patch_level_active_data_manager.getSize(); ln++) {
            if (!d_patch_level_active_data_manager[ln].isNull()) {
               os << "\n  Active patches on level " << ln << " :\n    ";
               const int asize = d_patch_level_active_data_manager[ln]->
                           getPatchLevel()->getNumberOfPatches();
               int ip = 0;
               bool found_first = false;
               while (!found_first && (ip < asize)) {
                  if (d_patch_level_active_data_manager[ln]->
                      getPatchDataActive(var_data_index, ip)) {
                     found_first = true;
                     os << ip;
                  }
                  ip++;
               }
               for (; ip < asize; ip++) {
                  if (d_patch_level_active_data_manager[ln]->
                      getPatchDataActive(var_data_index, ip)) {
                     os << " , " << ip;
                  }
               }
            }
         }
         os << "\n ";
      }
   }
   os << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
      << endl << flush;
   os << "Active patch data information by level and patch:";
   for (i = 0; i < d_patch_level_active_data_manager.getSize(); i++) {
      if (!d_patch_level_active_data_manager[i].isNull()) {
         d_patch_level_active_data_manager[i]->printClassData(os);
      }
   }
   
   hier::VariableDatabase<DIM>::printClassData(os);
}

}
}
#endif
