//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/variables/LocallyActiveVariableDatabase.C $
// Package:     SAMRAI hierarchy 
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2043 $
// Modified:    $LastChangedDate: 2008-03-12 09:14:32 -0700 (Wed, 12 Mar 2008) $
// Description: Singleton database for variables defined on subset of hierarchy patches.
//

#ifndef included_hier_LocallyActiveVariableDatabase_C
#define included_hier_LocallyActiveVariableDatabase_C

#include "LocallyActiveVariableDatabase.h"

#include "LocallyActiveDataPatchLevelManager.h"

#include "tbox/ShutdownRegistry.h"
#include "tbox/Utilities.h"

#include <typeinfo>

namespace SAMRAI {
    namespace hier {

template<int DIM> LocallyActiveVariableDatabase<DIM>*
LocallyActiveVariableDatabase<DIM>::
   s_locally_active_variable_database_instance = 0;

template<int DIM> bool 
   LocallyActiveVariableDatabase<DIM>::s_registered_callback = false;

template<int DIM> int 
LocallyActiveVariableDatabase<DIM>::s_patchlevel_array_alloc_size = 10;

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
   if (!s_registered_callback) {
      tbox::ShutdownRegistry::registerShutdownRoutine(freeDatabase,
                             tbox::ShutdownRegistry::priorityVariableDatabase);
      s_registered_callback = true;
   }
   return(s_locally_active_variable_database_instance);
}

template<int DIM> void 
LocallyActiveVariableDatabase<DIM>::freeDatabase()
{
   if (s_locally_active_variable_database_instance) {
      delete s_locally_active_variable_database_instance;
   }
   s_locally_active_variable_database_instance = 
      ((LocallyActiveVariableDatabase<DIM>*)0);
}

/*
*************************************************************************
*                                                                       *
* Protected LocallyActiveVariableDatabase constructor and destructor.   *
*                                                                       *
*************************************************************************
*/

template<int DIM>
LocallyActiveVariableDatabase<DIM>::LocallyActiveVariableDatabase()
{
   d_variable_database = VariableDatabase<DIM>::getDatabase();
   d_locally_active_context = 
      d_variable_database->getContext("LOCALLY_ACTIVE");
   d_num_registered_variables = 0;
}

template<int DIM>
LocallyActiveVariableDatabase<DIM>::~LocallyActiveVariableDatabase()
{
   d_patch_level_active_data_manager.resizeArray(0);
}

/*
*************************************************************************
*                                                                       *
* Accessory function to return patch descriptor.                        *
*                                                                       *
*************************************************************************
*/

template<int DIM> 
tbox::Pointer< PatchDescriptor<DIM> >
LocallyActiveVariableDatabase<DIM>::getPatchDescriptor() const
{
   return( d_variable_database->getPatchDescriptor() );
}

/*
*************************************************************************
*                                                                       *
* Accessory function to return number of registered variables           *
* (which is the same as the number of registered patch data ids).       *
*                                                                       *
*************************************************************************
*/

template<int DIM> 
int 
LocallyActiveVariableDatabase<DIM>::getNumberOfRegisteredVariables() const
{
   return(d_num_registered_variables);
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
   TBOX_ASSERT(!variable.isNull());
   TBOX_ASSERT(ghosts.min() >= 0);
#endif

   int desc_id = 
      d_variable_database->registerVariableAndContext(
                              variable,
                              d_locally_active_context,
                              ghosts);

   if ( desc_id >= 0 ) {
      if ( !d_locally_active_patch_data_ids.isSet(desc_id) ) {
         d_num_registered_variables++;
      }
      d_locally_active_patch_data_ids.setFlag(desc_id);
   }

   return(desc_id);

}

/*
*************************************************************************
*                                                                       *
* Return variable in database with given name.  If no such variable     *
* resides in locally-active database, return a null pointer.            *
*                                                                       *
*************************************************************************
*/

template<int DIM>
tbox::Pointer< hier::Variable<DIM> >
LocallyActiveVariableDatabase<DIM>::getVariable(const std::string& name) const
{
   tbox::Pointer< hier::Variable<DIM> > ret_variable(0);

   tbox::Pointer< hier::Variable<DIM> > test_variable = 
      d_variable_database->getVariable(name);

   if ( !test_variable.isNull() ) {
      int data_id = 
         d_variable_database->mapVariableAndContextToIndex(
                                 test_variable,
                                 d_locally_active_context);
      if (data_id >= 0 &&
           d_locally_active_patch_data_ids.isSet(data_id) ) {
         ret_variable = test_variable;
      }
   }

   return(ret_variable);
}

/*
*************************************************************************
*                                                                       *
* Return true if variable with given name exists in locally-active      *
* database.  Otherwise, return false.                                   *
*                                                                       *
*************************************************************************
*/

template<int DIM>
bool LocallyActiveVariableDatabase<DIM>::checkVariableExists(
   const std::string& name) const
{
   bool ret_value = false;

   tbox::Pointer< hier::Variable<DIM> > test_variable =
      d_variable_database->getVariable(name);

   if ( !test_variable.isNull() ) {
      int data_id =
         d_variable_database->mapVariableAndContextToIndex(
                                 test_variable,
                                 d_locally_active_context);
      if ( data_id >= 0 &&
           d_locally_active_patch_data_ids.isSet(data_id) ) {
         ret_value = true;
      }
   }

   return(ret_value);
}

/*
*************************************************************************
*                                                                       *
* Return true if variable exists in locally-active database.            *
* Otherwise, return false.                                              *
*                                                                       *
*************************************************************************
*/

template<int DIM>
bool LocallyActiveVariableDatabase<DIM>::checkVariableExists(
   const tbox::Pointer< hier::Variable<DIM> > variable) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!variable.isNull());
#endif
   
   bool ret_value = false;

   int data_id =
       d_variable_database->mapVariableAndContextToIndex(
                               variable,
                               d_locally_active_context);
   if ( data_id >= 0 &&
        d_locally_active_patch_data_ids.isSet(data_id) ) {
      ret_value = true;
   }

   return(ret_value);
 
}

/*
*************************************************************************
*                                                                       *
* Return true if the given variable is mapped to the given patch data   *
* index.  Otherwise, return false.                                      *
*                                                                       *
*************************************************************************
*/

template<int DIM> 
bool LocallyActiveVariableDatabase<DIM>::checkVariablePatchDataIndex(
   const tbox::Pointer< hier::Variable<DIM> > variable,
   int data_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!variable.isNull());
   TBOX_ASSERT(data_id >= 0 &&
          data_id < getPatchDescriptor()->getMaxNumberRegisteredComponents());
#endif

   bool ret_value = false;

   if ( data_id >= 0 ) {
      int test_data_id = 
         d_variable_database->mapVariableAndContextToIndex(
                                 variable,
                                 d_locally_active_context);
      
      if ( test_data_id >= 0 &&
           d_locally_active_patch_data_ids.isSet(data_id) ) {
        ret_value = ( data_id == test_data_id );
      }
   }

   return(ret_value);

}

/*
*************************************************************************
*                                                                       *
* Return true if the type of the variable matches the type of the       *
* patch data at the given patch data index.  Otherwise, return false.   *
*                                                                       *
*************************************************************************
*/
 
template<int DIM> bool
LocallyActiveVariableDatabase<DIM>::checkVariablePatchDataIndexType(
   const tbox::Pointer< hier::Variable<DIM> > variable,
   int data_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!variable.isNull());
   TBOX_ASSERT(data_id >= 0 &&
          data_id < getPatchDescriptor()->getMaxNumberRegisteredComponents());
#endif

   bool ret_value = false;
 
   if ( d_locally_active_patch_data_ids.isSet(data_id) ) {

       ret_value = 
          d_variable_database->checkVariablePatchDataIndexType(
                                  variable,
                                  data_id);
 
   }
 
   return(ret_value);

}

/*
*************************************************************************
*                                                                       *
* Return patch data index mapped to variable in locally-active          *
* variable database.                                                    *
*                                                                       *
*************************************************************************
*/

template<int DIM>
int LocallyActiveVariableDatabase<DIM>::mapVariableToIndex(
   const tbox::Pointer< hier::Variable<DIM> > variable) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!variable.isNull());
#endif

   int return_data_id = d_variable_database->idUndefined();

   int test_data_id = 
      d_variable_database->mapVariableAndContextToIndex(
                              variable,
                              d_locally_active_context);
   
   if ( test_data_id >= 0 && 
        d_locally_active_patch_data_ids.isSet(test_data_id) ) {
      return_data_id = test_data_id;
   }

   return(return_data_id);

}

/*
*************************************************************************
*                                                                       *
* Map patch data index to variable in locally-active variable database. *
*                                                                       *
*************************************************************************
*/

template<int DIM>
bool LocallyActiveVariableDatabase<DIM>::mapIndexToVariable(
   int data_id,
   tbox::Pointer< hier::Variable<DIM> >& variable) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(data_id >= 0);
#endif

   bool ret_value = false;

   if ( d_locally_active_patch_data_ids.isSet(data_id) ) {
      ret_value = d_variable_database->mapIndexToVariable(data_id,
                                                          variable);
   }

   return(ret_value);

}

/*
*************************************************************************
*                                                                       *
* Return locally-active data patch level manager object associated      *
* with given patch level.                                               *
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
                             s_patchlevel_array_alloc_size;
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
               << std::endl);
         }
      }

   }

   return(ret_plman);
}

/*
*************************************************************************
*                                                                       *
* Reset locally-active data patch level manager object associated       *
* with given patch level.                                               *
*                                                                       *
*************************************************************************
*/

template<int DIM>
void
LocallyActiveVariableDatabase<DIM>::resetLocallyActiveDataPatchLevelManager(
   const tbox::Pointer< hier::PatchLevel<DIM> > level)
{
   if ( validLevel(level) ) {

      const int level_number = level->getLevelNumber();
      if (d_patch_level_active_data_manager.size() <= level_number) {
         const int newsize = d_patch_level_active_data_manager.getSize() +
                             s_patchlevel_array_alloc_size;
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
void LocallyActiveVariableDatabase<DIM>::printClassData(std::ostream& os) const
{
   os << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
      << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
   os << "Printing LocallyActiveVariableDatabase<DIM> information...";
   os << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
   os << "Variable Context shared by all variables in locally-active database:" << std::endl;
   if (!d_locally_active_context.isNull()) {
      os << "   Context name = " << d_locally_active_context->getName() << std::endl;
      os << "   Context identifier = " << d_locally_active_context->getIndex() << std::endl;
   } else {
         os << " : NO SHARED VARIABLE CONTEXT IN DATABASE";
   }
   os << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
      << std::endl << std::flush;
   os << "Active patch data information by variable:";
   const int descriptor_size = 
      d_locally_active_patch_data_ids.getSize();
   for (int id = 0; id < descriptor_size; ++id) {
      if ( d_locally_active_patch_data_ids.isSet(id) ) {
         tbox::Pointer< hier::Variable<DIM> > variable;
         (void) mapIndexToVariable( id, variable);
         os << "\nVariable name = " << variable->getName();
         for (int ln = 0; ln < d_patch_level_active_data_manager.getSize(); ++ln) {
            if (!d_patch_level_active_data_manager[ln].isNull()) {
               os << "\n  Active patches on level " << ln << " :\n    ";
               const int asize = d_patch_level_active_data_manager[ln]->
                                 getPatchLevel()->getNumberOfPatches();
               int ip = 0;
               bool found_first = false;
               while (!found_first && (ip < asize)) {
                  if ( d_patch_level_active_data_manager[ln]->
                       getPatchDataActive( PatchDataId(id), 
                                           PatchNumber(ip) ) ) {
                     found_first = true;
                     os << ip;
                  }
                  ip++;
               }
               for (; ip < asize; ip++) {
                  if ( d_patch_level_active_data_manager[ln]->
                       getPatchDataActive( PatchDataId(id), 
                                           PatchNumber(ip) ) ) {
                     os << " , " << ip;
                  }
               }
            }  // if patch level manager exists
         }  // iterate over levels
         os << "\n ";
      }  // if patch data index is locally-active
   }  // iterate over patch descriptor ids
   os << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
      << std::endl << std::flush;
   os << "Active patch data information by level and patch:";
   for (int im = 0; im < d_patch_level_active_data_manager.getSize(); ++im) {
      if (!d_patch_level_active_data_manager[im].isNull()) {
         d_patch_level_active_data_manager[im]->printClassData(os);
      }
   }
   
   d_variable_database->printClassData(os);
}

/*
*************************************************************************
*                                                                       *
* Private member functions for managing locally-active data patch       *
* level manager objects.                                                *
*                                                                       *
*************************************************************************
*/

template<int DIM>
bool
LocallyActiveVariableDatabase<DIM>::validLevel(
   const tbox::Pointer< hier::PatchLevel<DIM> > level) const
{
   return(!level.isNull() &&
          level->getLevelNumber() >= 0 &&
          level->inHierarchy());
}

template<int DIM>
const hier::LocallyActiveDataPatchLevelManager<DIM>*
LocallyActiveVariableDatabase<DIM>::getLevelManager(
   const hier::PatchLevel<DIM>& pl) const
{
   const int ln = pl.getLevelNumber();
   if ( ln >= 0 &&
        d_patch_level_active_data_manager.size() > ln &&
        !d_patch_level_active_data_manager[ln].isNull() &&
        d_patch_level_active_data_manager[ln]->checkLevel(pl) ) {
      return(d_patch_level_active_data_manager[ln]);
   } else {
      return( (hier::LocallyActiveDataPatchLevelManager<DIM>*)0 );
   }
}

template<int DIM>
const hier::LocallyActiveDataPatchLevelManager<DIM>*
LocallyActiveVariableDatabase<DIM>::getLevelManager(
   const hier::PatchLevel<DIM>* pl) const
{
   const int ln = pl->getLevelNumber();
   if ( ln >= 0 &&
        d_patch_level_active_data_manager.size() > ln &&
        !d_patch_level_active_data_manager[ln].isNull() &&
        d_patch_level_active_data_manager[ln]->checkLevel(pl) ) {
      return(d_patch_level_active_data_manager[ln]);
   } else {
      return( (hier::LocallyActiveDataPatchLevelManager<DIM>*)0 );
   }
}


}
}
#endif
