//
// File:        VariableDatabase.C
// Package:     SAMRAI hierarchy 
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 551 $
// Modified:    $Date: 2005-08-17 11:15:27 -0700 (Wed, 17 Aug 2005) $
// Description: Manager class for variables used in a SAMRAI application.
//

#ifndef included_hier_VariableDatabase_C
#define included_hier_VariableDatabase_C

#include "VariableDatabase.h"

#include <typeinfo>
using namespace std;

#include "tbox/ShutdownRegistry.h"
#include "tbox/Utilities.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

#ifdef DEBUG_NO_INLINE
#include "VariableDatabase.I"
#endif

#include <stdio.h>

namespace SAMRAI {
    namespace hier {

#define CONTEXT_ARRAY_SCRATCH_SPACE (10)
#define VARIABLE_ARRAY_SCRATCH_SPACE (100)
#define DESCRIPTOR_ARRAY_SCRATCH_SPACE (200)

#ifndef NULL
#defined NULL (0)
#endif

template<int DIM> VariableDatabase<DIM>*
VariableDatabase<DIM>::s_variable_database_instance = NULL;
template<int DIM> bool VariableDatabase<DIM>::s_registered_callback = false;

/*
*************************************************************************
*                                                                       *
* Static database member functions.                                     *
*                                                                       *
*************************************************************************
*/

template<int DIM> VariableDatabase<DIM>* VariableDatabase<DIM>::getDatabase()
{
   if (!s_variable_database_instance) {
      s_variable_database_instance = new VariableDatabase<DIM>();
   }
   if (!s_registered_callback) {
      tbox::ShutdownRegistry::registerShutdownRoutine(freeDatabase,
                             tbox::ShutdownRegistry::priorityVariableDatabase);
      s_registered_callback = true;
   }
   return(s_variable_database_instance);
}

template<int DIM> void VariableDatabase<DIM>::freeDatabase()
{
   if (s_variable_database_instance) delete s_variable_database_instance;
   s_variable_database_instance = ((VariableDatabase<DIM>*) NULL);
}

/*
*************************************************************************
*                                                                       *
* Protected VariableDatabase constructor, destructor, and function to   *
* register Singleton subclass instance for inheritance.                 *
*                                                                       *
*************************************************************************
*/

template<int DIM>  VariableDatabase<DIM>::VariableDatabase()
{
   d_patch_descriptor = new hier::PatchDescriptor<DIM>();

   d_max_variable_id = idUndefined();
   d_max_context_id = idUndefined();
   d_max_descriptor_id = idUndefined();
}

template<int DIM>  VariableDatabase<DIM>::~VariableDatabase()
{
}

template<int DIM> void VariableDatabase<DIM>::registerSingletonSubclassInstance(
   VariableDatabase<DIM>* subclass_instance)
{
   if (!s_variable_database_instance) {
      s_variable_database_instance = subclass_instance;
      if (!s_registered_callback) {
         tbox::ShutdownRegistry::registerShutdownRoutine(freeDatabase,
                             tbox::ShutdownRegistry::priorityVariableDatabase);
         s_registered_callback = true;
      }
   } else {
      TBOX_ERROR("hier::VariableDatabase<DIM> internal error...\n"
                 << "Attemptng to set Singleton instance to subclass instance,"
                 << "\n but Singleton instance already set." << endl);
   }
}

/*
*************************************************************************
*                                                                       *
* Return the context in the database with the given name, or add a      *
* context to the database with that name if no such context exists.     *
*                                                                       *
*************************************************************************
*/

template<int DIM> tbox::Pointer<hier::VariableContext> 
VariableDatabase<DIM>::getContext(const string& name)
{
   tbox::Pointer<hier::VariableContext> context(NULL);

   if (!name.empty()) {

      int ctxt_id = getContextId_Private(name);

      if (ctxt_id == idUndefined()) {
         context = new hier::VariableContext(name);
         addContext_Private(context);
      } else {
         context = d_contexts[ctxt_id];
      }

   } 

   return(context);
}

/*
*************************************************************************
*                                                                       *
* Return true if context with given name exists in database;            *
* otherwise return false.                                               *
*                                                                       *
*************************************************************************
*/

template<int DIM> bool VariableDatabase<DIM>::checkContextExists(const string& name) const
{
   int ctxt_id = getContextId_Private(name); 

   return( (ctxt_id != idUndefined()) );
}

/*
*************************************************************************
*                                                                       *
* Add user-defined variable to database if it doesn't already exist in  *
* the database.                                                         *
*                                                                       *
*************************************************************************
*/

template<int DIM> void VariableDatabase<DIM>::addVariable(
   const tbox::Pointer< hier::Variable<DIM> > variable)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!variable.isNull());
#endif

   const bool user_variable = true;
   bool variable_added = addVariable_Private(variable, user_variable);

   if (!variable_added) {
      TBOX_ERROR("hier::VariableDatabase<DIM>::addVariable() error...\n"
         << "Attempt to add variable with duplicate name " << variable->getName()
         << " to database is not allowed.\n"
         << "Another variable with this name already exists in database." 
         << endl); 
   }

}

/*
*************************************************************************
*                                                                       *
* Return variable in database with given name.  If no such variable     *
* resides in database, return a null pointer.                           *
*                                                                       *
*************************************************************************
*/

template<int DIM> tbox::Pointer< hier::Variable<DIM> >
VariableDatabase<DIM>::getVariable(const string& name) const
{
   tbox::Pointer< hier::Variable<DIM> > variable(NULL);

   int var_id = getVariableId(name);

   if (var_id != idUndefined()) {
      variable = d_variables[var_id];
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

template<int DIM> bool VariableDatabase<DIM>::checkVariableExists(const string& name) const
{
   int var_id = getVariableId(name);

   return( (var_id != idUndefined()) );
}

/*
*************************************************************************
*                                                                       *
* Create new patch data index index by cloning factory for variable     *
* at the old index and return index of new factory.   Note that the     *
* function checkVariablePatchDataIndex() checks type of variable        *
* against given patch data index.   If these types match, then add      *
* variable and new patch data index to database.  If the types do not   *
* match, the program will abort with an error message in the private    *
* routine checkVariablePatchDataIndex().                                *
*                                                                       *
*************************************************************************
*/

template<int DIM> int VariableDatabase<DIM>::registerClonedPatchDataIndex(
   const tbox::Pointer< hier::Variable<DIM> > variable,
   int old_id)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!variable.isNull());
#endif

   int new_id = idUndefined();

   if ( checkVariablePatchDataIndex(variable, old_id) ) {

      string old_name = d_patch_descriptor->mapIndexToName(old_id);
      string new_name;

      char old_id_string[4+1/*for the null terminator*/];
      sprintf(old_id_string,"%04d",old_id);

      if (old_name.find("-clone_of_id=") == string::npos) {
         new_name = old_name + "-clone_of_id=" + old_id_string;
      } else {
         string::size_type last_dash = old_name.rfind("=");
         new_name = old_name.substr(0, last_dash+1) + old_id_string;
      }

      new_id = d_patch_descriptor->definePatchDataComponent(
               new_name,
               d_patch_descriptor->getPatchDataFactory(old_id)->
                                   cloneFactory() );

      const bool user_variable = true;
      addVariablePatchDataIndexPairToDatabase_Private(variable, new_id, user_variable);

   } else {

      TBOX_ERROR("hier::VariableDatabase<DIM>::registerClonedPatchDataIndex()"
         << "  error...\n"
         << "Variable with name " << variable->getName()
         << "\n does not match type at descriptor index = " << old_id
         << "\n That type is " << typeid(
            *(d_patch_descriptor->getPatchDataFactory(old_id))).name()
         << endl);
   }

   return(new_id);
}

/*
*************************************************************************
*                                                                       *
* Add patch data index and variable pair to the database.  Note         *
* that the function checkVariablePatchDataIndex() checks type of        *
* variable against given patch data index.  If the types do not match,  *
* the program will abort with an error message in the private routine   *
* checkVariablePatchDataIndex().   If the input index is undefined,     *
* we clone the default variable factory and add this new index to the   *
* database.  In any case, the index of the index-variable pair that     *
* is added to the database is returned.                                 * 
*                                                                       *
*************************************************************************
*/

template<int DIM> int VariableDatabase<DIM>::registerPatchDataIndex(
   const tbox::Pointer< hier::Variable<DIM> > variable,
   int data_id)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!variable.isNull());
#endif

   int new_id = data_id;

   if (new_id == idUndefined()) {
   
       new_id = d_patch_descriptor->definePatchDataComponent(
                  variable->getName(),
                  variable->getPatchDataFactory()->cloneFactory() );

       const bool user_variable = true;
       addVariablePatchDataIndexPairToDatabase_Private(variable, new_id, user_variable);

   } else {

      if ( checkVariablePatchDataIndex(variable, new_id) ) {

         const bool user_variable = true;
         addVariablePatchDataIndexPairToDatabase_Private(variable, new_id, user_variable);
   
      } else {

         TBOX_ERROR("hier::VariableDatabase<DIM>::registerPatchDataIndex()"
               << "  error...\n"
               << "Variable with name " << variable->getName()
               << "\n does not match type at patch data index = " << new_id
               << "\n That type is " << typeid(
                  *(d_patch_descriptor->getPatchDataFactory(data_id))).name() 
               << endl);
   
      }

   }

   return(new_id); 
}

/*
*************************************************************************
*                                                                       *
* Remove the given patch data index from the database.  Also, clear     *
* the index from the patch descriptor if the index is in the database.  *
*                                                                       *
*************************************************************************
*/

template<int DIM> void VariableDatabase<DIM>::removePatchDataIndex(int data_id)
{

   if ( (data_id >= 0) && (data_id <= d_max_descriptor_id) ) {

      tbox::Pointer< hier::Variable<DIM> > variable = d_index2variable_map[data_id]; 

      if (!variable.isNull()) {

         tbox::Array<int>& indx_array = 
            d_variable_context2index_map[variable->getInstanceIdentifier()];
         int array_size = indx_array.getSize();
         for (int i = 0; i < array_size; i++) {
            if (indx_array[i] == data_id) {
               indx_array[i] = idUndefined();
               break;
            }
         }

         d_patch_descriptor->removePatchDataComponent(data_id);

         d_index2variable_map[data_id].setNull();
         if (data_id == d_max_descriptor_id) {
            for (int id = d_max_descriptor_id; id >= 0; id--) {
               if (d_index2variable_map[id].isNull()) {
                  d_max_descriptor_id--;  
               } else {
                  break;
               }
            }
         }

      }

   }

}

/*
*************************************************************************
*                                                                       *
* Return true if the type of the variable matches the type of the       *
* patch data at the given patch data index.  Otherwise, return false.   *
*                                                                       *
*************************************************************************
*/

template<int DIM> bool VariableDatabase<DIM>::checkVariablePatchDataIndex(
   const tbox::Pointer< hier::Variable<DIM> > variable,
   int data_id)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!variable.isNull());
   assert(data_id >= 0 && data_id < d_patch_descriptor->getMaxNumberRegisteredComponents());
#endif

   bool ret_value = false;
   
   if ( !(d_patch_descriptor->getPatchDataFactory(data_id).isNull()) ) {

      tbox::Pointer< hier::PatchDataFactory<DIM> > dfact =
         d_patch_descriptor->getPatchDataFactory(data_id);

      if ( !dfact.isNull() && 
           (typeid(*(variable->getPatchDataFactory())) == typeid(*dfact)) ) {
         ret_value = true;
      }

   }

   return(ret_value);
}

/*
*************************************************************************
*                                                                       *
* Register variable-context pair with the database and return           *
* patch daya index corresponding to this pair and given ghost width.    *
* The steps are:                                                        *
*                                                                       *
* (1) Check whether variable-context pair maps to a valid patch data    *
*     index in the database.  If it does, then check whether the        *
*     index is null in the patch descriptor.  If it is, then we will    *
*     create a new patch data index.  If the index is not null in       *
*     the patch descriptor, the we check to see if the ghost width of   *
*     that patch data index matches that in the argument list.  If the  *
*     ghost width does not match, then we report an error and abort.    *
*                                                                       *
* (2) If we find a matching patch data index in step 1, we are done.    *
*     We return the index.                                              *
*                                                                       *
* (3) If we need to create a new patch data index, do the following:    * 
*                                                                       *
*     (3a) Create a new patch data factory, add it to the patch         *
*          descriptor, and record the index.                            *
*                                                                       *
*     (3b) We add the context to the database, if not already there.    *
*                                                                       *
*     (3c) We add the variable, and index to variable map to the        *
*          database, if not already there.                              *
*                                                                       *
*     (3d) We add the variable-context to index map to the database.    *
*                                                                       *
* (4) In the end, we return the patch data index for the                *
*     variable-context pair.                                            *
*                                                                       *
*************************************************************************
*/

template<int DIM> int VariableDatabase<DIM>::registerVariableAndContext(
   const tbox::Pointer< hier::Variable<DIM> > variable,
   const tbox::Pointer<hier::VariableContext> context,
   const hier::IntVector<DIM>& ghosts) 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!variable.isNull());
   assert(!context.isNull());
   assert(ghosts.min() >= 0);
#endif

   static string separator = "##";

   int desc_id = idUndefined();

   bool make_new_factory = true;
   int context_id = context->getIndex();
   int variable_id = variable->getInstanceIdentifier();

   if (variable_id <= d_max_variable_id) {
      tbox::Array<int>& test_indx_array = d_variable_context2index_map[variable_id];
      if (context_id < test_indx_array.getSize()) {
         desc_id = test_indx_array[context_id];
         if (desc_id != idUndefined()) {
            tbox::Pointer< hier::PatchDataFactory<DIM> > factory = 
               d_patch_descriptor->getPatchDataFactory(desc_id);
            if ( !factory.isNull() &&
                 (factory->getDefaultGhostCellWidth() != ghosts) ) {
               TBOX_ERROR("hier::VariableDatabase<DIM>::registerVariableAndContext"
                  << " error ...\n" << "Attempting to to register variable " 
                  << variable->getName() 
                  << " and context " << context->getName()
                  << " with ghost width = " << ghosts 
                  << "\n This variable-context pair is currently "
                  << "registered with a different ghost width. " << endl);
            } else {
               if (!factory.isNull()) {
                  make_new_factory = false;
               }
            }
         }     // if (desc_id != idUndefined())
      }     // if (context_id < test_indx_array.getSize())
   }     // if (variable_id <= d_max_variable_id)
     
   if (make_new_factory) { 

      tbox::Pointer< hier::PatchDataFactory<DIM> > new_factory =
         variable->getPatchDataFactory()->cloneFactory();
      new_factory->setDefaultGhostCellWidth(ghosts);
   
      string tmp(variable->getName());
      tmp += separator;
      tmp += context->getName();
      desc_id = d_patch_descriptor->definePatchDataComponent(tmp, new_factory);

      addContext_Private(context);

      const bool user_variable = true;
      addVariablePatchDataIndexPairToDatabase_Private(variable, desc_id, user_variable);

      tbox::Array<int>& var_indx_array = d_variable_context2index_map[variable_id];
      int oldsize = var_indx_array.getSize();
      int newsize = context_id + 1;
      if ( oldsize < newsize ) {
         var_indx_array.resizeArray(newsize);
         for (int i = oldsize; i < newsize; i++) {
            var_indx_array[i] = idUndefined();
         }
      }
      var_indx_array[context_id] = desc_id;

   }  // if (make_new_factory) 

   return(desc_id);

}

/*
*************************************************************************
*                                                                       *
* Return patch data index that is mapped to given variable-context      *
* pair.  If variable-context pair does not exist in database, return    *
* an undefined patch data index of idUndefined().                       *
*                                                                       *
*************************************************************************
*/

template<int DIM> int VariableDatabase<DIM>::mapVariableAndContextToIndex(
   const tbox::Pointer< hier::Variable<DIM> > variable,
   const tbox::Pointer<hier::VariableContext> context) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(variable.isNull()));
   assert(!(context.isNull()));
#endif

   int index = idUndefined();

   int var_id  = variable->getInstanceIdentifier();
   int ctxt_id = context->getIndex();

   if ( (var_id <= d_max_variable_id) &&
        (ctxt_id < d_variable_context2index_map[var_id].getSize()) ) {

      index = d_variable_context2index_map[var_id][ctxt_id];

   }

   return(index);

}

/*
*************************************************************************
*                                                                       *
* Return true if given patch data index is mapped to some variable      *
* in the database and set the variable pointer to that variable.        *
* Otherwise, return false and set the variable pointer to null.         *
*                                                                       *
*************************************************************************
*/

template<int DIM> bool VariableDatabase<DIM>::mapIndexToVariable(
   const int index,
   tbox::Pointer< hier::Variable<DIM> >& variable) const
{
   variable.setNull();

   if ( (index >= 0) && (index <= d_max_descriptor_id) ) {
      variable = d_index2variable_map[index]; 
   }

   return(!variable.isNull());
}

/*
*************************************************************************
*                                                                       *
* Return true if specified index is mapped to some variable-context     *
* pair in the database and set the variable and context pointers        *
* appropriately.  Otherwise, return false and set the pointers to null. *
*                                                                       *
*************************************************************************
*/

template<int DIM> bool VariableDatabase<DIM>::mapIndexToVariableAndContext(
   const int index,
   tbox::Pointer< hier::Variable<DIM> >& variable,
   tbox::Pointer<hier::VariableContext>& context) const
{
   bool found = false;

   variable.setNull();
   context.setNull();   

   if ( (index >= 0) && (index <= d_max_descriptor_id) ) {

      variable = d_index2variable_map[index];

      if (!variable.isNull()) {

         const tbox::Array<int>& var_indx_array = 
            d_variable_context2index_map[variable->getInstanceIdentifier()];
         int arr_size = var_indx_array.getSize();
         for (int i = 0; i < arr_size; i++) {
            if (var_indx_array[i] == index) {
               found = true;
               context = d_contexts[i];
               break;      
            }
         }

      }
    
   } 

   return(found);

}


/*
*************************************************************************
*                                                                       *
* Print all context, variable, and patch data index data                *
* contained in database to given output stream.                         *
*                                                                       *
*************************************************************************
*/

template<int DIM> void VariableDatabase<DIM>::printClassData(
   ostream& os,
   bool print_only_user_defined_variables) const
{
   int i;
   os << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
      << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
   os << "Printing hier::VariableDatabase<DIM> information...";
   os << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
   os << "Variable Contexts registered with database:";
   for (i = 0; i <= d_max_context_id; i++) {
      os << "\nContext id = " << i;
      if (!d_contexts[i].isNull()) {
         os << " : Context name = " << d_contexts[i]->getName();
      } else {
         os << " : NOT IN DATABASE"; 
      }
   }
   os << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
      << endl << flush;
   os << "Variables registered with database:";
   for (i = 0; i <= d_max_variable_id; i++) {
      os << "\nVariable instance = " << i;
      if (!d_variables[i].isNull()) {
         os << "\n";
         if ( !print_only_user_defined_variables ||
              (print_only_user_defined_variables &&
               d_is_user_variable[i]) ) {
            os << "   Variable name = " << d_variables[i]->getName();
            os << "\n   Variable type = " << typeid(*(d_variables[i])).name();
         } else {
            os << "   internal SAMRAI variable"; 
         }
      } else {
         os << " : NOT IN DATABASE"; 
      }
   }
   os << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
      << endl << flush;
   os << "Variable-Context pairs mapping to Patch Data Indices in database:"; 
   for (i = 0; i <= d_max_variable_id; i++) {
      if (!d_variables[i].isNull()) {
         if ( !print_only_user_defined_variables ||
              (print_only_user_defined_variables &&
               d_is_user_variable[i]) ) {
            os << "\nVariable name = " << d_variables[i]->getName();
            int nctxts = d_variable_context2index_map[i].getSize();
            if (nctxts > 0) {
               for (int j = 0; j < nctxts; j++) {
                  if (d_variable_context2index_map[i][j] != idUndefined()) {
                     os << "\n   context id = " << j << ", name = " 
                                                  << d_contexts[j]->getName()
                                                  << " :  patch data id = "
                                                  << d_variable_context2index_map[i][j];
                  } else {
                     os << "\n   context id = " << j << " UNDEFINED for this variable";
                  } 
               }
            } else {
               os << "\n   --- No contexts defined ---";
            }
         }
      }
   }
   os << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
      << endl << flush;
   os << "Mapping from Patch Data Indices to Variables:";
   for (i = 0; i <= d_max_descriptor_id; i++) {
      os << "\nPatch data id = " << i << " -- ";
      if (d_index2variable_map[i].isNull()) { 
         os << "UNDEFINED in database";
      } else {
         if (!print_only_user_defined_variables ||
              (print_only_user_defined_variables &&
               d_is_user_variable[i]) ) {
            os << "data factory name = " 
               << d_patch_descriptor->mapIndexToName(i); 
         } else {
            os << "internal SAMRAI patch data"; 
         }
      }
   }
   os << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
      << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
      << endl << flush;
   os << "Printing contents of patch descriptor for comparison to database..." << endl;
   d_patch_descriptor->printClassData(os);
   os << flush;
   os << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
      << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
      << endl << flush;
}

/*
*************************************************************************
*                                                                       *
* Add patch data index and variable pair to the database.  Note         *
* that the function checkVariablePatchDataIndex() checks type of        *
* variable against given patch data index.  If the types do not match,  *
* the program will abort with an error message in the private routine   *
* checkVariablePatchDataIndex().   If the input index is undefined,     *
* we clone the default variable factory and add this new index to the   *
* database.  In any case, the index of the index-variable pair that     *
* is added to the database is returned.                                 * 
*                                                                       *
*************************************************************************
*/

template<int DIM> int 
VariableDatabase<DIM>::makeInternalSAMRAIWorkVariablePatchDataIndex(
   const tbox::Pointer< hier::Variable<DIM> > variable,
   const hier::IntVector<DIM>& ghosts) 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!variable.isNull());
   assert(ghosts.min() >= 0);
#endif
   static string samrai_internal = "##SAMRAI-Internal-Variable";

   int new_id = idUndefined();

   int var_id = variable->getInstanceIdentifier();
   if ( (var_id <= d_max_variable_id) && 
         !d_variables[var_id].isNull() && 
         d_is_user_variable[var_id] ) {
       TBOX_ERROR("hier::VariableDatabase<DIM>::registerInternalSAMRAIWorkVariable error...\n"
          << "Attempt to register internal work variable named " << variable->getName()
          << " with database,\n" 
          << "But, that variable is already registered with the database as a user-defined variable."
          << endl);
   }

   tbox::Pointer< hier::PatchDataFactory<DIM> > new_factory =
      variable->getPatchDataFactory()->cloneFactory();
   new_factory->setDefaultGhostCellWidth(ghosts);

   string tmp(variable->getName());
   tmp += samrai_internal;
   new_id = d_patch_descriptor->definePatchDataComponent(tmp, new_factory);

   const bool user_variable = false;
   addVariablePatchDataIndexPairToDatabase_Private(variable, new_id, user_variable);

   return(new_id); 
}

/*
*************************************************************************
*                                                                       *
* Remove the given patch data index from the database if it has been    *
* generated as an internal SAMRAI work variable index.  Also, clear     *
* the index from the patch descriptor if the index is in the database.  *
*                                                                       *
*************************************************************************
*/

template<int DIM> void
VariableDatabase<DIM>::removeInternalSAMRAIWorkVariablePatchDataIndex(
   int data_id)
{
   if ( (data_id >= 0) && (data_id <= d_max_descriptor_id) ) {

      tbox::Pointer< hier::Variable<DIM> > variable = d_index2variable_map[data_id]; 

      if ( !variable.isNull() && 
           !d_is_user_variable[variable->getInstanceIdentifier()]) {
         removePatchDataIndex(data_id);
      }

   }
}

/*
*************************************************************************
*                                                                       *
* Protected member function to get variable id by string name.          *
*                                                                       *
*************************************************************************
*/

template<int DIM> int 
VariableDatabase<DIM>::getVariableId(const string& name) const
{
   int ret_id = idUndefined();

   if (!name.empty()) {
      for (int i = 0; i <= d_max_variable_id; i++) {
         if ( !d_variables[i].isNull() &&
              (d_variables[i]->getName() == name) ) {
            ret_id = i;
            break;
         }
      }
   }

   return(ret_id);
}

/*
*************************************************************************
*                                                                       *
* Private member functions to add contexts to database and to look up   *
* context by string name.                                               *
*                                                                       *
*************************************************************************
*/

template<int DIM> int 
VariableDatabase<DIM>::getContextId_Private(const string& name) const
{
   int ret_id = idUndefined();

   if (!name.empty()) {
      for (int i = 0; i <= d_max_context_id; i++) {
         if ( !d_contexts[i].isNull() &&
              (d_contexts[i]->getName() == name) ) {
            ret_id = i;
            break;
         }
      }
   }

   return(ret_id);
}

template<int DIM> void VariableDatabase<DIM>::addContext_Private(
   const tbox::Pointer<hier::VariableContext> context)
{
   int new_id = context->getIndex();
   int oldsize = d_contexts.getSize();
   int newsize = new_id + 1;
   if (oldsize < newsize) {
      newsize = tbox::Utilities::imax(oldsize+CONTEXT_ARRAY_SCRATCH_SPACE,
                                    newsize);
      d_contexts.resizeArray(newsize);
   }
   d_contexts[new_id] = context;
   d_max_context_id = tbox::Utilities::imax(d_max_context_id, new_id);
}

/*
*************************************************************************
*                                                                       *
* Private member functions to add mapping from data index to variable   *
* to the database.  Note that no error checking is done.                * 
*                                                                       *
*************************************************************************
*/

template<int DIM> void 
VariableDatabase<DIM>::addVariablePatchDataIndexPairToDatabase_Private(
   const tbox::Pointer< hier::Variable<DIM> > variable,
   int data_id,
   bool user_variable)
{
   bool variable_added = addVariable_Private(variable, user_variable);

   if (!variable_added) {
      TBOX_ERROR("Internal hier::VariableDatabase<DIM> error...\n"
         << "Attempt to add variable with duplicate name " << variable->getName()
         << " to database is not allowed.\n"
         << "Another variable with this name already exists in database."
         << endl);
   }

   int oldsize = d_index2variable_map.getSize();
   if (data_id >= oldsize) {
      d_index2variable_map.resizeArray(
         tbox::Utilities::imax(oldsize+DESCRIPTOR_ARRAY_SCRATCH_SPACE,
                               data_id+1) );
   }

   d_index2variable_map[data_id] = variable;
   d_max_descriptor_id = tbox::Utilities::imax(d_max_descriptor_id, data_id);
}

/*
*************************************************************************
*                                                                       *
* Add variable to database if it doesn't already exist in the database. *
* If variable already exists in the database, do nothing.  Note that    *
* we check ensure that no two distinct user-defined variables can exist *
* in the database with the same name.                                   *
*                                                                       *
*************************************************************************
*/

template<int DIM> bool 
VariableDatabase<DIM>::addVariable_Private(
   const tbox::Pointer< hier::Variable<DIM> > variable,
   bool user_variable)
{
   bool ret_value = true;

   int var_id = variable->getInstanceIdentifier();
   bool var_found = false;
   bool grow_array = false;

   if (var_id < d_variables.getSize()) {
      var_found = !d_variables[var_id].isNull();
   } else {
      grow_array = true;
   }

   if (!var_found) {

      if (user_variable) {
         if ( getVariableId(variable->getName()) != idUndefined() ) {
            ret_value = false;
         }
      }

      if (ret_value) {

         if (grow_array) {
            const int newsize =
               tbox::Utilities::imax(d_variables.getSize() +
                                     VARIABLE_ARRAY_SCRATCH_SPACE,
                                     var_id + 1);
            d_variables.resizeArray(newsize);
            d_variable_context2index_map.resizeArray(newsize);

            const int oldsize = d_is_user_variable.getSize();
            d_is_user_variable.resizeArray(newsize);
            for (int i = oldsize; i < newsize; i++) {
               d_is_user_variable[i] = false;
            }
         }

         d_variables[var_id] = variable;
         if (user_variable) {
            d_is_user_variable[var_id] = true;
         }
         d_max_variable_id = tbox::Utilities::imax(d_max_variable_id, var_id);

      }

   } // if !var_found

   return(ret_value);

}


}
}
#endif
