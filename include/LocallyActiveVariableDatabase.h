//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/hierarchy/variables/LocallyActiveVariableDatabase.h $
// Package:     SAMRAI hierarchy	
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1811 $
// Modified:	$LastChangedDate: 2007-12-20 01:19:26 -0800 (Thu, 20 Dec 2007) $
// Description:	Singleton database for variables defined on subset of hierarchy patches.
//

#ifndef included_hier_LocallyActiveVariableDatabase
#define included_hier_LocallyActiveVariableDatabase

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_hier_ComponentSelector
#include "ComponentSelector.h"
#endif
#ifndef included_hier_IntVector
#include "IntVector.h"
#endif
#ifndef included_hier_PatchDescriptor
#include "PatchDescriptor.h"
#endif
#ifndef included_hier_PatchLevel
#include "PatchLevel.h"
#endif
#ifndef included_hier_ProcessorMapping
#include "ProcessorMapping.h"
#endif
#ifndef included_hier_Variable
#include "Variable.h"
#endif
#ifndef included_hier_VariableContext
#include "VariableContext.h"
#endif
#ifndef included_hier_VariableDatabase
#include "VariableDatabase.h"
#endif
#ifndef included_tbox_Array
#include "tbox/Array.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif
#ifndef included_String
#include <string>
#define included_String
#endif
#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif

namespace SAMRAI {
    namespace hier {

// forward declarations for types used in this file
template<int DIM> class LocallyActiveDataPatchLevelIterator;
template<int DIM> class LocallyActiveDataPatchLevelManager;

/*!
 * @brief Class LocallyActiveVariableDatabase is a Singleton class that 
 * provides functionality for using the VariableDatabase to manage variables 
 * and data that live on different sets of patches in an AMR patch hierarchy, 
 * so-called "locally-active" data, can be managed.   
 *
 * This class uses the VariableDatabase directly and provides most of the 
 * operations available in that class through the interface declared here.
 * However, support for locally-active data is more limited than in the 
 * standard case.  For example, this class does not support registration of
 * locally-active variables with VariableContexts and each variable can be 
 * associated with only one patch data index. That is, there is a one-to-one 
 * mapping between variables and patch data indices and associating multiple 
 * VariableContexts with a variable is not allowed.  To make usage of this 
 * class reasonably transparent, all member functions are similar to the
 * analogous operations in the VariableDatabase.
 *
 * To avoid potentially improper, or unexpected, behavior, all locally-active
 * variables should be registered with and managed using operations provided
 * by this class.  Standard variable registration and variable management 
 * operations should use the VariableDatabase class.
 *
 * This class provides "manager" objects to define which patches are active
 * on each level for each locally-active variable; see the declaration of
 * LocallyActiveDataPatchLevelManager below.  
 *
 * @see hier::VariableDatabase
 * @see hier::PatchDescriptor
 * @see hier::LocallyActiveDataPatchLevelManager
 */

template<int DIM>
class LocallyActiveVariableDatabase 
{
   friend class LocallyActiveDataPatchLevelIterator<DIM>;

public:
   /*!
    * Return a pointer to the singleton instance of this locally-active variable 
    * database.  All access to the LocallyActiveVariableDatabase<DIM> object is 
    * through the getDatabase() function.  For example, to access the variable with 
    * string name "my_variable", use the following call: 
    * LocallyActiveVariableDatabase<DIM>::getVariable()->getVariable("my_variable").
    * 
    * Note that when the database is accessed for the first time, the
    * Singleton instance is registered with the ShutdownRegistry
    * class which destroys such objects at program completion.  Thus,
    * users of this class do not explicitly allocate or deallocate the
    * Singleton instance.
    *
    * @return  Bare pointer to variable database instance.
    */
   static LocallyActiveVariableDatabase<DIM>* getDatabase();

   /*!
    * Deallocate the LocallyActiveVariableDatabase<DIM> instance.  
    * It is not necessary to call this function at program termination, 
    * since it is automatically called by the ShutdownRegistry class.
    */
   static void freeDatabase();

   /*!
    * Return number of variables registered with the 
    * locally-active variable database.
    */
   virtual
   int getNumberRegisteredVariables() const;

   /*!
    * Return pointer to the patch descriptor managed by the database
    * (and shared by all patches in an SAMR hierarchy).
    *
    * @return  tbox::Pointer to patch descriptor instance.
    */
   tbox::Pointer< hier::PatchDescriptor<DIM> > getPatchDescriptor() const;

   /*!
    * Rgister the given variable and ghost cell width with the database of 
    * locally-active variables.  This function is similar to the variable 
    * registration member functions in the VariableDatabase<DIM> class,
    * but is more restrictive here since each variable can be registered with
    * only one patch data index.  This function imposes the same restrictions
    * on uniqueness of variable names as the VariableDatabase<DIM> base
    * class.
    *
    * Typically, this function will generate a new patch descriptor index for the
    * variable and ghost cell width and add the variable-ghost cell width pair
    * and index to the database.  If the variable-ghost cell width pair is already
    * mapped to some patch data identifier in the database, then that index
    * will be returned and the function will do nothing.   However, if
    * the variable-ghost cell width pair is already mapped to some patch data 
    * identifier with a different ghost cell width, the program will abort with a 
    * descriptive error message.  
    *
    * @return integer patch descriptor index corresponding to storage for variable.
    *
    * @param variable const pointer to variable to add to database.
    * @param ghosts   const reference to IntVector indicating ghost cell width of
    *                 variable patch data.
    *
    * When assertion checking is active, an assertion will result when
    * the variable pointer is null, or if the ghost width vector has a
    * negative entry.
    */
   int registerVariable(const tbox::Pointer< hier::Variable<DIM> > variable,
                        const hier::IntVector<DIM>& ghosts); 

   /*!
    * Get variable in locally-active database with given name string identifier.
    *
    * @return  tbox::Pointer to variable in database with given name.  If no such
    *          variable exists, a null pointer is returned.  In particular, if
    *          the variable exists in the VariableDatabase instance, but is not
    *          registered as locally-active, then a null pointer is returned.
    *
    * @param name  Const reference to name string identifying the variable.
    */
   tbox::Pointer< hier::Variable<DIM> > getVariable(const std::string& name) const;

   /*!
    * Check if variable with given name has been registered with the 
    * locally-active variable database.
    *
    * @return boolean true if a variable whose name matches the argument string
    *         exists in the database; otherwise, return false.  In particular, if
    *         the variable exists in the VariableDatabase instance, but is not
    *         registered as locally-active, then a value of false is returned.
    * 
    * @param name string name of variable to retrieve.
    */
   bool checkVariableExists(const std::string& name) const;

   /*!
    * Check if given variable has been registered with the locally-active 
    * variable database.
    *
    * @return boolean true if argument variable exists in the database; 
    *         otherwise, return false.
    *
    * @param variable smart pointer to variable to check whether it is in database.
    *        When assertion checking is active, an assertion will result when
    *        the variable pointer is null.
    */ 
   bool checkVariableExists(const tbox::Pointer< hier::Variable<DIM> > variable) const; 

   /*!
    * Check whether the given variable matches the patch data index in 
    * the locally-active variable database.  
    *
    * @return  Boolean true if the variable/patch data index pair reside in
    *          the locally-active variable database; false otherwise.
    *
    * @param  variable  tbox::Pointer to variable.  When assertion checking is
    *                   active, an unrecoverable assertion will result if
    *                   the variable pointer is null.
    * @param  data_id   Integer patch data identifier.  When assertion checking
    *                   is active, an unrecoverable assertion will result if
    *                   the value is an invalid identifier (either < 0 or    
    *                   larger than the maximum allowed patch data id).
    */
   bool
   checkVariablePatchDataIndex(const tbox::Pointer< Variable<DIM> > variable,
                               int data_id) const;

   /*!    
    * Check whether the given variable matches the patch data type
    * associated with the given patch data index in the locally-active
    * variable database.    
    *    
    * @return  Boolean true if the type of the variable matches the type of 
    *          the patch data at the given patch data index; false otherwise.
    *    
    * @param  variable  tbox::Pointer to variable.  When assertion checking
    *                   is active, an unrecoverable assertion will result if
    *                   the variable pointer is null.
    * @param  data_id   Integer patch data index.  When assertion checking
    *                   is active, an unrecoverable assertion will result
    *                   the value is an invalid identifier (either < 0 or
    *                   larger than the maximum allowed patch data id).
    */
   virtual
   bool checkVariablePatchDataIndexType(
        const tbox::Pointer< hier::Variable<DIM> > variable,
        int data_id) const;


   /*!
    * Map variable in locally-active variable database to patch data index.
    *
    * @return Integer patch data identifier for variable in database.
    *         If variable is not in the database (i.e., it has not been 
    *         registered), then an invalid patch data index (i.e., < 0) 
    *         is returned.
    *
    * @param  variable  tbox::Pointer to variable.  When assertion checking 
    *                   is active, an unrecoverable assertion results 
    *                   if the variable pointer is null.
    */
   int mapVariableToIndex(
      const tbox::Pointer< hier::Variable<DIM> > variable) const;

   /*!
    * Map patch data index to variable, if index is associated with
    * a variable registered with the locally-active variable database.
    *
    * @return  Boolean true if patch data index maps to variable in the
    *          locally-active database; otherwise false.
    *
    * @param   data_id  Integer patch data identifier.  When assertion checking
    *                   is active, an unrecoverable assertion will if the index
    *                   is invalid (i.e., < 0).
    * @param   variable  tbox::Pointer to variable that maps to patch data 
    *                    index in database.  If there is no patch data index
    *                    in the database matching the index input value, 
    *                    then the variable pointer is set to null.
    */
   bool mapIndexToVariable(
      int data_id,
      tbox::Pointer< hier::Variable<DIM> >& variable) const; 

   /*!
    * Return level manager in the database that maintains active patch data 
    * information for the input level, if possible.  Note that this database 
    * class only maintains level manager objects corresponding to patch levels 
    * that live in a patch hierarchy.  Also, the number of the patch level for 
    * each corresponding manager must be unique.  Thus, the following results 
    * are possible when this method is called:
    *
    * -#  The given patch level is not in a patch hierarchy or does not have a 
    *     valid patch hierarchy level number (i.e., it is < 0).  In this case, 
    *     an null level manager pointer is returned.
    * 
    * -#  The input patch level has a valid hierarchy level number and is in a 
    *     patch hierarchy.   If the database does not own a patch level manager
    *     corresponding to the level number of the input level, a new 
    *     patch level manager corresponding to the input level number is 
    *     added to the database and a pointer to it is returned.  If the database 
    *     owns a patch level manager and its level matches the input level, that 
    *     manager is returned.  If the database owns a patch level manager 
    *     corresponding to the level number of the input level, but its level 
    *     does not match the input level, then an unrecoverable error results.
    *
    * @return pointer to hier::LocallyActiveDataPatchLevelManager<DIM> object 
    *         for given patch level.
    *
    * @param level const patch level pointer.
    */
   tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> >
   getLocallyActiveDataPatchLevelManager(
      const tbox::Pointer< hier::PatchLevel<DIM> > level);

   /*!
    * Reset the level manager in the database that maintains active patch data
    * information for the given level, if possible.  Note that this database
    * class only maintains level manager objects corresponding to patch levels
    * that live in a patch hierarchy.  Thus, if the input level does not live
    * in a hierarchy, this method does nothing.  Also, the number of the patch 
    * level for each corresponding manager must be unique.  Thus, if the database
    * owns a level manager associated with a patch level whose number matches the 
    * input level, the existing manager is removed from the database and is replaced
    * with a new manager object for the given level.  However, if one has a smart
    * pointer to the pre-existing level manager, one can maintain that while resetting
    * the one in the database as long as the pointer reference count remains greater
    * than zero.
    *
    * @param level const patch level pointer.
    */
   void resetLocallyActiveDataPatchLevelManager(
      const tbox::Pointer< hier::PatchLevel<DIM> > level);

   /*!
    * Print all variable, context, and patch descriptor information
    * contained in the database to the specified output stream.  
    * Active patch information for the patch hierarchy is also printed 
    * by dumping the contents of the locally-active patch data level 
    * manager objects.
    *
    * @param os  Optional output stream.  If not given, tbox::plog is used.
    */
   virtual void printClassData(std::ostream& os = tbox::plog) const;

protected:
   /**
    * The constructor for LocallyActiveVariableDatabase<DIM> is protected. 
    * Consistent with the definition of a Singleton class, only the 
    * database object has access to the constructor for the class. 
    *
    * The constructor initializes the state of database contents.
    */
   LocallyActiveVariableDatabase();

   /**
    * The destructor for LocallyActiveVariableDatabase<DIM> is protected. 
    * See the comments for the constructor.
    *
    * The destructor deallocates database contents.
    */
   virtual ~LocallyActiveVariableDatabase();

private:

   /*
    * Private member function for checking database contents.
    * 
    * validLevel() returns true if argument level is known to the 
    * locally-active variable database; otherwise, returns false.
    */
   bool validLevel(const tbox::Pointer< hier::PatchLevel<DIM> > level) const;

   /*
    * Private member functions to access locally-active data information.
    * Note that these routines are used by the locally-active data 
    * patch level iterator initialization functions.
    */
   const hier::LocallyActiveDataPatchLevelManager<DIM>* 
      getLevelManager(const hier::PatchLevel<DIM>& pl) const;
   const hier::LocallyActiveDataPatchLevelManager<DIM>*
      getLevelManager(const hier::PatchLevel<DIM>* pl) const;

   /*
    * Static data members used to control access to and destruction of
    * singleton variable database instance.
    */
   static LocallyActiveVariableDatabase<DIM>* 
      s_locally_active_variable_database_instance;

   static bool s_registered_callback;

   /*
    * Static data member used to control allocation of arrays.
    */
   static int s_patchlevel_array_alloc_size;

   /*
    * Cached pointer to variable database.
    */
   hier::VariableDatabase<DIM>* d_variable_database;

   /*
    * Count of number of variables registered with 
    * locally-active variable database.
    */
   int d_num_registered_variables;

   /*
    * Single variable context shared by all locally-active variables.
    */
   tbox::Pointer<hier::VariableContext> d_locally_active_context;

   /*
    * Set of patch data identifiers associated with registered 
    * active variables.
    */
   hier::ComponentSelector d_locally_active_patch_data_ids;

   /* 
    * Array of manager objects that maintain active variable information for
    * levels in patch hierarchy. Indexing is d_patch_active_data[level #].
    */
   // Array of patch level manager pointers is indexed as 
   // d_patch_level_active_data_manager[ <level number> ]
   tbox::Array< tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > > 
      d_patch_level_active_data_manager; 

};

}
}

#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "LocallyActiveVariableDatabase.C"
#endif
