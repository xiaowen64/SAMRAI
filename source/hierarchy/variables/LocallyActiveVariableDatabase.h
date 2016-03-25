//
// File:	LocallyActiveVariableDatabase.h
// Package:     SAMRAI hierarchy	
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 731 $
// Modified:	$Date: 2005-11-15 14:05:53 -0800 (Tue, 15 Nov 2005) $
// Description:	Singleton database for variables defined on subset of hierarchy patches.
//

#ifndef included_hier_LocallyActiveVariableDatabase
#define included_hier_LocallyActiveVariableDatabase

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_hier_PatchLevel
#include "PatchLevel.h"
#endif
#ifndef included_hier_ProcessorMapping
#include "ProcessorMapping.h"
#endif
#ifndef included_hier_VariableDatabase
#include "VariableDatabase.h"
#endif

#ifndef included_iostream
#define included_iostream
#include <iostream>
using namespace std;
#endif

namespace SAMRAI {
    namespace hier {

// forward declarations for types used in this file
template<int DIM> class LocallyActiveDataPatchLevelIterator;
template<int DIM> class LocallyActiveDataPatchLevelManager;

/*!
 * @brief Class LocallyActiveVariableDatabase is a Singleton class for managing
 * variables and data that live on different sets of patches in an AMR patch
 * hierarchy; i.e., they are "locally-active".   Although this class is derived 
 * from the VariableDatabase class, this class does not extend the full range of 
 * functionality provided by the VariableDatabase base class.  Instead, this class 
 * uses some of the functionality of the VariableDatabase base class so that it is 
 * not  replicated here.  For example, this class assumes that each variable 
 * is associated with only one patch data index and all variables registered with 
 * this class share the same variable context.  However, the usage of this class is
 * similar to the VariableDatabase class.  To make usage of this class reasonably
 * transparent all functions defined in the base class that make sense for this 
 * class are redeclared here.
 *
 * @see hier::VariableContext
 * @see hier::VariableDatabase
 * @see hier::PatchDescriptor
 * @see hier::LocallyActiveDataPatchLevelManager
 */

template<int DIM>
class LocallyActiveVariableDatabase : 
   public hier::VariableDatabase<DIM>
{
   friend class LocallyActiveDataPatchLevelIterator<DIM>;
public:
   /*!
    * Return a pointer to the singleton instance of the locally-active variable 
    * database.  All access to the LocallyActiveVariableDatabase<DIM> object is 
    * through the getDatabase() function.  For example, to access the variable with 
    * string name "my_variable", use the following call: 
    * LocallyActiveVariableDatabase<DIM>::getVariable()->getVariable("my_variable").
    * 
    * Note that there is no freeDatabase() static member function as one might 
    * expect for a Singleton class.  This is unnecessary since the deallocation 
    * at program termination is handled by the VariableDatabase base class through
    * the virtual destructor.
    */
   static LocallyActiveVariableDatabase<DIM>* getDatabase();

   /*!
    * Return pointer to the patch descriptor managed by the database
    * (and shared by all patches in an SAMR hierarchy).
    *
    * @return  tbox::Pointer to patch descriptor instance.
    */
   tbox::Pointer< hier::PatchDescriptor<DIM> > getPatchDescriptor() const;

   /*!
    * Return pointer to variable context shared by all locally-active 
    * variables registered with this database.
    */
   tbox::Pointer<hier::VariableContext> getSharedContext() const;

   /*!
    * Add the given variable and ghost cell width to the database of 
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
    * the variable-ghost cell width pair is already mapped to some patch data identifier
    * with a different ghost cell width, the program will abort with a descriptive
    * error message.  
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
    *          variable exists, a null pointer is returned.
    *
    * @param name  Const reference to name string identifying the variable.
    */
   tbox::Pointer< hier::Variable<DIM> > getVariable(const string& name) const;

   /*!
    * Check if variable with given name exists in the locally-active database.
    *
    * @return boolean true if a variable whose name matches the argument string
    *         exists in the database; otherwise, return false. 
    * 
    * @param name string name of variable to retrieve.
    */
   bool checkVariableExists(const string& name) const;

   /*!
    * Check if given variable exists in the locally-active database.
    *
    * @return boolean true if argument variable exists in the database; 
    *         otherwise, return false.
    *
    * @param const smart pointer to variable to check whether it is in database.
    *
    * When assertion checking is active, an assertion will result when
    * the variable pointer is null.
    */ 
   bool checkVariableExists(const tbox::Pointer< hier::Variable<DIM> > variable) const; 

   /*!
    * Check whether the given variable matches the type of the patch data
    * at the given descriptor index.  Note this check can be performed
    * regardless of whether the variable and data index are in the variable
    * database and so this function does not provide information about the
    * contents of the database.
    *
    * @return  Boolean true if the type of the variable matches the type of
    *          the patch data at the given descriptor index; false otherwise.
    *
    * @param  variable  tbox::Pointer to variable.  When assertion checking is
    *                   active, an unrecoverable exception will result if
    *                   the variable pointer is null.
    * @param  data_id   Integer patch data identifier.  When assertion checking
    *                   is active, an unrecoverable exception will result if
    *                   the value is an invalid identifier (< 0).
    */
   bool
   checkVariablePatchDataIndex(const tbox::Pointer< Variable<DIM> > variable,
                               int data_id);

   /*!
    * Map variable in database to patch data identifier.  If variable is not 
    * in the database (i.e., it has not been added), then an invalid descriptor 
    * index (i.e., < 0) is returned.
    *
    * @return Integer patch data identifier for variable in database.
    *
    * @param  variable  tbox::Pointer to variable.  When assertion checking is active,
    *                   an unrecoverable exception will result if the variable
    *                   pointer is null.
    */
   int mapVariableToIndex(const tbox::Pointer< hier::Variable<DIM> > variable) const;

   /*!
    * Map patch data identifier to variable associated with the data, if
    * possible, and set the variable pointer to the variable in the database.
    *
    * @return  Boolean true if patch data identifier maps to variable in the
    *          database; otherwise false.
    *
    * @param   index  Integer patch data identifier.  When assertion checking
    *                 is active, an unrecoverable exception will if the index
    *                 is invalid (i.e., < 0).
    * @param   variable  tbox::Pointer to variable that maps to patch data identifier
    *                    in database.  If there is no identifier in the database
    *                    matching the index input value, then the variable pointer
    *                    is set to null.
    */
   bool mapIndexToVariable(
      const int index,
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
    *
    * @param os  Optional output stream.  If not given, tbox::plog is used.
    */
   void printClassData(ostream& os = tbox::plog) const;

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
    * Static pointer to singleton active variable database object.
    */
   static LocallyActiveVariableDatabase<DIM>* 
      s_locally_active_variable_database_instance;

   /*
    * Single variable context shared by all locally-active variables.
    */
   tbox::Pointer<hier::VariableContext> d_locally_active_context;

   /* 
    * Locally active variables information held in database:
    * 
    * d_num_variables_registered        number of locally-active vars in database.
    * d_max_locally_active_variable_id  max var instance among vars in database.
    * d_locally_active_variables        array of pointers to vars in database
    *                                   indexed by variable instance identifier.
    * d_locally_active_variable_descriptor_indices  
    *                                   array of integer patch data indices of
    *                                   storage locations for vars in database
    *                                   indexed by variable instance identifier. 
    */
   int d_num_variables_registered;
   int d_max_locally_active_variable_id;
   tbox::Array< tbox::Pointer< hier::Variable<DIM> > > d_locally_active_variables;
   tbox::Array<int> d_locally_active_variable_descriptor_indices;

   /* 
    * Array of manager objects that maintain active variable information for
    * levels in patch hierarchy. Indexing is d_patch_active_data[level #].
    */
   tbox::Array< tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > > 
      d_patch_level_active_data_manager; 

};

}
}

#ifndef DEBUG_NO_INLINE
#include "LocallyActiveVariableDatabase.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "LocallyActiveVariableDatabase.C"
#endif
