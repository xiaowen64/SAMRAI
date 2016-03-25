//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/variables/VariableDatabase.h $
// Package:     SAMRAI hierarchy	
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Singleton database class for managing variables and contexts.
//

#ifndef included_hier_VariableDatabase
#define included_hier_VariableDatabase

#include "SAMRAI_config.h"

#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif

#include "ComponentSelector.h"
#include "IntVector.h"
#include "PatchDescriptor.h"
#include "Variable.h"
#include "VariableContext.h"
#include "tbox/Array.h"
#include "tbox/List.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#ifndef included_String
#include <std::string>
#define included_String
#endif
#include "tbox/DescribedClass.h"


namespace SAMRAI {
    namespace hier {

// forward declarations for types used in this file
template<int DIM> class LocallyActiveVariableDatabase;

/*!
 * @brief Class VariableDatabase<DIM> is a Singleton class that manages 
 * mapping information between variables and patch data objects in ways
 * that are configurable to the needs of an application built using SAMRAI.
 * In SAMRAI, a single patch descriptor object is shared by all patches in 
 * an SAMR hierarchy.  This database generates and maintain mappings 
 * (for lookup) between integer patch data indices and either 
 * variable-context pairs or just variables when contexts are not used.  
 * Users define which variables are included in the database, which variable  
 * contexts are in the database, and which contexts are associated with which 
 * variables via registration operations described below.
 *
 * Since the database is a Singleton class, it serves as a globally accessible 
 * container for bookkeeping information related to variable storage on an 
 * SAMR patch hierarchy.  For example, a time-dependent problem may require 
 * distinct "OLD" and "NEW" storage for some quantities.  A context is a name 
 * describing some storage for a collection of variables in a given problem.  
 * Each variable and context pair maps to a unique integer patch data index.  
 * This integer indices are the same for every patch in the hierarchy and
 * so data objects are accessed identically for all patches.
 *
 * Typically, classes that contains numerical routines for patches access
 * the appropriate data in terms of variables.  Classes that provide 
 * solution algorithms manage variable storage in terms of storage "contexts",
 * such as "OLD" and "NEW".   Typical usage of the database for managing 
 * variable storage bookkeeping in this way is as follows:
 * - #   Add variables to the database via the addVariable() function.  
 *       Generally, this is not be needed if the registerVariableAndContext() 
 *       function is used.  However, it is included to deal with complexity 
 *       associated with multiple integrators that share variable data.  In
 *       other words, it circumvents potential problems associated with 
 *       multiple objects owning a single variable.
 * - #   Get variable contexts via the getContext() function.
 * - #   Define variable context pairs in the database using the function
 *       registerVariableAndContext().  This routine returns valid patch
 *       data indices that can be used to manage storage on the
 *       SAMR patch hierarchy.
 *
 * Alternatively, the database can be used to maintain a map between a 
 * single patch data index and a variable that created the patch data factory 
 * at that index.  Such a mapping can be constructed using the functions
 * registerVariablePatchDataIndex() and registerClonedVariablePatchDataIndex().
 * These operations may be helpful when managing variables and contexts is 
 * inconvenient.  However, in contrast to the functions above, these functions do 
 * not allow using a single variable to define multiple patch data indices. 
 * One may use this functionality as follows:
 * - #   Add variable and index pair to the database via the 
 *       addVariablePatchDataIndex() function.   Or, clone an existing
 *       variable index pair via registerClonedVariablePatchDataIndex(). 
 *       Either of these functions will check to make sure that the type 
 *       of the variable matches the patch data index.
 * - #   Get variable for the index via the mapIndexToVariable() function.
 * 
 * Note that the database can be used to allocate and free patch data 
 * indices using the registerClonedVariablePatchDataIndex() and 
 * removeVariablePatchDataIndex() functions.  The registerVariableAndContext()
 * function will also create a new patch data index when needed.  In general,
 * variables and contexts persist once create and cannot be removed from the
 * database.
 *
 * Finally, this database is used in SAMRAI to manage the reading/writing of
 * patch data objects from/to restart files.
 *
 * \verbatim
 * Important notes:
 * 
 *    (1) Having two different variables or two different contexts in 
 *        the database with the same name string identifier is not allowed.
 *        This eliminates the need to resolve ambiguities when querying the
 *        database for variables or contexts based on string identifiers.
 *        The error checking in this class implementation will catch
 *        any attempt to improperly register variables and contexts and
 *        report a descriptive log message.
 * 
 *    (2) It is STRONGLY recommended that variable contexts be generated by 
 *        using the getContext() function in this class.  While contexts 
 *        may be created by other means and passed into the member functions 
 *        of this database class, doing so may produce unexpected results 
 *        due to potentially non-unique integer and name string context 
 *        identifiers.
 *
 * \endverbatim
 *
 * @see hier::VariableContext
 * @see hier::Patch
 * @see hier::PatchDescriptor
 */

template<int DIM> class VariableDatabase 
{
public:
   friend class LocallyActiveVariableDatabase<DIM>;

   /*!
    * Return a pointer to the singleton instance of the variable database.  
    * All access to the VariableDatabase<DIM> object is through the 
    * getDatabase() function.  For example, to access the context with 
    * string name "my_context", use the following call: 
    * VariableDatabase<DIM>::getDatabase()->getContext("my_context").
    *
    * Note that when the database is accessed for the first time, the
    * Singleton instance is registered with the ShutdownRegistry
    * class which destroys such objects at program completion.  Thus,
    * users of this class do not explicitly allocate or deallocate the
    * Singleton instance.
    * 
    * @return  Bare pointer to variable database instance.
    */
   static VariableDatabase<DIM>* getDatabase();

   /*!
    * Deallocate the Singleton VariableDatabase<DIM> instance.  
    * It is not necessary to call this function at program termination, 
    * since it is automatically called by the ShutdownRegistry class.
    */
   static void freeDatabase();

   /*!
    * Return number of patch data indices registered with the database.
    */
   virtual
   int getNumberOfRegisteredPatchDataIndices() const;

   /*!
    * Return number of variable contexts registered with the database.
    */
   virtual
   int getNumberOfRegisteredVariableContexts() const;

   /*!
    * Return pointer to the patch descriptor managed by the database 
    * (and shared by all patches in an SAMR hierarchy).
    *
    * @return  tbox::Pointer to patch descriptor instance.
    */
   virtual 
   tbox::Pointer< hier::PatchDescriptor<DIM> > getPatchDescriptor() const;

   /*!
    * Return a variable context object with the given name string identifier.
    * If a context whose name matches the argument string exists in the
    * database, it is returned.  Otherwise, a new context is created 
    * and returned.  Note that it is impossible to add two distinct 
    * variable context objects to the database with the same name.
    *
    * If the string identifier is empty, the context pointer 
    * return value will be null.  Thus, you cannot create a variable 
    * context with an empty string identifier using this function.
    *
    * @return  tbox::Pointer to variable context.
    *
    * @param name  Const reference to name string identifying the context.
    */
   virtual 
   tbox::Pointer<hier::VariableContext> getContext(const std::string& name);

   /*!
    * Check whether context with given name exists in the database.  
    *
    * @return  Boolean true if context with name exists in database; 
    *          otherwise, false.
    *
    * @param name  Const reference to name string identifying the context.
    */
   virtual
   bool checkContextExists(const std::string& name) const;

   /*!
    * Add the given variable to the database.  This function checks 
    * whether the variable already exists in the database, or whether
    * a different variable exists with the same name string identifier.
    * If the variable does not already exist in the database and its name
    * string is distinct from all variables in the database, the variable
    * is added to the database.  Having two different variables in the 
    * database with the same name is not allowed.  If this is attempted,
    * an error message will be logged and the program will abort.  If the 
    * variable already exists in the database, the routine essentially
    * does nothing.  
    *
    * @param variable tbox::Pointer to variable.  When assertion checking 
    *                 is active, an assertion will result if the 
    *                 variable pointer is null.
    */
   virtual
   void addVariable(const tbox::Pointer< hier::Variable<DIM> > variable); 

   /*! 
    * Get variable in database with given name string identifier.  
    *
    * @return  tbox::Pointer to variable in database with given name.  
    *          If no such variable exists, a null pointer is returned.
    *
    * @param name  Const reference to name string identifying the variable. 
    */
   virtual
   tbox::Pointer< hier::Variable<DIM> > getVariable(const std::string& name) const;

   /*!
    * Check whether variable with given name exists in the database.
    *
    * @return  Boolean true if variable with name exists in database;
    *          otherwise, false.
    *
    * @param name  Const reference to name string identifying the variable.
    */
   virtual
   bool checkVariableExists(const std::string& name) const; 

   /*!
    * Create a new patch data index by cloning the data factory
    * for the given variable at the given index and return the patch data
    * index of the new data location.  The new index and variable pair is 
    * added to the variable database.  A variable-patch data index pair 
    * generated using this function cannot be looked up using a 
    * VariableContext.  Note that this function does not deallocate  
    * any patch data storage associated with the new patch data index.
    *
    * @return New integer patch data index. If new patch data not added, 
    *         return value is an invalid patch data index (< 0).
    *
    * @param  variable tbox::Pointer to variable.  If the variable is 
    *         unknown to the database, then an invalid patch data index 
    *         (< 0) will be returned. When assertion checking is active,  
    *         an assertion will result when the variable pointer is null.
    * @param  old_id Integer patch data index currently associated with variable.
    *         If this value is not a valid patch data index (e.g., < 0) or 
    *         does not map to patch data matching the type of the given 
    *         variable, the program will abort with an error message.
    */
   virtual 
   int registerClonedPatchDataIndex(
      const tbox::Pointer< hier::Variable<DIM> > variable,
      int old_id);

   /*!
    * Add given patch data index and variable pair to the database
    * if not already there.  This registration function is primarily
    * intended for variables (i.e., not internal SAMRAI variables) that
    * are not associated with a VariableContext and for which a patch data
    * index is already known. 
    *
    * If the index is unspecified (default case), the default variable 
    * factory is cloned and the variable and new index are added to the 
    * database; in this case, the patch data will have the default number 
    * of ghost cells associated with the given variable (defined
    * by the patch data factory it generates).  Also, a variable-patch data 
    * index pair generated with this function cannot be looked up using a 
    * VariableContext.  Note that this function does not allocate
    * any patch data storage associated with the integer index.
    *
    * NOTE: This function must not be used by SAMRAI developers for 
    * creating patch data indices for internal SAMRAI variables.  The
    * routine registerInternalSAMRAIVariable() must be used for that 
    * case. 
    *
    * @return New integer patch data index.  If new patch data index not 
    *         added, return value is an invalid patch data index (< 0).
    *
    * @param  variable tbox::Pointer to variable.  When assertion checking is 
    *                  active, an assertion will result when the variable 
    *                  pointer is null.
    * @param  data_id  Optional integer patch data index to be added 
    *                  (along with variable) to the database.  If the value 
    *                  is unspecified (default case), the default variable 
    *                  patch data factory is used to generate a new factory. 
    *                  If the value is provided and does not map to patch 
    *                  data matching the type of th given variable, the 
    *                  program will abort with an error message.
    */
   virtual
   int registerPatchDataIndex(
       const tbox::Pointer< hier::Variable<DIM> > variable,
       int data_id = -1);

   /*!
    * Remove the given patch index from the variable database if it exists in 
    * the database.  Also, remove the given index from the patch descriptor 
    * and remove any mapping between the index and a variable from the 
    * variable database.  Note that this function does not deallocate 
    * any patch data storage associated with the integer index.
    * 
    * @param  data_id  Integer patch data index to be removed from
    *                  the database.  When assertion checking is active, 
    *                  an assertion will result when the patch data index 
    *                  is invalid (i.e., < 0).
    */
   virtual
   void removePatchDataIndex(int data_id);

   /*!
    * Check whether the given variable is mapped to the given patch data 
    * index in the database.
    *
    * @return  Boolean true if the variable is mapped the given patch 
    * data index; false otherwise.  
    * 
    * @param  variable  tbox::Pointer to variable.  When assertion checking 
    *                   is active, an unrecoverable assertion will result if 
    *                   the variable pointer is null.
    * @param  data_id   Integer patch data index.  When assertion checking 
    *                   is active, an unrecoverable assertion will result if
    *                   the value is an invalid identifier (either < 0 or
    *                   larger than the maximum allowed patch data id).
    */
   virtual
   bool checkVariablePatchDataIndex(
        const tbox::Pointer< hier::Variable<DIM> > variable,
        int data_id) const;

   /*!
    * Check whether the given variable matches the patch data type 
    * associated with the given patch data index in the database.
    *
    * @return  Boolean true if the type of the variable matches the type of
    *          the patch data at the given patch data index; false otherwise.  
    *
    * @param  variable  tbox::Pointer to variable.  When assertion checking
    *                   is active, an unrecoverable assertion will result if
    *                   the variable pointer is null.
    * @param  data_id   Integer patch data index.  When assertion checking
    *                   is active, an unrecoverable assertion will result if
    *                   the value is an invalid identifier (either < 0 or
    *                   larger than the maximum allowed patch data id).
    */
   virtual
   bool checkVariablePatchDataIndexType(
        const tbox::Pointer< hier::Variable<DIM> > variable,
        int data_id) const;

   /*!
    * Register variable and context pair along with the number of ghost 
    * cells required for storage with the variable database.  Typically, 
    * this function will generate a new patch data index for the
    * variable and ghost cell width and add the variable-context pair
    * and index to the database.  If the variable-context pair is already
    * mapped to some patch data index in the database, then that index
    * will be returned and the function will do nothing.   However, if
    * the variable-context pair is already mapped to some patch data index 
    * with a different ghost cell width, the program will abort with a 
    * descriptive error message.   When the number of ghost cells is not 
    * provided (default case), a default value of zero ghosts will be used.  
    *
    * If either the variable or the context is unknown to the database 
    * prior to calling this routine, both items will be added to the
    * database, if possible.  The constraints for the getContext() and
    * addVariable() routines apply.  
    * 
    * @return Integer patch data index of variable-context pair in database.
    *
    * @param  variable  tbox::Pointer to variable.  When assertion checking 
    *                   is active, an unrecoverable assertion will result if 
    *                   the variable pointer is null.
    * @param context    tbox::Pointer to variable context.  When assertion 
    *                   checking is active, an unrecoverable assertion 
    *                   will result if the context pointer is null.
    * @param ghosts     Optional ghost cell width for patch data associated 
    *                   with variable-context pair.  If the ghost cell width 
    *                   is given, all entries of the vector must be >= 0.  
    *                   When assertion checking is active, an unrecoverable 
    *                   assertion will result if the ghost width vector 
    *                   contains a negative entry.
    */
   virtual
   int registerVariableAndContext(
      const tbox::Pointer< hier::Variable<DIM> > variable, 
      const tbox::Pointer<hier::VariableContext> context,
      const hier::IntVector<DIM>& ghosts = hier::IntVector<DIM>(0));

   /*!
    * Map variable-context pair in database to patch data index. 
    * If there is no such pair in the database (either the variable does 
    * not exist, the context does not exist, or the pair has not been 
    * registered), then an invalid patch data index (i.e., < 0) is returned. 
    * Note that for this function to operate as expected, the database mapping
    * information must have been generated using the 
    * registerVariableAndContext() function.  If the variable was 
    * registered without a variable context, then the patch data index 
    * associated with the variable will not be returned.  See the other
    * map...() functions declared in this class. 
    *
    * @return Integer patch data index of variable-context pair in database.
    *         If the variable-context pair was not regisetered with the 
    *         database, then an invalid data index (< 0) will be returned.
    *   
    * @param  variable  tbox::Pointer to variable.  When assertion checking 
    *                   is active, an unrecoverable assertion will result 
    *                   if the variable pointer is null.
    * @param  context   tbox::Pointer to variable context.  When assertion 
    *                   checking is active, an unrecoverable assertion 
    *                   will result if the variable context pointer is null.
    */
   virtual 
   int mapVariableAndContextToIndex(
      const tbox::Pointer< hier::Variable<DIM> > variable,
      const tbox::Pointer<hier::VariableContext> context) const;

   /*!
    * Map patch data index to variable associated with the data, if
    * possible, and set the variable pointer to the variable in the database.
    * 
    * @return  Boolean true if patch data index maps to variable in the
    *          database; otherwise false.
    *
    * @param   index  Integer patch data index.  When assertion checking 
    *                 is active, an unrecoverable assertion will if the index 
    *                 is invalid (i.e., < 0).
    * @param   variable  tbox::Pointer to variable that maps to patch data 
    *                    index in database.  If there is no index in the 
    *                    database matching the index input value, then the 
    *                    variable pointer is set to null.
    */
   virtual
   bool mapIndexToVariable(
      const int index,
      tbox::Pointer< hier::Variable<DIM> >& variable) const;

   /*!
    * Map patch data index to variable-context pair associated with
    * the data, if possible, and set the variable and context pointers to
    * the corresponding database entries.  Note that for this function 
    * to operate as expected, the database mapping information must 
    * have been generated using the registerVariableAndContext() function.
    * If the variable was registered without a variable context, then 
    * the variable and variable context returned may not be what is
    * expected by the user; e.g., they may be associated with internal
    * SAMRAI variables.
    *
    * @return  Boolean true if patch data index maps to variable-context
    *          pair in the database; otherwise false.
    * @param   index patch data index
    * @param   variable tbox::Pointer to variable set to matching variable 
    *          in database.  If no match is found, it is set to null.
    * @param   context  tbox::Pointer to variable context set to matching 
    *          variable context in database. If no match is found, it is 
    *          set to null.
    */
   virtual 
   bool mapIndexToVariableAndContext(
      const int index,
      tbox::Pointer< hier::Variable<DIM> >& variable,  
      tbox::Pointer<hier::VariableContext>& context) const;  

   /*!
    * Return copy of component selector that holds information about which 
    * patch data entries are written to restart.
    * 
    * @return Component selector describing patch data items registered 
    *         for restart.  That is, the flags set in the component
    *         selector will correspond to the patch data indices 
    *         that have been registered for restart.
    */
   virtual 
   hier::ComponentSelector getPatchDataRestartTable() const;

   /*!
    * Check whether given patch data index is registered with the database
    * for restart.
    * 
    * @return Boolean true if the patch data with the given index 
    *         is registered for restart; otherwise false.
    * 
    * @param  index  Integer patch data index to check.
    */
   virtual 
   bool isPatchDataRegisteredForRestart(int index) const;

   /*!
    * Register the given patch data index for restart.
    * 
    * @param  index  Integer patch data index to set.
    */
   virtual 
   void registerPatchDataForRestart(int index);

   /*!
    * Unregister the given patch data index for restart.
    *
    * @param  index  Integer patch data index to unset.
    */
   virtual 
   void unregisterPatchDataForRestart(int index);


   /*!
    * Print variable, context, and patch descriptor information 
    * contained in the database to the specified output stream.
    * 
    * @param os  Optional output stream.  If not given, tbox::plog is used.
    * @param print_only_user_defined_variables Optional boolean value 
    *        indicating whether to print information for all variables 
    *        in database or only those that are associated with user-
    *        defined quantities; i.e., not internal SAMRAI variables.
    *        The default is true, indicating that only user-defined 
    *        variable information will be printed.
    */
   virtual void printClassData(std::ostream& os = tbox::plog,
                       bool print_only_user_defined_variables = true) const;

   /*!
    * Register internal SAMRAI variable and number of ghost cells 
    * with the variable database.
    *
    * This function will generate a new patch data index for the variable 
    * and ghost cell width unless the variable is already mapped to some 
    * patch data index in the database with a different ghost cell width
    * or as a user-defined variable.  If the variable is unknown to the 
    * database prior to calling this routine, it will be added to the 
    * database.  
    * 
    * Note that this routine is intended for managing internal SAMRAI 
    * work variables that are typically unseen by users.  It should not be 
    * called by users for registering variables, or within SAMRAI for 
    * registering any user-defined variables with the variable database.  
    * This function enforces the same constrants on variable registration 
    * that are applied for registering user-defined variables; e.g., using
    * the routine registerVariableAndContext().  Thus, it is the 
    * responsibility of SAMRAI developers to avoid naming and ghost cell 
    * width conflicts with user-defined variables or other internal 
    * SAMRAI variables.
    * 
    * @return Integer patch data index of variable-ghost cell width pair 
    * in database.
    *
    * @param  variable  tbox::Pointer to variable.  When assertion checking 
    *                   is active, an unrecoverable assertion will result 
    *                   if the variable pointer is null.
    * @param ghosts     Ghost cell width for patch data associated with the
    *                   variable.  All entries of the vector must be >= 0.  
    *                   When assertion checking is active, an unrecoverable 
    *                   assertion results if the vector contains a negative 
    *                   entry.
    */
   virtual
   int registerInternalSAMRAIVariable(
      const tbox::Pointer< hier::Variable<DIM> > variable, 
      const hier::IntVector<DIM>& ghosts);

   /*!
    * Remove the given index from the variable database if it exists in the
    * database and is associated with an internal SAMRAI variable registered
    * with the function registerInternalSAMRAIVariable().  Also, remove 
    * the given index from the patch descriptor and remove any mapping 
    * between the index and a variable from the variable database.  Note 
    * that this function does not deallocate any patch data storage 
    * associated with the integer index.
    * 
    * Note that this routine is intended for managing internal SAMRAI 
    * work variables that are typically unseen by users.  It should not be 
    * called by users for removing patch data indices, or within SAMRAI for 
    * removing any patch data indices associated with user-defined variables.
    * 
    * Note that the given index will not be removed if is not associated with 
    * an internal SAMRAI variable in the variable database; i.e., a 
    * user-defined variable.
    *
    * @param  data_id  Integer patch data identifier to be removed from
    *                  the database.
    */
   virtual
   void removeInternalSAMRAIVariablePatchDataIndex(
      int data_id);

protected:
   /**
    * The constructor for VariableDatabase<DIM> is protected. 
    * Consistent with the definition of a Singleton class, only the 
    * database object has access to the constructor for the class. 
    *
    * The constructor initializes the state of database contents.
    */
   VariableDatabase();

   /**
    * The destructor for VariableDatabase<DIM> is protected. See the
    * comments for the constructor.
    *
    * The destructor deallocates database contents.
    */
   virtual ~VariableDatabase();

   /**
    * Return integer value used to indicate undefined variable or 
    * context identifier.  This routine is protected to allow subclasses
    * to be consistent with this database class.
    */
   int idUndefined() const;

   /**
    * Return integer identifier for first variable found matching given 
    * string name identifier, or return an undefined integer id if no such 
    * variable exists in the database.
    */ 
   int getVariableId(const std::string& name) const;

   /**
    * Initialize Singleton instance with instance of subclass.  This function
    * is used to make the singleton object unique when inheriting from this
    * base class.
    */ 
   void registerSingletonSubclassInstance(
      hier::VariableDatabase<DIM>* subclass_instance);

private:

   /*
    * Private member functions to search for indices of contexts 
    * given a string identifier, and to add a new context
    * to database.
    */

   int getContextId_Private(const std::string& name) const;  
   void addContext_Private(const tbox::Pointer<hier::VariableContext> context); 

   /*
    * Private member function to add variable to database (either 
    * user-defined or internal SAMRAI variable depending on boolean 
    * argument).  Boolean return value is true if variable is either 
    * added to variable database or already found in variable database. 
    * Boolean return value is false only when a user variable has 
    * the same name string identifier of another variable already in 
    * the database.  In this case, the variable is not added to the database.
    */

   bool addVariable_Private(
      const tbox::Pointer< hier::Variable<DIM> > variable,
      bool user_variable);

   /*
    * Private member function to add variable-patch data index mapping
    * to the database.  Boolean indicates whether variable is a user-defined
    * variable or an internal SAMRAI work variable.
    */

   void addVariablePatchDataIndexPairToDatabase_Private(
      const tbox::Pointer< hier::Variable<DIM> > variable,
      int data_id,
      bool user_variable);

   /*
    * Private member function to add variable-context pair to the database.  
    * Boolean indicates whether variable is a user-defined variable or 
    * an internal SAMRAI work variable.
    */

   int registerVariableAndContext_Private(
      const tbox::Pointer< hier::Variable<DIM> > variable,
      const tbox::Pointer<hier::VariableContext> context,
      const hier::IntVector<DIM>& ghosts,
      bool user_variable);

   /*
    * Static data members used to control access to and destruction of
    * singleton variable database instance.
    */

   static VariableDatabase<DIM>* s_variable_database_instance;
   static bool s_registered_callback;

   /*
    * Static data members used to control allocation of arrays.
    */

   static int s_context_array_alloc_size;
   static int s_variable_array_alloc_size;
   static int s_descriptor_array_alloc_size;

   /*
    * Data members used to store variable, context, patch data id information.
    */

   tbox::Pointer< hier::PatchDescriptor<DIM> > d_patch_descriptor;

   tbox::Pointer<hier::VariableContext> d_internal_SAMRAI_context;

   int d_num_registered_patch_data_ids;

   // Array of VariableContext pointers is indexed as 
   // d_contexts[ <context id> ]
   int d_max_context_id;
   tbox::Array< tbox::Pointer<hier::VariableContext> > d_contexts;

   // Array of Variable pointers is indexed as d_variables[ <variable id> ]
   int d_max_variable_id;
   tbox::Array< tbox::Pointer< hier::Variable<DIM> > > d_variables;

   // Array of VariableContext to patch descriptor indices is indexed as 
   // d_variable_context2index_map[ <context id> ]
   tbox::Array< tbox::Array<int> > d_variable_context2index_map;

   // Array of patch descriptor indices to Variables is indexed as 
   // d_index2variable_map[ <descriptor id> ]
   int d_max_descriptor_id;
   tbox::Array< tbox::Pointer< hier::Variable<DIM> > > d_index2variable_map;

   // Array of user variable booleans is indexed as 
   // d_is_user_variable[ <variable id> ]
   tbox::Array<bool> d_is_user_variable;

   /*
    * tbox::Array of bits that determine which pieces of patch data
    * need to be written to the restart database.  The
    * bit in position j corresponds to the patch data stored in
    * the j-th index of the patch descriptor object.
    */

   hier::ComponentSelector d_patchdata_restart_table;

};

}
}
#ifndef DEBUG_NO_INLINE
#include "VariableDatabase.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "VariableDatabase.C"
#endif
