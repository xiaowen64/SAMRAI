//
// File:	LocallyActiveDataPatchLevelManager.h
// Package:     SAMRAI hierarchy	
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 468 $
// Modified:	$Date: 2005-06-27 15:48:50 -0700 (Mon, 27 Jun 2005) $
// Description:	Class for managing locally-active data on a single patch level.
//

#ifndef included_hier_LocallyActiveDataPatchLevelManager
#define included_hier_LocallyActiveDataPatchLevelManager

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_hier_ComponentSelector
#include "ComponentSelector.h"
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
#ifndef included_tbox_DescribedClass
#include "tbox/DescribedClass.h"
#endif

#ifndef NULL
#define NULL 0
#endif

namespace SAMRAI {
    namespace hier {

template<int DIM> class LocallyActiveDataPatchLevelIterator;
template<int DIM> class LocallyActiveVariableDatabase;

/*!
 * @brief Class LocallyActiveDataPatchLevelManager is a utility class for
 * managing data on a patch level where each data item may be defined
 * on a different set of patches; i.e., the data is "locally-active".  
 * A separate object of this class is needed for each patch level on 
 * which locally-active data is defined.  Typical usage involves constructing 
 * an instance of this class with a patch level and then defining the active
 * patches for each patch data integer identifier. Then, this class supports
 * various patch level operations asscociated with the locally-active data,
 * such as allocation and deallocation of data and iteration over patches for
 * which a particular data id is active.
 *
 * @see hier::PatchLevel
 * @see hier::ComponentSelector
 * @see hier::LocallyActiveDataPatchLevelIterator
 */

template<int DIM>
class LocallyActiveDataPatchLevelManager 
  : public tbox::DescribedClass
{
   friend class LocallyActiveDataPatchLevelIterator<DIM>;
public:
   /*!
    * Default constructor for LocallyActiveDataPatchLevelManager class. 
    * The object state is invalid, hence the object cannot do anything
    * useful, until it is set using the initialize() member function.
    */
   LocallyActiveDataPatchLevelManager();

   /*!
    * Construct a new LocallyActiveDataPatchLevelManager object and
    * initialize it based on the given patch level reference. 
    * 
    * @param level  const reference to patch level.
    */
   LocallyActiveDataPatchLevelManager(
      const hier::PatchLevel<DIM>& level);

   /*!
    * Construct a new LocallyActiveDataPatchLevelManager object and
    * initialize it based on the given patch level pointer. 
    *
    * @param level  const pointer to patch level.
    *
    * When assertion checking is active, an assertion will result when
    * the patch level pointer is null.
    */
   LocallyActiveDataPatchLevelManager(
      const tbox::Pointer< hier::PatchLevel<DIM> > level);

   /*!
    * Destructor for LocallyActiveDataPatchLevelManager class frees
    * internal storage.
    */
   ~LocallyActiveDataPatchLevelManager();

   /*!
    * An iterator over patches on the patch level.  The iterator will enumerate 
    * the patches that live on the local processor and on which a given patch 
    * data index is active (see constructor for LocallyActiveDataPatchLevelIterator).
    *
    * Use iterator LocallyActiveDataPatchLevelManager::Iterator instead of 
    * LocallyActiveDataPatchLevelIterator<DIM>, since the iterator may 
    * be defined as a nested class in the future.
    */
   typedef LocallyActiveDataPatchLevelIterator<DIM> Iterator;

   /*!
    * Return an iterator that will enumerate the patches on the local 
    * processor and on which the given patch data index is active.
    * 
    * @return iterator.
    *
    * @param patch_data_index integer patch data index.
    *
    * When assertion checking is active, an assertion will result when the
    * manager has not been initialized with a level, or the data index is 
    * invalid (< 0).
    */
   Iterator getIterator(int patch_data_index) const;

   /*!
    * Return an iterator that will enumerate the patches on the local
    * processor and on which data for the given variable is active.
    * Note that we assume that a variable is associated with only one
    * patch data index.
    *
    * @return iterator.
    *
    * @param variable pointer to variable.
    *
    * When assertion checking is active, an assertion will result when the
    * manager has not been initialized with a level, or if the variable
    * is not registered with the locally-active variable database.
    *
    * @param variable const smart pointer to variable.
    */
   Iterator getIterator(
      const tbox::Pointer< hier::Variable<DIM> > variable) const;

   /*!
    * Return pointer to patch level associated with this manager object.
    */
   tbox::Pointer< hier::PatchLevel<DIM> > getPatchLevel() const;

   /*!
    * Return true if argument level is same as that with which this 
    * LocallyActiveDataPatchLevelManager object was initialized; 
    * otherwise return false.
    *
    * @param level  const reference to level.
    */
   bool checkLevel(const hier::PatchLevel<DIM>& level) const;

   /*!
    * Return true if argument level is same as that with which this
    * LocallyActiveDataPatchLevelManager object was initialized; 
    * otherwise return false.
    *
    * @param level  const smart pointer to level.
    */
   bool checkLevel(const tbox::Pointer< hier::PatchLevel<DIM> > level) const;

   /*!
    * Return true if argument level is same as that with which this
    * LocallyActiveDataPatchLevelManager object was initialized; 
    * otherwise return false.
    *
    * @param level  const pointer to level.
    */
   bool checkLevel(const hier::PatchLevel<DIM>* level) const;

   /*!
    * Reset the state of the LocallyActiveDataPatchLevelManager object
    * to that associated with the given level.  If the object was previously 
    * initialized based on a different level, that information is destroyed.
    * and replaced with information from the argument level.  Note that, at 
    * that point, it is impossible to recover the manager state associated with 
    * the previous level via this object.
    *
    * @param level  const reference to patch level.
    */
   void reset(const hier::PatchLevel<DIM>& level);

   /*!
    * Initialize the state of the LocallyActiveDataPatchLevelManager object
    * to that associated with the given level.  If the object was previously
    * initialized based on a different level, that information is destroyed.
    * and replaced with information from the argument level.  Note that, at
    * that point, it is impossible to recover the manager state associated with
    * the previous level via this object.
    *
    * @param level  const smart pointer to patch level.
    *
    * When assertion checking is active, an assertion will result when
    * the patch level pointer is null.
    */
   void reset(const tbox::Pointer< hier::PatchLevel<DIM> > level);

   /*!
    * Check whether given patch data index is active on given patch and
    * return boolean true if data is active on patch; false otherwise.
    *
    * @param patch_data_index integer patch data index.
    * @param patch_number     integer patch number.
    *
    * When assertion checking is active, an assertion will result when the
    * manager has not been initialized with a level, when the patch number 
    * is invalid for the level, or the data index is invalid (< 0).
    */
   bool getPatchDataActive(int patch_data_index,
                           int patch_number) const;

   /*!
    * Return const reference to component selector indicating active/inactive 
    * patch data indices for given patch.
    *
    * @param patch_number     integer patch number.
    *
    * When assertion checking is active, an assertion will result when the
    * manager has not been initialized with a level, or hen the patch number
    * is invalid for the level.
    */
   const hier::ComponentSelector& getAllPatchDataActive(int patch_number) const;

   /*!
    * Set specified patch data active on given patch.  Note that this function
    * does not allocate the corresponding patch data.
    *
    * @param patch_data_index integer patch data index.
    * @param patch_number     integer patch number.
    *
    * When assertion checking is active, an assertion will result when the
    * manager has not been initialized with a level, when the patch number
    * is invalid for the level, or the data index is invalid (< 0).
    */ 
   void setPatchDataActive(int patch_data_index,
                           int patch_number);

   /*!
    * Set patch data active/inactive for given patch based on component
    * selector information.  Note that this function does not allocate/deallocate
    * the corresponding patch data.
    *
    * @param active_indices const reference to component selector containing
    *                       active/inactive patch data index information.
    * @param patch_number   integer patch number.
    *
    * When assertion checking is active, an assertion will result when the 
    * manager has not been initialized with a level or when the patch
    * number is invalid for the level.
    */
   void setPatchDataActive(const hier::ComponentSelector& active_indices,
                           int patch_number);

   /*!
    * Set all patch data active for given patch.  Note that this function does 
    * not allocate the corresponding patch data.
    *
    * @param patch_number   integer patch number.
    *
    * When assertion checking is active, an assertion will result when the 
    * manager has not been initialized with a level or when the patch
    * number is invalid for the level.
    */
   void setAllPatchDataActive(int patch_number);

   /*!
    * Set patch data active/inactive for all patches based on component
    * selector information. Note that this function does not allocate/deallocate
    * the corresponding patch data.
    *
    * @param active_indices const reference to component selector containing
    *                       active/inactive patch data index information.
    *
    * When assertion checking is active, an assertion will result when the
    * manager has not been initialized with a level.
    */
   void setPatchDataActive(const hier::ComponentSelector& active_indices);

   /*!
    * Set specified patch data inactive on given patch.  Note that this 
    * function does not deallocate the corresponding patch data.
    *
    * @param patch_data_index integer patch data index.
    * @param patch_number     integer patch number.
    *
    * When assertion checking is active, an assertion will result when the
    * manager has not been initialized with a level, when the patch number
    * is invalid for the level, or the data index is invalid (< 0).
    */
   void setPatchDataInactive(int patch_data_index,
                             int patch_number);

   /*!
    * Set all patch data inactive for given patch.  Note that this function does
    * not deallocate the corresponding patch data.
    *
    * @param patch_number   integer patch number.
    *
    * When assertion checking is active, an assertion will result when the
    * manager has not been initialized with a level or when the patch
    * number is invalid for the level.
    */
   void setAllPatchDataInactive(int patch_number);

   /*!
    * Set all patch data inactive for all patches.  Note that this function does 
    * not deallocate any patch data.
    *
    * When assertion checking is active, an assertion will result when the
    * manager has not been initialized with a level.
    */
   void setAllPatchDataInactive();

   /*!
    * Clear all information from locally-active data patch level manager object, 
    * setting object state to that created by the default constructor.
    *
    * Note that this function does not deallocate any patch data.
    */
   void clearAllActiveDataInfo();

   /*!
    * Check if data corresponding to given index is allocated on all active
    * patches on the level.
    *
    * @return bool true if data is allocated on all active patches; false otherwise. 
    *              Note that if no patch is active for data, true is returned.
    *
    * @param patch_data_index integer patch data index.
    *
    * When assertion checking is active, an assertion will result when the
    * manager has not been initialized with a level or when the data index 
    * is invalid (< 0).
    */
   bool checkAllocated(int patch_data_index) const;

   /*!
    * Allocate data for given patch data index on all level patches on which
    * the data is active.  Each allocated patch data object will be stamped with
    * the given time value.
    *
    * @param patch_data_index integer patch data index.
    * @param timestamp        optional double data timestamp.
    * @param pool             optional pointer to memory arena for data.
    *
    * When assertion checking is active, an assertion will result when the
    * manager has not been initialized with a level or when the data index
    * is invalid (< 0).
    */
   void allocatePatchData(int patch_data_index,
                          double timestamp = 0.0,
                          tbox::Pointer<tbox::Arena> pool = NULL) const;

   /*!
    * Allocate all active patch data on all level patches associated with this
    * manager object.  Each allocated patch data object will be stamped with 
    * the given time value.  
    *
    * @param timestamp    optional double data timestamp.
    * @param pool         optional pointer to memory arena for data.
    *
    * When assertion checking is active, an assertion will result when the
    * manager has not been initialized with a level.
    */
   void allocateAllPatchData(double timestamp = 0.0,
                             tbox::Pointer<tbox::Arena> pool = NULL) const;

   /*!
    * Allocate all active patch data for given patch.  Each allocated patch data 
    * object will be stamped with the given time value.
    *
    * @param patch_number integer patch number.
    * @param timestamp    optional double data timestamp.
    * @param pool         optional pointer to memory arena for data.
    *
    * When assertion checking is active, an assertion will result when the
    * manager has not been initialized with a level or when the patch
    * number is invalid for the level.
    */
   void allocateAllPatchData(int patch_number,
                             double timestamp = 0.0,
                             tbox::Pointer<tbox::Arena> pool = NULL) const;

   /*!
    * Deallocate data for given patch data index on all level patches on which
    * the data is active.  Note that the state of this manager object remains
    * intact after this operation.
    *
    * @param patch_data_index integer patch data index.
    *
    * When assertion checking is active, an assertion will result when the
    * manager has not been initialized with a level or when the data index
    * is invalid (< 0).
    */
   void deallocatePatchData(int patch_data_index) const;

   /*!
    * Deallocate all active patch data on all level patches associated with this
    * manager object.  Note that the state of this manager object remains
    * intact after this operation.
    *
    * When assertion checking is active, an assertion will result when the
    * manager has not been initialized with a level.
    */
   void deallocateAllPatchData() const;

   /*!
    * Deallocate all active patch data for given patch.  Note that the state of 
    * this manager object remains intact after this operation.
    *
    * @param patch_number integer patch number.
    *
    * When assertion checking is active, an assertion will result when the
    * manager has not been initialized with a level or when the patch
    * number is invalid for the level.
    */
   void deallocateAllPatchData(int patch_number) const;

   /*!
    * Print all active patch data information contained in the locally-active
    * box set to the specified output stream.
    *
    * @param os optional reference to output stream (default is plog).
    */
   void printClassData(ostream& os = tbox::plog) const;

private:
   /*
    * Cached patch level information.
    */ 
   tbox::Pointer< hier::PatchLevel<DIM> > d_patch_level;
   int d_number_patches;

   /*
    * Array of component selectors (indexed by patch number), each of which 
    * holds active patch data index information for corresponding patch.
    */
   tbox::Array< tbox::Pointer< hier::ComponentSelector > > d_active_patch_data;

};

/*!
 * Class hier::LocallyActiveDataPatchLevelIterator iterates over the locally-owned
 * (i.e., on processor) patches of a patch level on which data for a given variable
 * or patch data index is active (i.e., can exist).  This is consistent with the 
 * standard owner-computes rule and mimicks the behavior of the SAMRAI 
 * hier::PatchLevelIterator class.  The iterator should be declared in user code 
 * as LocallyActiveDataPatchLevelManager<DIM>::Iterator since the implementation 
 * may change to be a nested class in the future.  Also, since this class is not
 * set up for general usage, it is recommended that an iterator object be obtained
 * using one of the getIterator() functions in the LocallyActiveDataPatchLevelManager
 * class.
 *
 * @see hier::PatchLevel
 * @see hier::LocallyActiveDataPatchLevelManager
 */

template <int DIM>
class LocallyActiveDataPatchLevelIterator
{
public:
   /*!
    * Default constructor for the locally-active data patch iterator.  This
    * iterator must be initialized before it can be used to iterate over patches
    * on a level.  
    * 
    * @see initialize()
    */
   LocallyActiveDataPatchLevelIterator();

   /*!
    * Constructor for the locally-active data patch iterator.  The iterator 
    * will enumerate the local patches in the patch level belonging
    * to the local processor on which data for the given variable is active.
    *
    * @param variable  smart pointer to variable.
    * @param pl        const reference to patch level.
    */
   LocallyActiveDataPatchLevelIterator(
      const tbox::Pointer< hier::Variable<DIM> > variable,
      const hier::PatchLevel<DIM>& pl);

   /*!
    * Constructor for the locally-active data patch iterator.  The iterator 
    * will enumerate the local patches in the patch level belonging
    * to the local processor on which data for the given variable is active.
    *
    * @param variable  smart pointer to variable.
    * @param pl        const pointer to patch level.
    */
   LocallyActiveDataPatchLevelIterator(
      const tbox::Pointer< hier::Variable<DIM> > variable,
      const hier::PatchLevel<DIM>* pl);

   /*!
    * Constructor for the locally-active data patch iterator.  The iterator
    * will enumerate the local patches in the patch level belonging
    * to the local processor on which the give patch data index is active.
    * 
    * Note that this is a very special constructor used by the 
    * LocallyActiveDataPatchLevelManager class.
    *
    * @param variable             smart pointer to variable.
    * @param pl                   const pointer to patch level.
    * @param active_data_indices  const array of pointers to component selectors
    *                             describing active data on patch level.
    */
   LocallyActiveDataPatchLevelIterator(
      int patch_data_index,
      const hier::PatchLevel<DIM>* pl,
      const tbox::Array< tbox::Pointer<hier::ComponentSelector> >* active_data_indices);

   /*!
    * Const copy constructor for the locally-active data patch iterator.
    */
   LocallyActiveDataPatchLevelIterator(
      const LocallyActiveDataPatchLevelIterator<DIM>& iterator);

   /*!
    * Initializer for the locally-active data patch iterator.  The iterator
    * will enumerate the local patches in the patch level belonging to the
    * local processor on which data for the given variable is active.
    *
    * @param variable             const smart pointer to variable.
    * @param pl                   const reference to patch level.
    */
   void initialize(const tbox::Pointer< hier::Variable<DIM> > variable,
                   const hier::PatchLevel<DIM>& pl);

   /*!
    * Initializer for the locally-active data patch iterator.  The iterator
    * will enumerate the local patches in the patch level belonging to the
    * local processor on which data for the given variable is active.
    *
    * @param variable             const smart pointer to variable.
    * @param pl                   const pointer to patch level.
    */
   void initialize(const tbox::Pointer< hier::Variable<DIM> > variable,
                   const hier::PatchLevel<DIM>* pl);

   /*!
    * Assignment operator for the iterator sets calling object to state
    * of the argument iterator.
    */
   LocallyActiveDataPatchLevelIterator<DIM>&
      operator=(const LocallyActiveDataPatchLevelIterator<DIM>& iterator);

   /*!
    * Destructor for the iterator releases all internal storage.
    */
   ~LocallyActiveDataPatchLevelIterator<DIM>();

   /*!
    * Extract the integer patch index corresponding to the current patch in 
    * the patch level.
    */
   int operator*() const;

   /*!
    * Extract the integer patch index corresponding to the current patch in 
    * the patch level.
    */
   int operator()() const;

   /*!
    * Return true if the iterator points to a valid patch on the level; i.e.,
    * patch exists on the level and variable with which iterator is initialized 
    * is active on that patch.
    */
   operator bool() const;

#ifndef LACKS_BOOL_VOID_RESOLUTION
   /*!
    * Return non-NULL if the iterator points to a valid patch on the level;
    * i.e., patch exists on the level and variable with which iterator is 
    * initialized is active on that patch.
    */
   operator const void*() const;
#endif

   /*!
    * Return true if the iterator points to a valid patch in the level (patch exists
    * on the level and variable with which iterator is initialized is active on that
    * patch); false otherwise.  This operator mimics the !p operation applied to 
    * a pointer p.
    */
   bool operator!() const;

   /*!
    * Increment the iterator to point to the next local patch on which the variable
    * with which iterator is initialized is active on the level.
    */
   void operator++(int);

   /*!
    * Test whether two iterators point to the same patch index.
    */
   bool operator==(const LocallyActiveDataPatchLevelIterator<DIM>& iterator) const;

   /*!
    * Test whether two iterators point to different patch indices.
    */
   bool operator!=(const LocallyActiveDataPatchLevelIterator& iterator) const;

private:
   /*
    * Private member function used to skip over inactive patches.
    */
   bool notActivePatch(int patch_number) const;

   int d_patch;
   int d_data_index;
   const tbox::Array< tbox::Pointer<hier::ComponentSelector> >* d_active_patch_info;
   int d_number_patches;
   const hier::ProcessorMapping* d_mapping;

   static hier::LocallyActiveVariableDatabase<DIM>* s_variable_database;
};

}
}
#ifndef DEBUG_NO_INLINE
#include "LocallyActiveDataPatchLevelManager.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "LocallyActiveDataPatchLevelManager.C"
#endif
