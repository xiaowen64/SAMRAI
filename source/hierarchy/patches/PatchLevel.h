//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/patches/PatchLevel.h $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	A collection of patches at one level of the AMR hierarchy
//

#ifndef included_hier_PatchLevel
#define included_hier_PatchLevel

#include "SAMRAI_config.h"
#include "tbox/Arena.h"
#include "tbox/Array.h"
#include "BasePatchLevel.h"
#include "BoxArray.h"
#include "BinaryTree.h"
#include "tbox/Database.h" 
#include "ComponentSelector.h"
#include "GridGeometry.h"
#include "Patch.h"
#include "PatchFactory.h"
#include "tbox/Pointer.h"
#include "ProcessorMapping.h"
#include "BoxGraph.h"
#include "BoxTop.h"
#include "BoxTree.h"
#include "tbox/DescribedClass.h"


#ifndef NULL
#define NULL 0
#endif

namespace SAMRAI {
    namespace hier {

template<int DIM> class PatchLevelIterator;

/**
 * Class PatchLevel<DIM> is a container class for the patches defined
 * at one level of the AMR hierarchy.  The patches in a patch level are
 * distributed across the processors of a parallel machine, so not all
 * patches reside on the local processor.  Each patch is assigned to one
 * and only one processor.
 *
 * To iterate over the local patches in a patch level, use the patch
 * level iterator class (PatchLevel<DIM>::Iterator).  This iterator
 * will use the owner-computes rule to pick off local patches for
 * computation.
 *
 * @see hier::BasePatchLevel
 * @see hier::Patch
 * @see hier::PatchDescriptor
 * @see hier::PatchFactory
 * @see hier::PatchLevelIterator
 */

template<int DIM> class PatchLevel  : public hier::BasePatchLevel<DIM>
{
public:

   /**
    * Default constructor.  PatchLevel must be initialized before it can
    * be used.
    */
   PatchLevel();

   /**
    * Construct a new patch level given an array of boxes and a processor
    * mapping.  The ratio information establishes the ratio between the
    * index space of the new level and some reference level.  Typically,
    * the reference level is level zero (i.e., the coarsest level) in some
    * patch hierarchy.  Also, the ratio is used by the grid geometry object
    * to initialize the geometry information of the level and the patches
    * on the level.  All patch data owned by the local processor are allocated
    * using the specified descriptor.  A patch factory can be passed into
    * the level constructor to specify the appropriate patch factory; if
    * no factory is specified, then the standard patch factory will be used.
    *
    * When assertion checking is active, an unrecoverable assertion results
    * if either the grid geometry pointer or patch descriptor pointer is
    * null, or if the number of boxes in the array does not match the
    * mapping array.
    *
    * defer_boundary_box_creation, if set to true, will suppress the
    * construction of the boundary boxes.
    */
   PatchLevel(
      const BoxArray<DIM>& boxes,
      const ProcessorMapping& mapping,
      const IntVector<DIM>& ratio_to_level_zero,
      const tbox::Pointer< GridGeometry<DIM> > grid_geometry,
      const tbox::Pointer< PatchDescriptor<DIM> > descriptor,
      tbox::Pointer< PatchFactory<DIM> > factory = NULL,
      bool defer_boundary_box_creation = false);

   /**
    * Construct a new patch level from the specified PatchLevel database.  
    * The box, mapping, and ratio to level zero data which are normally
    * passed in during the construction of a new patch level are
    * retrieved from the specified database.  The component_selector
    * argument specifies which patch data components should be allocated
    * and read in from the level_database.  By default, all bits in the
    * component selector are set to false so that no patch data are 
    * allocated.
    *
    * When assertion checking is turned on, the level_database,
    * grid_geometry, and descriptor are checked to make sure that
    * they are not null.
    *
    * defer_boundary_box_creation, if set to true, will suppress the
    * construction of the boundary boxes.
    */

   PatchLevel(   
      tbox::Pointer<tbox::Database> level_database,
      tbox::Pointer< GridGeometry<DIM> > grid_geometry,
      tbox::Pointer< PatchDescriptor<DIM> > descriptor,
      tbox::Pointer< PatchFactory<DIM> > factory,
      ComponentSelector component_selector = 
         ComponentSelector(false),
      bool defer_boundary_box_creation = false);

   /**
    * The virtual destructor for patch level deallocates all patches.
    */
   virtual ~PatchLevel<DIM>();

   /**
    * Return the number of this level in a hierarchy, or the number of
    * a hierarchy level matching the index space of this level. If this
    * level does not align with the index space of a level in the hierarchy,
    * then this value is -1.  When the level is in a hierarchy, the return
    * value os the number of the level in the hierarchy.  See member 
    * function inHierarchy() below.
    */
   int getLevelNumber() const;

   /**
    * Set the number of this level to the level in the hierarchy 
    * aligning with the index space of this level.  The default value
    * is -1 meaning the level index space doesn't align with that of 
    * any hierarchy level.
    */
   void setLevelNumber(const int l);

   /**
    * Return the number of the next coarser level in a hierarchy for
    * the purposes of data interpolation from coarser levels.  If the
    * level is in a hierarchy, then this value is getLevelNumber() - 1.
    * See member function inHierarchy() below.
    */
   int getNextCoarserHierarchyLevelNumber() const;

   /**
    * Set the number of of the next coarser level in a hierarchy for
    * the purposes of data interpolation from coarser levels.  The
    * default is -1 meaning the level doesn't relate to any hierarchy.
    */
   void setNextCoarserHierarchyLevelNumber(const int l);

   /**
    * Return true when this level resides in a hierarchy and false otherwise.
    */
   bool inHierarchy() const; 
   
   /**    
    * Set to true if this level resides in a hierarchy; false otherwise.
    * The default setting is false.  
    */
   void setLevelInHierarchy(bool in_hierarchy);

   /**
    * Return the total number of patches on the level.
    */
   int getNumberOfPatches() const;

   /**
    * Return a pointer to the specified patch on the level.  The patches
    * are numbered starting at zero.
    */
   tbox::Pointer< Patch<DIM> > getPatch(const int p) const;

   /**
    * Return pointer to the patch descriptor for the level.
    */
   tbox::Pointer< PatchDescriptor<DIM> > getPatchDescriptor() const;

   /**
    * Return the factory object used to created patches in the level.
    */
   tbox::Pointer< PatchFactory<DIM> > getPatchFactory() const;

   /**
    * Return the grid geometry description.
    */
   tbox::Pointer< GridGeometry<DIM> > getGridGeometry() const;

   /*
    * Set data members of this patch level by refining information on
    * given coarse level by the given ratio between the two levels.  The
    * fine level will cover the same physical space as the coarse level and
    * will have the same number of patches with the same mapping of those
    * patches to processors.  However, the index space of the level will
    * be refined by the specified ratio.
    *
    * If the fine grid geometry is null (default case), then it is assumed
    * that this level is to use the same grid geometry as the given coarse
    * level and the ratio to level zero is set relative to the given coarse
    * level.  Otherwise, we use the given grid geometry (assumed to be a proper
    * refinement of the grid geometry used on the given coarse level) and copy
    * ratio to level zero from given coarse level.  In other words, the function
    * can be used to produce two different results.  First, when passed a null
    * grid geometry pointer, the refined patch level can be used for data exchange
    * operations with the AMR hierarchy in which the coarse level resides -- both
    * levels are defined with respect to the index space of the grid geometry
    * object which they share.  Thus, the refined patch level can be used in data
    * exchanges with the AMR hierarchy of the coarse level automatically.  Second,
    * when passed a non-null fine grid geometry pointer, the level is defined relative
    * to that geometry and the refined patch level cannot be used in data exchanges 
    * with the AMR hierarchy of the coarse level automatically in general.  This mode 
    * is used to construct a refined copy of an entire patch hierarchy, typically.
    */
   void setRefinedPatchLevel(
      const tbox::Pointer<hier::PatchLevel<DIM> > coarse_level,
      const hier::IntVector<DIM>& refine_ratio,
      const tbox::Pointer<hier::GridGeometry<DIM> > fine_grid_geometry = NULL,
      bool defer_boundary_box_creation = false);

   /**
    * Set data members of this patch level by coarsening information on
    * given fine level by the given ratio between the two levels.  The
    * coarse level will cover the same physical space as the fine level and
    * will have the same number of patches with the same mapping of those
    * patches to processors.  However, the index space of the level will
    * be coarsened by the specified ratio.
    *
    * If the coarse grid geometry is null (default case), then it is assumed
    * that this level is to use the same grid geometry as the given fine
    * level and the ratio to level zero is set relative to the given fine
    * level.  Otherwise, we use the given grid geometry (assumed to be a proper
    * coarsening of the grid geometry used on the given fine level) and copy
    * ratio to level zero from given fine level.  In other words, the function
    * can be used to produce two different results.  First, when passed a null
    * grid geometry pointer, the coarsened patch level can be used for data exchange
    * operations with the AMR hierarchy in which the fine level resides -- both
    * levels are defined with respect to the index space of the grid geometry
    * object which they share.  Thus, the coarsened patch level can be used in data
    * exchanges with the AMR hierarchy of the fine level automatically.  Second,
    * when passed a non-null coarse grid geometry pointer, the level is defined 
    * relative to that geometry and the coarsened patch level cannot be used in 
    * data exchanges with the AMR hierarchy of the fine level automatically in 
    * general.  This mode is used to construct a coarsened copy of an entire patch 
    * hierarchy, typically.
    */
   void setCoarsenedPatchLevel(
      const tbox::Pointer<hier::PatchLevel<DIM> > fine_level,
      const hier::IntVector<DIM>& coarsen_ratio,
      const tbox::Pointer<hier::GridGeometry<DIM> > coarse_grid_geom = NULL,
      bool defer_boundary_box_creation = false);

   /**
    * Create and store the boundary boxes for this level, if they do not
    * already exist.  If they have already been construction, this function
    * does nothing.  If the level is constructed with boundary box
    * creation deferred, this method must be called before any attempt at
    * filling data at physical boundaries.  This function is called from
    * xfer::RefineSchedule<DIM> prior to any physical boundary operations.
    */
   void setBoundaryBoxes();

   /**
    * Return a const reference to the box array that defines 
    * the extent of the index space on the level.
    */
   const BoxArray<DIM>& getPhysicalDomain() const;

   /**
    * Return a const reference to the box array that defines 
    * the patches on the level.
    */
   const BoxArray<DIM>& getBoxes() const;

   /**
    * Return a const reference to the shift array for the patches on the level;
    * the shift array contains a list of shift vectors for each patch
    * when the domain has some periodic direction.
    *
    */
   const tbox::Array< tbox::List< IntVector<DIM> > >& getShiftsForLevel() const;

   /**
    * Return a const reference to the mapping of patches to processors.
    */
   const ProcessorMapping& getProcessorMapping() const;

   /**
    * Return a const reference to the vector ratio between the index
    * space of this patch level and that of a reference level in AMR 
    * hierarchy (typically, level zero).  Specifically, this is the 
    * ratio passed to the constructor.
    */
   const IntVector<DIM>& getRatio() const;

   /**
    * Return the vector ratio between this level and the next coarser
    * level in the patch hierarchy.  This vector is set with the
    * setRatioToCoarserLevel() function.  If the level is not in a
    * hierarchy, a default ratio of zero is returned.
    */
   const IntVector<DIM>& getRatioToCoarserLevel() const;

   /**
    * Set the vector ratio between this level and the next coarser
    * level in the patch hierarchy.  This is required only when level
    * resides in a hierarchy.
    */
   void setRatioToCoarserLevel(const IntVector<DIM>& ratio);

   /**
    * Return the processor that owns the specified patch.  The patches
    * are numbered starting at zero.
    */
   int getMappingForPatch(const int p) const;

   /**
    * Return the box for the specified patch.  The patches are numbered
    * starting at zero.
    */
   const Box<DIM>& getBoxForPatch(const int p) const;

   /**
    * Return the list of valid periodic shifts for the specified patch.
    * The patches are numbered starting at zero.
    */
   const tbox::List< IntVector<DIM> >& getShiftsForPatch(const int p) const;

   /**
    * Return true if patch with given number is adjacent to a non-periodic
    * physical domain boundary.  Otherwise, return false.
    */
   bool patchTouchesRegularBoundary(const int p) const;

   /**
    * Return true if patch with given number is adjacent to a periodic 
    * physical domain boundary.  Otherwise, return false.
    */
   bool patchTouchesPeriodicBoundary(const int p) const;

   /**
    * Allocate the specified component on all patches.  If no memory
    * arena is specified, then the standard memory arena will be used.
    */
   void allocatePatchData(
      const int id,
      const double timestamp = 0.0,
      tbox::Pointer<tbox::Arena> pool = NULL);
 
   /**
    * Allocate the specified components on all patches.  If no memory
    * arena is specified, then the standard memory arena will be used.
    */
   void allocatePatchData(
      const ComponentSelector& components,
      const double timestamp = 0.0,
      tbox::Pointer<tbox::Arena> pool = NULL);
 
   /**
    * Check whether the specified patch data index has been allocated.
    * This function will return true if (1) there are no patches in this
    * patch level or (2) all of the patches have allocated the patch data
    * component.
    */
   bool checkAllocated(const int id) const;

   /**
    * Deallocate the specified component on all patches.  This component 
    * will need to be reallocated before its next use.
    */
   void deallocatePatchData(const int id);
 
   /**
    * Deallocate the specified components on all patches.  These components 
    * will need to be reallocated before their next use.
    */
   void deallocatePatchData(const ComponentSelector& components); 

   /**
    * Set the simulation time for the specified patch component.
    */
   void setTime(const double timestamp, const int id);

   /**
    * Set the simulation time for the specified patch components.
    */
   void setTime(
      const double timestamp,
      const ComponentSelector& components);

   /**
    * Set the simulation time for all allocated patch components.
    */
   void setTime(const double timestamp);

   /**
    * The iterator for the patches on a patch level.  The iterator will
    * return the patches that live on the local processor.  Use iterator
    * PatchLevel<DIM>::Iterator instead of PatchLevelIterator<DIM>,
    * since the iterator may be defined as a nested class in the future.
    */
   typedef PatchLevelIterator<DIM> Iterator;

   /**
    * Uses the PatchLevel database to set the state of the PatchLevel 
    * and to create all patches on the local processor.  
    *
    * Assertions: check that database is a non-null Pointer,
    * that the data being retrieved from the database are of
    * the type expected, that the number of patches is positive, 
    * that the number of patches and size of processor
    * mapping array are the same, and that the number of patches and
    * the number of boxes on the level are equal.
    */
   void getFromDatabase(tbox::Pointer<tbox::Database> database,
                        ComponentSelector component_selector);

   /**
    * Writes the data from the PatchLevel to the database.  
    * Also tells all local patches to write out their state to 
    * the database.  The patchdata_write_table specifies 
    * which patchdata are to be written to the database.
    *
    * Assertions: check that database is a non-null Pointer.
    */
   void putToDatabase(tbox::Pointer<tbox::Database> database,
                      const ComponentSelector& patchdata_write_table);


   /*!
    * @brief Print a patch level to varying details.
    *
    * If depth>0, print function will be called for each patch in the level.
   */
   int recursivePrint( std::ostream &os ,
                       const std::string &border=std::string() ,
                       unsigned short depth=0 );

   /**
    * Returns pointer to a BoxGraph graph object that was constructed
    * using only boxes on this level.  
    */
   tbox::Pointer< BoxGraph<DIM> > getBoxGraph();

   /**
    * Returns pointer to a BoxTop object that was constructed
    * using boxes on this level.
    */
   tbox::Pointer< BoxTop<DIM> > getBoxTop();

   /**
    * Returns pointer to a BoxTree object that was constructed
    * using boxes on this level.
    */
   tbox::Pointer< BoxTree<DIM> > getBoxTree();

   /**
    * Returns pointer to BinaryTree 
    */
   tbox::Pointer< BinaryTree<DIM> > getBinaryTree();   

private:

   /*
    * Private member function that reduces patch boundary information
    * across all (distributed) patches on level.
    */
   void setPatchTouchesBoundaryArrays();

   /*
    * Data mambers that describe contents of every patch level.
    */
   BoxArray<DIM> d_boxes;                  // boxes for all level patches
   ProcessorMapping d_mapping;             // patch mapping to processors
   IntVector<DIM> d_ratio_to_level_zero;   // ratio to reference level
                                           // typically, coarsest level
                                           // in a patch hierarchy
   tbox::Pointer< GridGeometry<DIM> > d_geometry; // grid geometry description
                                                  // used to set patch geometry
   tbox::Pointer< PatchDescriptor<DIM> > d_descriptor;
                                                // patch data info shared by
                                                // all patches in hierarchy
   tbox::Pointer< PatchFactory<DIM> > d_factory;// factory for creating patches
   int d_number_patches;                        // number of patches on level
   BoxArray<DIM> d_physical_domain;             // extent of level index space

   /*
    * The ratio to coarser level applies only when the level resides
    * in a hierarchy.  The level number is that of the hierarchy level
    * that aligns with the index space of the level; if level aligns with
    * no such level then the value is -1 (default value).  The next coarser
    * level number is the next coarser level in the hierarchy for the
    * purposes of filling data from coarser levels.   It is -1 by default
    * but is usually a valid level number more often than level number.
    * The boolean is true when the level is in a hierarchy, false otherwise.
    */
   IntVector<DIM> d_ratio_to_coarser_level;    // ratio to coarser level
   int d_level_number;                          // level number in hierarchy
                                                // aligning with index space 
                                                // of level
   int d_next_coarser_level_number;             // next coarser level number 
                                                // for filling data
   bool d_in_hierarchy;                         // is level in a hierarchy?

    

   /*
    * tbox::Array of distributed patches on level; local patches will be
    * non-null.
    */
   tbox::Array< tbox::Pointer< Patch<DIM> > > d_patches;

   /*
    * Arrays of global boundary information for patches on level.
    */
   tbox::Array<bool> d_patch_touches_regular_boundary;
                                                // boolean value for each patch
                                                // patch indicating whether
                                                // it abuts non-periodic
                                                // physical boundary 
   tbox::Array<bool> d_patch_touches_periodic_boundary;
                                                // boolean value for each
                                                // patch indicating whether
                                                // it abuts periodic physical
                                                // boundary 
   tbox::Array< tbox::List< IntVector<DIM> > > d_shifts;  
                                                // list of shift vectors
                                                // for each when domain
                                                // has some periodic direction

   bool d_boundary_boxes_created;

   /*
    * Used to for various functions, when operating on two
    * levels that are the same, e.g., in 
    * xfer_RefineSchedule<DIM>::generateCommunicationSchedule.
    */
   tbox::Pointer< BoxGraph<DIM> > d_box_graph;
   tbox::Pointer< BoxTop<DIM> > d_box_top;
   tbox::Pointer< BoxTree<DIM> > d_box_tree;
 
   /*
    * Manages communications in clustering algorithms
    */
   tbox::Pointer< BinaryTree<DIM> > d_binary_tree;
};

/**
 * Class PatchLevelIterator<DIM> iterates over the locally owned patches
 * of a patch level.  This is consistent with the standard owner-computes
 * rule.  The iterator should be accessed via PatchLevel<DIM>::Iterator
 * since the implementation may change to be a nested class in the future.
 *
 * @see hier::PatchLevel
 */

template <int DIM> class PatchLevelIterator
{
public:
   /**
    * Default constructor for the patch iterator.  This iterator must be
    * initialized before it can be used to iterate over the patches.
    *
    * @see initialize()
    */
   PatchLevelIterator();

   /**
    * Constructor for the patch level iterator.  The iterator will enumerate
    * the local patches in the patch level belonging to the local processor.
    */
   PatchLevelIterator(const PatchLevel<DIM>& pl);

   /**
    * Constructor for the patch level iterator.  The iterator will enumerate
    * the local patches in the patch level belonging to the local processor.
    */
   PatchLevelIterator(const PatchLevel<DIM>* pl);

   /**
    * Const copy constructor for the iterator.
    */
   PatchLevelIterator(const PatchLevelIterator<DIM>& iterator);

   /**
    * Initializer for the patch level iterator.  The iterator will enumerate
    * the local patches in the patch level belonging to the local processor.
    */
   void initialize(const PatchLevel<DIM>& pl);

   /**
    * Initializer for the patch level iterator.  The iterator will enumerate
    * the local patches in the patch level belonging to the local processor.
    */
   void initialize(const PatchLevel<DIM>* pl);

   /**
    * Assignment operator for the iterator.
    */
   PatchLevelIterator<DIM>&
   operator=(const PatchLevelIterator<DIM>& iterator);

   /**
    * Destructor for the patch level iterator.
    */
   ~PatchLevelIterator<DIM>();

   /**
    * Extract the integer patch index corresponding to the current patch 
    * in the patch level.
    */
   int operator*() const;

   /**
    * Extract the integer patch index corresponding to the current patch 
    * in the patch level.
    */
   int operator()() const;

   /**
    * Return true if the iterator points to a valid patch on the level.
    */
   operator bool() const;

#ifndef LACKS_BOOL_VOID_RESOLUTION
   /**
    * Return non-NULL if the iterator points to a valid patch on the level.
    */
   operator const void*() const;
#endif

   /**
    * Return whether the iterator points to a valid patch in the level.
    * This operator mimics the !p operation applied to a pointer p.
    */
   bool operator!() const;

   /**
    * Increment the iterator to point to the next patch in the level.
    */
   void operator++(int);

   /**
    * Test whether two patch level iterators point to the same patch index.
    */
   bool operator==(const PatchLevelIterator<DIM>& iterator) const;

   /**
    * Test whether two patch level iterators point to different patch indices. 
    */
   bool operator!=(const PatchLevelIterator<DIM>& iterator) const;

private:
   int d_patch;
   int d_number_patches;
   const tbox::Array<int> *d_local_box_indices;
};

}
}

#ifndef DEBUG_NO_INLINE
#include "PatchLevel.I"
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "PatchLevel.C"
#endif

#endif
