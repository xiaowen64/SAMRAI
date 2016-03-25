//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/mesh/gridding/GriddingAlgorithm.h $
// Package:     SAMRAI mesh
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2148 $
// Modified:    $LastChangedDate: 2008-04-24 10:04:24 -0700 (Thu, 24 Apr 2008) $
// Description: AMR hierarchy generation and regridding routines.
//

#ifndef included_mesh_GriddingAlgorithm
#define included_mesh_GriddingAlgorithm

#include "SAMRAI_config.h"
#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif
#include "BaseGriddingAlgorithm.h"
#include "TagAndInitializeStrategy.h"
#include "BoxArray.h"
#include "BoxIOUtility.h"
#include "BoxList.h"
#include "IntVector.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "ProcessorMapping.h"
#include "BoxGeneratorStrategy.h"
#include "LoadBalanceStrategy.h"
#include "CellVariable.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Serializable.h"
#ifndef included_String
#include <string>
#define included_String
#endif
#include "tbox/Timer.h"
#include "RefineAlgorithm.h"
#include "RefineSchedule.h"



namespace SAMRAI {
   namespace mesh {

/*!
 * @brief Class GriddingAlgorithm<DIM> manages gridding operations in
 * SAMRAI.  Specifically, it provides AMR patch hierarchy generation and 
 * regridding routines that may be used with a variety of AMR solution 
 * algorithms and application-specific numerical routines.
 *
 * The three main functions provided by this class are:
 *    - @b    makeCoarsestLevel()   
 *       This routine constructs or repartitions
 *       the coarsest hierarchy level (level 0).
 *
 *    - @b    makeFinerLevel()   
 *       This routine will attempt to add a new 
 *       finest level to the hierarchy if the
 *       maximum number of levels allows it and
 *       cells on the current finest level are
 *       selected for refinement.
 *
 *    - @b    regridAllFinerLevels()   
 *       This routine will regrid all levels finer
 *       than some specified level based on cells 
 *       that are selected for refinement on each 
 *       level finer than and including the given
 *       level.  This routine may add a new finest
 *       hierarchy level if the maximum number of
 *       levels allows it and cells on the current
 *       finest level are selected for refinement.
 *       Levels may also be removed from the 
 *       hierarchy if no cells are tagged.
 * 
 * 
 * These basic AMR operations are used to generate of individual levels in 
 * the AMR patch hierarchy at the beginning of a simulation, and regridding 
 * collections of levels during an adaptive calculation.  More details are
 * found in the comments accompanying each member function below.
 *
 * The operations that identify cells for refinement on a single level and 
 * initialize data and solution algorithm-specific information that depend
 * on the AMR hierarchy configuration are provided by the data member of
 * type TagAndInitializeStrategy<DIM>.  Operations that cluster
 * tagged cells into a collection of box regions are provided by the 
 * BoxGeneratorStrategy<DIM> data member.  Routines that load balancing 
 * patches on each level are provided by the LoadBalanceStrategy<DIM> data 
 * member.  The collaboration between this class and each of those objects
 * follows the Strategy design pattern.  Each instantiation of this gridding
 * algorithm class is configured with concrete implementations of those
 * routines by passing appropriate objects into this constructor. 
 *
 * Initialization of an GriddingAlgorithm<DIM> object is performed via a
 * combination of default parameters and values read from input.  Data 
 * read from input is summarized as follows:
 *
 * Required input keys and data types: 
 *
 *    - \b    max_levels   
 *       Integer value specifying maximum number of levels
 *       allowed in the AMR patch hierarchy.
 * 
 *    - \b    largest_patch_size   
 *       An array of integer vectors (each has length = DIM) that specify
 *       the dimensions of largest patch allowed on each level in the 
 *       hierarchy.  The index of the vector in the array corresponds to the
 *       number of the level to which it applies.  If more values are given 
 *       than the maximum number of levels, extra entries will be ignored.
 *       If fewer values are given, then the last element in the array will 
 *       be used on each level without a specified input value.  For example, 
 *       if only a single value is specified, then that value will be used on 
 *       all levels.  See sample input below for more information.
 * 
 *    - \b    ratio_to_coarser   
 *       Set of max_levels - 1 integer vectors, each of which indicates 
 *       the ratio of the index space of a patch level to that of the next 
 *       coarser level.  The input for each level must correspond to the 
 *       format ``level_n = vector'', where n is the level number and each 
 *       vector must have length DIM.  If more values are given 
 *       than max_levels - 1 , extra entries will be ignored.  If fewer 
 *       values are given, then the last element in the array will be used 
 *       for finer levels.
 *
 * Optional input keys, data types, and defaults:
 *
 *    - \b    efficiency_tolerance   
 *       An array of double values, each of which specifies the minimum 
 *       percentage of tagged cells in each box used to construct patches
 *       on a finer level.  If the ratio of the number of tagged cells in 
 *       a box to total cells in the box is below the tolerance value, the 
 *       box may be split into smaller boxes and pieces removed until the 
 *       ratio becomes greater than or equal to the the tolerance.  The index 
 *       of the value in the array corresponds to the number of the level to 
 *       which the tolerance value applies.  If more values are given 
 *       than max_levels - 1 , extra entries will be ignored.  If fewer 
 *       values are given, then the last element in the array will be used 
 *       on each level without a specified input value.  For example, 
 *       if only a single value is specified, then that value will be used on 
 *       all levels.  If no input values are given, a default of 0.8 is used.
 *       See sample input below for input file format.
 *
 *    - \b    combine_efficiency   
 *       An array of double values, each of which serves as a threshold 
 *       for the ratio of the total number of cells in two boxes into which 
 *       a box may be split and the number of cells in the original box.  
 *       If that ratio is greater than combine efficiency, the box will not 
 *       be split.  This avoids splitting up portions of the domain into 
 *       potentially more costly smaller pieces if there appears to be little
 *       to be gained by splitting up the boxes.  The index of the value in 
 *       the array corresponds to the number of the level to which the 
 *       efficiency value applies.  If more values are given than 
 *       max_levels - 1 , extra entries will be ignored.  If fewer values 
 *       are given, then the last element in the array will be used on each 
 *       level without a specified input value.  For example, if only a single 
 *       value is specified, then that value will be used on all levels.  If 
 *       no input values are given, a default of 0.8 is used.  See 
 *       sample input below for input file format. 
 *
 *    - \b    smallest_patch_size   
 *       An array of integer vectors (each has length = DIM) that specify
 *       the dimensions of smallest patch allowed on each level in the
 *       hierarchy.  The smallest patch allowed must be at least as great as 
 *       the maximum ghost cell width for all variables in the problem.  If
 *       some smaller patch size is given in input, then it will be overridden 
 *       with a value consistent with the maximum ghost width.   The index of 
 *       the vector in the array corresponds to the number of the level to 
 *       which it applies.  If more values are given than the maximum number 
 *       of levels, extra entries will be ignored.  If fewer values are given,
 *       then the last element in the array will be used on each level without 
 *       a specified input value.  For example, if only a single value is 
 *       specified, then that value will be used on all levels.  If no input 
 *       is given, a default of the maximum ghost cell width over all 
 *       variables is used.  See sample input below for input file format.
 *  
 *    - \b    proper_nesting_buffer   
 *       An integer array specifying the number of coarse cells by which the 
 *       next finer level is nested within the interior of the domain of the 
 *       next coarser level when creating a new level.  If more values are 
 *       given than max_levels - 1 , extra entries will be ignored.  If fewer
 *       values are given, then the last element in the array will be used
 *       on each level without a specified input value.  For example,
 *       if only a single value is specified, then that value will be used on
 *       all levels.  If no values are given, a default of 1 is used for each 
 *       nesting buffer value.  See sample input below for input file format.
 *
 *    - \b    allow_patches_smaller_than_ghostwidth
 *       In order to enforce minimum patch size restrictions, boxes may be
 *       grown in gridding operations.  This may in turn lead to overlaps
 *       created between boxes. If overlaps are undesirable and you are 
 *       willing to relax the minimum size constraints, set this parameter
 *       true.  By default, it is false.
 *
 *    - \b    resolve_nonnesting_tags
 *       How to resolve user-specified tags that violates proper nesting.
 *       Set to one of these characters:
 *       @b 'i' for ignore (violating tags will be quietly disregarded).
 *       @b 'w' for warn (violating tags will cause a warning and be
 *       disregarded).
 *       @b 'e' for error (violating tags will cause an unrecoverable
 *       exception).
 *       The default is 'w'.  Ignoring is the most efficient because
 *       no checks are required (but this puts the burden of dealing
 *       with unrefined tagged cells on the user.
 *
 *    - \b    write_regrid_boxes   
 *       Output sequence of refine boxes to file. 
 *
 *    - \b    read_regrid_boxes   
 *       Read sequence of refine boxes from file. 
 *
 *    - \b    regrid_boxes_filename  
 *       Filename used for writing or reading refine boxes.
 *       The read and write boxes option is usable only under very specific
 *       circumstances. Please contact the SAMRAI developers at 
 *       samrai\@llnl.gov) for more information.
 *
 *
 * Note that when continuing from restart, the input values in the 
 * input file override all values read in from the restart database.
 *
 * The following represents sample input data for a three-dimensional problem:
 *
 * \verbatim
 *
 *    // Required input: maximum bumber of levels in patch hierarchy
 *    max_levels = 4
 *
 *    // Required input: vector ratio between each finer level and next coarser
 *    ratio_to_coarser {
 *       level_1 = 2, 2, 2
 *       level_2 = 2, 2, 2
 *       level_3 = 4, 4, 4
 *    }
 *
 *    // Required input: int vector for largest patch size on each level.
 *    largest_patch_size {
 *       level_0 = 40, 40, 40
 *       level_1 = 30, 30, 30
 *       // all finer levels will use same values as level_1...
 *    }
 *
 *    // Optional input:  buffer of one cell used on each level  
 *    proper_nesting_buffer = 1
 *    grow_after_nesting = FALSE
 *
 *    // Optional input: int vector for smallest patch size on each level.
 *    smallest_patch_size {
 *       level_0 = 16, 16, 16
 *       // all finer levels will use same values as level_0...
 *    }
 *
 *    // Optional input: different efficiency tolerance for each coarser level
 *    efficiency_tolerance = 0.80e0, 0.85e0, 0.90e0
 *
 *    // Optional input: combine efficiency is same for all levels.
 *    combine_efficiency = 0.95e0
 *
 *    write_regrid_boxes = TRUE
 *    regrid_boxes_filename = "regrid_boxes_32proc"
 *
 * \endverbatim
 *
 * @see mesh::TagAndInitializeStrategy
 * @see mesh::LoadBalanceStrategy
 * @see mesh::BoxGeneratorStrategy
 */

template<int DIM> class GriddingAlgorithm 
:
public BaseGriddingAlgorithm<DIM>
{
public:
   /*!
    * The constructor for GriddingAlgorithm<DIM> configures the gridding
    * algorithm with the concrete strategy objects in the argument list.
    * Gridding parameters are initialized from values provided in the 
    * specified input and in the restart database corresponding to the
    * specified object_name argument.  The constructor also registers 
    * this object for restart using the specified object name when the 
    * boolean argument is true.  Whether object will write its state to 
    * restart files during program execution is determined by this argument.
    * Note that it has a default state of true.
    *
    * If assertion checking is turned on, an unrecoverable assertion will
    * result if any of the input database, level strategy,
    * box generator, or load balancer pointers is null.  Exceptions
    * may also be thrown if any checks for consistency among input 
    * parameters fail.
    */
   GriddingAlgorithm(
      const std::string& object_name,
      tbox::Pointer<tbox::Database> input_db,
      tbox::Pointer< TagAndInitializeStrategy<DIM> > level_strategy,
      tbox::Pointer< BoxGeneratorStrategy<DIM> > generator,
      tbox::Pointer< LoadBalanceStrategy<DIM> > balancer,
      bool register_for_restart = true);

   /*!
    * Virtual destructor for GriddingAlgorithm<DIM>.
    */
   virtual ~GriddingAlgorithm<DIM>();

   /*!
    * This routine will attempt to construct the coarsest level in an AMR 
    * patch hierarchy (i.e., level 0).  If level 0 does not already exist, 
    * then the domain specification is checked against the constraints of
    * the grid generation procedures.  The level gridding strategy data 
    * member defines these constraints.  Recall that the domain specification
    * is maintained by the grid geometry object associated with the hierarchy.
    * Generally, an unrecoverable exception will result if the constraints
    * are not satisfied.
    * 
    * If level 0 already exists in the hierarchy, then the routine will 
    * generate a new level by re-applying the load balancing procedure to 
    * the existing level.  Data will be moved from the old level to the 
    * new level and the pre-existing level 0 will be discarded.  Note that 
    * this routine is different than the routine makeFinerLevel() below, 
    * which is used to construct levels 1 and finer.  In particular, this
    * routine does not select cells for refinement, whereas the other 
    * routine does.
    *
    * Important note: If assertion checking is turned on, then an 
    * unrecoverable assertion will result if either the patch hierarchy 
    * or its grid geometry is NULL.
    *
    * The two optional arguments are only to be used for a special case where
    * the user wishes to manually specify a box decomposition and load
    * balance for the coarsest level of the hierarchy.  The BoxArray argument
    * must be a decomposition of the the coarsest level, and must exactly
    * fill the index space of the physical domain of the hierarchy.  The
    * ProcessorMapping must be constructed to map each box in the BoxArray
    * to a processor.  The size of the mapping must be equal to the length
    * of the box array, or an assertion failure will result.
    *
    * @param hierarchy The hierarchy on which coarse level is constructed.
    * @param level_time Simulation time when level is constructed
    * @param override_boxes box array representing a decomposition of level
    *                       zero of the hierarchy
    * @param override_mapping processor mapping that maps each box in the
    *                         above array to a processor.
    */
   virtual void makeCoarsestLevel(
      tbox::Pointer< hier::BasePatchHierarchy<DIM> > hierarchy,
      const double level_time,
      const hier::BoxArray<DIM>& override_boxes = 0,
      const hier::ProcessorMapping& override_mapping = 0);

   /*!
    * This routine attempts to create a new level in an AMR patch hierarchy 
    * finer than the finest level currently residing in the hierarchy.  
    * It will select cells for refinement on the finest level and construct 
    * a new finest level, if necessary.  If no cells are selected for 
    * refinement, no new level will be added to the hierarchy.   The boolean 
    * argument initial_time indicates whether the routine is called at the 
    * initial simulation time.  If true, this routine is used to build 
    * individual levels during the construction of the AMR hierarchy at the 
    * initial simulation time.  If false, the routine is being used to add 
    * new levels to the hierarchy at some later point.  In either case, the 
    * time value is the current simulation time.  Note that this routine 
    * cannot be used to construct the coarsest level in the hierarchy
    * (i.e., level 0).  The routine makeCoarsestLevel() above must be used
    * for that purpose.
    *
    * The tag buffer indicates the number of cells by which cells selected
    * for refinement will be buffered before new finer level boxes are 
    * constructed.  The buffer is important to keep phenomena of interest 
    * on refined regions of the mesh until adaptive regridding occurs next.
    * Thus, the buffer size should take into account how the simulation may 
    * evolve before regridding occurs (e.g., number of timesteps taken).
    *
    * Important note: If assertion checking is activated, several checks
    * are applied to the functions arguments.  If any check is violated,
    * an unrecoverable assertion will result.  In particular, the hierarchy 
    * pointer must be non-NULL and the given level number must match that of 
    * the finest level currently residing in the hierarchy.  Also, the 
    * the tag buffer must be positive.
    */
   virtual void makeFinerLevel(
      tbox::Pointer< hier::BasePatchHierarchy<DIM> > hierarchy,
      const double level_time,
      const bool initial_time,
      const int tag_buffer,
      const double regrid_start_time = 0.);

   /*!
    * This routine attempts to reconfigure the patches on each level in
    * an AMR patch hierarchy which is finer than the specified level.  The 
    * given level number is that of the coarsest level on which cells will be
    * will be selected for refinement.  In other words, that level is the 
    * finest level that will not be subject to a change in its patch 
    * configuration during the regridding process.  Generally, this routine
    * should be used to alter a pre-existing AMR patch hierarchy based on
    * the need to adapt the computational mesh around some phenomenon of 
    * interest.  The routine makeFinerLevel() above should be used to 
    * construct an initial hierarchy configuration or to add more than one
    * new level into the hierarchy.  Also, this routine will not reconfigure
    * the patches on level 0 (i.e., the coarsest in any hierarchy).  The
    * routine makeCoarsestLevel() above is provided for that purpose.
    *
    * Note that the current algorithm permits at most one new finest level 
    * to be added to the hierarchy with each invocation of the regridding 
    * process.  This constraint, though seemingly restrictive makes the 
    * process of maintaining properly nested levels much easier.
    *
    * The tag buffer array indicates the number of cells by which cells 
    * selected for refinement on a level will be buffered before new finer
    * level boxes are constructed.  The buffer is important to keep phenomena 
    * of interest on refined regions of the mesh until adaptive regridding 
    * occurs next.  Thus, the buffer size should take into account how the 
    * simulation may evolve before regridding occurs (e.g., number of 
    * timesteps taken on each level).
    *
    * The boolean argument is used for regridding in time-dependent 
    * problems.  When true, it indicates that the specified level is
    * the coarsest level to synchronize at the current regrid time 
    * before this regridding method is called.  This is a pretty
    * idiosyncratic argument but allows some flexibility in the way 
    * memory is managed during time-dependent regridding operations.
    *
    * Important note: If assertion checking is activated, several checks
    * are applied to the functions arguments.  If any check is violated,
    * an unrecoverable assertion will result.  In particular, the hierarchy
    * pointer must be non-NULL and the given level number must match that of
    * of some level in the hierarchy.  Also, the tag buffer array must 
    * contain a positive value for each level in the hierarchy.
    */
   virtual void regridAllFinerLevels(
      tbox::Pointer< hier::BasePatchHierarchy<DIM> > hierarchy,
      const int level_number,
      const double regrid_time,
      const tbox::Array<int>& tag_buffer,
      tbox::Array<double> regrid_start_time = tbox::Array<double>(),
      const bool level_is_coarsest_to_sync = true);

   /*!
    * Return true if error estimation process uses time integration;
    * otherwise, return false.
    */
   virtual bool errorEstimationUsesTimeIntegration() const;

   /*!
    * Return the error coarsen ratio.  This is needed for cases where 
    * an error estimation scheme uses time integration (e.g. Richardson
    * extrapolation) to determine how many time levels to maintain to
    * properly apply the estimtion scheme. In general, an even refine 
    * ratio (e.g. 2, 4, 8) will maintain two time levels, while an odd 
    * refine ratio (e.g. 3) will maintain three. 
    */
   virtual int getErrorCoarsenRatio() const;

   /*!
    * Return true if level associated with the specified level number can 
    * be refined; i.e., the level number is less than that of the finest
    * level allowed in the hierarchy.  Otherwise, false is returned.
    */
   virtual bool levelCanBeRefined(const int level_number) const;

   /*!
    * Return pointer to level gridding strategy data member.
    */
   virtual
   tbox::Pointer< TagAndInitializeStrategy<DIM> > 
      getTagAndInitializeStrategy() const;

   /*!
    * Return pointer to load balance strategy data member.
    */
   virtual
   tbox::Pointer< LoadBalanceStrategy<DIM> >
      getLoadBalanceStrategy() const;

   /*!
    * Return maximum number of levels allowed in hierarchy.
    */
   virtual int getMaxLevels() const;

   /*!
    * Return const reference to ratio between specified level and next coarser.
    */
   virtual const hier::IntVector<DIM>& 
      getRatioToCoarserLevel(const int level_number) const;

   /*!
    * Return efficiency tolerance for clustering tags on level.
    */
   virtual double getEfficiencyTolerance(const int level_number) const;

   /*!
    * Return combine efficiency for clustering tags on level.
    */
   virtual double getCombineEfficiency(const int level_number) const;

   /*!
    * Return proper nesting buffer width for level.
    */
   virtual int getProperNestingBuffer(const int level_number) const;

   /*!
    * Return const reference to smallest patch size for level.
    */
   virtual const hier::IntVector<DIM>&
      getSmallestPatchSize(const int level_number) const;

   /*!
    * Return const reference to largest patch size for level.
    */
   virtual const hier::IntVector<DIM>&
      getLargestPatchSize(const int level_number) const;

   /*!
    * Print out all members of the class instance to given output stream.
    */
   virtual void printClassData(std::ostream& os) const;

   /*!
    * Write object state out to the given database.
    * 
    * When assertion checking is active, the database pointer must be non-null.
    */
   virtual void putToDatabase(tbox::Pointer<tbox::Database> db);

private:
   /*
    * Read input data from specified database and initialize class members.
    * The argument is_from_restart should be set to true if the simulation 
    * is from restart.  Otherwise, it should be set to false.
    * 
    * When assertion checking is active, the database pointer must be non-null.
    */
   void getFromInput(tbox::Pointer<tbox::Database> db, 
                     bool is_from_restart);

   /*
    * Read object state from the restart file and initialize class data 
    * members.  The database from which the restart data is read is
    * determined by the object_name specified in the constructor. 
    *
    * Unrecoverable Errors:
    *
    *    -The database corresponding to object_name is not found
    *     in the restart file.
    *
    *    -The class version number and restart version number do not 
    *     match.
    *
    */
   void getFromRestart();

   /*
    * Recursively regrid the specified hierarchy level and all finer levels
    * in the hierarchy.  This private member function is invoked by the 
    * regridAllFinerLevels() routine.  
    */
   void regridFinerLevel(
      tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
      const int level_number,
      const double regrid_time,
      const int finest_level_not_regridded,
      const bool level_is_coarsest_to_sync,
      const tbox::Array<int>& tag_buffer,
      const tbox::Array<double>& regrid_start_time = NULL);

   /*
    * Compute the proper nesting domain data for gridding
    * based at level base_ln in the given hierarchy.
    *
    * Sets d_domain_complement, d_proper_nesting_complement
    * and d_proper_nesting_boxes.
    */
   void findProperNestingData(
        const tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
        const int base_ln);

   /*
    * Set integer tags to specified value on intersection between patch
    * level and given box array.  The index value corresponds to the patch
    * descriptor entry of the cell-centered integer tag array.  The boolean
    * flag indicates whether tags are to be set on the regions corresponding
    * to the interior of the level only, if the tag arrays contain ghosts.
    */
   void setTagsOnLevel(
      const int tag_value,
      const tbox::Pointer< hier::PatchLevel<DIM> > level,
      const int index,
      const hier::BoxArray<DIM>& boxes,
      const bool interior_only = true) const;

   /*
    * Buffer each integer tag on patch level matching given tag value
    * with a border of matching tags.
    */
   void bufferTagsOnLevel(
      const int tag_value,
      const tbox::Pointer< hier::PatchLevel<DIM> > level,
      const int buffer_size) const;

   /*
    * Extend tags that are within the ghost cell width to the physical
    * domain boundary to the boundary.
    */
   void extendTagsToBoundary(
      const int tag_value,
      const tbox::Pointer< hier::PatchLevel<DIM> > level,
      const tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy) const;

   /*
    * Set the new level boxes using information stored in a file.
    * If cell tagging is not performed, the new level boxes may
    * be specified either from user input (by specifying "REFINE_BOXES" 
    * input) or from output from a previous run stored in an
    * HDF file.  This method sets the "new_level_boxes" based on
    * the information in the file.  It also sets the boolean
    * parameter "remove_old_fine_level" which indicates whether
    * the level box configuration has changed and, consequently, 
    * whether we need to remove the old level.
    */
   void readLevelBoxes(
      hier::BoxArray<DIM>& new_level_boxes,
      hier::ProcessorMapping& mapping,
      const tbox::Pointer<hier::PatchHierarchy<DIM> > hierarchy,
      const int level_number,
      const double regrid_time,
      bool& remove_old_fine_level);

   /*
    * Given a hierarchy and a level number, determine an array of boxes from
    * which a refinement of the level may be constructed.   It is assumed 
    * that the integer tags that identify cells for refinement have been set 
    * on the level before this routine is called.  Note that this routine 
    * also returns the processor mapping for the fine boxes.  The processor
    * mapping indicates the assignments of the new patches to processors
    * when the new fine patch level is constructed.  The processor mapping
    * and final box array are set by the LoadBalanceStrategy<DIM> data 
    * member.
    */
   void findRefinementBoxes(
      hier::BoxArray<DIM>& fine_boxes,
      hier::ProcessorMapping& mapping,
      const tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
      const int coarse_level_number ) const;

   /*
    * Set patch size and ghost cell information needed to create new
    * refinement levels.  This method applies to levels that are being 
    * used to build new finer levels (i.e. level_number is a coarser 
    * level in the hierarchy) and to levels which are simply being 
    * reconstructed (i.e. the same level in the hierarchy).  The boolean 
    * value "for_building_finer"  controls the logic for the two cases - 
    * in the former case, it is true while in the latter case it is false.   
    *
    * When a finer level is being constructed, the maximum number of ghost 
    * cells needed in any variable is used to compute the smallest patch 
    * size allowed and the extent to which patches may be extended to touch 
    * the physical boundary.  This avoids problems in setting ghost cell 
    * data that may occur when ghost cell regions intersect the physical 
    * boundary in strange ways.
    *
    * This routine sets the smallest and largest patch sizes for the specified
    * level, the smallest box to refine on the next coarser level, and the
    * number of cells that a patch may be extended to the boundary if it
    * sufficiently close to the boundary (extend_ghosts).
    */
   virtual void getGriddingParameters(
      hier::IntVector<DIM>& smallest_patch,
      hier::IntVector<DIM>& smallest_box_to_refine,
      hier::IntVector<DIM>& largest_patch,
      hier::IntVector<DIM>& extend_ghosts,
      const tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
      const int level_number,
      const bool for_building_finer) const;

   /*!
     @brief Check for user tags that violate proper nesting.
   */
   void checkNonnestingTags(
      const tbox::Pointer<hier::PatchLevel<DIM> > &level,
      const hier::BoxTree<DIM> &proper_nesting_complement );

   /*
    * Temporary method to sort box, used in performance comparison
    * of clustering algorithms.
    */
   static int qsortBoxCompare(const void *v, const void *w);

   /*
    * Static members for managing shared tag data among multiple
    * GriddingAlgorithm objects.
    */
   static int s_instance_counter;
   static int s_tag_indx;
   static int s_buf_tag_indx;

   /*
    * The object name is used for error reporting and accessing
    * restart file information.  Whether the object writes restart
    * file data depends on the value of this boolean which is
    * set in the constructor.
    */
   std::string d_object_name;
   bool d_registered_for_restart;

   /*
    * Data members that manage application-specific level initialization
    * and cell tagging, clustering of tagged cells into boxes, and load 
    * balancing of patches to processors, respectively.
    */
   tbox::Pointer< TagAndInitializeStrategy<DIM> >  d_tag_init_strategy;
   tbox::Pointer< BoxGeneratorStrategy<DIM> >   d_box_generator;
   tbox::Pointer< LoadBalanceStrategy<DIM> >    d_load_balancer;

   /*
    * Cell-centered integer variables use to tag cells for refinement. 
    * The descriptor index d_tag_indx is used to obtain tag information 
    * from user-defined routines on patch interiors.  The descriptor index 
    * d_buf_tag_indx is used to buffer tags on patches that may be 
    * distributed across processors.  The refine algorithm and schedule are 
    * used for interprocessor communication.
    */
   tbox::Pointer< pdat::CellVariable<DIM,int> > d_tag;
   tbox::Pointer< pdat::CellVariable<DIM,int> > d_buf_tag;
   int d_tag_indx;
   int d_buf_tag_indx;

   tbox::Pointer< xfer::RefineAlgorithm<DIM> > d_bdry_fill_tags;
   tbox::Array< tbox::Pointer< xfer::RefineSchedule<DIM> > > d_bdry_sched_tags;
    
   /*
    * True and false integer tag values set in constructor and used to 
    * set and compare integer tag values on the patch hierarchy.  Generally,
    * these variables are easier to read in the code than integer constants,
    * such as 0 and 1.
    */
   int d_true_tag;
   int d_false_tag;

   /*
    * Arrays of box lists (one for each level on which cell tagging occurs) 
    * used to ensure that finer levels reside within coarser level interiors.
    */
   tbox::Array< hier::BoxList<DIM> > d_proper_nesting_boxes;

   /*
    * Arrays of box lists (one for each level on which cell tagging occurs) 
    * used to ensure that finer levels reside within coarser level interiors.
    * These boxes mark the limits to which a new or regridded level may
    * extend and still be properly nested in the next coarser level.
    */
   tbox::Array< tbox::Pointer< hier::BoxTree<DIM> > >
      d_proper_nesting_complement;

   /*
    * The remaining data members are parameters that govern the 
    * construction and reconfiguration of the patch hierarchy:
    *
    * Maximum number of levels allowed in the AMR patch hierarchy.
    * This value MUST BE SPECIFIED IN THE INPUT FILE.
    */
   int d_max_levels;

   /*
    * Each ratio to coarser is a vector of integers > 1 that specifies the
    * refinement of the index space (for the patch level with the same number)
    * in terms of the index space on the next coarser level.  The coarsening
    * ratio for level 0 (i.e., the coarsest level in the hierarchy) will
    * always be the vector of ones.
    *
    * The coarsening ratios for each finer level in the hierarchy (up to 
    * the maximum number of levels) MUST BE SPECIFIED IN THE INPUT FILE.
    * These values will be checked against the constraints of the level 
    * error estimation routines.
    */
   tbox::Array< hier::IntVector<DIM> > d_ratio_to_coarser;

   /*
    * Parameters for box generation routines that govern the splitting
    * of boxes containing tagged cells into smaller boxes:
    *
    * The efficiency tolerance is a threshold value for the percentage of
    * tagged cells in each box.  If this percentage is below the tolerance,
    * the box will continue to be split into smaller boxes.  
    *
    * The combine efficiency is a threshold value for the sum of the volumes
    * of two boxes into which a box may be potentially split.  If that sum 
    * is greater than the combine efficiency multiplied by the volume of 
    * the original box, the box will not be split.
    *
    * For each of these parameters, an array of values may be given.  Each 
    * value in the array will be used for cell clustering on the level whose 
    * number matches the index of the value in the array.   If more values 
    * are given than the maximum number of levels, extra values will 
    * be ignored.  If fewer values are given, then the last element in the 
    * array will be used on each level without a specified input value.  
    * For example, if only a single value is specified, then that value 
    * will be used for all levels.
    *
    * These values are optional input parameters.  If not given, a default
    * value of 0.8 is set for each parameter.
    */
   tbox::Array<double> d_efficiency_tolerance;
   tbox::Array<double> d_combine_efficiency;

   /*
    * Largest and smallest patches allowed in the AMR patch hierarchy.
    *
    * The smallest patch value may be specified in the input file.  However,
    * the value must always be at least as great as the maximum ghost cell
    * width for all variables in the problem.  Thus, if some smaller patch
    * size is given in input, then it will be overridden with a value 
    * consistent with the maximum ghost width.  If the value is not specified
    * in the input file, a default value of the maximum ghost width is used.
    *
    * The largest patch value MUST BE SPECIFIED IN THE INPUT FILE.  
    * 
    * For each of these parameters, an array of values may be given.  Each 
    * vector in the array will be used on the level whose number matches 
    * the index of the vector in the array.   If more vectors are given than 
    * the maximum number of levels, extra values will be ignored.  If fewer 
    * vectors are given, then the last element in the array will be used on 
    * each level without a specified input value.  For example, if only a 
    * single vector is specified, then that vector will be used for all levels.
    */
   tbox::Array< hier::IntVector<DIM> > d_smallest_patch_size;
   tbox::Array< hier::IntVector<DIM> > d_largest_patch_size;

   /*
    * The proper nesting buffer specifies the number of coarse cells 
    * by which the next finer level is nested within the interior of 
    * the domain of the next coarser level.  These buffer values 
    * are used to compute the proper nesting boxes on each level.  They
    * can be used to guarantee that:
    * 
    * (1) Adjacent cells on the composite grid differ by no more than 
    *     one refinement level.
    * 
    * (2) Data on the interior of fine patches can be filled by 
    *     interpolation from the next coarser level if the buffer 
    *     width is at least as large as the maximum stencil width 
    *     of all interpolation operators.
    *
    * The proper nesting buffer size for each level may be specified in the
    * input file.  If not, a nesting buffer of 1 is used.  This value should
    * be suitable for most problems.  It is not recommended that the buffer 
    * be set to 0 since this may cause features of the solution to move off
    * of refined mesh regions before subsequent regridding occurs.  If 
    * additional buffering is required (e.g., if the stencil width for some 
    * interpolation operator is greater than one), then it may be necessary 
    * to increase this value.
    */
   tbox::Array<int> d_proper_nesting_buffer;

   /*
    * Boxes generated from tagging should satisfy the constraint that none
    * are smaller than the specified minimum box size, which is generally
    * the max ghost width.  But intersecting boxes with the proper nesting 
    * buffer may create boxes smaller than the specified minimum size.  The 
    * default behavior is to grow these boxes to meet the constraint and, 
    * in the process, some overlapping boxes may be constructed.  
    *
    * To avoid this behavior, the "allow_patches_smaller_than_ghostwidth" 
    * option may be set true.  Setting it true relaxes the min box size 
    * constraint so there may be boxes that are less than the specified min 
    * box size, but it will avoid generating overlapping boxes. It is 
    * by default set false.
    */
   bool d_allow_patches_smaller_than_ghostwidth;

   bool d_barrier_before_clustering;
   bool d_sort_boxes_after_clustering;

   /*
    * The following data members are used for writing and reading refinement
    * boxes for the case where one wants to specify refinement rather than
    * have the TagAndInitializeStrategy do it. The two booleans 
    * specify whether a user wants to dump or read a set of refine boxes.  
    * These may be set on input and are false by default. 
    */
   bool d_write_dumped_level_boxes;
   bool d_read_dumped_level_boxes;
   std::string d_regrid_boxes_filename;
   hier::BoxIOUtility<DIM>* d_regrid_box_utility;
   tbox::Array<int> d_regrid_box_counter; 

   bool d_extend_tags_to_bdry;

   /*!
     @brief How to resolve user tags that violate nesting requirement.

     Values can be:
     - 'i' Ignore (violating tags will be quietly disregarded)
     - 'w' Warn (violating tags will cause a warning and be disregarded)
     - 'e' Error (violating tags will cause an unrecoverable exception)

     This defaults to 'w' and set by input parameter
     "resolve_nonnesting_tags".
   */
   char d_resolve_nonnesting_tags;

   /*
    * Timers interspersed throughout the class (those timers used 
    * only locally within a function are defined in that method).
    */
   tbox::Pointer<tbox::Timer> t_find_domain_complement;
   tbox::Pointer<tbox::Timer> t_intersect_boxes_find_refinement;
   tbox::Pointer<tbox::Timer> t_load_balance_boxes; 
   tbox::Pointer<tbox::Timer> t_bdry_fill_tags_create;
   tbox::Pointer<tbox::Timer> t_make_coarsest;
   tbox::Pointer<tbox::Timer> t_make_finer;
   tbox::Pointer<tbox::Timer> t_remove_intersections_make_finer;
   tbox::Pointer<tbox::Timer> t_regrid_all_finer;
   tbox::Pointer<tbox::Timer> t_remove_intersections_regrid_all;
   tbox::Pointer<tbox::Timer> t_remove_intersections_find_proper;
   tbox::Pointer<tbox::Timer> t_intersect_boxes_find_proper;
   tbox::Pointer<tbox::Timer> t_set_tags;
   tbox::Pointer<tbox::Timer> t_buffer_tags;
   tbox::Pointer<tbox::Timer> t_bdry_fill_tags_comm;
   tbox::Pointer<tbox::Timer> t_find_refinement;
   tbox::Pointer<tbox::Timer> t_find_boxes_containing_tags;
   tbox::Pointer<tbox::Timer> t_find_nesting_restriction;
   tbox::Pointer<tbox::Timer> t_apply_nesting_restriction;
   tbox::Pointer<tbox::Timer> t_extend_domain_to_boundary;
   tbox::Pointer<tbox::Timer> t_coalesce_boxes;
   tbox::Pointer<tbox::Timer> t_grow_boxes_within_domain;
   tbox::Pointer<tbox::Timer> t_before_load_balance;

   // The following are not yet implemented:
   GriddingAlgorithm(const GriddingAlgorithm<DIM>&);
   void operator=(const GriddingAlgorithm<DIM>&);

};

}
}
#ifndef DEBUG_NO_INLINE
#include "GriddingAlgorithm.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "GriddingAlgorithm.C"
#endif

