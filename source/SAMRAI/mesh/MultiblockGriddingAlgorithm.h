/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   AMR hierarchy generation and regridding routines. 
 *
 ************************************************************************/

#ifndef included_mesh_MultiblockGriddingAlgorithm
#define included_mesh_MultiblockGriddingAlgorithm

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Connector.h"
#include "SAMRAI/hier/MappedBoxLevel.h"
#include "SAMRAI/xfer/MultiblockRefineAlgorithm.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/mesh/BaseGriddingAlgorithm.h"
#include "SAMRAI/mesh/BoxGeneratorStrategy.h"
#include "SAMRAI/mesh/LoadBalanceStrategy.h"
#include "SAMRAI/mesh/MultiblockGriddingTagger.h"

#include <iostream>
#include <string>
#include <vector>

#define MGA_RECORD_STATS
// #undef MGA_RECORD_STATS

#ifdef MGA_RECORD_STATS
#include "SAMRAI/tbox/Statistic.h"
#include "SAMRAI/tbox/Statistician.h"
#endif

namespace SAMRAI {
namespace mesh {

/*!
 * @brief Class MultiblockGriddingAlgorithm manages gridding operations in
 * SAMRAI.  Specifically, it provides AMR patch hierarchy generation and
 * regridding routines that may be used with a variety of AMR solution
 * algorithms and application-specific numerical routines.
 *
 * The three main functions provided by this class are:
 *   - @b    makeCoarsestLevel()
 *      This routine constructs or repartitions
 *      the coarsest hierarchy level (level 0).
 *
 *   - @b    makeFinerLevel()
 *      This routine will attempt to add a new
 *      finest level to the hierarchy if the
 *      maximum number of levels allows it and
 *      cells on the current finest level are
 *      selected for refinement.
 *
 *   - @b    regridAllFinerLevels()
 *      This routine will regrid all levels finer
 *      than some specified level based on cells
 *      that are selected for refinement on each
 *      level finer than and including the given
 *      level.  This routine may add a new finest
 *      hierarchy level if the maximum number of
 *      levels allows it and cells on the current
 *      finest level are selected for refinement.
 *      Levels may also be removed from the
 *      hierarchy if no cells are tagged.
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
 * type TagAndInitializeStrategy.  Operations that cluster
 * tagged cells into a collection of box regions are provided by the
 * BoxGeneratorStrategy data member.  Routines that load balancing
 * patches on each level are provided by the LoadBalanceStrategy data
 * member.  The collaboration between this class and each of those objects
 * follows the Strategy design pattern.  Each instantiation of this gridding
 * algorithm class is configured with concrete implementations of those
 * routines by passing appropriate objects into this constructor.
 *
 * Initialization of an MultiblockGriddingAlgorithm object is performed via a
 * combination of default parameters and values read from input.  Data
 * read from input is summarized as follows:
 *
 * Optional input keys, data types, and defaults:
 *
 *   - \b    efficiency_tolerance
 *      An array of double values, each of which specifies the minimum
 *      percentage of tagged cells in each box used to construct patches
 *      on a finer level.  If the ratio of the number of tagged cells in
 *      a box to total cells in the box is below the tolerance value, the
 *      box may be split into smaller boxes and pieces removed until the
 *      ratio becomes greater than or equal to the the tolerance.  The index
 *      of the value in the array corresponds to the number of the level to
 *      which the tolerance value applies.  If more values are given
 *      than max_levels - 1 , extra entries will be ignored.  If fewer
 *      values are given, then the last element in the array will be used
 *      on each level without a specified input value.  For example,
 *      if only a single value is specified, then that value will be used on
 *      all levels.  If no input values are given, a default of 0.8 is used.
 *      See sample input below for input file format.
 *
 *   - \b    combine_efficiency
 *      An array of double values, each of which serves as a threshold
 *      for the ratio of the total number of cells in two boxes into which
 *      a box may be split and the number of cells in the original box.
 *      If that ratio is greater than combine efficiency, the box will not
 *      be split.  This avoids splitting up portions of the domain into
 *      potentially more costly smaller pieces if there appears to be little
 *      to be gained by splitting up the boxes.  The index of the value in
 *      the array corresponds to the number of the level to which the
 *      efficiency value applies.  If more values are given than
 *      max_levels - 1 , extra entries will be ignored.  If fewer values
 *      are given, then the last element in the array will be used on each
 *      level without a specified input value.  For example, if only a single
 *      value is specified, then that value will be used on all levels.  If
 *      no input values are given, a default of 0.8 is used.  See
 *      sample input below for input file format.
 *
 *   - \b    check_nonrefined_tags
 *      How to resolve user-specified tags that violates proper nesting.
 *      Set to one of these characters:
 *      @b "IGNORE" - violating tags will be quietly disregarded.
 *      @b "WARN" - violating tags will cause a warning and be
 *      disregarded.
 *      @b "ERROR" - violating tags will cause an unrecoverable
 *      assertion.
 *      The default is "WARN".  It is fastest to ignore non-nesting tags
 *      because no checking has to be done.
 *
 *   - \b    check_overlapping_patches
 *      Specify whether to check for overlapping patches on a newly created
 *      level, and what to do if any are found.
 *      Set to one of these characters:
 *      @b "IGNORE" - there is no check for overlapping patches,
 *      and they will be quietly disregarded.
 *      @b "WARN" - overlapping patches will cause a warning and be
 *      disregarded.
 *      @b "ERROR" - violating tags will cause an unrecoverable
 *      assertion.
 *      The default is "IGNORE".  The check for overlapping patches may be
 *      expensive, so the use of "WARN" and "ERROR" is recommended only for debugging
 *      purposes.  To prevent the creation of any levels with overlapping
 *      patches, see the input flag
 *      "allow_patches_smaller_than_minimum_size_to_prevent_overlaps"
 *
 *   - \b   sequentialize_patch_indices
 *      Whether to globally sequentialize patch indices.  This is not
 *      scalable, but is required for writing correct VisIt files.
 *      Due to the current VisIt requirement, this is currently true by default.
 *      It will evetually be set back to false after we remove the VisIt
 *      requirement.
 *
 *
 * Note that when continuing from restart, the input values in the
 * input file override all values read in from the restart database.
 *
 * The following represents sample input data for a three-dimensional
 * problem:
 *
 * \verbatim
 *
 *   // Optional input: different efficiency tolerance for each coarser level
 *   efficiency_tolerance = 0.80e0, 0.85e0, 0.90e0
 *
 *   // Optional input: combine efficiency is same for all levels.
 *   combine_efficiency = 0.95e0
 *
 * \endverbatim
 *
 * @see mesh::TagAndInitializeStrategy
 * @see mesh::LoadBalanceStrategy
 * @see mesh::BoxGeneratorStrategy
 */

class MultiblockGriddingAlgorithm:
   public BaseGriddingAlgorithm
{
public:
   /*!
    * The constructor for MultiblockGriddingAlgorithm configures the
    * gridding algorithm with the concrete strategy objects in the
    * argument list.
    *
    * Gridding parameters are initialized from values provided in the
    * specified input and in the restart database corresponding to the
    * specified object_name argument.  The constructor also registers
    * this object for restart using the specified object name when the
    * boolean argument is true.  Whether object will write its state
    * to restart files during program execution is determined by this
    * argument.  Note that it has a default state of true.
    *
    * If assertion checking is turned on, an unrecoverable assertion will
    * result if any of the input database, level strategy,
    * box generator, or load balancer pointers is null.  Assertions
    * may also be thrown if any checks for consistency among input
    * parameters fail.
    *
    * @param mb_hierarchy The hierarchy that this
    * MultiblockGriddingAlgorithm will work on.  The pointer cached.
    * All hierarchy operations will be on this hierarchy.
    *
    * @param[in] object_name For registering the object in the restart
    * database.
    *
    * @param[in] input_db
    *
    * @param[in] level_strategy
    *
    * @param[in] generator
    *
    * @param[in] balancer Load balancer
    *
    * @param balancer0 Special load balancer to use when a single
    * process owns all the unbalanced load (such as during
    * initialization).  If omitted, will use @c balancer instead.
    *
    * @param[in] mb_tagger_strategy
    *
    * @param[in] register_for_restart
    */
   MultiblockGriddingAlgorithm(
      const tbox::Pointer<hier::PatchHierarchy> &mb_hierarchy,
      const std::string& object_name,
      tbox::Pointer<tbox::Database> input_db,
      tbox::Pointer<TagAndInitializeStrategy> tag_init_strategy,
      tbox::Pointer<BoxGeneratorStrategy> generator,
      tbox::Pointer<LoadBalanceStrategy> balancer,
      tbox::Pointer<LoadBalanceStrategy> balancer0,
      MultiblockGriddingTagger* mb_tagger_strategy =
         (MultiblockGriddingTagger *)NULL,
      bool register_for_restart = true);

   /*!
    * @brief Destructor
    *
    * Virtual destructor for MultiblockGriddingAlgorithm.
    */
   virtual ~MultiblockGriddingAlgorithm();

   /*!
    * @brief Create or recreate the coarsest level.
    *
    * This is an implementation of interface defined in BaseGriddingAlgorithm.
    *
    * This routine will attempt to construct the coarsest level in the AMR
    * patch hierarchy (i.e., level 0).  If level 0 does not already exist,
    * then the domain specification is checked against the constraints of
    * the grid generation procedures.  The level gridding strategy data
    * member defines these constraints.  Recall that the domain specification
    * is maintained by the grid geometry object associated with the hierarchy.
    * Generally, an unrecoverable assertion will result if the constraints
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
    * @param[in] level_time Simulation time.
    *
    * @param[in] override_mapped_box_level For specifying the
    * configuration of the new level (instead of applying the load
    * balancer).  If not using the override, use the other alternate
    * version of this method, defined in BaseGriddingAlgorithm.
    */
   void
   makeCoarsestLevel(
      const double level_time,
      const hier::MappedBoxLevel& override_mapped_box_level);

   using BaseGriddingAlgorithm::makeCoarsestLevel;

   /*!
    * @brief Attempts to create a new level in the hierarchy finer
    * than the finest level currently residing in the hierarchy.
    *
    * This is an implementation of interface defined in BaseGriddingAlgorithm.
    *
    * This routine attempts to create a new level in the AMR patch hierarchy
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
    * an unrecoverable assertion will result.  In particular,
    * the tag buffer must be non-negative.
    *
    * @param[in] level_time See text.
    *
    * @param[in] initial_time See text.
    *
    * @param[in] tag_buffer See text.
    *
    * @param regrid_start_time[in] The simulation time when the
    * regridding operation began (this parameter is ignored except
    * when using Richardson extrapolation)
    */
   void
   makeFinerLevel(
      const double level_time,
      const bool initial_time,
      const int tag_buffer,
      const double regrid_start_time = 0.);

   /*!
    * @brief Attempts to reconfigure the patches on each level in the
    * AMR patch hierarchy which is finer than the specified level.
    *
    * The given level number is that of the coarsest level on which
    * cells will be will be selected for refinement.  In other words,
    * that level is the finest level that will not be subject to a
    * change in its patch configuration during the regridding process.
    * Generally, this routine should be used to alter a pre-existing
    * AMR patch hierarchy based on the need to adapt the computational
    * mesh around some phenomenon of interest.  The routine
    * makeFinerLevel() above should be used to construct the initial
    * hierarchy configuration or to add more than one new level into
    * the hierarchy.  Also, this routine will not reconfigure the
    * patches on level 0 (i.e., the coarsest in any hierarchy).  The
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
    * an unrecoverable assertion will result.  In particular,
    * and the given level number must match that of an existing level
    * in the hierarchy.  Also, the tag buffer array must
    * contain a non-negative value for each level in the hierarchy.
    *
    * @param[in] level_number Coarsest level on which cells will be
    * tagged for refinement
    *
    * @param[in] regrid_time Simulaition time when regridding occurs
    *
    * @param[in] tag_buffer Size of buffer on each level around tagged
    * cells that will be covered by the next finer level
    *
    * @param[in] regrid_start_time The simulation time when the
    * regridding operation began on each level (this parameter is
    * ignored except when using Richardson extrapolation)
    *
    * @param[in] level_is_coarsest_to_sync Level is the coarsest to sync
    */
   void
   regridAllFinerLevels(
      const int level_number,
      const double regrid_time,
      const tbox::Array<int>& tag_buffer,
      tbox::Array<double> regrid_start_time = tbox::Array<double>(),
      const bool level_is_coarsest_to_sync = true);

   /*!
    * @brief Return true if error estimation process uses time integration;
    * otherwise, return false.
    *
    * @return true if error estimation process uses time integration;
    * otherwise, return false.
    */
   bool
   errorEstimationUsesTimeIntegration() const;

   /*!
    * @brief Return pointer to level gridding strategy data member.
    *
    * @return pointer to level gridding strategy data member.
    */
   virtual
   tbox::Pointer<TagAndInitializeStrategy>
   getTagAndInitializeStrategy() const;

   /*!
    * Return pointer to load balance strategy data member.
    */
   virtual
   tbox::Pointer<LoadBalanceStrategy>
   getLoadBalanceStrategy() const;

   /*!
    * @brief Return efficiency tolerance for clustering tags on level.
    *
    * @return efficiency tolerance for clustering tags on level.
    */
   double
   getEfficiencyTolerance(
      const int level_number) const;

   /*!
    * @brief Return combine efficiency for clustering tags on level.
    *
    * @return combine efficiency for clustering tags on level.
    */
   double
   getCombineEfficiency(
      const int level_number) const;

   /*!
    * Print out all members of the class instance to given output stream.
    */
   void
   printClassData(
      std::ostream& os) const;

   /*!
    * Write object state out to the given database.
    *
    * When assertion checking is active, the database pointer must be non-null.
    */
   void
   putToDatabase(
      tbox::Pointer<tbox::Database> db);

   /*
    * Write out statistics recorded on numbers of cells and patches generated.
    */
   void
   printStatistics(
      std::ostream& s = tbox::plog) const;

   /*!
    * @brief Get the name of this object.
    *
    * @return The name of this object.
    */
   const std::string&
   getObjectName() const;

private:
   /*
    * Static integer constant describing this class's version number.
    */
   static const int ALGS_GRIDDING_ALGORITHM_VERSION;

   //! @brief Shorthand typedef.
   typedef hier::Connector::NeighborSet NeighborSet;

   /*
    * @brief Read input data from specified database and initialize class members.
    *
    * When assertion checking is active, the database pointer must be non-null.
    *
    * @param[in] db
    *
    * @param[in] is_from_restart Should be set to true if the
    * simulation is from restart.  Otherwise, it should be set to
    * false.
    */
   void
   getFromInput(
      tbox::Pointer<tbox::Database> db,
      bool is_from_restart);

   /*
    * @brief Read object state from the restart file and initialize
    * class data members.
    *
    * The database from which the restart data is read is determined
    * by the object_name specified in the constructor.
    *
    * Unrecoverable Errors:
    *
    *   -The database corresponding to object_name is not found
    *    in the restart file.
    *
    *   -The class version number and restart version number do not
    *    match.
    */
   void
   getFromRestart();

   /*
    * @brief Recursively regrid the hierarchy level and all finer
    * levels in the hierarchy.
    *
    * This private member function is invoked by the
    * regridAllFinerLevels() routine.
    */
   void
   regridFinerLevel(
      const int level_number,
      const double regrid_time,
      const int finest_level_not_regridded,
      const bool level_is_coarsest_to_sync,
      const tbox::Array<int>& tag_buffer,
      const tbox::Array<double>& regrid_start_time = tbox::Array<double>(0));

   /*
    * @brief Set integer tags to specified value on intersection
    * between patch level and the MappedBoxLevel provided by the
    * Connector.
    *
    * The index value corresponds to the patch descriptor entry of the
    * cell-centered integer tag array.  The boolean flag indicates
    * whether tags are to be set on the regions corresponding to the
    * interior of the level only, if the tag arrays contain ghosts.
    *
    * @param[in] tag_value
    *
    * @param[in] level
    *
    * @param[in] index
    *
    * @param[in] level_to_fill_mapped_box_level Connector from the level with
    * the tags to the MappedBoxLevel describing where to fill.
    *
    * @param[in] interior_only
    *
    * @param fill_box_growth How much to grow fill boxes before using them
    *       to tag.  Must be in index space of level holding tags..
    */
   void
   fillTagsFromMappedBoxLevel(
      const int tag_value,
      const tbox::Pointer<hier::PatchLevel> level,
      const int index,
      const hier::Connector& level_to_fill_mapped_box_level,
      const bool interior_only, // Default v2.x.x  = true
      const hier::IntVector& fill_box_growth // Default v2.x.x = hier::IntVector::getZero(tbox::Dimension(DIM))
      ) const;

   /*
    * @brief Make a map that, when applied to an improperly nested MappedBoxLevel,
    * removes the nonnesting parts.
    */
   void
   makeProperNestingMap(
      const hier::MappedBoxLevel& unnested_mapped_box_level,
      const hier::Connector& hierarchy_to_unnested,
      const hier::Connector& unnested_to_hierarchy,
      const int unnested_ln,
      const hier::BlockId& block_id,
      hier::MappedBoxLevel& nested_mapped_box_level,
      hier::Connector& unnested_to_nested) const;

   /*
    * @brief Make a map that, when applied to an MappedBoxLevel that extends too far
    * outside of a nominal mapped_box_level, removes those parts too far outside.
    */
   void
   makeOverflowNestingMap(
      const hier::MappedBoxLevel& unnested_mapped_box_level,
      const hier::Connector& unnested_to_nominal,
      hier::MappedBoxLevel& nested_mapped_box_level,
      hier::Connector& unnested_to_nested,
      const hier::BlockId& block_id,
      const hier::IntVector& allowed_overflow) const;

   /*!
    * @brief Make a map from a MappedBoxLevel to parts of that MappedBoxLevel that
    * violate proper nesting.
    *
    * @param[in] candidate MappedBoxLevel being examined for nesting violation.
    *
    * @param[out] violator MappedBoxLevel containing violating parts of candidate.
    *
    * @param[in] candidate_ln Level number corresponding to candidate's refinement ratio.
    *
    * @param[in] candidate_to_hierarchy Connector to mapped_box_level number
    *       candidate_ln in the hierarchy.
    *
    * @param[in] hierarchy_to_candidate Connector from mapped_box_level number
    *       candidate_ln in the hierarchy.
    */
   void
   computeNestingViolator(
      const hier::MappedBoxLevel& candidate,
      hier::MappedBoxLevel& violator,
      hier::Connector& candidate_to_violator,
      const hier::Connector& candidate_to_hierarchy,
      const hier::Connector& hierarchy_to_candidate,
      const int candidate_ln,
      const hier::BlockId& block_id) const;

   void
   extendMappedBoxesToDomainBoundary(
      hier::MappedBoxLevel& new_mapped_box_level,
      hier::Connector& tag_to_new,
      hier::Connector& new_to_tag,
      const hier::Connector& tag_to_tag,
      const hier::BoxList& physical_domain_array,
      const hier::IntVector& extend_ghosts,
      const hier::PatchHierarchy& hierarchy,
      const int new_ln) const;

   /*!
    * @brief Precompute data used to define proper nesting.
    *
    * Data is associated with level number ln, to be used for
    * constructing level number ln+1.
    */
   void
   computeNestingData(
      const int ln,
      const hier::BlockId& block_id);

   /*!
    * @brief Make a mapping Connector that can be used to grow boxes
    * within domain by the minimum amount needed to make all boxes in a
    * MappedBoxLevel satisfy the min_size requirement.
    *
    * The map is local.
    */
   void
   makeBoxGrowingMap(
      const hier::MappedBoxLevel& orig_mapped_box_level,
      const hier::Connector& hierarchy_to_orig,
      const hier::Connector& orig_to_hierarchy,
      const int ln,
      const hier::BlockId& block_id,
      hier::MappedBoxLevel& grown_mapped_box_level,
      hier::Connector& orig_to_grown,
      const hier::IntVector& min_size) const;

   void
   refineNewMappedBoxLevel(
      hier::MappedBoxLevel& new_mapped_box_level,
      hier::Connector& tag_to_new,
      hier::Connector& new_to_tag,
      const hier::IntVector& ratio) const;

   /*
    * @brief Sort the nodes in the new MappedBoxLevel.
    *
    * Sorting by the corners of the MappedBoxes removes randomness in
    * the MappedBoxLevel.  Sequentializing the global indices numbers
    * them sequentially, like Patch numbering in the traditional SAMR
    * approach.
    */
   void
   sortNodes(
      hier::MappedBoxLevel& new_mapped_box_level,
      hier::Connector& tag_to_new,
      hier::Connector& new_to_tag,
      bool sort_by_corners,
      bool sequentialize_global_indices) const;

   /*
    * @brief Version of sortNodes that does no Connector operations.
    *
    * If sequentialize_global_indices is true, then the LocalIds of the 
    * MappedBoxes in the level will be globally and uniquely sequentialized
    * across all processors.  If false, the LocalIds will be sequentialized
    * locally, but not uniquely across all processors.
    *
    * @param[in,out] new_mapped_box_level  level to be sorted
    * @param[in] sequentialize_global_indices
    */
   void
   sortNodes(
      hier::MappedBoxLevel& new_mapped_box_level,
      bool sequentialize_global_indices) const;

   /*
    * @brief Buffer each integer tag on patch level matching given tag
    * value with a border of matching tags.
    */
   void
   bufferTagsOnLevel(
      const int tag_value,
      const tbox::Pointer<hier::PatchLevel> level,
      const int buffer_size) const;

   /*
    * @brief Set the new level boxes using information stored in a file.
    *
    * If cell tagging is not performed, the new level boxes may
    * be specified either from user input (by specifying "REFINE_BOXES"
    * input) or from output from a previous run stored in an
    * HDF file.  This method sets the "new_level_boxes" based on
    * the information in the file.  It also sets the boolean
    * parameter "remove_old_fine_level" which indicates whether
    * the level box configuration has changed and, consequently,
    * whether we need to remove the old level.
    */
   void
   readLevelBoxes(
      hier::BoxList& new_level_boxes,
      hier::MappedBoxLevel& new_mapped_box_level,
      hier::Connector& coarser_to_new,
      hier::Connector& new_to_coarser,
      const hier::BlockId& block_id,
      const int level_number,
      const double regrid_time,
      bool& remove_old_fine_level);

   /*
    * @brief Given a level number, determine an array of boxes from
    * which a refinement of the level may be constructed.
    *
    * It is assumed that the integer tags that identify cells for
    * refinement have been set on the level before this routine is
    * called.  At the end of this function, new_mapped_box_level will
    * represent the box extents of a new fine level on the given
    * block, but in the index space of the coarser level (that is, the
    * level that was tagged).
    */
   void
   findRefinementBoxes(
      hier::MappedBoxLevel& new_mapped_box_level,
      hier::Connector& tag_to_new,
      hier::Connector& new_to_tag,
      const int coarse_level_number,
      const hier::BlockId& block_id) const;

   /*
    * @brief Set patch size and ghost cell information needed to create new
    * refinement levels.
    *
    * This method applies to levels that are being used to build new
    * finer levels (i.e. level_number is a coarser level in the
    * hierarchy) and to levels which are simply being reconstructed
    * (i.e. the same level in the hierarchy).  The boolean value
    * "for_building_finer" controls the logic for the two cases - in
    * the former case, it is true while in the latter case it is
    * false.
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
   void
   getGriddingParameters(
      hier::IntVector& smallest_patch,
      hier::IntVector& smallest_box_to_refine,
      hier::IntVector& largest_patch,
      hier::IntVector& extend_ghosts,
      const hier::PatchHierarchy& hierarchy,
      const int level_number,
      const bool for_building_finer) const;

   /*
    * @brief Given a new_mapped_box_level that is outputted from
    * findRefinementBoxes, load balance the boxes and refine them to the
    * resolution of level tag_ln+1.
    *
    * If the rank_group has been set to
    * a non-default state, it will restrict the load balance to a certain
    * set of processors.  The Connector agruments should be references to
    * the same objects that were passed into findRefinementBoxes.
    */
   void
   loadBalanceAndRefineBoxes(
      hier::MappedBoxLevel& new_mapped_box_level,
      hier::Connector& tag_to_new,
      hier::Connector& new_to_tag,
      const int tag_ln,
      const tbox::RankGroup& rank_group,
      const hier::BlockId& block_id) const;

   /*
    * @brief Set an array of RankGroup objects in order to assign processors to
    * specific blocks when load balancing a new level.
    *
    * Each entry in
    * work_on_block contains an integer representing the amount of work
    * expected to be on a specific block within the new level.  Each entry
    * of rank_groups will contain a RankGroup that determines which processors
    * will receive patches for a specific block.
    *
    * If the number of active blocks on the new level (defined as the number
    * of blocks that have nonzero work according to the work_on_block array),
    * is greater than or equal to the number of ranks in the SAMRAI_MPI
    * object, then all entries of the rank_groups array will be set to a
    * default state, meaning that all processors will be eligible to receive
    * patches on all blocks.
    */
   void
   setRankGroups(
      tbox::Array<tbox::RankGroup>& rank_groups,
      const tbox::Array<int>& work_on_block,
      const tbox::SAMRAI_MPI& mpi) const;

   /*!
    * @brief Return proper nesting buffer width for level.
    *
    * Level number level_number+1 must nest inside level level_number
    * by the width returned (except where the levels touch the domain
    * boundary).
    */
   int
   getProperNestingBuffer(
      const int level_number) const;

   /*!
    * @brief Return const reference to smallest patch size for level.
    */
   const hier::IntVector&
   getSmallestPatchSize(
      const int level_number) const;

   /*!
    * @brief Return const reference to largest patch size for level.
    */
   const hier::IntVector&
   getLargestPatchSize(
      const int level_number) const;

   /*!
    * @brief Return pointer to load balance strategy specialized for the case
    * where one processor owns all the initial loads.
    *
    * @return pointer to load balance strategy specialized for the case
    * where one processor owns all the initial loads.
    */
   const hier::IntVector&
   getRatioToCoarserLevel(
      const int level_number) const;

   /*!
    * @brief Check for user tags that violate proper nesting.
    */
   void
   checkNonrefinedTags(
      const hier::PatchLevel& level,
      int tag_ln,
      const hier::BlockId& block_id) const;

   /*!
    * @brief Check for overlapping patches within the level.
    */
   void
   checkOverlappingPatches(
      const hier::Connector& mapped_box_level_to_self) const;

   void
   warnIfDomainTooSmallInPeriodicDir(
      const hier::IntVector& smallest_patch_size,
      const hier::IntVector& domain_bounding_box_size) const;

   /*!
    * @brief Initialize a single-block MappedBoxLevel from a multiblock
    * MappedBoxLevel.
    */
   void makeSingleBlockMappedBoxLevel(
      hier::MappedBoxLevel &singleblock_mapped_box_level,
      const hier::MappedBoxLevel &multiblock_mapped_box_level,
      const hier::BlockId &block_id) const;

   /*
    * @brief Record statistics on how many patches and cells were generated.
    */
   void
   recordStatistics(
      double current_time);

   /*!
    * @brief Initialize static objects and register shutdown routine.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   startupCallback();

   /*!
    * @brief Method registered with ShutdownRegister to cleanup statics.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   shutdownCallback();

   /*!
    * @brief Initialize static objects and register shutdown routine.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   initializeCallback();

   /*!
    * @brief Method registered with ShutdownRegister to cleanup statics.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   finalizeCallback();

   /*!
    * @brief The hierarchy that this GriddingAlgorithm works on.
    */
   const tbox::Pointer<hier::PatchHierarchy> d_mb_hierarchy;

   /*!
    * @brief Single-block subsets of the multiblock domain.
    */
   std::vector<hier::MappedBoxLevel> d_singleblock_domain_mapped_box_level;

   /*
    * Static members for managing shared tag data among multiple
    * GriddingAlgorithm objects.
    */
   static tbox::Array<int>* s_tag_indx;
   static tbox::Array<int>* s_buf_tag_indx;

   const tbox::Dimension d_dim;

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
   tbox::Pointer<TagAndInitializeStrategy> d_tag_init_strategy;
   tbox::Pointer<BoxGeneratorStrategy> d_box_generator;
   tbox::Pointer<LoadBalanceStrategy> d_load_balancer;
   tbox::Pointer<LoadBalanceStrategy> d_load_balancer0;
   MultiblockGriddingTagger* d_mb_tagger_strategy;
   bool d_internal_tagger_strategy;

   /*
    * Cell-centered integer variables use to tag cells for refinement.
    * The descriptor index d_tag_indx is used to obtain tag information
    * from user-defined routines on patch interiors.  The descriptor index
    * d_buf_tag_indx is used to buffer tags on patches that may be
    * distributed across processors.  The refine algorithm and schedule are
    * used for interprocessor communication.
    */
   tbox::Pointer<pdat::CellVariable<int> > d_tag;
   tbox::Pointer<pdat::CellVariable<int> > d_buf_tag;
   int d_tag_indx;
   int d_buf_tag_indx;

   tbox::Pointer<xfer::MultiblockRefineAlgorithm> d_mblk_bdry_fill_tags;
   tbox::Array<tbox::Pointer<xfer::MultiblockRefineSchedule> >
   d_bdry_sched_tags;

   /*
    * True and false integer tag values set in constructor and used to
    * set and compare integer tag values on the patch hierarchy.  Generally,
    * these variables are easier to read in the code than integer constants,
    * such as 0 and 1.
    */
   int d_true_tag;
   int d_false_tag;

   /*!
    * @brief Finest level not changed during regrid.
    *
    * This member is temporary, used to coordinate with private methods.
    */
   int d_base_ln;

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
    * @brief When regridding level ln+1, the new level ln must not flow into
    * d_proper_nesting_complement[ln].
    *
    * Has length d_mb_hierarchy->getMaxNumberOfLevels().  The objects
    * are initialized only during gridding/regridding.
    */
   std::vector<std::vector<SAMRAI::hier::MappedBoxLevel> >
   d_proper_nesting_complement;

   /*
    * @brief Connectors from the hierarchy to d_proper_nesting_complement.
    */
   std::vector<std::vector<SAMRAI::hier::Connector> > d_to_nesting_complement;
   /*
    * @brief Connectors from d_proper_nesting_complement to the hierarchy.
    */
   std::vector<std::vector<SAMRAI::hier::Connector> > d_from_nesting_complement;
   tbox::Array<hier::MappedBoxTree> d_domain_complement_tree;

   /*!
    * @brief How to resolve user tags that violate nesting requirements.
    *
    * If a tag violates the nesting requirements, its location in index space
    * will not be refined when creating a finer level.  This flag allows the
    * user to determine what to do when this occurs
    *
    * Values can be:
    * - 'i' Ignore (violating tags will be quietly disregarded)
    * - 'w' Warn (violating tags will cause a warning and be disregarded)
    * - 'e' Error (violating tags will cause an unrecoverable assertion)
    *
    * This defaults to 'w' and set by input parameter
    * "check_nonrefined_tags".
    */
   char d_check_nonrefined_tags;

   /*!
    * @brief Whether or not to check for overlapping patches on a level.
    *
    * This determines whether to check if a new level has any patches that
    * overlap in indes space.
    *
    * Values can be:
    * - 'i' Ignore (overlapping patches will be quietly disregarded)
    * - 'w' Warn (overlapping patches will cause a warning and be disregarded)
    * - 'e' Error (overlapping patches will cause an unrecoverable assertion)
    *
    * This defaults to 'i' and set by input parameter
    * "check_overlapping_patches".
    */
   char d_check_overlapping_patches;

   /*!
    * @brief Whether to globally sequentialize patch indices on every level.
    */
   bool d_sequentialize_patch_indices;

   /*
    * Switches for massaging boxes after clustering.
    * Should be on for most AMR applications.
    * Turning off is mainly for debugging purposes.
    */
   bool d_enforce_proper_nesting;
   bool d_extend_to_domain_boundary;
   bool d_load_balance;

   /*
    * Timers interspersed throughout the class.
    */
   static tbox::Pointer<tbox::Timer> t_find_domain_complement;
   static tbox::Pointer<tbox::Timer> t_load_balance;
   static tbox::Pointer<tbox::Timer> t_load_balance0;
   static tbox::Pointer<tbox::Timer> t_load_balance_setup;
   static tbox::Pointer<tbox::Timer> t_bdry_fill_tags_create;
   static tbox::Pointer<tbox::Timer> t_make_coarsest;
   static tbox::Pointer<tbox::Timer> t_make_finer;
   static tbox::Pointer<tbox::Timer> t_make_finer_setup;
   static tbox::Pointer<tbox::Timer> t_make_finer_tagging;
   static tbox::Pointer<tbox::Timer> t_make_finer_create;
   static tbox::Pointer<tbox::Timer> t_regrid_all_finer;
   static tbox::Pointer<tbox::Timer> t_regrid_finer_create;
   static tbox::Pointer<tbox::Timer> t_bridge_links;
   static tbox::Pointer<tbox::Timer> t_fill_tags_from_mapped_box_level;
   static tbox::Pointer<tbox::Timer> t_tag_cells_for_refinement;
   static tbox::Pointer<tbox::Timer> t_buffer_tags;
   static tbox::Pointer<tbox::Timer> t_bdry_fill_tags_comm;
   static tbox::Pointer<tbox::Timer> t_second_finer_tagging;
   static tbox::Pointer<tbox::Timer> t_find_refinement;
   static tbox::Pointer<tbox::Timer> t_bridge_new_to_new;
   static tbox::Pointer<tbox::Timer> t_find_new_to_new;
   static tbox::Pointer<tbox::Timer> t_bridge_new_to_coarser;
   static tbox::Pointer<tbox::Timer> t_bridge_new_to_finer;
   static tbox::Pointer<tbox::Timer> t_bridge_new_to_old;
   static tbox::Pointer<tbox::Timer> t_find_boxes_containing_tags;
   static tbox::Pointer<tbox::Timer> t_enforce_nesting;
   static tbox::Pointer<tbox::Timer> t_make_nesting_map;
   static tbox::Pointer<tbox::Timer> t_make_nesting_map_compute;
   static tbox::Pointer<tbox::Timer> t_make_nesting_map_convert;
   static tbox::Pointer<tbox::Timer> t_use_nesting_map;
   static tbox::Pointer<tbox::Timer> t_make_overflow_map;
   static tbox::Pointer<tbox::Timer> t_make_overflow_map_compute;
   static tbox::Pointer<tbox::Timer> t_make_overflow_map_convert;
   static tbox::Pointer<tbox::Timer> t_use_overflow_map;
   static tbox::Pointer<tbox::Timer> t_compute_external_parts;
   static tbox::Pointer<tbox::Timer> t_compute_nesting_violator;
   static tbox::Pointer<tbox::Timer> t_extend_to_domain_boundary;
   static tbox::Pointer<tbox::Timer> t_extend_within_domain;
   static tbox::Pointer<tbox::Timer> t_grow_boxes_within_domain;
   static tbox::Pointer<tbox::Timer> t_sort_nodes;
   static tbox::Pointer<tbox::Timer> t_modify_connector;
   static tbox::Pointer<tbox::Timer> t_misc1;
   static tbox::Pointer<tbox::Timer> t_misc2;
   static tbox::Pointer<tbox::Timer> t_misc3;
   static tbox::Pointer<tbox::Timer> t_misc4;
   static tbox::Pointer<tbox::Timer> t_misc5;
   static tbox::Pointer<tbox::Timer> t_make_domain;
   static tbox::Pointer<tbox::Timer> t_get_balance;
   static tbox::Pointer<tbox::Timer> t_use_balance;
   static tbox::Pointer<tbox::Timer> t_make_new;
   static tbox::Pointer<tbox::Timer> t_process_error;
   static tbox::Pointer<tbox::Timer> t_limit_overflow;
   static tbox::Pointer<tbox::Timer> t_reset_hier;
   static tbox::Pointer<tbox::Timer> t_box_massage;

#ifdef MGA_RECORD_STATS
   /*
    * Statistics on number of cells and patches generated.
    */
   tbox::Array<tbox::Pointer<tbox::Statistic> > d_boxes_stat;
   tbox::Array<tbox::Pointer<tbox::Statistic> > d_cells_stat;
   tbox::Array<tbox::Pointer<tbox::Statistic> > d_timestamp_stat;
#endif

   // The following are not yet implemented:
   MultiblockGriddingAlgorithm(
      const MultiblockGriddingAlgorithm&);
   void
   operator = (
      const MultiblockGriddingAlgorithm&);

   // Verbose flags.
   static char s_check_overflow_nesting;
   static char s_check_proper_nesting;
   static char s_check_connectors;
   static char s_print_mapped_box_level_hierarchy;
   static char s_print_steps;

   /* 
    * Static initialization and cleanup handler.
    */
   static tbox::StartupShutdownManager::Handler
   s_initialize_handler;

   /* 
    * 
    */
   static tbox::StartupShutdownManager::Handler
   s_startup_shutdown_handler;

};

}
}
#ifdef SAMRAI_INLINE
#include "SAMRAI/mesh/MultiblockGriddingAlgorithm.I"
#endif
#endif
