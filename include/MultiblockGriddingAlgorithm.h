//
// File:        MultiblockGriddingAlgorithm.h.sed
// Package:     SAMRAI multiblock
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 459 $
// Modified:    $Date: 2005-06-16 14:58:16 -0700 (Thu, 16 Jun 2005) $
// Description: AMR hierarchy generation and regridding routines.
//

#ifndef included_mblk_MultiblockGriddingAlgorithm
#define included_mblk_MultiblockGriddingAlgorithm

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_iostream
#define included_iostream
#include <iostream>
using namespace std;
#endif
#ifndef included_mesh_BaseGriddingAlgorithm
#include "BaseGriddingAlgorithm.h"
#endif
#ifndef included_mesh_BoxGeneratorStrategy
#include "BoxGeneratorStrategy.h"
#endif
#ifndef included_hier_BoxIOUtility
#include "BoxIOUtility.h"
#endif
#ifndef included_pdat_CellVariable
#include "CellVariable.h"
#endif
#ifndef included_mesh_LoadBalanceStrategy
#include "LoadBalanceStrategy.h"
#endif
#ifndef included_mblk_MultiblockGriddingTagger
#include "MultiblockGriddingTagger.h"
#endif
#ifndef included_mblk_MultiblockPatchHierarchy
#include "MultiblockPatchHierarchy.h"
#endif
#ifndef included_mblk_MultiblockRefineAlgorithm
#include "MultiblockRefineAlgorithm.h"
#endif
#ifndef include_hier_PatchLevel
#include "PatchLevel.h"
#endif
#ifndef included_mesh_TagAndInitializeStrategy
#include "TagAndInitializeStrategy.h"
#endif

namespace SAMRAI {
    namespace mblk {

/*!
 * @brief Class MultiblockGriddingAlgorithm<DIM> manages gridding operations in
 * SAMRAI for problems on multiblock domains.  Specifically, it provides AMR
 * patch hierarchy generation and regridding routines that may be used with
 * a variety of AMR solution algorithms and application-specific numerical
 * routines.
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
 * type mesh::TagAndInitializeStrategy<DIM>.  Operations that cluster
 * tagged cells into a collection of box regions are provided by the
 * mesh::BoxGeneratorStrategy<DIM> data member.  Routines that load balancing
 * patches on each level are provided by the mesh::LoadBalanceStrategy<DIM> data
 * member.  The collaboration between this class and each of those objects
 * follows the Strategy design pattern.  Each instantiation of this gridding
 * algorithm class is configured with concrete implementations of those
 * routines by passing appropriate objects into this constructor.
 *
 * Initialization of an MultiblockGriddingAlgorithm<DIM> object is performed
 * via a combination of default parameters and values read from input.  Data
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
 *       vector must have length DIM.  See sample input below.
 *
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
 * @see mesh::GriddingAlgorithm
 * @see mesh::TagAndInitializeStrategy
 * @see mesh::LoadBalanceStrategy
 * @see mesh::BoxGeneratorStrategy
 * @see mesh::MultiblockGriddingTagger
 */

template<int DIM>
class MultiblockGriddingAlgorithm 
: public mesh::BaseGriddingAlgorithm<DIM>
{

public:

   /*!
    * @brief Construct the multiblock gridding algorithm.
    *
    * The constructor for MultiblockGriddingAlgorithm<DIM> configures the
    * gridding algorithm with the concrete strategy objects in the argument
    * list.  Gridding parameters are initialized from values provided in the
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
    *
    * @param object_name     Name of object, to be used by RestartManager
    * @param input_db        Input database
    * @param multiblock      Multiblock patch hierarchy for the application
    * @param level_strategy  Pointer to object that manages cell tagging
    * @param generator       Pointer to object that generates boxes for a level
    * @param balancer        Pointer to load-balancing object
    * @param mb_tagger_strategy Optional pointer to multiblock tagging strategy 
    *                           used to communicate tags across block
    *                           boundaries.  If not provided, a default will
    *                           be used.  Should be NULL unless user has
    *                           implemented a class inheriting from
    *                           MultiblockGriddingTagger. 
    * @param register_for_restart Optional boolean flag indicating whether
    *                             object will be written to restart files.  The
    *                             default is true indicating that
    *                             the object will be written to resatrt files.
    */
   MultiblockGriddingAlgorithm(
      const string& object_name,
      tbox::Pointer<tbox::Database> input_db,
      tbox::Pointer< MultiblockPatchHierarchy<DIM> > multiblock,
      tbox::Pointer< mesh::TagAndInitializeStrategy<DIM> > level_strategy,
      tbox::Pointer< mesh::BoxGeneratorStrategy<DIM> > generator,
      tbox::Pointer< mesh::LoadBalanceStrategy<DIM> > balancer,
      MultiblockGriddingTagger<DIM>* mb_tagger_strategy = (MultiblockGriddingTagger<DIM>*)NULL, 
      bool register_for_restart = true);

   /*!
    * Virtual destructor for MultiblockGriddingAlgorithm<DIM>.
    */
   virtual ~MultiblockGriddingAlgorithm<DIM>();

   /*!
    * @brief Create level 0 for a hierarchy
    *
    * This routine will attempt to construct the coarsest level in an AMR
    * multiblock patch hierarchy (i.e., level 0).  If level 0 does not already
    * exist, then the domain specification is checked against the constraints
    * of the grid generation procedures.  The level gridding strategy data
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
    * @param multiblock Multiblock patch hierarchy on which level is created
    * @param level_time Simulation time for data on the new level
    */
   virtual void makeCoarsestLevel(
      tbox::Pointer< hier::BasePatchHierarchy<DIM> > multiblock,
      const double level_time);

   /*!
    * @brief Create new finer level in hierarchy
    *
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
    * are applied to the function arguments.  If any check is violated,
    * an unrecoverable assertion will result.  In particular, the hierarchy
    * pointer must be non-NULL and the given level number must match that of
    * the finest level currently residing in the hierarchy.  Also, the
    * the tag buffer must be positive.
    *
    * @param multiblock        Multiblock patch hierarchy on which level is
    *                          created
    * @param level_time        Simulation time at which this level is created
    * @param initial_time      Boolean flag that should be true only if
    *                          level_time is the initial time of the simulation
    * @param tag_buffer        Number of cells to buffer cells tagged for
    *                          refinement
    * @param regrid_start_time Time of the previous regrid
    *
    */
   virtual void makeFinerLevel(
      tbox::Pointer< hier::BasePatchHierarchy<DIM> > multiblock,
      const double level_time,
      const bool initial_time,
      const int tag_buffer,
      const double regrid_start_time);

   /*!
    * @brief Regrid all levels finer that a specified level
    *
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
    * The boolean argument level_is_coarsest_to_sync is used for regridding in
    * time-dependent problems.  When true, it indicates that the specified
    * level is the coarsest level to synchronize at the current regrid time
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
    *
    * @param multiblock                Multiblock patch hierarchy where the
    *                                  levels exist
    * @param level_number              Level number of a coarse level for which
    *                                  all finer levels will be regridded
    * @param regrid_time               Simulation time when regridding occurs 
    * @param tag_buffer                Array that stores the tag buffer around
    *                                  tagged cells that will be applied at
    *                                  each level. 
    * @param regrid_start_time         Array that stores the time that each
    *                                  level of the hierarchy was last
    *                                  regridded
    * @param level_is_coarsest_to_sync See comments above
    */
   virtual void regridAllFinerLevels(
      tbox::Pointer< hier::BasePatchHierarchy<DIM> > multiblock,
      const int level_number,
      const double regrid_time,
      const tbox::Array<int>& tag_buffer,
      tbox::Array<double> regrid_start_time = tbox::Array<double>(),
      const bool level_is_coarsest_to_sync = true);

   /*! @brief Return true if error estimation process uses time integration;
    *  otherwise, return false.
    */
   virtual bool
   errorEstimationUsesTimeIntegration () const;

   /*!
    * @brief Return the error coarsen ratio.
    *
    * This is needed for cases where an error estimation schem uses time
    * integration (e.g. Richardson extrapolation) to determine how many time
    * levels to maintain to properly apply the estimtion scheme. In general, an
    * even refine ratio (e.g. 2, 4, 8) will maintain two time levels, while an
    * odd refine ratio (e.g. 3) will maintain three.
    */
   virtual int
   getErrorCoarsenRatio () const;

   /*!
    * @brief Return true if level associated with the specified level number can    * be refined; i.e., the level number is less than that of the finest
    * level allowed in the hierarchy.  Otherwise, false is returned.
    */
   virtual bool levelCanBeRefined(const int level_number) const;

   /*!
    * @brief Return pointer to level gridding strategy data member.
    */
   virtual
   tbox::Pointer< mesh::TagAndInitializeStrategy<DIM> >
      getTagAndInitializeStrategy () const;

   /*!
    * @brief Return maximum number of levels allowed in hierarchy.
    */
   virtual int getMaxLevels() const;

   /*!
    * @brief Return const reference to ratio between specified level and next
    * coarser.
    */
   virtual
   const hier::IntVector<DIM>& getRatioToCoarserLevel(int level_number) const;

   /*!
    * @brief Return efficiency tolerance for clustering tags on level.
    *
    */
   virtual double getEfficiencyTolerance(int level_number) const;

   /*!
    * @brief Return combine efficiency for clustering tags on level.
    */
   virtual double getCombineEfficiency(int level_number) const;

   /*!
    * @brief Return proper nesting buffer width for level.
    */
   virtual int getProperNestingBuffer(int level_number) const;

   /*!
    * @brief Return const reference to smallest patch size for level.
    */
   virtual const hier::IntVector<DIM>&
      getSmallestPatchSize(const int level_number) const;

   /*!
    * @brief Return const reference to largest patch size for level.
    */
   virtual const hier::IntVector<DIM>&
      getLargestPatchSize(const int level_number) const;

   /*!
    * @brief Write object state out to the given database.
    *
    * When assertion checking is active, the database pointer must be non-null.
    */
   virtual void putToDatabase(tbox::Pointer<tbox::Database> db);

private:

   /*!
    * @brief Read input data from database and initialize class members.
    *
    * When assertion checking is active, the database pointer must be non-null.
    *
    * @param db              input database
    * @param is_from_restart set to true if simulation needs to be initialized
    *                         from restart
    */
   void getFromInput(tbox::Pointer<tbox::Database> db, bool is_from_restart);

   /*!
    * @brief Initialize members from restart file
    *
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

   /*!
    * @brief Determine nesting boxes for a level
    *
    * Recursively determine proper nesting boxes for specified level and
    * each finer level in the hierarchy.  The proper nesting boxes for a
    * level are an array of boxes whose union represents the allowable
    * extent of the next finer level within the interior of the level.
    *
    * @param hierarchy    Patch hierarchy where the level exists
    * @param level_number Level to be nested
    * @param block_number Block where this level exists
    * @param complement   Output box list containing nesting boxes
    */
   void findProperNestingBoxes(
        const tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
        const int level_number,
        const int block_number,
        hier::BoxList<DIM>& complement);

   /*!
    * @brief Set tags on specified boxes
    *
    * Set integer tags to specified value on intersection between patch
    * level and given box array.  The index value corresponds to the patch
    * descriptor entry of the cell-centered integer tag array.  The boolean
    * flag indicates whether tags are to be set on the regions corresponding
    * to the interior of the level only, if the tag arrays contain ghosts.
    *
    * @param tag_value Value of the tags that will be set
    * @param level     Level on which tags will be set
    * @param index     Patch data index corresponding to data that stores tags
    * @param boxes     Tags will be set on all cells contained in the
    *                  intersection of these boxes and the patch level
    * @param interior_only If true, tags will be set only on the interior
    *                      of the tag patch data, else they will be set also
    *                      on ghost regions.
    */
   void setTagsOnLevel(
      const int tag_value,
      const tbox::Pointer< hier::PatchLevel<DIM> > level,
      const int index,
      const hier::BoxArray<DIM>& boxes,
      const bool interior_only = true) const;

   /*!
    * @brief Buffer each integer tag on patch level matching given tag value
    * with a border of matching tags.
    *
    * @param tag_value   If a cell has been tagged with this value, a buffer
    *                    of surrounding cells will also be tagged with this
    *                    value.
    * @param level       Multiblock patch level on which tags will be buffered.
    * @param buffer_size Size of buffer to be created. 
    */
   void bufferTagsOnLevel(
      const int tag_value,
      const tbox::Pointer< MultiblockPatchLevel<DIM> > level,
      const int buffer_size) const;

   /*!
    * @brief Find boxes where a fine level will be created.
    *
    * Given a hierarchy and a level number, determine an array of boxes from
    * which a refinement of the level may be constructed.   It is assumed
    * that the integer tags that identify cells for refinement have been set
    * on the level before this routine is called.  Note that this routine
    * also returns the processor mapping for the fine boxes.  The processor
    * mapping indicates the assignments of the new patches to processors
    * when the new fine patch level is constructed.  The processor mapping
    * and final box array are set by the mesh::LoadBalanceStrategy<DIM> data
    * member.
    *
    * @param fine_boxes   The boxes on which a new fine level will be made
    * @param mapping      Mapping of new patches to processors
    * @param hierarchy    Patch hierarchy being operated on
    * @param coarse_level_number Boxes are being created for level finer than
    *                            this level.
    * @param nesting_boxes Fine boxes will be nested within these boxes
    */
   void findRefinementBoxes(
      hier::BoxArray<DIM>& fine_boxes,
      hier::ProcessorMapping& mapping,
      const tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
      const int coarse_level_number,
      const hier::BoxList<DIM>& nesting_boxes) const;

   /*!
    * @brief Regrid an existing level and all finer levels
    *
    * Recursively regrid the specified hierarchy level and all finer levels
    * in the hierarchy.  This private member function is invoked by the
    * regridAllFinerLevels() routine.
    *
    * @param hierarchy     Hierarchy containing the levels
    * @param level_number  Level number to be regridded
    * @param regrid_time   Simulation time when called
    * @param finest_level_not_regridded Level number of the finest level that
    *                                   has not yet been regridded
    * @param level_is_coarsest_to_sync Set to true if the level specified
    *                                  by level_number is the coarsest
    *                                  level that needs to be synchronized
    * @param tag_buffer    Array that stored the size of tag buffers for each
    *                      level
    * @param regrid_start_time Array that stores the times of the last time
    *                          each level was regridded
    */
   void regridFinerLevel(
      tbox::Pointer< MultiblockPatchHierarchy<DIM> > hierarchy,
      const int level_number,
      const double regrid_time,
      const int finest_level_not_regridded,
      const bool level_is_coarsest_to_sync,
      const tbox::Array<int>& tag_buffer,
      const tbox::Array<double>& regrid_start_time = NULL);

   /*!
    * @brief Get parameters that control gridding process
    *
    * Set patch size and ghost cell information needed to create new
    * refinement levels.  The maximum number of ghost cells needed in any
    * variable is used to compute the smallest patch size allowed and the
    * extent to which patches may be extended to touch the physical boundary.
    * This avoids problems in setting ghost cell data that may occur when
    * ghost cell regions intersect the physical boundary in strange ways.
    *
    * This routine sets the smallest and largest patch sizes for the specified
    * level, the smallest box to refine on the next coarser level, and the
    * number of cells that a patch may be extended to the boundary if it
    * sufficiently close to the boundary (extend_ghosts).
    *
    * @param smallest_patch         Smallest size allowed for a patch
    * @param smallest_box_to_refine Smallest box that can be refined on the
    *                               next coarser level and still fall
    *                               within the smallest and largest patch
    *                               parameters.
    * @param largest_patch          Largest size allowed for a patch
    * @param extend_ghosts          If distance between a patch and a domain
    *                               boundary is less than this parameter,
    *                               the patch will be extended to the boundary
    * @param hierarchy              Patch hierarchy containing patch levels
    * @param level_number           Number of coarser level from which a finer
    *                               level is being computed
    */
   virtual void getGriddingParameters(
      hier::IntVector<DIM>& smallest_patch,
      hier::IntVector<DIM>& smallest_box_to_refine,
      hier::IntVector<DIM>& largest_patch,
      hier::IntVector<DIM>& extend_ghosts,
      const tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
      const int level_number) const;

   /*
    * Static members for managing shared tag data among multiple
    * MultiblockGriddingAlgorithm objects.
    */
   static int s_tag_indx;
   static int s_buf_tag_indx;
   static int s_buf_tag_src_indx;

   string d_object_name;
   bool d_registered_for_restart;

   int d_max_levels;

   tbox::Array<hier::IntVector<DIM> > d_ratio_to_coarser;

   /*
    * Data members that manage application-specific level initialization
    * and cell tagging, clustering of tagged cells into boxes, and load
    * balancing of patches to processors, respectively.
    */
   tbox::Pointer< mesh::TagAndInitializeStrategy<DIM> >  d_tag_init_strategy;
   tbox::Pointer< mesh::BoxGeneratorStrategy<DIM> >   d_box_generator;
   tbox::Pointer< mesh::LoadBalanceStrategy<DIM> >    d_load_balancer;
   MultiblockGriddingTagger<DIM>* d_mb_tagger_strategy;
   bool d_internal_tagger_strategy;

   /*
    * Cell-centered integer variables use to tag cells for refinement.
    * The descriptor index d_tag_indx is used to obtain tag information
    * from user-defined routines on patch interiors.  The descriptor index
    * d_buf_tag_indx is used to buffer tags on patches that may be
    * distributed across processors.   The descriptor index
    * d_buf_tag_src_indx is used to transfer tag data between patches on different
    * patch hierarchies in the multiblock configuration. 
    */
   tbox::Pointer< pdat::CellVariable<DIM,int> > d_tag;
   tbox::Pointer< pdat::CellVariable<DIM,int> > d_buf_tag;
   tbox::Pointer< pdat::CellVariable<DIM,int> > d_buf_src_tag;
   int d_tag_indx;
   int d_buf_tag_indx;
   int d_buf_tag_src_indx;

   tbox::Array<hier::IntVector<DIM> > d_smallest_patch_size;
   tbox::Array<hier::IntVector<DIM> > d_largest_patch_size;

   /*
    * The following data members are used for writing and reading refinement
    * boxes for the case where one wants to specify refinement rather than
    * have the TagAndInitializeStrategy do it. The two booleans
    * specify whether a user wants to dump or read a set of refine boxes.
    * These may be set on input and are false by default.
    */
   bool d_write_dumped_level_boxes;
   bool d_read_dumped_level_boxes;
   string d_regrid_boxes_filename;
   tbox::Array< tbox::Array<int> > d_regrid_box_counter;
   hier::BoxIOUtility<DIM>* d_regrid_box_utility;

   tbox::Array<int> d_proper_nesting_buffer;
   bool d_allow_patches_smaller_than_ghostwidth;
   tbox::Array<double> d_efficiency_tolerance;
   tbox::Array<double> d_combine_efficiency;

   tbox::Pointer< MultiblockRefineAlgorithm<DIM> > d_bdry_fill_tags;
   tbox::Array< tbox::Pointer< MultiblockRefineSchedule<DIM> > > d_bdry_sched_tags;

   tbox::Array< tbox::Array<hier::BoxList<DIM> > > d_proper_nesting_boxes;

   /*
    * True and false integer tag values set in constructor and used to
    * set and compare integer tag values on the patch hierarchy.  Generally,
    * these variables are easier to read in the code than integer constants,
    * such as 0 and 1.
    */
   int d_true_tag;
   int d_false_tag;

   // The following are not yet implemented:
   MultiblockGriddingAlgorithm(const MultiblockGriddingAlgorithm<DIM>&);
   void operator=(const MultiblockGriddingAlgorithm<DIM>&);

   /*
    * Timers interspersed throughout the class (those timers used
    * only locally within a function are defined in that method).
    */
   tbox::Pointer<tbox::Timer> t_find_proper_nesting;
   tbox::Pointer<tbox::Timer> t_intersect_boxes_find_refinement;
   tbox::Pointer<tbox::Timer> t_load_balance_boxes;
   tbox::Pointer<tbox::Timer> t_bdry_fill_tags_create;

};

}
}

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "MultiblockGriddingAlgorithm.C"
#endif

#endif
