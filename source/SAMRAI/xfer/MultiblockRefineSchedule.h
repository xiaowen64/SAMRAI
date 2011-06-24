/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   class to manage multiblocks 
 *
 ************************************************************************/

#ifndef included_xfer_MultiblockRefineSchedule
#define included_xfer_MultiblockRefineSchedule

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/xfer/PatchLevelFillPattern.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"

namespace SAMRAI {
namespace xfer {

/*!
 * @brief Class MultiblockRefineSchedule is an extension of the
 * concept of xfer::RefineSchedule to be used in applications that require
 * a multiblock domain.
 *
 * This class contains two constructors called from MultiblockRefineAlgorithm
 * for different refinement patterns.  In the fillData() routine, it first
 * uses xfer::RefineSchedule to fill data within the interiors of each
 * block of the multiblock domain, then communicates or copies data to
 * fill boundary conditions at the boundaries between blocks.
 *
 * @see PatchHierarchy
 * @see MultiblockRefineAlgorithm
 * @see RefineSchedule
 */

class MultiblockRefinePatchStrategy;
class MultiblockRefineAlgorithm;

class MultiblockRefineSchedule
{

public:
   /*!
    * Constructor to create a MultiblockRefineSchedule that copies data
    * from the interiors of source patch data components on the source level
    * into the interiors and ghost cells of destination patch data components
    * on the destination level.  Only data on the intersection of the
    * source and destination patch components will be copied.  The source
    * and destination patch levels must reside in the same Multiblock domain.
    * In general, this constructor is called by a
    * MultiblockRefineAlgorithm object.
    *
    * @param fill_pattern Indicates which parts of the destination level
    *                     to fill.  See RefineSchedule for valid values.
    * @param dst_level        Pointer to destination level.
    * @param src_level        Pointer to source level.
    * @param multiblock       Multiblock patch hierarchy object containing all
    *                         of the levels that hold the data being
    *                         communicated
    * @param refine_alg       Pointer to an xfer::RefineAlgorithm that
    *                         will be used to create refine schedules that
    *                         will do the data transfers and communication.
    *                         In general, this is a data member of the
    *                         calling MultiblockRefineAlgorithm
    *                         object.
    * @param transaction_factory Transaction factory.
    * @param strategy    Pointer to a multiblock patch strategy object
    *                         that provides user-defined boundary filling
    *                         operations for patch boundaries that touch a
    *                         multiblock singularity point, as well as user-
    *                         defined physical boundary filling operations.
    *                         This pointer may be null, in which case no
    *                         boundary filling operations will occur.
    *                         case no boundary filling operations will occur.
    * @param use_time_refinement Let the destination level be filled using
    *                            time refinement operations.  This defaults
    *                            to false because it should only be used
    *                            in recursive calls within this class
    * @param enable_singularity_patches Enable the use of the singularity
    *                                   patch functionality for filling
    *                                   ghost data at reduced or enhanced
    *                                   connectivity patch boundaries.
    */
   explicit MultiblockRefineSchedule(
      tbox::Pointer<PatchLevelFillPattern> pl_fill_pattern,
      tbox::Pointer<hier::PatchLevel> dst_level,
      tbox::Pointer<hier::PatchLevel> src_level,
      tbox::Pointer<hier::PatchHierarchy> hierarchy,
      tbox::Pointer<xfer::RefineAlgorithm> refine_alg,
      tbox::Pointer<xfer::RefineTransactionFactory> transaction_factory,
      MultiblockRefinePatchStrategy* strategy = NULL,
      bool use_time_refinement = false,
      bool enable_singularity_patches = true);

   /*!
    * Constructor to create a MultiblockRefineShedule that moves data
    * from the interiors of source patch data components on the source level
    * and coarser levels in the patch hierarchy into the interiors and ghost
    * cells of destination patch data components on the destination level.
    * Only data on the intersection of the source and destination patch
    * components will be copied.  If portions of the destination level
    * remain unfilled, then the algorithm recursively fills those unfilled
    * portions from coarser levels in the AMR hierarchy.  The source
    * and destination patch levels must reside in the same index space.
    * However, the levels do not have to be in the same AMR patch hierarchy.
    * In general, this constructor is called by a
    * MultiblockRefineAlgorithm object.
    *
    * @param fill_pattern Indicates which parts of the destination level
    *                     to fill.  See RefineSchedule for valid values.
    * @param dst_level        Pointer to destination level.
    * @param src_level        Pointer to source level.  This pointer may be
    *                         null, in which case the destination level will
    *                         be filled only using data interpolated from
    *                         coarser levels.
    * @param next_coarser_level Integer number of next coarser level in
    *                           relative to the destination level.  Note that
    *                           when the destination level has number zero
    *                           (i.e., the coarsest level), this value should
    *                           be < 0.
    * @param multiblock       Multiblock object containing all of the levels
    *                         that hold the data being communicated
    * @param refine_alg       Pointer to an xfer::RefineAlgorithm that
    *                         will be used to create refine schedules that
    *                         will do the data transfers and communication.
    *                         In general, this is a data member of the
    *                         calling MultiblockRefineAlgorithm
    *                         object.
    * @param transaction_factory Transaction factory.
    * @param strategy    Pointer to a multiblock patch strategy object
    *                         that provides user-defined boundary filling
    *                         operations for patch boundaries that touch a
    *                         multiblock singularity point, as well as user-
    *                         defined physical boundary filling operations.
    *                         This pointer may be null, in which case no
    *                         boundary filling operations will occur.
    *                         case no boundary filling operations will occur.
    * @param use_time_refinement Let the destination level be filled using
    *                            time refinement operations.  This defaults
    *                            to false because it should only be used
    *                            in recursive calls within this class
    * @param enable_singularity_patches Enable the use of the singularity
    *                                   patch functionality for filling
    *                                   ghost data at reduced or enhanced
    *                                   connectivity patch boundaries.
    */
   explicit MultiblockRefineSchedule(
      tbox::Pointer<PatchLevelFillPattern> pl_fill_pattern,
      tbox::Pointer<hier::PatchLevel> dst_level,
      tbox::Pointer<hier::PatchLevel> src_level,
      const int next_coarser_level,
      tbox::Pointer<hier::PatchHierarchy> hierarchy,
      tbox::Pointer<xfer::RefineAlgorithm> refine_alg,
      tbox::Pointer<xfer::RefineTransactionFactory> transaction_factory,
      MultiblockRefinePatchStrategy* strategy = NULL,
      bool use_time_refinement = false,
      bool enable_singularity_patches = true);

   /*!
    * Virtual destructor
    */
   virtual ~MultiblockRefineSchedule();

   /*!
    * @brief Reset this refine schedule to perform data transfers
    * asssociated with refine class items in function argument.
    *
    * In general, this function is
    * called by a MultiblockRefineAlgorithm object.
    *
    * @param refine_algorithm  Pointer to a RefineAlgorithm that holds
    *                          the refine items that will be used by the
    *                          schedule after it is reset
    */
   void
   reset(
      const tbox::Pointer<xfer::RefineAlgorithm> refine_algorithm);

   /*!
    * @brief Execute the stored communication schedule and perform the
    * data movement.
    *
    * @param fill_time Double simulation time when the fill take place.
    * @param do_physical_boundary_fill Boolean flag used to bypass the
    *                                  physical boundary data filling
    *                                  operations on the destination level.
    *                                  The default value is true indicating
    *                                  that boundary data will be filled
    *                                  (assuming a non-null refine patch
    *                                  strategy pointer was passed to the
    *                                  constructor.  Note that even when the
    *                                  value is false, boundary filling methods
    *                                  may be called on levels coarser
    *                                  than the destination level if such data
    *                                  is needed for proper interpolation.
    */
   void
   fillData(
      double fill_time,
      bool do_physical_boundary_fill = true) const;

   /*!
    * @brief Get the equivalence classes associated with the algorithm that
    * created this schedule
    */
   const
   tbox::Pointer<xfer::RefineClasses>&
   getEquivalenceClasses() const;

private:
   /*!
    * @brief Initialize a vector to contain source data components
    *
    * The component selector argument will be filled with the data component
    * id's that have been registered as source items for the algorithm that
    * created this schedule
    */

   void
   initializeSourceVector(
      hier::ComponentSelector& allocate_vector) const;

   /*!
    * @brief Allocate scratch space on the specified level and
    * return the allocated patch data indices in the component
    * selector for later deallocation.
    *
    * @param allocate_vector  The patch data indices associated with
    *                         scratch space will be stored here to be
    *                         used later for deallocation
    * @param level            Level on which to allocate scratch space
    * @param fill_time        Simulation time
    */
   void
   allocateScratchSpace(
      hier::ComponentSelector& allocate_vector,
      tbox::Pointer<hier::PatchLevel> level,
      double fill_time) const;

   /*!
    * @brief Private data filling method, called from the public method,
    * or recursively.
    *
    * @param fill_time Double simulation time when the fill take place.
    * @param do_physical_boundary_fill Boolean flag used to bypass the
    *                                  physical boundary data filling
    *                                  operations on the destination level.
    *                                  The default value is true indicating
    *                                  that boundary data will be filled
    *                                  (assuming a non-null refine patch
    *                                  strategy pointer was passed to the
    *                                  constructor.  Note that even when the
    *                                  value is false, boundary filling methods
    *                                  may be called on levels coarser
    *                                  than the destination level if such data
    *                                  is needed for proper interpolation.
    * @param filling_coarse_scratch Boolean flag indicating whether this call
    *                               is filling scratch data on a temporary
    *                               coarse level.
    * @param filling_coarse_scr_recursive Boolean flag indicating if the
    *                                     method is being called recursively
    *                                     to fill scratch data on a coarse
    *                                     level.
    *
    */
   void
   fillData(
      double fill_time,
      bool do_physical_boundary_fill,
      bool filling_coarse_scratch,
      bool filling_crse_scr_recursive = false) const;

   /*!
    * @brief Private method to create the RefineSchedules that transfer data
    * at block boundaries
    *
    * @param dst_level The level where data is being filled
    * @param src_level The source level
    * @param refine_strategy The refine strategy that will be given to
    *                        the constructed RefineSchedules
    * @param level_number The level number being filled.  If dst_level is not
    *                     on a hierarchy, then MULTIBLOCK_FAKE_LEVEL_NUMBER
    *                     should be used
    * @param use_time_refinement Boolean to determine if time refinement is
    *                            to be used
    */
   void
   createInterblockSchedules(
      tbox::Pointer<hier::PatchLevel> dst_level,
      tbox::Pointer<hier::PatchLevel> src_level,
      xfer::RefinePatchStrategy* refine_strategy,
      int level_number,
      bool use_time_refinement);

   /*!
    * @brief Private method to create the BoxOverlaps used in
    * refineScratchData
    * TODO: src_level not needed
    * @param dst_level The level where data is being filled
    * @param src_level The source level
    */
   void
   createOverlapsForRefine(
      tbox::Pointer<hier::PatchLevel>& dst_level,
      tbox::Pointer<hier::PatchLevel>& src_level);

   /*!
    * @brief Private method to create the BoxOverlaps used in
    * copying neighbor data in fillData
    */
   void
   createOverlapsForCopy(
      tbox::Pointer<hier::PatchLevel>& dst_level,
      tbox::Pointer<hier::PatchLevel>& src_level);

   /*!
    * @brief Create BoxOverlaps that will be used in refineScratchData.
    *
    * @param overlaps Container to store the overlaps
    * @param dst_patch_level Destination pactch level on the block being filled
    * @param neighbor Neighbor object representing source block
    * @param block_number Block number of the destination
    * @param fill_fine_gcw Array index identifying the neighbor being used
    * @param equiv_classes Equivalence classes holding all of the items that
    *                      will be copied
    */
   void
   createCopyOverlaps(
      tbox::List<tbox::Pointer<hier::BoxOverlap> >& overlaps,
      const tbox::Pointer<hier::PatchLevel>& dst_patch_level,
      const tbox::Pointer<hier::PatchLevel>& src_patch_level,
      const hier::GridGeometry::Neighbor& neighbor,
      const int block_number,
      const int neighbor_counter,
      const tbox::Pointer<xfer::RefineClasses>& equiv_classes);

   /*
    * Private function that manages the copying of data from src_level of
    * one block to dst_level of another block, using the given rotation and
    * shift.  The levels used in these function are created internally within
    * this Multiblock class.
    */
   void
   copyBetweenBlocks(
      tbox::Pointer<hier::PatchLevel> dst_level,
      const tbox::Pointer<hier::PatchLevel> src_level,
      const hier::IntVector& shift,
      const hier::Transformation::RotationIdentifier rotate,
      const tbox::Pointer<xfer::RefineClasses> refine_classes,
      const hier::BlockId& dst_block_id,
      const hier::BlockId& src_block_id) const;

   /*
    * Private function that manages the filling of data from src_level of
    * one block to dst_level of another block, allowing for a user-defined
    * filling rather than a copy, using the given rotation and
    * shift.  The levels used in these function are created internally within
    * this Multiblock class.
    */
   void
   fillBetweenBlocks(
      tbox::Pointer<hier::PatchLevel> dst_level,
      const tbox::Pointer<hier::PatchLevel> src_level,
      const hier::IntVector& shift,
      const hier::Transformation::RotationIdentifier rotate,
      const tbox::Pointer<xfer::RefineClasses> refine_classes,
      const hier::BlockId& dst_block_id,
      const hier::BlockId& src_block_id) const;

   /*
    * Use patch strategy to fill ghost regions of level
    * at the singularity point(s) of the block indicated by block_number.
    */
   void
   fillSingularityBoundary(
      tbox::Array<tbox::List<tbox::Pointer<hier::Patch> > >&
      singularity_patches,
      const int block_number,
      const double fill_time) const;

   /*!
    * @brief Private method to determine if data from other blocks is needed
    * to fill data on a particular block.
    *
    * @param dst_boxes BoxList describing all of the boxes to be filled on a
    *                  single block.
    * @param src_boxes BoxList describing all available source boxes on the
    *                  same block as dst_boxes that can be used to fill data
    *                  in dst_boxes.
    * @param domain_outside_block BoxList representing the physical domain of
    *                             all neighboring blocks, translated to the
    *                             index space of the destination block.
    */
   bool
   needOtherSourceBlocks(
      hier::BoxList& dst_boxes,
      hier::BoxList& src_boxes,
      hier::BoxList& domain_outside_block) const;

   /*!
    * @brief Private method to create a MultiblockRefineSchedule to fill
    * data on a temporary coarse level.
    *
    * @param fine_level A fine level that we are attempting to fill.
    * @param next_coarser_level The level number of the next coarser level
    */
   void
   createCoarseSchedule(
      const tbox::Pointer<hier::PatchLevel>& fine_level,
      const int next_coarser_level);

   /*!
    * @brief Private method to create a MultiblockRefineSchedule to fill
    * data on a temporary coarse level that represents space on a neighboring
    * block.
    *
    * When filling data on a PatchLevel that exists on a particular block
    * and touches a block boundary, there will be ghost regions along the
    * block boundary to be filled.  We first attempt to find data from the
    * neighboring block at the same level of resolution.  If none exists, then
    * data from a coarser level on the neighboring block must be interpolated.
    * This method creates a MultiblockRefineSchedule that will fill data on
    * that neighboring coarser level before that interpolation is executed.
    *
    * @param fine_level A fine level that we are attempting to fill.
    * @param next_coarser_level The level number of the next coarser level.
    * @param ratio_to_coarser Refinement ratio between fine_level and the
    *                         next coarser level.
    * @param hierarchy The PatchHierarchy where fine_level exists.  This
    *                  represents one block in a multiblock PatchHierarchy.
    */
   void
   createNeighborCoarseSchedule(
      const tbox::Pointer<hier::PatchLevel>& fine_level,
      const int next_coarser_level,
      const hier::IntVector& ratio_to_coarser,
      const tbox::Pointer<hier::PatchHierarchy>& hierarchy);

   /*!
    * @brief Private method to find boxes that cannot be filled by
    * copies or communications within a single level of resolution and must
    * be filled by interpolation.
    *
    * @param unfilled_boxes The output BoxList specifying the unfilled boxes.
    * @param block_number The number of the block where we are filling data.
    * @param coarse_hierarchy_level The next coarsest PatchLevel
    * @param pseudo_domain A BoxList representing the entire multiblock domain
    *                      as if all blocks were in the index space of the
    *                      destination block.
    * @param gcw Maximum ghost width to be filled.
    */
   void
   findUnfilledBoxes(
      hier::BoxList& unfilled_boxes,
      const int block_number,
      tbox::Pointer<hier::PatchLevel> coarse_hierarchy_level,
      const hier::BoxList& pseudo_domain,
      const hier::IntVector& gcw);

   /*!
    * @brief Construct the RefineAlgorithm to be used to create RefineSchedules
    * that fill scratch data on temporary coarse levels.
    *
    * The temporary coarse levels need to be filled with scratch data
    * before being used for interpolation.  createCoarseSchedule() and
    * createNeighborCoarseSchedule() create schedules to fill these levels.
    * The destination data components for these schedules must be the scratch
    * data components for the overall MultiblockRefineSchedule that is being
    * called by the user.  This method constructs a MultiblockRefineAlgorithm
    * that makes the scratch data components also be destination data
    * components, and this algorithm is then used to construct schedules
    * in createCoarseSchedule() and createNeighborCoarseSchedule().
    */
   void
   constructScratchRefineAlgorithm();

   /*!
    * @brief Interpolate data from a temporary coarse level onto unfilled
    * boxes on a finer level.
    *
    * @param coarse_level A temporary coarse level that provides the coarse
    *                     data for the interpolation.
    * @param fine_level A fine level onto which data is interpolated.
    * @param overlaps Container of pre-computed BoxOverlaps
    * @param crs_to_dst Connector to be used for the interpolation
    * @param unfilled_boxes Boxes representing regions where interpolated data
    *                       is needed.
    * @param block_number Block number of the block where fine_level exists.
    * @param fill_fine_gcw Boolean that tells whether to fill ghost regions
    *                      on fine_level
    * @param filling_neighbor True if filling ghost regions on a neighboring
    *                         block
    */
   void
   refineScratchData(
      tbox::Pointer<hier::PatchLevel> coarse_level,
      tbox::Pointer<hier::PatchLevel> fine_level,
      const std::vector<std::vector< tbox::Pointer<hier::BoxOverlap> > >& overlaps,
      const hier::Connector& crs_to_dst,
      const hier::BoxList& unfilled_boxes,
      const int block_number,
      const bool fill_fine_gcw,
      const bool filling_neighbor) const;

   /*!
    * @brief Copy data that has been filled by interpolation operations from
    * scratch space to destination space.
    *
    * @param level The patch level on a particular block where data has been
    *              filled.
    * @param unfilled_boxes The boxes where scratch data has been filled and
    *                       now needs to be copied to the destination space.
    * @param refine_classes The refine equivalence classes being used for this
    *                       refinement operation.
    */
   void
   copyScratchToDestination(
      tbox::Pointer<hier::PatchLevel> level,
      const hier::BoxList& unfilled_boxes,
      const hier::BlockId& block_id,
      const tbox::Pointer<xfer::RefineClasses> refine_classes) const;


   /*!
    * @brief Create BoxOverlaps that will be used in refineScratchData.
    *
    * @param overlaps Container to store the overlaps
    * @param coarse_mb_level A temporary coarse level that provides the coarse
    *                        data for the interpolation.
    * @param fine_level A fine level onto which data is interpolated.
    * @param crs_to_dst Connector to be used for the interpolation
    * @param unfilled_boxes Boxes representing regions where interpolated data
    *                       is needed.
    * @param block_number Block number of the block where fine_level exists.
    * @param fill_fine_gcw Boolean that tells whether to fill ghost regions
    *                      on fine_level
    * @param filling_neighbor True if filling ghost regions on a neighboring
    *                         block
    */
   void createRefineOverlaps(
      std::vector<std::vector< tbox::Pointer<hier::BoxOverlap> > >& overlaps,
      tbox::Pointer<hier::PatchLevel> coarse_mb_level,
      tbox::Pointer<hier::PatchLevel> fine_level,
      const hier::Connector& crs_to_dst,
      const hier::BoxList& unfilled_boxes,
      const int block_number,
      const bool fill_fine_gcw,
      const bool filling_neighbor);


   /*!
    * @brief Calculate overlap between a destination and source patch for
    * a single member of the refine equivalence classes.
    *
    * @param dst_patch The destination patch
    * @param src_patch The source patch
    * @param refine_classes The refine equivalence classes being used for this
    *                       refinement operation.
    * @param restrict_box   A box that creates a restriction on the overlap.
    *                       None of the overlap will lie outside this box.
    */
   tbox::Pointer<hier::BoxOverlap>
   calculateOverlap(
      const hier::Patch& dst_patch,
      const hier::Patch& src_patch,
      const xfer::RefineClasses::Data& refine_classes,
      const hier::Box& restrict_box) const;

   /*!
    * @brief Calculate overlap that will be used to fill a singularity
    * patch for a single member of the refine equivalence classes.
    *
    * @param patch          The Patch from a SingularityPatch object
    * @param refine_classes The refine equivalence classes being used for this
    *                       refinement operation.
    * @param restrict_box   A box that creates a restriction on the overlap.
    *                       None of the overlap will lie outside this box.
    */
   tbox::Pointer<hier::BoxOverlap>
   calculateSingularityPatchOverlap(
      const hier::Patch& patch,
      const xfer::RefineClasses::Data& refine_classes,
      const hier::Box& restrict_box) const;

   /*!
    * @brief Get width of ghost cell region to fill passed to user supplied
    * physical boundary condition routine.
    */
   hier::IntVector
   getBoundaryFillGhostWidth() const;

   /*!
    * @brief Set up ComponentSelector that contains all destination data
    * components for a schedule
    */
   void
   initializeDestinationVector(
      hier::ComponentSelector& dst_vector) const;

   tbox::Pointer<hier::PatchHierarchy> d_multiblock_hierarchy;

   /*
    * Refine algorithms to be used on a single block
    */
   tbox::Pointer<xfer::RefineAlgorithm> d_single_block_refine_alg;
   tbox::Pointer<xfer::RefineAlgorithm> d_single_block_scratch_refine_alg;

   /*
    * Fill pattern identifier and transaction factory used for that pattern
    */
   tbox::Pointer<PatchLevelFillPattern> d_fill_pattern;
   tbox::Pointer<xfer::RefineTransactionFactory> d_transaction_factory;

   /*
    * Schedule arrays indexed by block number
    */
   tbox::Pointer<xfer::RefineSchedule> d_single_block_fill_local;
   tbox::Pointer<MultiblockRefineSchedule> d_multiblock_coarse_schedule;

   /*
    * Component selector to store destination data components.
    */
   tbox::Pointer<hier::ComponentSelector> d_coarse_selector;

   /*
    * Boolean array that tells whether destination level for a particular
    * block can be filled without getting data from neighboring blocks.
    */
   bool d_local_fill_only;

   /*
    * Array to store unfilled boxes on each block that need to be filled
    * by interpolation.
    */
   tbox::Array<hier::BoxList> d_unfilled_boxes;

   /*
    * Array of temporary coarse levels, each of which will be used to
    * interpolated data onto finer destination level.
    */
   tbox::Pointer<hier::PatchLevel> d_multiblock_coarse_scratch_level;

   /*
    * Connectors used to refine data from d_multiblock_coarse_scratch_level
    */
   hier::Connector d_crs_to_dst_connector;

   /*
    * The d_neighbor_* arrays are used to manage the filling of data from block
    * neighbors.
    *
    * Arrays are indexed [block_number][neighbor_number].
    */

   /*
    * The unfilled boxes in ghost regions on a neighboring block, requiring
    * interpolation from a coarser level.
    */
   tbox::Array<tbox::Array<hier::BoxList> > d_neighbor_unfilled_boxes;

   /*
    * The temporary coarser levels to be used for interpolation
    */
   tbox::Pointer<hier::PatchLevel>
   d_neighbor_multiblock_coarse_level;

   /*
    * Connectors used to refine data from d_neighbor_multiblock_coarse_level
    */
   hier::Connector d_neighbor_crs_to_dst_connector;

   /*
    * The schedules to be used to fill temporarly coarser levels
    */
   tbox::Pointer<MultiblockRefineSchedule>
   d_neighbor_multiblock_coarse_schedule;

   /*
    * Levels that represent ghost regions lying in a neighboring block.  They
    * are filled by copies (if d_neighbor_copy_only is true) or interpolation.
    * These levels exist in the neighboring block's index space.
    */
   tbox::Pointer<hier::PatchLevel>
   d_neighbor_ghost_level;

   /*
    * Levels that represent the same regions as d_neighbor_ghost level, but
    * in the destination block's index space.  Data will be translated from
    * d_neighbor_ghost_level to d_finalize_ghost_level.
    */
   tbox::Pointer<hier::PatchLevel> d_finalize_ghost_level;

   /*
    * Probably going away for DLBG
    */
   tbox::Array<tbox::Array<std::vector<int> > >
   d_finalize_ghost_patch_numbers;
   tbox::Array<tbox::Array<std::vector<int> > >
   d_finalize_ghost_num_src_patches;

   /*
    * RefineSchedules to fill data on ghost regions that lie in a neighboring
    * block, when d_copy_neighbor_only is true and no spatial interpolation
    * is needed.
    */
   tbox::Pointer<xfer::RefineSchedule>
   d_neighbor_single_block_refine_schedule;

   tbox::Array<tbox::Array<tbox::List< tbox::Pointer<hier::BoxOverlap> > > >
   d_neighbor_copy_overlaps;

   /*
    * When a schedule is created that does not use a coarser level for
    * interpolation, this is used to make sure overlaps on block boundaries
    * only allow copying data that originated on the source level.
    */ 
   bool d_only_copy_source_intersections;

   /*
    * Containers to hold BoxOverlaps for refineScratchData
    */
   std::vector<std::vector<std::vector< tbox::Pointer<hier::BoxOverlap> > > >
      d_local_refine_overlaps;
   std::vector<std::vector<std::vector<std::vector< tbox::Pointer<hier::BoxOverlap> > > > >
      d_neighbor_refine_overlaps;

   /*
    * True if the transaction factory is a standard factory.
    */
   bool d_using_standard_transaction;

   /*
    * Destination level for the MultiblockRefineSchedule
    */
   tbox::Pointer<hier::PatchLevel> d_multiblock_dst_level;

   /*
    * Strategy class to handle filling of physical boundaries and
    * boundaries at singularities.  Also may have implementations of
    * interpolation operations.
    */
   MultiblockRefinePatchStrategy* d_multiblock_strategy;

   /*
    * ghost cell width used in creating the schedule
    */
   hier::IntVector d_data_fill_gcw;

   /*
    * Tells if singularity patch functionality enabled.
    */
   bool d_enable_singularity_patches;

   /*
    * Constant to be used in cases where destination level is outside the
    * hierarchy and has no level number.
    */
   static const int MULTIBLOCK_FAKE_LEVEL_NUMBER;

};

}
}

#endif
