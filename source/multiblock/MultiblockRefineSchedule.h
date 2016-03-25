//
// File:        MultiblockRefineSchedule.h
// Package:     SAMRAI multiblock package
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 750 $
// Modified:    $Date: 2005-12-02 10:54:58 -0800 (Fri, 02 Dec 2005) $
// Description: class to manage multiblocks
//

#ifndef included_mblk_MultiblockRefineSchedule
#define included_mblk_MultiblockRefineSchedule

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_mblk_MultiblockPatchHierarchy
#include "MultiblockPatchHierarchy.h"
#endif
#ifndef included_mblk_MultiblockPatchLevel
#include "MultiblockPatchLevel.h"
#endif
#ifndef included_hier_PatchLevel
#include "PatchLevel.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif
#ifndef included_xfer_RefineAlgorithm
#include "RefineAlgorithm.h"
#endif
#ifndef included_xfer_RefineClasses
#include "RefineClasses.h"
#endif
#ifndef included_xfer_RefinePatchStrategy
#include "RefinePatchStrategy.h"
#endif

#ifndef MULTIBLOCK_FAKE_LEVEL_NUMBER
#define MULTIBLOCK_FAKE_LEVEL_NUMBER (-24)
#endif

namespace SAMRAI {
    namespace mblk {

/*!
 * @brief Class MultiblockRefineSchedule<DIM> is an extension of the
 * concept of xfer::RefineSchedule<DIM> to be used in applications that require
 * a multiblock domain.
 *
 * This class contains several constructors called from
 * MultiblockRefineAlgorithm.  In the fillData() routine, it first
 * uses xfer::RefineSchedule<DIM> to fill data within the interiors of each
 * block of the multiblock domain, then communicates or copies data to
 * fill boundary conditions at the boundaries between blocks.
 *
 * @see MultiblockPatchHierarchy
 * @see MultiblockRefineAlgorithm
 * @see xfer::RefineSchedule
 */

template<int DIM> class MultiblockRefinePatchStrategy;
template<int DIM> class MultiblockRefineAlgorithm;

template<int DIM> class MultiblockRefineSchedule
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
    * MultiblockRefineAlgorithm<DIM> object.
    *
    * @param dst_level        Pointer to destination level.
    * @param src_level        Pointer to source level.
    * @param multiblock       Multiblock patch hierarchy object containing all
    *                         of the levels that hold the data being
    *                         communicated
    * @param refine_alg       Pointer to an xfer::RefineAlgorithm<DIM> that
    *                         will be used to create refine schedules that
    *                         will do the data transfers and communication.
    *                         In general, this is a data member of the
    *                         calling MultiblockRefineAlgorithm<DIM>
    *                         object.
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
    */
   MultiblockRefineSchedule(
      tbox::Pointer< MultiblockPatchLevel<DIM> > dst_level,
      tbox::Pointer< MultiblockPatchLevel<DIM> > src_level,
      tbox::Pointer< MultiblockPatchHierarchy<DIM> > multiblock,
      tbox::Pointer< xfer::RefineAlgorithm<DIM> > refine_alg,
      MultiblockRefinePatchStrategy<DIM>* strategy,
      bool use_time_refinement = false);

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
    * MultiblockRefineAlgorithm<DIM> object.
    *
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
    * @param refine_alg       Pointer to an xfer::RefineAlgorithm<DIM> that
    *                         will be used to create refine schedules that
    *                         will do the data transfers and communication.
    *                         In general, this is a data member of the
    *                         calling MultiblockRefineAlgorithm<DIM>
    *                         object.
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
    */
   MultiblockRefineSchedule(
      tbox::Pointer< MultiblockPatchLevel<DIM> > dst_level,
      tbox::Pointer< MultiblockPatchLevel<DIM> > src_level,
      const int next_coarser_level,
      tbox::Pointer< MultiblockPatchHierarchy<DIM> > multiblock,
      tbox::Pointer< xfer::RefineAlgorithm<DIM> > refine_alg,
      MultiblockRefinePatchStrategy<DIM>* strategy,
      bool use_time_refinement = false);

   /*!
    * Virtual destructor
    */
   virtual ~MultiblockRefineSchedule<DIM>();

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
    *                                  createSchedule() function.  Note that
    *                                  even when the value is false, boundary
    *                                  routines may be called on levels coarser
    *                                  than the destination level if such data
    *                                  is needed for proper interpolation.
    * @param filling_coarse_scratch Boolean flag that should only be true
    *                               if fillData is being called recursively
    *                               to fill scratch data on a coarser
    *                               multiblock level.
    * @param filling crse_scr_recursive Boolean flag that should only be true
    *                                   when fillData is called in a second or
    *                                   greater level of recursion 
    */
   void fillData(double fill_time,
                 bool do_physical_boundary_fill = true) const;

   /*!
    * @brief Initialize a vector to contain source data components
    *
    * The component selector argument will be filled with the data component
    * id's that have been registered as source items for the algorithm that
    * created this schedule
    */
   void initializeSourceVector(hier::ComponentSelector& allocate_vector) const;

   /*!
    * @brief Get the equivalence classes associated with the algorithm that
    * created this schedule
    */
   const
   tbox::Pointer< xfer::RefineClasses<DIM> >& getEquivalenceClasses() const;

   /*!
    * @brief Allocate scratch space on the specified level and
    * return the allocated patch data indices in the component
    * selector for later deallocation.
    *
    * @param level            Level on which to allocate scratch space
    * @param fill_time        Simulation time
    * @param allocate_vector  The patch data indices associated with
    *                         scratch space will be stored here to be
    *                         used later for deallocation 
    */
   void allocateScratchSpace(
      tbox::Pointer< hier::PatchLevel<DIM> > level,
      double fill_time,
      hier::ComponentSelector& allocate_vector) const;

   /*!
    * @brief struct SingularityPatch allows a temporary patch that contains
    * data near a singularity to be paired with a block id number.
    *
    * @param d_id     The number of block where the data on the patch exists
    * @param d_patch  A temporary patch that will hold data from the block
    *                 signified by d_id.  Data from this patch will be used
    *                 to fill ghost data at a singularity.
    */ 
   struct SingularityPatch {
      int d_id;
      tbox::Pointer< hier::Patch<DIM> > d_patch;
   };


private:

   void fillData(double fill_time,
                 bool do_physical_boundary_fill,
                 bool filling_coarse_scratch,
                 bool filling_crse_scr_recursive = false) const;

   void createInterblockSchedules(
      tbox::Pointer< MultiblockPatchLevel<DIM> > dst_level,
      tbox::Pointer< MultiblockPatchLevel<DIM> > src_level,
      xfer::RefinePatchStrategy<DIM>* refine_strategy,
      int level_number = MULTIBLOCK_FAKE_LEVEL_NUMBER,
      bool use_time_refinement = false);

   /*
    * Private function that manages the copying of data from src_level of
    * one block to dst_level of another block, using the given rotation and
    * shift.  The levels used in these function are created internally within
    * this Multiblock class.
    */
   void copyBetweenBlocks(
      tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
      const tbox::Pointer< hier::PatchLevel<DIM> > src_level,
      const hier::IntVector<DIM>& shift,
      const typename MultiblockPatchHierarchy<DIM>::RotationIdentifier rotate,
      const tbox::Pointer< xfer::RefineClasses<DIM> > refine_classes) const;

   /*
    * Use patch strategy to fill ghost regions of level
    * at the singularity point(s) of the block indicated by block_number.
    */
   void fillSingularityBoundary(
      tbox::Pointer< hier::PatchLevel<DIM> >& level,
      tbox::Array< tbox::List<SingularityPatch> >&
         singularity_patches,
      const int block_number,
      const double fill_time) const;

   bool needOtherSourceBlocks(hier::BoxList<DIM>& dst_boxes,
                              hier::BoxList<DIM>& src_boxes,
                              hier::BoxList<DIM>& domain_outside_block) const;

   void createCoarseSchedule(
      tbox::Pointer< hier::PatchLevel<DIM> >& fine_level,
      int next_coarser_level,
      const hier::IntVector<DIM>& ratio_to_coarser,
      tbox::Pointer< hier::PatchHierarchy<DIM> >& hierarchy,
      int block_number);

   void createNeighborCoarseSchedule(
      tbox::Pointer< hier::PatchLevel<DIM> >& fine_level,
      int next_coarser_level,
      const hier::IntVector<DIM>& ratio_to_coarser,
      tbox::Pointer< hier::PatchHierarchy<DIM> >& hierarchy,
      int neighbor_block_number,
      int dst_block_number,
      int neighbor_counter);

   void findUnfilledBoxes(
      hier::BoxList<DIM>& unfilled_boxes,
      const int block_number,
      tbox::Pointer< MultiblockPatchLevel<DIM> > coarse_hierarchy_level,
      const hier::BoxList<DIM>& pseudo_domain,
      const hier::IntVector<DIM>& gcw);

   void constructScratchRefineAlgorithm();

   void refineScratchData(
      tbox::Pointer< MultiblockPatchLevel<DIM> > coarse_level,
      tbox::Pointer< hier::PatchLevel<DIM> > fine_level,
      const hier::BoxList<DIM>& unfilled_boxes,
      const int block_number,
      const bool use_fine_gcw) const;

   void copyScratchToDestination(
      tbox::Pointer< hier::PatchLevel<DIM> > level,
      const hier::BoxList<DIM>& unfilled_boxes,
      const tbox::Pointer< xfer::RefineClasses<DIM> > refine_classes) const;

   tbox::Pointer< hier::BoxOverlap<DIM> > calculateOverlap(
      const hier::Patch<DIM>& dst_patch,
      const hier::Patch<DIM>& src_patch,
      const xfer::RefineClasses<DIM>& refine_classes,
      const int refine_class_id) const;

   hier::IntVector<DIM> getBoundaryFillGhostWidth() const;

   void initializeDestinationVector(hier::ComponentSelector& dst_vector) const;

   tbox::Pointer< MultiblockPatchHierarchy<DIM> > d_multiblock_hierarchy;
   tbox::Pointer< xfer::RefineAlgorithm<DIM> > d_single_block_refine_alg;
   tbox::Pointer< xfer::RefineAlgorithm<DIM> > d_single_block_scratch_refine_alg;

   tbox::Array< tbox::Pointer< xfer::RefineSchedule<DIM> > >
      d_single_block_fill_local;
   tbox::Array< tbox::Pointer< MultiblockRefineSchedule<DIM> > >
      d_multiblock_coarse_schedule;
   tbox::Array< tbox::Pointer<  hier::ComponentSelector > > d_coarse_selector;

   tbox::Array<bool> d_local_fill_only;
   tbox::Array<hier::BoxList<DIM> > d_unfilled_boxes;
   tbox::Array< tbox::Pointer< MultiblockPatchLevel<DIM> > >
      d_multiblock_coarse_scratch_level;

   tbox::Array< tbox::Array<bool> > d_neighbor_copy_only;
   tbox::Array< tbox::Array<hier::BoxList<DIM> > > d_neighbor_unfilled_boxes;
   tbox::Array< tbox::Array< tbox::Pointer< MultiblockPatchLevel<DIM> > > >
      d_neighbor_multiblock_coarse_level;
   tbox::Array< tbox::Array< tbox::Pointer< MultiblockRefineSchedule<DIM> > > >
      d_neighbor_multiblock_coarse_schedule;

   tbox::Array< tbox::Array< tbox::Pointer< hier::PatchLevel<DIM> > > >
   d_neighbor_ghost_level;

   tbox::Array< tbox::Array< tbox::Pointer< hier::PatchLevel<DIM> > > >
   d_finalize_ghost_level;

   tbox::Array< tbox::Array< tbox::Array<int> > >
   d_finalize_ghost_patch_numbers;
   tbox::Array< tbox::Array< tbox::Array<int> > >
   d_finalize_ghost_num_src_patches;

   tbox::Array< tbox::Array< tbox::Pointer< xfer::RefineSchedule<DIM> > > >
   d_neighbor_single_block_refine_schedule;

   tbox::Pointer< MultiblockPatchLevel<DIM> > d_multiblock_dst_level;

   MultiblockRefinePatchStrategy<DIM>* d_multiblock_strategy;

};

}
}

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "MultiblockRefineSchedule.C"
#endif

#endif
