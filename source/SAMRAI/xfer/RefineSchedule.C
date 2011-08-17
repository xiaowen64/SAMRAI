/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Refine schedule for data transfer between AMR levels 
 *
 ************************************************************************/

#ifndef included_xfer_RefineSchedule_C
#define included_xfer_RefineSchedule_C

#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/xfer/BoxGeometryVariableFillPattern.h"
#include "SAMRAI/xfer/PatchLevelFullFillPattern.h"
#include "SAMRAI/xfer/RefineCopyTransaction.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"
#include "SAMRAI/xfer/RefineScheduleConnectorWidthRequestor.h"
#include "SAMRAI/xfer/RefineTimeTransaction.h"
#include "SAMRAI/hier/BoxGeometry.h"
#include "SAMRAI/hier/BoxOverlap.h"
#include "SAMRAI/hier/BoxUtilities.h"
#include "SAMRAI/hier/ComponentSelector.h"
#include "SAMRAI/hier/MappedBoxLevelConnectorUtils.h"
#include "SAMRAI/hier/MappingConnectorAlgorithm.h"
#include "SAMRAI/hier/MappedBoxSet.h"
#include "SAMRAI/hier/OverlapConnectorAlgorithm.h"
#include "SAMRAI/hier/PeriodicShiftCatalog.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/hier/PatchGeometry.h"
#include "SAMRAI/tbox/AsyncCommPeer.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Utilities.h"

#include <cassert>

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace xfer {

const int RefineSchedule::BIG_GHOST_CELL_WIDTH = 10;

static const std::string logbord;
  static const std::string errbord("E-> ");

bool RefineSchedule::s_extra_debug = false;
bool RefineSchedule::s_barrier_and_time = false;

tbox::Pointer<tbox::Timer> RefineSchedule::t_fill_data;
tbox::Pointer<tbox::Timer> RefineSchedule::t_recursive_fill;
tbox::Pointer<tbox::Timer> RefineSchedule::t_refine_scratch_data;
tbox::Pointer<tbox::Timer> RefineSchedule::t_finish_sched_const;
tbox::Pointer<tbox::Timer> RefineSchedule::t_finish_sched_const_recurse;
tbox::Pointer<tbox::Timer> RefineSchedule::t_gen_comm_sched;
tbox::Pointer<tbox::Timer> RefineSchedule::t_bridge_connector;
tbox::Pointer<tbox::Timer> RefineSchedule::t_modify_connector;
tbox::Pointer<tbox::Timer> RefineSchedule::t_make_seq_map;
tbox::Pointer<tbox::Timer> RefineSchedule::t_shear;
tbox::Pointer<tbox::Timer> RefineSchedule::t_misc1;
tbox::Pointer<tbox::Timer> RefineSchedule::t_barrier_and_time;
tbox::Pointer<tbox::Timer> RefineSchedule::t_get_global_mapped_box_count;
tbox::Pointer<tbox::Timer> RefineSchedule::t_coarse_shear;
tbox::Pointer<tbox::Timer> RefineSchedule::t_setup_supp_mapped_box_level;
tbox::Pointer<tbox::Timer> RefineSchedule::t_misc2;
tbox::Pointer<tbox::Timer> RefineSchedule::t_bridge_supp_hiercoarse;
tbox::Pointer<tbox::Timer> RefineSchedule::t_bridge_dst_hiercoarse;
tbox::Pointer<tbox::Timer> RefineSchedule::t_make_supp_level;
tbox::Pointer<tbox::Timer> RefineSchedule::t_make_supp_to_unfilled;
tbox::Pointer<tbox::Timer> RefineSchedule::t_invert_edges;
tbox::Pointer<tbox::Timer> RefineSchedule::t_construct_send_trans;
tbox::Pointer<tbox::Timer> RefineSchedule::t_construct_recv_trans;

tbox::StartupShutdownManager::Handler
RefineSchedule::s_initialize_finalize_handler(
   RefineSchedule::initializeCallback,
   0,
   0,
   RefineSchedule::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);

/*
 **************************************************************************
 *
 * Create a refine schedule that copies data from the source level into
 * the destination level on the components represented by the refine
 * classes.  Ony data on the intersection of the two levels will be
 * copied. It is assumed that the index spaces of the source and
 * destination levels are "consistent"; i.e., they represent the same
 * grid resolution.  The levels do not have to be part of the same
 * AMR patch hierarchy, however.
 *
 **************************************************************************
 */

RefineSchedule::RefineSchedule(
   tbox::Pointer<PatchLevelFillPattern> dst_level_fill_pattern,
   tbox::Pointer<hier::PatchLevel> dst_level,
   tbox::Pointer<hier::PatchLevel> src_level,
   const tbox::Pointer<xfer::RefineClasses> refine_classes,
   tbox::Pointer<xfer::RefineTransactionFactory> transaction_factory,
   xfer::RefinePatchStrategy* patch_strategy,
   bool use_time_refinement):
   d_max_stencil_width(dst_level->getDim()),
   d_max_scratch_gcw(dst_level->getDim()),
   d_boundary_fill_ghost_width(dst_level->getDim()),
   d_periodic_shift(dst_level->getDim()),
   d_unfilled_mapped_box_level(),
   d_unfilled_encon_box_level(),
   d_src_masks(dst_level->getDim()),
   d_dst_level_fill_pattern(dst_level_fill_pattern),
   d_constructing_internal_schedule(false)
{
   TBOX_ASSERT(!dst_level.isNull());
   TBOX_ASSERT(!src_level.isNull());
   TBOX_ASSERT(!refine_classes.isNull());
#ifdef DEBUG_CHECK_DIM_ASSERTIONS
   TBOX_DIM_ASSERT_CHECK_ARGS2(*dst_level, *src_level);
   if (patch_strategy) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*dst_level, *patch_strategy);
   }
#endif

   const tbox::Dimension& dim(dst_level->getDim());

   d_transaction_factory = transaction_factory;

   /*
    * Initial values; some will change in setup operations.
    */

   d_dst_level = dst_level;
   d_src_level = src_level;

   d_refine_patch_strategy = patch_strategy;

   d_number_refine_items = 0;
   d_refine_items = (const xfer::RefineClasses::Data **)NULL;

   setRefineItems(refine_classes);
   initialCheckRefineClassItems();

   d_force_boundary_fill = false;

   d_domain_is_one_box.resizeArray(
      dst_level->getGridGeometry()->getNumberBlocks(), false);
   d_num_periodic_directions = 0;

   d_coarse_priority_level_schedule = new tbox::Schedule();
   d_fine_priority_level_schedule = new tbox::Schedule();
   d_coarse_priority_level_schedule->setTimerPrefix("xfer::RefineSchedule");
   d_fine_priority_level_schedule->setTimerPrefix("xfer::RefineSchedule");

   d_supp_schedule.setNull();
   d_supp_level.setNull();

   d_max_fill_boxes = 0;

   /*
    * Initialize destination level, ghost cell widths,
    * and domain information data members.
    */

   bool recursive_schedule = false;
   initializeDomainAndGhostInformation(recursive_schedule);

   hier::IntVector min_connector_width = d_max_scratch_gcw;
   min_connector_width.max(d_boundary_fill_ghost_width);

   const Connector &dst_to_src =
      dst_level->getMappedBoxLevel()->getPersistentOverlapConnectors().findConnector(
         *src_level->getMappedBoxLevel(),
         min_connector_width);

   const Connector &src_to_dst =
      src_level->getMappedBoxLevel()->getPersistentOverlapConnectors().findConnector(
         *dst_level->getMappedBoxLevel(),
         Connector::convertHeadWidthToBase(src_level->getMappedBoxLevel()->
            getRefinementRatio(),
            dst_level->getMappedBoxLevel()->getRefinementRatio(),
            min_connector_width));

   TBOX_ASSERT(dst_to_src.getBase() == *dst_level->getMappedBoxLevel());
   TBOX_ASSERT(src_to_dst.getHead() == *dst_level->getMappedBoxLevel());
   TBOX_ASSERT(dst_to_src.getConnectorWidth() >= d_max_scratch_gcw);
   TBOX_ASSERT(dst_to_src.getConnectorWidth() >= d_boundary_fill_ghost_width);

   if ( s_extra_debug ) {
      /*
       * This check may be redundant because
       * PersistentOverlapConnectors should already guarantee
       * completeness.
       */
      hier::OverlapConnectorAlgorithm oca;
      oca.assertOverlapCorrectness(dst_to_src);
      oca.assertOverlapCorrectness(src_to_dst);
   }

   /*
    * Create fill_mapped_box_level, representing all parts of the
    * destination level, including ghost regions if desired, that this
    * schedule will attempt to fill. 
    */

   MappedBoxLevel fill_mapped_box_level(dim);
   Connector dst_to_fill;
   FillSet dst_to_fill_on_src_proc;
   setDefaultFillMappedBoxLevel(
      fill_mapped_box_level,
      dst_to_fill,
      dst_to_fill_on_src_proc,
      *dst_level->getMappedBoxLevel(),
      &dst_to_src,
      &src_to_dst,
      d_boundary_fill_ghost_width);


   /*
    * Generation the communication transactions that will move data from
    * the source to the destination.  generateCommunicationSchedule will
    * initialize the "unused" objects with information about the parts of
    * fill_mapped_box_level that cannot be filled from the source level.
    * They are unused because this RefineSchedule constructor creates
    * schedules that do not do anything to fill the parts of the
    * destination that can't be filled directly from the source.
    */ 
   tbox::Pointer<MappedBoxLevel> unused_unfilled_mapped_box_level;
   tbox::Pointer<Connector> unused_dst_to_unfilled;
   tbox::Pointer<MappedBoxLevel> unused_unfilled_encon_box_level;
   tbox::Pointer<Connector> unused_encon_to_unfilled_encon;
   bool create_transactions = true;
   generateCommunicationSchedule(
      unused_unfilled_mapped_box_level,
      unused_dst_to_unfilled,
      unused_unfilled_encon_box_level,
      unused_encon_to_unfilled_encon,
      dst_to_src,
      src_to_dst,
      dst_to_fill,
      dst_to_fill_on_src_proc,
      use_time_refinement,
      create_transactions);

   if (!d_supp_level.isNull()) {
      computeRefineOverlaps(d_refine_overlaps,
                            d_dst_level,
                            d_supp_level,
                            d_supp_to_dst,
                            d_supp_to_unfilled);
   }

   if (!d_supp_encon_level.isNull()) {
      computeRefineOverlaps(d_encon_refine_overlaps,
                            d_encon_level,
                            d_supp_encon_level,
                            d_supp_encon_to_encon,
                            d_supp_encon_to_unfilled_encon);
   }

}

/*
 **************************************************************************
 *
 * Create a refine schedule that copies data from the source level into
 * the destination level on the components represented by the refine
 * classes.  If portions of the destination level remain unfilled, then
 * the algorithm recursively fills those unfilled portions from coarser
 * levels in the AMR hierarchy.  It is assumed that the index spaces of
 * the source and destination levels are "consistent"; i.e., they
 * represent the same grid resolution.  Also, the next coarser level
 * integer argument must be the number of level in the specified
 * hierarchy representing the next coarser level of mesh resolution to
 * the destination level.
 *
 * IMPORTANT NOTES: The source level may be NULL, in which case the
 * destination level will be filled only using data interpolated from
 * coarser levels in the AMR hierarchy.  The hierarchy may be NULL only
 * if the next coarser level is -1 (that is, there is no coarser level).
 *
 **************************************************************************
 */

RefineSchedule::RefineSchedule(
   tbox::Pointer<PatchLevelFillPattern> dst_level_fill_pattern,
   tbox::Pointer<hier::PatchLevel> dst_level,
   tbox::Pointer<hier::PatchLevel> src_level,
   int next_coarser_ln,
   tbox::Pointer<hier::PatchHierarchy> hierarchy,
   const tbox::Pointer<xfer::RefineClasses> refine_classes,
   tbox::Pointer<xfer::RefineTransactionFactory> transaction_factory,
   xfer::RefinePatchStrategy* patch_strategy,
   bool use_time_refinement):
   d_max_stencil_width(dst_level->getDim()),
   d_max_scratch_gcw(dst_level->getDim()),
   d_boundary_fill_ghost_width(dst_level->getDim()),
   d_periodic_shift(dst_level->getDim()),
   d_unfilled_mapped_box_level(),
   d_unfilled_encon_box_level(),
   d_src_masks(dst_level->getDim()),
   d_dst_level_fill_pattern(dst_level_fill_pattern),
   d_constructing_internal_schedule(false)
{
   TBOX_ASSERT(!dst_level.isNull());
   TBOX_ASSERT((next_coarser_ln == -1) || !hierarchy.isNull());
   TBOX_ASSERT(!refine_classes.isNull());
#ifdef DEBUG_CHECK_DIM_ASSERTIONS
   if (!src_level.isNull()) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*dst_level, *src_level);
   }
   if (!hierarchy.isNull()) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*dst_level, *hierarchy);
   }
   if (patch_strategy) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*dst_level, *patch_strategy);
   }
#endif

   const tbox::Dimension& dim(dst_level->getDim());

   d_transaction_factory = transaction_factory;

   /*
    * Initial values; some will change in setup operations.
    */

   d_dst_level = dst_level;
   d_src_level = src_level;

   d_refine_patch_strategy = patch_strategy;

   d_number_refine_items = 0;
   d_refine_items = (const xfer::RefineClasses::Data **)NULL;

   setRefineItems(refine_classes);
   initialCheckRefineClassItems();

   d_force_boundary_fill = false;

   d_domain_is_one_box.resizeArray(
      dst_level->getGridGeometry()->getNumberBlocks(), false);
   d_num_periodic_directions = 0;

   d_coarse_priority_level_schedule.setNull();
   d_fine_priority_level_schedule.setNull();

   d_supp_schedule.setNull();
   d_supp_level.setNull();

   d_max_fill_boxes = 0;

   /*
    * Initialize destination level, ghost cell widths,
    * and domain information data members.
    */

   bool recursive_schedule = false;
   initializeDomainAndGhostInformation(recursive_schedule);

   const Connector dummy_connector;

   const Connector* dst_to_src = &dummy_connector;
   const Connector* src_to_dst = &dummy_connector;

   if (!src_level.isNull()) {
      hier::IntVector min_connector_width = d_max_scratch_gcw;
      min_connector_width.max(d_boundary_fill_ghost_width);

      dst_to_src =
         &dst_level->getMappedBoxLevel()->getPersistentOverlapConnectors().findConnector(
            *src_level->getMappedBoxLevel(),
            min_connector_width);

      src_to_dst =
         &src_level->getMappedBoxLevel()->getPersistentOverlapConnectors().findConnector(
            *dst_level->getMappedBoxLevel(),
            Connector::convertHeadWidthToBase(src_level->getMappedBoxLevel()->
                                              getRefinementRatio(),
                                              dst_level->getMappedBoxLevel()->getRefinementRatio(),
                                              min_connector_width));

      TBOX_ASSERT(dst_to_src->getBase() == *dst_level->getMappedBoxLevel());
      TBOX_ASSERT(src_to_dst->getHead() == *dst_level->getMappedBoxLevel());
      TBOX_ASSERT(dst_to_src->getConnectorWidth() >= d_max_scratch_gcw);
      TBOX_ASSERT(dst_to_src->getConnectorWidth() >= d_boundary_fill_ghost_width);
   }

   /*
    * Create fill_mapped_box_level, representing all parts of the
    * destination level, including ghost regions if desired, that this
    * schedule will fill.
    */

   MappedBoxLevel fill_mapped_box_level(dim);
   Connector dst_to_fill;
   FillSet dst_to_fill_on_src_proc;

   setDefaultFillMappedBoxLevel(
      fill_mapped_box_level,
      dst_to_fill,
      dst_to_fill_on_src_proc,
      *dst_level->getMappedBoxLevel(),
      dst_to_src,
      src_to_dst,
      d_boundary_fill_ghost_width);

   const bool skip_first_generate_schedule =
      !d_dst_level_fill_pattern->doesSourceLevelCommunicateToDestination();

   const hier::IntVector dummy_intvector(dim, -1);
   const bool dst_is_supplemental_level = false;

   /*
    * finishScheduleConstruction sets up all transactions to communicate
    * data from source to destination, and sets up recursive schedules to
    * fill whatever cannot be filled by the source.
    */
   finishScheduleConstruction(
      next_coarser_ln,
      hierarchy,
      *dst_to_src,
      *src_to_dst,
      dst_is_supplemental_level,
      dummy_intvector,
      fill_mapped_box_level,
      dst_to_fill,
      dst_to_fill_on_src_proc,
      use_time_refinement,
      skip_first_generate_schedule);

   /*
    * Compute the BoxOverlap objects that will be used to refine the
    * data from coarser levels onto the destination.
    */
   if (!d_supp_schedule.isNull()) {
      computeRefineOverlaps(d_refine_overlaps,
                            d_dst_level,
                            d_supp_level,
                            d_supp_to_dst,
                            d_supp_to_unfilled);
   }

   if (!d_supp_encon_schedule.isNull()) {
      computeRefineOverlaps(d_encon_refine_overlaps,             
                            d_encon_level,
                            d_supp_encon_level,
                            d_supp_encon_to_encon,
                            d_supp_encon_to_unfilled_encon);
   }

}

/*
 **************************************************************************
 *
 * This private constructor is used to create internal schedules that
 * fill internal levels that are used as coarse levels in refinement
 * operations. 
 *
 **************************************************************************
 */

RefineSchedule::RefineSchedule(
   tbox::Pointer<hier::PatchLevel> dst_level,
   tbox::Pointer<hier::PatchLevel> src_level,
   int next_coarser_ln,
   tbox::Pointer<hier::PatchHierarchy> hierarchy,
   const hier::IntVector& src_growth_to_nest_dst,
   const hier::Connector &dst_to_src,
   const hier::Connector &src_to_dst,
   const tbox::Pointer<xfer::RefineClasses> refine_classes,
   tbox::Pointer<xfer::RefineTransactionFactory> transaction_factory,
   xfer::RefinePatchStrategy* patch_strategy):
   d_max_stencil_width(dst_level->getDim()),
   d_max_scratch_gcw(dst_level->getDim()),
   d_boundary_fill_ghost_width(dst_level->getDim()),
   d_periodic_shift(dst_level->getDim()),
   d_unfilled_mapped_box_level(),
   d_unfilled_encon_box_level(),
   d_src_masks(dst_level->getDim()),
   d_constructing_internal_schedule(true)
{
   TBOX_ASSERT(!dst_level.isNull());
   TBOX_ASSERT((next_coarser_ln == -1) || !hierarchy.isNull());
   TBOX_ASSERT(!refine_classes.isNull());
#ifdef DEBUG_CHECK_DIM_ASSERTIONS
   if (!src_level.isNull()) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*dst_level, *src_level);
   }
   if (!hierarchy.isNull()) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*dst_level, *hierarchy);
   }
   if (patch_strategy) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*dst_level, *patch_strategy);
   }
#endif

   const tbox::Dimension& dim(dst_level->getDim());

   d_transaction_factory = transaction_factory;

   /*
    * Initial values; some will change in setup operations.
    * Note that we do not check refine items here, since this
    * constructor is private and called recursively (i.e., the
    * items have been checked already).
    */

   d_dst_level = dst_level;
   d_src_level = src_level;

   d_refine_patch_strategy = patch_strategy;

   d_number_refine_items = 0;
   d_refine_items = (const xfer::RefineClasses::Data **)NULL;

   setRefineItems(refine_classes);

   d_force_boundary_fill = false;

   d_domain_is_one_box.resizeArray(
      dst_level->getGridGeometry()->getNumberBlocks(), false);
   d_num_periodic_directions = 0;

   d_coarse_priority_level_schedule.setNull();
   d_fine_priority_level_schedule.setNull();

   d_supp_schedule.setNull();
   d_supp_level.setNull();

   d_max_fill_boxes = 0;

   d_dst_level_fill_pattern = new PatchLevelFullFillPattern();

   /*
    * Initialize destination level, ghost cell widths,
    * and domain information data members.
    */

   bool recursive_schedule = true;
   initializeDomainAndGhostInformation(recursive_schedule);

   /*
    * Note that we cannot assert dst<==>src are complete, because they
    * are supp<==>hiercoarse from the recursion.
    * finishScheduleConstruction should ensure that supp<==>hiercoarse
    * are complete enough to fill dst and connect the unfilled portion
    * to the next hierarchy coarser level.  finishScheduleConstruction
    * does not (cannot) guarantee that supp<==>hiercoarse is complete.
    */
   TBOX_ASSERT( !src_level.isNull() );
   TBOX_ASSERT( dst_to_src.isInitialized() );
   TBOX_ASSERT( src_to_dst.isInitialized() );
   TBOX_ASSERT(dst_to_src.getConnectorWidth() >= d_max_stencil_width);
   TBOX_ASSERT(dst_to_src.getBase() == *dst_level->getMappedBoxLevel());
   TBOX_ASSERT(src_to_dst.getHead() == *dst_level->getMappedBoxLevel());

   /*
    * Create fill_mapped_box_level, representing all parts of the
    * destination level, including ghost regions if desired, that this
    * schedule will fill.  Here, the destination is always a supplemental
    * level constructed by coarsening another RefineSchedule's unfilled
    * boxes.  As the destination will be used as a coarse level in a
    * refinement operation, the fill_mapped_box_level will be the boxes
    * of the destination level grown by the maximum interplation stencil
    * width.
    */

   MappedBoxLevel fill_mapped_box_level(dim);
   Connector dst_to_fill;
   FillSet dst_to_fill_on_src_proc;
   setDefaultFillMappedBoxLevel(
      fill_mapped_box_level,
      dst_to_fill,
      dst_to_fill_on_src_proc,
      *dst_level->getMappedBoxLevel(),
      &dst_to_src,
      &src_to_dst,
      d_max_stencil_width);

   bool use_time_refinement = true;
   const bool dst_is_supplemental_level = true;

   /*
    * finishSchedule construction sets up all transactions to communicate
    * data from source to destination, and sets up recursive schedules to
    * fill whatever cannot be filled by the source.
    */

   finishScheduleConstruction(
      next_coarser_ln,
      hierarchy,
      dst_to_src,
      src_to_dst,
      dst_is_supplemental_level,
      src_growth_to_nest_dst,
      fill_mapped_box_level,
      dst_to_fill,
      dst_to_fill_on_src_proc,
      use_time_refinement);

   /*
    * Compute the BoxOverlap objects that will be used to refine the
    * data from coarser levels onto the destination.
    */

   if (!d_supp_schedule.isNull()) {
      computeRefineOverlaps(d_refine_overlaps,
                            d_dst_level,
                            d_supp_level,
                            d_supp_to_dst,
                            d_supp_to_unfilled);
   }

   if (!d_supp_encon_schedule.isNull()) {
      computeRefineOverlaps(d_encon_refine_overlaps,             
                            d_encon_level,
                            d_supp_encon_level,
                            d_supp_encon_to_encon,
                            d_supp_encon_to_unfilled_encon);
   }

}

/*
 **************************************************************************
 *
 * The destructor for the refine schedule class implicitly deallocates
 * all of the data associated with the communication schedule.
 *
 **************************************************************************
 */

RefineSchedule::~RefineSchedule()
{
   clearRefineItems();

   d_coarse_priority_level_schedule.setNull();
   d_fine_priority_level_schedule.setNull();

   d_supp_schedule.setNull();
   d_supp_level.setNull();
}

/*
 *************************************************************************
 *
 * Reset schedule with new set of refine items.
 *
 ************************************************************************
 */

void RefineSchedule::reset(
   const tbox::Pointer<xfer::RefineClasses> refine_classes)
{
   TBOX_ASSERT(!refine_classes.isNull());

   setRefineItems(refine_classes);
   if (!d_supp_schedule.isNull()) {
      d_supp_schedule->reset(refine_classes);
   }
   if (!d_supp_encon_schedule.isNull()) {
      d_supp_encon_schedule->reset(refine_classes);
   }
}

/*
 **************************************************************************
 *
 * Return const pointer to equivalence classes used in schedule.
 *
 **************************************************************************
 */

const tbox::Pointer<xfer::RefineClasses>&
RefineSchedule::getEquivalenceClasses() const
{
   return d_refine_classes;
}

const hier::IntVector&
RefineSchedule::getBoundaryFillGhostWidth() const
{
   return d_boundary_fill_ghost_width;
}

/*
 ************************************************************************
 * Construct transactions for schedule and set up recursive schedules if
 * needed.
 *
 * Generate communication schedules to transfer data from src to
 * fillboxes associated with dst boxes.  What parts cannot be filled
 * becomes the "unfilled" boxes.  If no source, all fill boxes become
 * "unfilled" boxes.  We also construct unfilled boxes at enhanced
 * connectivity block boundaries.
 *
 * If there are any unfilled boxes, we coarsen them to create a
 * supplemental level and set up a recursive schedule for filling the
 * supplemental level.  The idea is to interpolate data from the
 * supplemental level to fill the unfilled boxes.
 ************************************************************************
 */

void RefineSchedule::finishScheduleConstruction(
   int next_coarser_ln,
   tbox::Pointer<hier::PatchHierarchy> hierarchy,
   const Connector& dst_to_src,
   const Connector& src_to_dst,
   const bool dst_is_supplemental_level,
   const hier::IntVector& src_growth_to_nest_dst,
   const MappedBoxLevel& fill_mapped_box_level,
   const Connector& dst_to_fill,
   const FillSet& dst_to_fill_on_src_proc,
   bool use_time_interpolation,
   bool skip_generate_schedule)
{
   t_finish_sched_const->start();
   TBOX_ASSERT((next_coarser_ln == -1) || !hierarchy.isNull());

   static int recursion_level = -1;
   ++recursion_level;
   if (s_extra_debug) {
      tbox::plog << "finishScheduleConstruction entered recursion_level="
                 << recursion_level << " next_coarser_ln=" << next_coarser_ln
                 << std::endl;
   }

   const tbox::Dimension& dim(hierarchy->getDim());

   const bool fully_periodic = d_num_periodic_directions == dim.getValue();

   hier::MappedBoxLevelConnectorUtils edge_utils;
   hier::OverlapConnectorAlgorithm oca;

   const MappedBoxLevel& dst_mapped_box_level = dst_to_fill.getBase();
   if (!d_src_level.isNull()) {
      // Should never have a source without connection from destination.
      TBOX_ASSERT(dst_to_src.isInitialized());
   }

   if (s_extra_debug) {
      tbox::plog << "finishScheduleConstruction in recursion_level="
                 << recursion_level << " next_coarser_ln=" << next_coarser_ln
                 << " before computing unfilled boxes."
                 << "\ndst_mapped_box_level:\n"
                 << dst_mapped_box_level.format("D->", 2)
                 << "\nfill_mapped_box_level:\n"
                 << fill_mapped_box_level.format("F->", 2)
                 << "\ndst_to_fill:\n"
                 << dst_to_fill.format("DF->", 2)
                 << std::endl;
   }

   const int nblocks = d_dst_level->getGridGeometry()->getNumberBlocks();

   d_coarse_priority_level_schedule = new tbox::Schedule();
   d_fine_priority_level_schedule = new tbox::Schedule();



   /*
    * Generate the schedule for filling the boxes in dst_to_fill.
    * Any portions of the fill boxes that cannot be filled from
    * the source is placed in d_unfilled_mapped_box_level.
    *
    * If the source is not given or skip_generate_schedule==true,
    * the schedule generation degenates to turning all the fill boxes
    * into unfilled boxes.
    */

   tbox::Pointer<Connector> dst_to_unfilled;
   const tbox::Pointer<hier::GridGeometry> &grid_geometry(
      d_dst_level->getGridGeometry());
   tbox::Pointer<Connector> encon_to_unfilled_encon;

   bool create_transactions = true;
   if (d_src_level.isNull() || skip_generate_schedule) {
      create_transactions = false;
   }

   generateCommunicationSchedule(
      d_unfilled_mapped_box_level,
      dst_to_unfilled,
      d_unfilled_encon_box_level,
      encon_to_unfilled_encon,
      dst_to_src,
      src_to_dst,
      dst_to_fill,
      dst_to_fill_on_src_proc,
      use_time_interpolation,
      create_transactions);

   if (s_extra_debug) {
      tbox::plog << "finishScheduleConstruction in recursion_level="
                 << recursion_level << " next_coarser_ln=" << next_coarser_ln
                 << " after computing unfilled boxes."
                 << "\nd_unfilled_mapped_box_level:\n"
                 << d_unfilled_mapped_box_level->format("UF->", 2)
                 << "\ndst_to_unfilled:\n"
                 << dst_to_unfilled->format("DUF->", 2)
                 << std::endl;
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   for ( int bn=0; bn<nblocks; ++bn ) {
      TBOX_ASSERT(fill_mapped_box_level.getLocalBoundingBox(bn).contains(
                     d_unfilled_mapped_box_level->getLocalBoundingBox(bn)));
   }
#endif



   /*
    * d_unfilled_mapped_box_level may now include boxes outside of the
    * physical domain.  These parts must be removed if they are not at
    * periodic boundaries.  They will be filled through a user
    * call-back method.
    */

   if (!fully_periodic) {
      finishScheduleConstruction_shearUnfilledBoxesOutsideNonperiodicBoundaries(
         *d_unfilled_mapped_box_level,
         *dst_to_unfilled,
         hierarchy );
   }


   t_get_global_mapped_box_count->barrierAndStart();

   const bool need_to_fill =
      (d_unfilled_mapped_box_level->getGlobalNumberOfBoxes() > 0);

   const bool need_to_fill_encon =
      grid_geometry->hasEnhancedConnectivity() &&
      (d_unfilled_encon_box_level->getGlobalNumberOfBoxes() > 0);

   t_get_global_mapped_box_count->stop();



   /*
    * If there remain boxes to be filled from coarser levels, then set
    * up data for recursive schedule generation:
    *
    * 1. Generate a supplemental MappedBoxLevel
    * (supp_mapped_box_level) by coarsening the unfilled boxes.
    *
    * 2. Connect supp_mapped_box_level to the next coarser level on
    * the hierarchy.
    *
    * 3: Construct the supplemental PatchLevel (d_supp_level) and
    * construct d_supp_schedule to fill d_supp_level.  The coarser
    * level on the hierarchy will be the source for filling
    * d_supp_level, which is why we need step 2..
    *
    * The idea is that once d_supp_level is filled, we can refine its
    * data to fill the current unfilled boxes.
    */

   if (need_to_fill) {

      t_finish_sched_const_recurse->start();

      if (s_extra_debug) {
         tbox::plog << "finishScheduleConstruction in recursion_level="
                    << recursion_level << " next_coarser_ln=" << next_coarser_ln
                    << " needs to recurse"
                    << std::endl;
      }


      /*
       * If there are no coarser levels in the hierarchy or the
       * hierarchy is null, then throw an error.  Something is messed
       * up someplace and code execution cannot proceed.
       */
      if (next_coarser_ln < 0) {
         TBOX_ERROR(
            "Internal error in RefineSchedule::finishScheduleConstruction..."
            << "\n In finishScheduleConstruction() -- "
            << "\n No coarser levels...will not fill from coarser."
            << "\n dst_mapped_box_level:\n" << dst_mapped_box_level.format("DEST->", 2)
            << "\n dst_to_fill:\n" << dst_to_fill.format("DF->", 2)
            << "\n d_unfilled_mapped_box_level:\n" << d_unfilled_mapped_box_level->format("UF->", 2)
            << "\n dst_to_unfilled:\n" << dst_to_unfilled->format("DU->", 2)
            << "\n dst_to_src:\n" << dst_to_src.format("DS->", 2)
            << std::endl);
      } else {
         if (hierarchy.isNull()) {
            TBOX_ERROR("Internal RefineSchedule error..."
               << "\n In finishScheduleConstruction() -- "
               << "\n Need to fill from coarser hierarchy level and \n"
               << "hierarchy is unavailable." << std::endl);
         }
      }


      /*
       * hiercoarse is the coarse level on the hierarchy.  It is to be
       * differentiated from the supplemental (supp) level, which is at
       * the same resolution and level number but is not on the hierarchy.
       */
      const tbox::Pointer<hier::PatchLevel> hiercoarse_level =
         hierarchy->getPatchLevel(next_coarser_ln);

      const hier::MappedBoxLevel &hiercoarse_mapped_box_level(
         *hiercoarse_level->getMappedBoxLevel());


      /*
       * Ratio to the next coarser level in the hierarchy.
       */
      const hier::IntVector dst_hiercoarse_ratio(
         d_dst_level->getRatioToLevelZero()
         / hiercoarse_level->getRatioToLevelZero());


      /*
       * Set up the supplemental MappedBoxLevel and also set up
       * d_dst_to_supp, d_supp_to_dst and d_supp_to_unfilled.  These
       * Connectors are easily generated using dst_to_unfilled.
       */

      hier::MappedBoxLevel supp_mapped_box_level(dim);
      finishScheduleConstruction_setupSupplementalMappedBoxLevel(
         supp_mapped_box_level,
         hiercoarse_mapped_box_level,
         *dst_to_unfilled);

      /*
       * Connect the supplemental MappedBoxLevel (the next recursion's
       * dst) to the hiercoarse MappedBoxLevel (the next recursion's
       * src).
       */

      Connector supp_to_hiercoarse;
      Connector hiercoarse_to_supp;

      finishScheduleConstruction_connectSuppToHiercoarse(
         supp_to_hiercoarse,
         hiercoarse_to_supp,
         supp_mapped_box_level,
         hierarchy,
         next_coarser_ln,
         dst_to_src,
         src_to_dst,
         dst_is_supplemental_level);



      if (d_num_periodic_directions > 0) {
         /*
          * Add periodic images to supp_mapped_box_level and
          * corresponding edges to hiercoarse_to_supp and
          * supp_to_hiercoarse.  (supp_mapped_box_level was built
          * from unfilled_mapped_box_level, which did not have
          * periodic images.  Since supp_mapped_box_level is the
          * next dst MappedBoxLevel, the next recursion may have to
          * fill it using periodic neighbors, supp_mapped_box_level
          * needs periodic images.)
          */
         TBOX_ASSERT(
            supp_mapped_box_level.getLocalNumberOfBoxes() ==
            supp_mapped_box_level.getMappedBoxes().size());
         const hier::Connector& hiercoarse_to_hiercoarse =
            hiercoarse_level->getMappedBoxLevel()->
            getPersistentOverlapConnectors().
            findConnector(*hiercoarse_level->getMappedBoxLevel(),
                          hiercoarse_to_supp.getConnectorWidth());
         edge_utils.addPeriodicImagesAndRelationships(
            supp_mapped_box_level,
            supp_to_hiercoarse,
            hiercoarse_to_supp,
            hierarchy->getDomainSearchTree(hier::BlockId(0)),
            hiercoarse_to_hiercoarse);

         if (s_extra_debug) {
            sanityCheckSupplementalAndHiercoarseLevels(
               supp_to_hiercoarse,
               hiercoarse_to_supp,
               hierarchy,
               next_coarser_ln);
         }

      }



      /*
       * Construct the supplemental PatchLevel and reset
       * supp<==>hiercoarse connectors to use the PatchLevel's
       * MappedBoxLevel.  Note that supp<==>hiercoarse is not
       * guaranteed to be complete, so we cannot put them into
       * PersistentOverlapConnectors.
       */

      t_make_supp_level->start();
      d_supp_level = new hier::PatchLevel(
            supp_mapped_box_level,
            hiercoarse_level->getGridGeometry(),
            hiercoarse_level->getPatchDescriptor());
      t_make_supp_level->stop();
      d_supp_level->setLevelNumber(next_coarser_ln);
      d_supp_level->setNextCoarserHierarchyLevelNumber(next_coarser_ln - 1);;

      if (hiercoarse_level->getGridGeometry()->getNumberBlocks() > 1) {
         hiercoarse_level->getGridGeometry()->
            adjustMultiblockPatchLevelBoundaries(*d_supp_level);
      }

      supp_to_hiercoarse.initialize( *d_supp_level->getMappedBoxLevel(),
                                     supp_to_hiercoarse.getHead(),
                                     supp_to_hiercoarse.getConnectorWidth(),
                                     supp_to_hiercoarse.getNeighborhoodSets() );
      hiercoarse_to_supp.initialize( hiercoarse_to_supp.getBase(),
                                     *d_supp_level->getMappedBoxLevel(),
                                     hiercoarse_to_supp.getConnectorWidth(),
                                     hiercoarse_to_supp.getNeighborhoodSets()  );


      /*
       * Compute how much hiercoarse has to grow to nest supp, a
       * required parameter in the private constructor.
       *
       * If dst is a supplemental level (generated by RefineSchedule),
       * we have the info to compute the growth.  If not, we make some
       * assumptions about where dst came from in order to determine
       * how its fill boxes nest in hiercoarse.
       */
      hier::IntVector hiercoarse_growth_to_nest_supp(dim);
      if (dst_is_supplemental_level) {
         /*
          * Assume that src barely nests in hiercoarse.  (In most
          * places, it nests by a margin equal to the nesting buffer,
          * but we don't count on that because the nesting buffer is
          * relaxed at physical boundaries.)  To nest dst, hiercoarse
          * has to grow as much as the src does, plus the ghost width
          * of the fill.
          *
          * FIXME: We may in fact be able to count on the nesting
          * buffer because extending boxes to physical boundaries do
          * not create any extra relationships.  However, we don't
          * currently have access to the size of the nesting buffer.
          */
         hiercoarse_growth_to_nest_supp =
            src_growth_to_nest_dst + dst_to_fill.getConnectorWidth();
      } else {
         /*
          * dst may be:
          * 1. The hierarchy level just finer than level number next_coarser_ln.
          * 2. A level that nests in level number next_coarser_ln:
          *    a. A new level generated by GriddingAlgorithm.
          *    b. The hierarchy level just finer than level number next_coarser_ln,
          *       coarsened for Richardson extrapolation.
          * In any case, dst should nest in hiercoarse.  Furthermore, it does
          * not grow when coarsened into the hiercoarse index space.
          * To nest dst and its fill boxes, hiercoarse just has to grow by
          * the ghost width of the fill.
          */
         hiercoarse_growth_to_nest_supp = dst_to_fill.getConnectorWidth();
      }
      hiercoarse_growth_to_nest_supp.ceiling(dst_hiercoarse_ratio);


      t_finish_sched_const_recurse->stop();


      /*
       * We now have all the data for building the supplemental
       * schedule using the private constructor.
       *
       * We need to make sure that the coarse schedule uses
       * BoxGeometryVariableFillPattern, so that it fills all needed
       * parts of d_supp_level
       */
      tbox::Pointer<BoxGeometryVariableFillPattern> bg_fill_pattern(
         new BoxGeometryVariableFillPattern());

      tbox::Pointer<RefineClasses> coarse_schedule_refine_classes(
         new RefineClasses());

      const int num_refine_items =
         d_refine_classes->getNumberOfRefineItems();

      for (int nd = 0; nd < num_refine_items; nd++) {
         RefineClasses::Data item = d_refine_classes->getRefineItem(nd);
         item.d_var_fill_pattern = bg_fill_pattern;
         coarse_schedule_refine_classes->insertEquivalenceClassItem(item);
      }

      t_finish_sched_const->stop();

      if (s_extra_debug) {
         const std::string dbgbord;
         tbox::plog << "finishScheduleConstruction in recursion_level="
                    << recursion_level << " next_coarser_ln=" << next_coarser_ln
                    << " creating coarse schedule to fill supp\n"
                    << "supp_mapped_box_level:\n"
                    << supp_mapped_box_level.format(dbgbord, 2);
      }
      d_supp_schedule = new RefineSchedule(d_supp_level,
            hiercoarse_level,
            next_coarser_ln - 1,
            hierarchy,
            hiercoarse_growth_to_nest_supp,
            supp_to_hiercoarse,
            hiercoarse_to_supp,
            coarse_schedule_refine_classes,
            d_transaction_factory,
            d_refine_patch_strategy);

   } else {
      t_finish_sched_const->stop();
   }

   if (need_to_fill_encon) {

      /*
       * Create schedule to fill unfilled boxes at enhanced connectivity
       */

      const tbox::Pointer<hier::PatchLevel> hiercoarse_level =
         hierarchy->getPatchLevel(next_coarser_ln);

      createEnconFillSchedule(
         hierarchy,
         hiercoarse_level,
         dst_is_supplemental_level,
         src_growth_to_nest_dst,
         *encon_to_unfilled_encon);
   }

   if (s_extra_debug) {
      tbox::plog << "finishScheduleConstruction exiting recursion_level="
                 << recursion_level << " next_coarser_ln=" << next_coarser_ln
                 << std::endl;
      --recursion_level;
   }
}

/*
 ***********************************************************************
 * Create schedule for filling unfilled boxes at enhanced connectivity
 *
 * d_supp_encon_level is created by coarsening d_unfilled_encon_level.
 * d_supp_encon_schedule is created to communicate data from the
 * hierarchy to fill d_supp_encon_level.
 ***********************************************************************
 */
void
RefineSchedule::createEnconFillSchedule(
   const tbox::Pointer<hier::PatchHierarchy>& hierarchy,
   const tbox::Pointer<hier::PatchLevel>& hiercoarse_level,
   const bool dst_is_supplemental_level,
   const hier::IntVector& src_growth_to_nest_dst,
   const hier::Connector& encon_to_unfilled_encon)
{
   TBOX_ASSERT(!hiercoarse_level.isNull());

   const int next_coarser_ln = hiercoarse_level->getLevelNumber();

   hier::MappedBoxLevel supp_encon_box_level(
      hiercoarse_level->getRatioToLevelZero(),
      hiercoarse_level->getGridGeometry(),
      d_unfilled_mapped_box_level->getMPI());

   hier::NeighborhoodSet encon_to_supp_encon_nbrhood_set;
   hier::NeighborhoodSet supp_encon_to_unfilled_nbrhood_set;

   const hier::IntVector dst_hiercoarse_ratio(
      d_dst_level->getRatioToLevelZero()
      / hiercoarse_level->getRatioToLevelZero());

   /*
    * Loop over encon_to_unfilled_encon_nbrhood_set, which mappes the boxes
    * of d_encon_level to the unfilled encon boxes.
    */
   const hier::NeighborhoodSet& encon_to_unfilled_encon_nbrhood_set =
      encon_to_unfilled_encon.getNeighborhoodSets();

   for (hier::NeighborhoodSet::const_iterator ei =
        encon_to_unfilled_encon_nbrhood_set.begin();
        ei != encon_to_unfilled_encon_nbrhood_set.end(); ++ei) {

      const hier::BoxId& encon_mapped_box_mbid = ei->first;
      const NeighborSet& encon_unfilled_parts = ei->second;

      /*
       * For each unfilled box, coarsen and add to supp_encon_box_level.
       */
      for (hier::MappedBoxSet::const_iterator ni =
         encon_unfilled_parts.begin();
         ni != encon_unfilled_parts.end(); ++ni) {

         const hier::Box& unfilled_mapped_box = *ni;
         hier::Box supp_box(unfilled_mapped_box);
         supp_box.coarsen(dst_hiercoarse_ratio);

         const hier::Box& supp_mapped_box =
            *supp_encon_box_level.addBox(supp_box, (*ni).getBlockId());

         /*
          * Set up neighbor relationships for supp_encon_box_level
          */
         NeighborSet& encon_to_supp_encon_nabrs =
            encon_to_supp_encon_nbrhood_set[encon_mapped_box_mbid];
         encon_to_supp_encon_nabrs.insert(supp_mapped_box);

         const hier::BoxId& supp_mapped_box_id =
            supp_mapped_box.getId();
         NeighborSet& supp_encon_to_unfilled_nabrs =
            supp_encon_to_unfilled_nbrhood_set[supp_mapped_box_id];
         supp_encon_to_unfilled_nabrs.insert(unfilled_mapped_box);
      }
   }

   const tbox::Dimension& dim = hiercoarse_level->getDim();

   /*
    * Initialize Connectors
    */
   d_encon_to_supp_encon.swapInitialize(
      *(d_encon_level->getMappedBoxLevel()),
      supp_encon_box_level,
      hier::IntVector::getZero(dim),
      encon_to_supp_encon_nbrhood_set);

   d_supp_encon_to_unfilled_encon.initialize(
      supp_encon_box_level,
      *d_unfilled_encon_box_level,
      hier::IntVector::getZero(dim),
      supp_encon_to_unfilled_nbrhood_set);

   if (d_encon_to_supp_encon.isInitialized()) {
      d_supp_encon_to_encon.initializeToLocalTranspose(
         d_encon_to_supp_encon);
   }

   /*
    * hiercoarse_level is the level in the hierarchy that is one level
    * coarser than d_dst_level.  It will be uses as the source of data
    * to fill the supplemental encon level.
    */
   const hier::MappedBoxLevel &hiercoarse_mapped_box_level(
      *hiercoarse_level->getMappedBoxLevel());

   /*
    * The computation of Connectors between supp_encon and hiercoarse is 
    * done the same way as the Connectors between supp and hiercoarse
    * above.
    *
    * TODO: Merge replicated code into a single private method?
    */
   Connector supp_encon_to_hiercoarse;
   Connector hiercoarse_to_supp_encon;
   {

      hier::OverlapConnectorAlgorithm oca;

      RefineScheduleConnectorWidthRequestor rscwri;
      std::vector<hier::IntVector> self_connector_widths;
      std::vector<hier::IntVector> fine_connector_widths;
      rscwri.computeRequiredConnectorWidths(self_connector_widths,
                                            fine_connector_widths,
                                            *hierarchy);

      const Connector* encon_to_hiercoarse = NULL;
      const Connector* hiercoarse_to_encon = NULL;
      Connector bridged_encon_to_hiercoarse;
      Connector bridged_hiercoarse_to_encon;

      hier::IntVector min_encon_to_hiercoarse_width(d_max_scratch_gcw);
      min_encon_to_hiercoarse_width.max(d_max_stencil_width);
      min_encon_to_hiercoarse_width.max(d_boundary_fill_ghost_width);
      hier::IntVector min_hiercoarse_to_encon_width =
         Connector::convertHeadWidthToBase(
            hiercoarse_mapped_box_level.getRefinementRatio(),
            d_dst_level->getMappedBoxLevel()->getRefinementRatio(),
            min_encon_to_hiercoarse_width);

      const bool has_encon_to_hiercoarse =
         d_encon_level->getMappedBoxLevel()->
            getPersistentOverlapConnectors().
               hasConnector(
                  hiercoarse_mapped_box_level,
                  min_encon_to_hiercoarse_width);
      const bool has_hiercoarse_to_encon =
         hiercoarse_mapped_box_level.getPersistentOverlapConnectors().
         hasConnector(
            *d_encon_level->getMappedBoxLevel(),
            min_hiercoarse_to_encon_width);

      if (has_encon_to_hiercoarse && has_hiercoarse_to_encon) {

         encon_to_hiercoarse =
            &d_encon_level->getMappedBoxLevel()->
               getPersistentOverlapConnectors()
                  .findConnector(
                     hiercoarse_mapped_box_level,
                     min_encon_to_hiercoarse_width);

         hiercoarse_to_encon =
            &hiercoarse_mapped_box_level.getPersistentOverlapConnectors()
            .findConnector(
               *d_encon_level->getMappedBoxLevel(),
               min_hiercoarse_to_encon_width);

      } else {
         /*
          * Connectors encon<==>hiercoarse are not provided.
          * We have to bridge through src for it.
          * (This requires src<==>hiercoarse.)
          */
         const hier::MappedBoxLevel& src_mapped_box_level =
            *(d_src_level->getMappedBoxLevel());

         if (*hierarchy->getMappedBoxLevel(next_coarser_ln+1) !=
             src_mapped_box_level) {
            TBOX_ERROR("Missing encon<==>hiercoarse connector and\n"
               << "src is not from hierarchy.");
         }

         hier::IntVector hiercoarse_to_src_width(
            fine_connector_widths[next_coarser_ln]);
         hier::IntVector src_to_hiercoarse_width(
            hiercoarse_to_src_width *
            d_src_level->getRatioToCoarserLevel());

         /*
          * Using hierarchy to get required Connector width assumes that
          * the src level has same refinement ratio as
          * next_coarser_ln+1, but for Richardson extrapolation,
          * that is not the case, so we have to adjust.
          */
         if (d_src_level->getMappedBoxLevel()->getRefinementRatio() <=
             hierarchy->getMappedBoxLevel(next_coarser_ln + 1)->
                getRefinementRatio()) {

            src_to_hiercoarse_width *=
               d_src_level->getMappedBoxLevel()->getRefinementRatio();
            src_to_hiercoarse_width /= hierarchy->getMappedBoxLevel(
               next_coarser_ln + 1)->getRefinementRatio();

         } else if (d_src_level->getMappedBoxLevel()->getRefinementRatio()
                    >=
                    hierarchy->getMappedBoxLevel(next_coarser_ln
                       + 1)->getRefinementRatio()) {

            src_to_hiercoarse_width *= hierarchy->getMappedBoxLevel(
                  next_coarser_ln + 1)->getRefinementRatio();
            src_to_hiercoarse_width /=
               d_src_level->getMappedBoxLevel()->getRefinementRatio();

         }

         const hier::Connector& src_to_hiercoarse =
            d_src_level->getMappedBoxLevel()->
               getPersistentOverlapConnectors().findConnector(
                  hiercoarse_mapped_box_level,
                  src_to_hiercoarse_width);

         const hier::Connector& hiercoarse_to_src =
            hiercoarse_mapped_box_level.getPersistentOverlapConnectors()
            .findConnector(
               *d_src_level->getMappedBoxLevel(),
               hiercoarse_to_src_width);

         oca.setSanityCheckMethodPostconditions(true);
         /*
          * Don't use the strict bridge theorem here because it
          * cannot guarantee sufficient width.  We know from how
          * dst nests in hiercoarse what output Connector width
          * can guarantee that all dst Boxes are seen by a
          * hiercoarse Box.
          */
         oca.bridge(
            bridged_encon_to_hiercoarse,
            bridged_hiercoarse_to_encon,
            d_encon_to_src,
            src_to_hiercoarse,
            hiercoarse_to_src,
            d_src_to_encon,
            fine_connector_widths[next_coarser_ln]);

         encon_to_hiercoarse = &bridged_encon_to_hiercoarse;
         hiercoarse_to_encon = &bridged_hiercoarse_to_encon;

      } // End block bridging for encon<==>hiercoarse.

      /*
       * Compute supp_encon<==>hiercoarse by bridging
       * supp_encon<==>encon<==>hiercoarse.
       */

      /*
       * Don't use the strict bridge theorem here because it
       * cannot guarantee sufficient width.  We know from how
       * dst nests in hiercoarse what output Connector width
       * can guarantee that all dst Boxes are seen by a
       * hiercoarse Box.
       */
      oca.bridge(
         supp_encon_to_hiercoarse,
         hiercoarse_to_supp_encon,
         d_supp_encon_to_encon,
         *encon_to_hiercoarse,
         *hiercoarse_to_encon,
         d_encon_to_supp_encon,
         fine_connector_widths[next_coarser_ln]);

   }

   /*
    * Create d_supp_encon_level
    */
   d_supp_encon_level = new hier::PatchLevel(
      supp_encon_box_level,
      d_dst_level->getGridGeometry(),
      hiercoarse_level->getPatchDescriptor());
   d_supp_encon_level->setLevelNumber(next_coarser_ln);
   d_supp_encon_level->getGridGeometry()->
      adjustMultiblockPatchLevelBoundaries(*d_supp_encon_level);

   /*
    * Reset supp_encon<==>hiercoarse connectors to use the PatchLevel's
    * MappedBoxLevel.
    */
   supp_encon_to_hiercoarse.initialize(
      *d_supp_encon_level->getMappedBoxLevel(),
      supp_encon_to_hiercoarse.getHead(),
      supp_encon_to_hiercoarse.getConnectorWidth(),
      supp_encon_to_hiercoarse.getNeighborhoodSets());
   hiercoarse_to_supp_encon.initialize(
      hiercoarse_to_supp_encon.getBase(),
      *d_supp_encon_level->getMappedBoxLevel(),
      hiercoarse_to_supp_encon.getConnectorWidth(),
      hiercoarse_to_supp_encon.getNeighborhoodSets()  );

   /*
    * Compute this nesting value the same as for supp
    */
   hier::IntVector hiercoarse_growth_to_nest_supp_encon(dim);
   if (dst_is_supplemental_level) {
      hiercoarse_growth_to_nest_supp_encon =
         src_growth_to_nest_dst + encon_to_unfilled_encon.getConnectorWidth();
      hiercoarse_growth_to_nest_supp_encon.ceiling(dst_hiercoarse_ratio);
   } else {
      hiercoarse_growth_to_nest_supp_encon =
         encon_to_unfilled_encon.getConnectorWidth();
      hiercoarse_growth_to_nest_supp_encon.ceiling(dst_hiercoarse_ratio);
   }

   /*
    * We need to make sure that the coarse schedule uses
    * BoxGeometryVariableFillPattern, so that it fills all needed parts of
    * d_supp_encon_level
    */
   tbox::Pointer<BoxGeometryVariableFillPattern> bg_fill_pattern(
      new BoxGeometryVariableFillPattern());

   tbox::Pointer<RefineClasses> coarse_schedule_refine_classes(
      new RefineClasses());

   const int num_refine_items =
      d_refine_classes->getNumberOfRefineItems();

   for (int nd = 0; nd < num_refine_items; nd++) {
      RefineClasses::Data item = d_refine_classes->getRefineItem(nd);
      item.d_var_fill_pattern = bg_fill_pattern;
      coarse_schedule_refine_classes->insertEquivalenceClassItem(item);
   }

   /*
    * Schedule to fill d_supp_encon_level
    */
   d_supp_encon_schedule = new RefineSchedule(d_supp_encon_level,
      hiercoarse_level,
      next_coarser_ln - 1,
      hierarchy,
      hiercoarse_growth_to_nest_supp_encon,
      supp_encon_to_hiercoarse,
      hiercoarse_to_supp_encon,
      coarse_schedule_refine_classes,
      d_transaction_factory,
      d_refine_patch_strategy);

}



/*
 ******************************************************************
 * Create level for unfilled boxes at enhanced connectivity.
 ******************************************************************
 */
void
RefineSchedule::createUnfilledEnconLevelWithNoSource(
   tbox::Pointer<hier::Connector>& encon_to_unfilled_encon,
   const hier::Connector& dst_to_fill)
{
   hier::MappedBoxSet unfilled_encon_boxes;
   hier::LocalId last_unfilled_local_id(-1);

   hier::NeighborhoodSet unfilled_encon_nbrhood_set;

   const hier::NeighborhoodSet& dst_to_fill_nbrhood_set =
      dst_to_fill.getNeighborhoodSets();

   /*
    * Loop over dst_to_fill_nbrhood_set, finding the set of fill
    * mapped boxes for each local destination mapped box.
    */
   for (hier::NeighborhoodSet::const_iterator
        df_iter = dst_to_fill_nbrhood_set.begin();
        df_iter != dst_to_fill_nbrhood_set.end(); ++df_iter) {

      const hier::BoxId& dst_mapped_box_id = df_iter->first;
      const int dst_blk = dst_mapped_box_id.getBlockId().getBlockValue();

      const tbox::Pointer<hier::GridGeometry> grid_geometry
         (d_dst_level->getGridGeometry());
 
      const hier::BoxList& sing_boxes =
         grid_geometry->getSingularityBoxList(dst_blk);

      if (sing_boxes.size()) {

         const tbox::List<hier::GridGeometry::Neighbor>& neighbors =
            grid_geometry->getNeighbors(dst_blk);

         const NeighborSet& fill_nabrs = df_iter->second;

         /*
          * Loop over the fill boxes for the current destination box
          */
         for (NeighborSet::iterator fn_iter = fill_nabrs.begin();
              fn_iter != fill_nabrs.end(); ++fn_iter) {

            const hier::Box& fill_mapped_box = *fn_iter;

            const hier::NeighborhoodSet& dst_to_encon_nbrhood_set =
               d_dst_to_encon.getNeighborhoodSets();

            /*
             * Determine there is anything in dst_to_encon_nbrhood_set
             * for dst_mapped_box_id.  If not, the destination box 
             * does not touch enhanced connectivity.
             */
            hier::NeighborhoodSet::const_iterator find_encon_nabrs =
               dst_to_encon_nbrhood_set.find(dst_mapped_box_id);

            if (find_encon_nabrs != dst_to_encon_nbrhood_set.end()) {

               const NeighborSet& dst_encon_nabrs =
                  find_encon_nabrs->second;

               /*
                * Loop over mapped boxes representing the
                * ghost regions of the destination at the enhanced
                * connectivity boundary.
                */
               for (NeighborSet::const_iterator
                    de_iter = dst_encon_nabrs.begin();
                    de_iter != dst_encon_nabrs.end(); ++de_iter) {

                  const hier::BoxId& encon_mb_id =
                     de_iter->getId();

                  const hier::BlockId& nbr_blk_id =
                     encon_mb_id.getBlockId();
                  const int nbr_blk = nbr_blk_id.getBlockValue();

                  /*
                   * The intersection, if any, of the fill_mapped_box 
                   * with the translated domain of the neighbor block
                   * will be an unfilled box added to
                   * unfilled_encon_boxes.
                   */
                  for (tbox::List<hier::GridGeometry::Neighbor>::Iterator
                       ni(neighbors); ni; ni++) {

                     if(ni().getBlockNumber() == nbr_blk) {

                        hier::BoxList unfilled_box_for_encon(
                           ni().getTransformedDomain());
                        unfilled_box_for_encon.refine(
                           d_dst_level->getRatioToLevelZero());
                        unfilled_box_for_encon.intersectBoxes(
                           fill_mapped_box);

                        if (!unfilled_box_for_encon.isEmpty()) {

                           NeighborSet& unfilled_encon_nabrs =
                              unfilled_encon_nbrhood_set[encon_mb_id];

                           hier::IntVector offset(ni().getShift());
                           offset *= d_dst_level->getRatioToLevelZero();
                           hier::Transformation transformation(
                              ni().getRotationIdentifier(), offset);

                           for (hier::BoxList::Iterator
                                bi(unfilled_box_for_encon); bi; bi++) {

                              hier::Box unfilled_box(bi());
                              transformation.inverseTransform(
                                 unfilled_box);

                              hier::Box unfilled_mapped_box(
                                 unfilled_box,
                                 ++last_unfilled_local_id,
                                 fill_mapped_box.getOwnerRank(),
                                 nbr_blk_id);

                              unfilled_encon_boxes.insert(
                                 unfilled_encon_boxes.end(),
                                 unfilled_mapped_box);

                              unfilled_encon_nabrs.insert(
                                 unfilled_encon_nabrs.end(),
                                 unfilled_mapped_box);

                           }

                        }

                        /*
                         * Break out of loop after finding the neighbor
                         * block.
                         */
                        break;
                     }
                  }
               }
            }
         }
      }
   }

   const MappedBoxLevel& dst_mapped_box_level = dst_to_fill.getBase();

   /*
    * Initialize MappedBoxLevel and Connector.
    */
   d_unfilled_encon_box_level = new MappedBoxLevel(
      unfilled_encon_boxes,
      dst_mapped_box_level.getRefinementRatio(),
      dst_mapped_box_level.getGridGeometry(),
      dst_mapped_box_level.getMPI());

   encon_to_unfilled_encon = new Connector(
      *(d_encon_level->getMappedBoxLevel()),
      *d_unfilled_encon_box_level,
      dst_to_fill.getConnectorWidth(),
      unfilled_encon_nbrhood_set,
      MappedBoxLevel::DISTRIBUTED);

}


/*
 **************************************************************************
 * Shear off parts of unfilled boxes that lie outside non-periodic
 * domain boundaries.
 *
 * Update the overlap Connector dst_to_unfilled.
 *
 * If domain is fully periodic, this is a no-op.
 **************************************************************************
 */

void RefineSchedule::finishScheduleConstruction_shearUnfilledBoxesOutsideNonperiodicBoundaries(
   hier::MappedBoxLevel &unfilled,
   hier::Connector &dst_to_unfilled,
   const tbox::Pointer<hier::PatchHierarchy> &hierarchy)
{
   const tbox::Dimension &dim(unfilled.getDim());

   const bool fully_periodic = d_num_periodic_directions == dim.getValue();

   if ( fully_periodic ) {
      return;
   }

   /*
    * Shearing is the removal of parts of unfilled_boxes that lie
    * along physical boundaries.  We should not fill these because
    * they would be filled by boundary conditions.
    *
    * We bypass shearing for fully_periodic domains, where it would
    * be a no-op anyway.
    */

   const MappedBoxLevel& periodic_domain_mapped_box_level =
      hierarchy->getDomainMappedBoxLevel();
   TBOX_ASSERT(
      periodic_domain_mapped_box_level.getParallelState() ==
      MappedBoxLevel::GLOBALIZED);

   t_shear->start();

   hier::MappedBoxLevelConnectorUtils edge_utils;
   hier::OverlapConnectorAlgorithm oca;

   // Shearing for the mapped_box_level.
   Connector unfilled_to_periodic_domain(
      *d_unfilled_mapped_box_level,
      periodic_domain_mapped_box_level,
      dst_to_unfilled.getConnectorWidth());
   oca.findOverlaps(unfilled_to_periodic_domain);

   Connector unfilled_to_sheared;
   MappedBoxLevel sheared_mapped_box_level(dim);
   edge_utils.computeInternalPartsForMultiblock(
      sheared_mapped_box_level,
      unfilled_to_sheared,
      unfilled_to_periodic_domain,
      hier::IntVector::getZero(dim) );

   t_modify_connector->start();
   hier::MappingConnectorAlgorithm mca;
   mca.modify(dst_to_unfilled,
              unfilled_to_sheared,
              d_unfilled_mapped_box_level.getPointer());
   t_modify_connector->stop();
   dst_to_unfilled.eraseEmptyNeighborSets();

   t_shear->stop();

   return;
}



/*
 **************************************************************************
 * Set up the supplemental MappedBoxLevel The supplemental level will
 * be used to fill data in d_unfilled_mapped_box_level.
 *
 * - Start with the d_unfilled_mapped_box_level.
 * - Coarsen unfilled boxes
 * - Shear off parts outside non-periodic boundaries, if needed.
 *
 * The supp_mapped_box_level is to be filled by the next coarser level
 * in the hierarchy.
 *
 * Build Connectors d_dst_to_supp, d_supp_to_dst and d_supp_to_unfilled.
 * The data to be filled on the supp level
 * will be the max stencil width.
 *
 * We set d_dst_to_supp's width big enough so each dst Box,
 * grown by this width, nests its potential supp MappedBoxes.  Note
 * that d_dst_to_supp is incomplete because each dst Box only
 * has edges to supp MappedBoxes it generated.
 **************************************************************************
 */
void RefineSchedule::finishScheduleConstruction_setupSupplementalMappedBoxLevel(
   hier::MappedBoxLevel &supp_mapped_box_level,
   const hier::MappedBoxLevel &hiercoarse_mapped_box_level,
   const hier::Connector &dst_to_unfilled)
{
   t_setup_supp_mapped_box_level->start();

   const tbox::Dimension &dim(dst_to_unfilled.getBase().getDim());

   const hier::IntVector dst_hiercoarse_ratio(
      d_dst_level->getRatioToLevelZero()
      / hiercoarse_mapped_box_level.getRefinementRatio() );

   const hier::NeighborhoodSet& dst_eto_unfilled =
      dst_to_unfilled.getNeighborhoodSets();

   const bool fully_periodic = d_num_periodic_directions == dim.getValue();

   const tbox::ConstPointer<hier::GridGeometry> &grid_geometry(
      dst_to_unfilled.getBase().getGridGeometry());

   const int nblocks(grid_geometry->getNumberBlocks());


   hier::IntVector big_grow_vector(dim, 0);
   if (d_num_periodic_directions > 0) {
      for (int dir = 0; dir < dim.getValue(); dir++) {
         if (d_periodic_shift(dir)) {
            big_grow_vector(dir) = BIG_GHOST_CELL_WIDTH;
         }
      }
   }

   tbox::Array<hier::BoxList> coarser_physical_domain(
      nblocks, hier::BoxList(dim));
   tbox::Array<hier::BoxList> coarser_shear_domain(
      nblocks, hier::BoxList(dim));
   tbox::Array<bool> do_coarse_shearing(nblocks);
   for (int b = 0; b < nblocks; b++) {
      // coarser_physical_domain[b] = hiercoarse_level->getPhysicalDomain(hier::BlockId(b));
      grid_geometry->computePhysicalDomain(
         coarser_physical_domain[b],
         hiercoarse_mapped_box_level.getRefinementRatio(),
         hier::BlockId(b));

      do_coarse_shearing[b] = (!fully_periodic && !d_domain_is_one_box[b]);

      coarser_shear_domain[b] = coarser_physical_domain[b];

      if (do_coarse_shearing[b]) {
         t_coarse_shear->start();

         if (d_num_periodic_directions > 0) {
            coarser_shear_domain[b].grow(big_grow_vector);
         }

         coarser_shear_domain[b].simplifyBoxes();

         t_coarse_shear->stop();
      }
   }



   supp_mapped_box_level.initialize(
      hiercoarse_mapped_box_level.getRefinementRatio(),
      grid_geometry,
      d_unfilled_mapped_box_level->getMPI());

   hier::NeighborhoodSet dst_eto_supp, supp_eto_unfilled;

   /*
    * This loop builds up supp_mapped_box_level, dst_eto_supp and
    *  supp_eto_unfilled.
    */
   for (hier::NeighborhoodSet::const_iterator ei = dst_eto_unfilled.begin();
        ei != dst_eto_unfilled.end(); ++ei) {

      const hier::BoxId& dst_mapped_box_mbid = ei->first;
      const NeighborSet& dst_unfilled_parts = ei->second;

      const int dst_blk = dst_mapped_box_mbid.getBlockId().getBlockValue();

      for (hier::MappedBoxSet::const_iterator ni = dst_unfilled_parts.begin();
           ni != dst_unfilled_parts.end(); ++ni) {

         const hier::Box& unfilled_mapped_box = *ni;
         hier::Box supp_box(unfilled_mapped_box);
         supp_box.coarsen(dst_hiercoarse_ratio);


         if (do_coarse_shearing[dst_blk] &&
             (d_dst_level->patchTouchesRegularBoundary(
                 dst_mapped_box_mbid))) {

            hier::BoxList sheared_supp_boxes(supp_box);
            sheared_supp_boxes.intersectBoxes(coarser_shear_domain[dst_blk]);
            sheared_supp_boxes.simplifyBoxes();

            (void)hier::BoxUtilities::extendBoxesToDomainBoundary(
               sheared_supp_boxes,
               coarser_physical_domain[dst_blk],
               d_max_stencil_width);
            /*
             * Connector widths must be big enough to make sure
             * we have complete sets after extending mapped_boxes to boundary!
             */
            if (sheared_supp_boxes.size() > 0) {

               NeighborSet& supp_nabrs = dst_eto_supp[dst_mapped_box_mbid];

               for (hier::BoxList::Iterator b(sheared_supp_boxes); b; b++) {
                  const hier::Box& supp_mapped_box =
                     *supp_mapped_box_level.addBox(*b, (*ni).getBlockId());
                  supp_nabrs.insert(supp_mapped_box);

                  /*
                   * Note that each supp_mapped_box must have at least one
                   * unfilled_nabr and may have multiple.
                   */
                  supp_eto_unfilled[supp_mapped_box.getId()].insert(
                     unfilled_mapped_box);

               }

            }

         } else {

            /*
             * If the supp_box is less than a ghost width
             * (d_max_stencil_width for the supplemental level) of a
             * physical boundary, extend it to the boundary.
             *
             * Note: If we end up extending the supp_box, we may
             * fail sanity checks further down because have not
             * accounted for this extension in the various
             * Connector widths.  However, it is not likely that
             * the extension would have any effects on computations
             * because, being so close to the physical boundary,
             * they should not create any new relationships.  This
             * is just a hunch and need to be rigously verified.
             *
             * For now, just warn if the box is grown.
             */
            const hier::Box save_supp_box(supp_box);
            (void)hier::BoxUtilities::extendBoxToDomainBoundary(
               supp_box,
               coarser_physical_domain[dst_blk],
               d_max_stencil_width);
            if (! supp_box.isSpatiallyEqual(save_supp_box)) {
               TBOX_WARNING("Supplemental box " << save_supp_box
                            << " was extended to " << supp_box
                            << " at a physical boundary.  This is"
                            << " probably ok but a rigorous proof"
                            << " that won't cause problem is currently"
                            << " lacking.  Expect some sanity checks"
                            << " to fail and a slim chance that the"
                            << " schedule generate will be bad.");
            }

            const hier::Box& supp_mapped_box = *supp_mapped_box_level.addBox(
               supp_box, (*ni).getBlockId());
            NeighborSet& supp_nabrs = dst_eto_supp[dst_mapped_box_mbid];
            supp_nabrs.insert(supp_mapped_box);
            supp_eto_unfilled[supp_mapped_box.getId()].insert(
               unfilled_mapped_box);

         }

      }
   }


   /*
    * Width of dst-->supp is
    *
    * - width of dst-->fill, but rounded up so it extends
    *   the growth of supp caused by coarsening unfilled.
    *
    * - extended by the stencil width, where supp has its ghost data.
    *
    * This width states that each dst box sees all of its
    * supplemental boxes, including the ghost cells in the
    * supplemental boxes.
    */
   const hier::IntVector dst_to_supp_width =
      (hier::IntVector::ceiling(dst_to_unfilled.getConnectorWidth(),
                                dst_hiercoarse_ratio) + d_max_stencil_width)
      * dst_hiercoarse_ratio;

   d_dst_to_supp.swapInitialize(
      dst_to_unfilled.getBase(),
      supp_mapped_box_level,
      dst_to_supp_width,
      dst_eto_supp);

   d_supp_to_unfilled.swapInitialize(
      supp_mapped_box_level,
      *d_unfilled_mapped_box_level,
      hier::IntVector::getZero(dim),
      supp_eto_unfilled);


   /*
    * Get the transpose of d_dst_to_supp, which is simple to compute
    * because we know the edges are all local.
    */
   t_misc2->start();
   d_supp_to_dst.initializeToLocalTranspose(d_dst_to_supp);
   t_misc2->stop();


   if (s_extra_debug) {
      /*
       * We have set up supp to nest in dst^dst_to_supp_width to
       * ensure dst sees all of supp and also supp's ghosts.  Note
       * that supp's relevant ghost data width is
       * d_max_stencil_width.
       *
       * The nesting assures that when bridging across dst<==>supp
       * for supp<==>hiercoarse, we get a complete overlap
       * Connectors.
       */
      const tbox::Dimension& dim(supp_mapped_box_level.getDim());
      const hier::IntVector &zero_vector(hier::IntVector::getZero(dim));
      hier::MappedBoxLevelConnectorUtils mblc_utils;
      mblc_utils.setSanityCheckMethodPreconditions(false);
      mblc_utils.setSanityCheckMethodPostconditions(false);
      tbox::plog << "\nsupp_mapped_box_level:\n" << supp_mapped_box_level.format("S-> ", 2)
                 << "\ndst_mapped_box_level:\n" << d_dst_to_supp.getBase().format("D-> ", 2)
                 << "\nd_dst_to_supp:\n" << d_dst_to_supp.format("DS-> ", 2)
                 << "\nd_supp_to_dst:\n" << d_supp_to_dst.format("SD-> ", 2)
                 << std::endl;
      bool locally_nests;
      if ( ! mblc_utils.baseNestsInHeadForMultiblock(
              &locally_nests,
              supp_mapped_box_level,
              d_dst_to_supp.getBase(),
              zero_vector,
              d_dst_to_supp.getConnectorWidth(),
              zero_vector,
              NULL) ) {
         TBOX_ERROR("RefineSchedule::finishScheduleConstruction: supp does\n"
                    <<"to nest in dst.\n");
      }

      TBOX_ASSERT(d_supp_to_dst.checkTransposeCorrectness(d_dst_to_supp) == 0);
      TBOX_ASSERT(d_dst_to_supp.checkTransposeCorrectness(d_supp_to_dst) == 0);
   }

   t_setup_supp_mapped_box_level->stop();

   return;
}




/*
 **************************************************************************
 * Build the Connectors between the supplementary MappedBoxLevel and the
 * coarser level on the hierarchy.  The coarser level on the hierarchy
 * at the same resolution as the supplementary and will be used as the
 * source for filling the supplementary level.
 **************************************************************************
 */
void RefineSchedule::finishScheduleConstruction_connectSuppToHiercoarse(
   hier::Connector &supp_to_hiercoarse,
   hier::Connector &hiercoarse_to_supp,
   const hier::MappedBoxLevel &supp_mapped_box_level,
   const tbox::Pointer<hier::PatchHierarchy> &hierarchy,
   const int next_coarser_ln,
   const hier::Connector &dst_to_src,
   const hier::Connector &src_to_dst,
   const bool dst_is_supplemental_level )
{
   const tbox::Dimension &dim(hierarchy->getDim());
   const hier::IntVector &zero_vector(hier::IntVector::getZero(dim));
   hier::OverlapConnectorAlgorithm oca;
   hier::MappedBoxLevelConnectorUtils edge_utils;


   /*
    * Compute the Connector widths required to properly
    * set up overlap Connectors.
    */

   std::vector<hier::IntVector> self_connector_widths;
   std::vector<hier::IntVector> fine_connector_widths;

   RefineScheduleConnectorWidthRequestor rscwri;
   rscwri.computeRequiredConnectorWidths(self_connector_widths,
                                         fine_connector_widths,
                                         *hierarchy);

   const hier::MappedBoxLevel &dst_mapped_box_level(
      *d_dst_level->getMappedBoxLevel());

   const tbox::Pointer<hier::PatchLevel> hiercoarse_level(
      hierarchy->getPatchLevel(next_coarser_ln));

   const hier::MappedBoxLevel &hiercoarse_mapped_box_level(
      *hiercoarse_level->getMappedBoxLevel());

   /*
    * Correspondance between consecutive recursions:
    *
    * This recursion:   supp<==>hiercoarse, supp_to_hiercoarse
    *
    * Next recursion:    dst<==>src,            dst_to_src
    *
    * We will bridge across dst for supp<==>hiercoarse.
    *
    * The bridge legs d_supp_to_dst and d_dst_to_supp are not
    * complete, as usually required to guarantee that the bridge is
    * complete.  However, supp nests in dst^gcw, so they need not be
    * complete.
    */

   /*
    * We now need supp<==>hiercoarse and can get it with the
    * bridge supp<==>dst<==>hiercoarse.
    *
    * To do this bridge, we require dst<==>hiercoarse.  Look
    * them up in the PersistentOverlapConnectors.  The Connector
    * width we need and expect to have been precomputed depends
    * on whether dst is a supplemental level from a previous
    * recursive RefineSchedule.
    *
    * If we cannot look up dst<==>hiercoarse, we will attempt to
    * bridge for it.
    */

   const Connector* dst_to_hiercoarse = NULL;
   const Connector* hiercoarse_to_dst = NULL;
   Connector bridged_dst_to_hiercoarse;
   Connector bridged_hiercoarse_to_dst;

   if (dst_is_supplemental_level) {
      /*
       * Here, dst is a supplemental (temporary) level created
       * from the unfilled boxes of another dst (which itself
       * may also be another supplemental level).  We cannot use
       * the required Connector widths for levels in the
       * hierarchy.  We must compute the correct required width
       * for this particular temporary level.  We know how
       * to do this because the temporary level is created by
       * RefineSchedule code.
       */

      hier::IntVector min_dst_to_hiercoarse_width(d_max_scratch_gcw);
      min_dst_to_hiercoarse_width.max(d_max_stencil_width);
      min_dst_to_hiercoarse_width.max(d_boundary_fill_ghost_width);
      hier::IntVector min_hiercoarse_to_dst_width =
         Connector::convertHeadWidthToBase(
            hiercoarse_mapped_box_level.getRefinementRatio(),
            d_dst_level->getMappedBoxLevel()->getRefinementRatio(),
            min_dst_to_hiercoarse_width);

      const bool has_dst_to_hiercoarse =
         d_dst_level->getMappedBoxLevel()->getPersistentOverlapConnectors().
         hasConnector(
            hiercoarse_mapped_box_level,
            min_dst_to_hiercoarse_width);
      const bool has_hiercoarse_to_dst =
         hiercoarse_mapped_box_level.getPersistentOverlapConnectors().
         hasConnector(
            *d_dst_level->getMappedBoxLevel(),
            min_hiercoarse_to_dst_width);

      if (has_dst_to_hiercoarse && has_hiercoarse_to_dst) {

         dst_to_hiercoarse =
            &d_dst_level->getMappedBoxLevel()->getPersistentOverlapConnectors()
            .findConnector(
               hiercoarse_mapped_box_level,
               min_dst_to_hiercoarse_width);

         hiercoarse_to_dst =
            &hiercoarse_mapped_box_level.getPersistentOverlapConnectors()
            .findConnector(
               *d_dst_level->getMappedBoxLevel(),
               min_hiercoarse_to_dst_width);

         TBOX_ASSERT(
            dst_to_hiercoarse->getBase() == *d_dst_level->getMappedBoxLevel());
         TBOX_ASSERT(
            dst_to_hiercoarse->getConnectorWidth() >= d_max_scratch_gcw);
         TBOX_ASSERT(
            dst_to_hiercoarse->getConnectorWidth() >= d_max_stencil_width);
         TBOX_ASSERT(
            dst_to_hiercoarse->getConnectorWidth() >=
            d_boundary_fill_ghost_width);
         TBOX_ASSERT(
            hiercoarse_to_dst->getHead() == *d_dst_level->getMappedBoxLevel());

      } else {
         /*
          * Connectors dst<==>hiercoarse are not provided.
          * We have to bridge through src for it.
          * (This requires dst<==>src<==>hiercoarse.)
          */
         const hier::MappedBoxLevel& src_mapped_box_level =
            dst_to_src.getHead();
         if (*hierarchy->getMappedBoxLevel(next_coarser_ln+1) !=
             src_mapped_box_level) {
            TBOX_ERROR("Missing dst<==>hiercoarse connector and\n"
                       << "src is not from hierarchy.");
         }
         if (s_extra_debug) {
            TBOX_WARNING(
               "RefineSchedule::finishScheduleConstruction bridging through src for dst<==>hiercoarse\n");
         }

         hier::IntVector hiercoarse_to_src_width =
            fine_connector_widths[next_coarser_ln];
         hier::IntVector src_to_hiercoarse_width =
            hiercoarse_to_src_width*d_src_level->getRatioToCoarserLevel();

         /*
          * Using lh to get required Connector width assumes that
          * the src level has same refinement ratio as
          * next_coarser_ln+1, but for Richardson extrapolation,
          * that is not the case, so we have to adjust.
          */
         if (d_src_level->getMappedBoxLevel()->getRefinementRatio() <=
             hierarchy->getMappedBoxLevel(next_coarser_ln + 1)->getRefinementRatio()) {
            src_to_hiercoarse_width *=
               d_src_level->getMappedBoxLevel()->getRefinementRatio();
            src_to_hiercoarse_width /= hierarchy->getMappedBoxLevel(
               next_coarser_ln + 1)->getRefinementRatio();
         } else if (d_src_level->getMappedBoxLevel()->getRefinementRatio()
                    >=
                    hierarchy->getMappedBoxLevel(next_coarser_ln
                                                 + 1)->getRefinementRatio()) {
            src_to_hiercoarse_width *= hierarchy->getMappedBoxLevel(
               next_coarser_ln + 1)->getRefinementRatio();
            src_to_hiercoarse_width /=
               d_src_level->getMappedBoxLevel()->getRefinementRatio();
         }

         const hier::Connector& src_to_hiercoarse =
            d_src_level->getMappedBoxLevel()->getPersistentOverlapConnectors()
            .findConnector(
               hiercoarse_mapped_box_level,
               src_to_hiercoarse_width);
         const hier::Connector& hiercoarse_to_src =
            hiercoarse_mapped_box_level.getPersistentOverlapConnectors()
            .findConnector(
               *d_src_level->getMappedBoxLevel(),
               hiercoarse_to_src_width);

         if (s_barrier_and_time) {
            t_bridge_dst_hiercoarse->barrierAndStart();
         }

         oca.setSanityCheckMethodPostconditions(true);
         /*
          * Don't use the strict bridge theorem here because it
          * cannot guarantee sufficient width.  We know from how
          * dst nests in hiercoarse what output Connector width
          * can guarantee that all dst MappedBoxes are seen by a
          * hiercoarse MappedBox.
          */
         oca.bridge(
            bridged_dst_to_hiercoarse,
            bridged_hiercoarse_to_dst,
            dst_to_src,
            src_to_hiercoarse,
            hiercoarse_to_src,
            src_to_dst,
            fine_connector_widths[next_coarser_ln]);

         if (s_barrier_and_time) {
            t_bridge_dst_hiercoarse->stop();
         }

         /*
          * hiercoarse has periodic images but dst does not.
          * Consequently, the bridge for dst<==>hiercoarse
          * may not catch all periodic neighbors.  Because we
          * don't need periodic neighbors for
          * dst<==>hiercoarse, we will remove them to make
          * dst<==>hiercoarse proper transposes.
          */
         bridged_dst_to_hiercoarse.removePeriodicRelationships();
         bridged_hiercoarse_to_dst.removePeriodicRelationships();

         if (s_extra_debug) {
            size_t err1 =
               bridged_dst_to_hiercoarse.checkTransposeCorrectness(
                  bridged_hiercoarse_to_dst);
            if (err1) tbox::perr
               << "bridged_dst_to_hiercoarse failed transpose correctness."
               << std::endl;
            size_t err2 =
               bridged_hiercoarse_to_dst.checkTransposeCorrectness(
                  bridged_dst_to_hiercoarse);
            if (err2) tbox::perr
               << "bridged_hiercoarse_to_dst failed transpose correctness."
               << std::endl;
            if (err1 || err2) {
               tbox::perr << "dst:\n"
                          << dst_mapped_box_level.format("DEST->", 2)
                          << "hiercoarse:\n"
                          << hierarchy->getMappedBoxLevel(next_coarser_ln)->format("HCRS->", 2)
                          << "bridged_dst_to_hiercoarse:\n"
                          << bridged_dst_to_hiercoarse.format("bDH->", 2)
                          << "bridged_hiercoarse_to_dst:\n"
                          << bridged_hiercoarse_to_dst.format("bHD->", 2)
                          << std::endl;
               dst_mapped_box_level.getMPI().Barrier();
               TBOX_ERROR(
                  "Bridged dst<==>hiercoarse (next_coarser_ln="
                  << next_coarser_ln << ") have problems as reported above.\n" );
            }
         }
         dst_to_hiercoarse = &bridged_dst_to_hiercoarse;
         hiercoarse_to_dst = &bridged_hiercoarse_to_dst;
      } // End block bridging for dst<==>hiercoarse.

   } else { /* !dst_is_supplemental_level */

      /*
       * dst may be the level next_coarser_ln+1 on the
       * hierarchy, or it could be a coarsened version (generated by Richardson extrapolation).
       * Make sure Connector width based on dst is
       * properly adjusted for the correct case.  The
       * refinement ratio of next_coarser_ln is not
       * affected by the choice of which case, so use it
       * as the basis for the transpose width.
       *
       * Note: This method always assumes that dst is
       * related to a level in the hierarchy.  src need
       * not be.  If that is not the case, we need to
       * rewrite a few things.  For scalability, we expect
       * to have precomputed Connectors dst<==>hiercoarse.
       * We simply don't have enough information in
       * RefineSchedule to compute them efficiently.
       */
      const hier::IntVector &hiercoarse_to_dst_width(fine_connector_widths[next_coarser_ln]);
      const hier::IntVector dst_to_hiercoarse_width(
         hier::Connector::convertHeadWidthToBase(
            d_dst_level->getMappedBoxLevel()->getRefinementRatio(),
            hiercoarse_mapped_box_level.getRefinementRatio(),
            hiercoarse_to_dst_width));
      dst_to_hiercoarse = &d_dst_level->getMappedBoxLevel()->
         getPersistentOverlapConnectors().
         findConnector(hiercoarse_mapped_box_level,
                       dst_to_hiercoarse_width);
      hiercoarse_to_dst = &hiercoarse_mapped_box_level.
         getPersistentOverlapConnectors().
         findConnector(*d_dst_level->getMappedBoxLevel(),
                       hiercoarse_to_dst_width);
      TBOX_ASSERT(
         dst_to_hiercoarse->getBase() == *d_dst_level->getMappedBoxLevel());
      TBOX_ASSERT(
         dst_to_hiercoarse->getConnectorWidth() >= d_max_scratch_gcw);
      TBOX_ASSERT(
         dst_to_hiercoarse->getConnectorWidth() >= d_max_stencil_width);
      TBOX_ASSERT(
         dst_to_hiercoarse->getConnectorWidth() >=
         d_boundary_fill_ghost_width);
      TBOX_ASSERT(
         hiercoarse_to_dst->getHead() == *d_dst_level->getMappedBoxLevel());
      if (s_extra_debug) {
         if (dst_to_hiercoarse != NULL || hiercoarse_to_dst != NULL) {
            size_t err1 = dst_to_hiercoarse ? oca.checkOverlapCorrectness(
               *dst_to_hiercoarse) : 0;
            if (err1) tbox::perr
               << "\ndst_to_hiercoarse failed overlap correctness." << std::endl;
            size_t err2 = hiercoarse_to_dst ? oca.checkOverlapCorrectness(
               *hiercoarse_to_dst) : 0;
            if (err2) tbox::perr
               << "\nhiercoarse_to_dst failed overlap correctness." << std::endl;
            size_t err3 = (dst_to_hiercoarse && hiercoarse_to_dst) ?
               dst_to_hiercoarse->checkTransposeCorrectness(*hiercoarse_to_dst) :
               0;
            if (err3) tbox::perr
               << "dst_to_hiercoarse failed transpose correctness."
               << std::endl;
            size_t err4 = (dst_to_hiercoarse && hiercoarse_to_dst) ?
               hiercoarse_to_dst->checkTransposeCorrectness(*dst_to_hiercoarse) :
               0;
            if (err4) tbox::perr
               << "hiercoarse_to_dst failed transpose correctness."
               << std::endl;
            if (err1 || err2 || err3 || err4) {
               TBOX_ERROR(
                  "dst<==>hiercoarse have problems as reported above."
                  << "dst:\n" << d_dst_level->getMappedBoxLevel()->format("DEST->", 2)
                  << "hiercoarse:\n" << hiercoarse_mapped_box_level.format("HCRS->", 2)
                  << "dst_to_hiercoarse:\n" << dst_to_hiercoarse->format("DH->", 2)
                  << "hiercoarse_to_dst:\n" << hiercoarse_to_dst->format("HD->", 2)
                  << std::endl);
            }
         }
      }
   }


   /*
    * Compute supp<==>hiercoarse by bridging
    * supp<==>dst<==>hiercoarse.
    *
    * Don't use the strict bridge theorem here because it
    * cannot guarantee sufficient width.  We know from how
    * dst nests in hiercoarse what output Connector width
    * can guarantee that all dst MappedBoxes are seen by a
    * hiercoarse MappedBox.
    */

   if (s_barrier_and_time) {
      t_bridge_supp_hiercoarse->barrierAndStart();
   }
   oca.bridge(
      supp_to_hiercoarse,
      hiercoarse_to_supp,
      d_supp_to_dst,
      *dst_to_hiercoarse,
      *hiercoarse_to_dst,
      d_dst_to_supp,
      fine_connector_widths[next_coarser_ln]);
   if (s_barrier_and_time) {
      t_bridge_supp_hiercoarse->stop();
   }

   if (s_extra_debug) {
      bool locally_nests = false;
      if ( !edge_utils.baseNestsInHead(
              &locally_nests,
              supp_to_hiercoarse.getBase(),
              supp_to_hiercoarse.getHead(),
              d_max_stencil_width,
              supp_to_hiercoarse.getConnectorWidth(),
              zero_vector,
              &hierarchy->getDomainSearchTree(hier::BlockId(0))) ) {
         TBOX_ERROR(
            "RefineSchedule::finishScheduleConstruction: library error:\n"
            <<"hiercoarse--->supp Connector may not see all of supp.\n"
            <<"This is because supp, grown by the max stencil width\n"
            <<"of " << d_max_stencil_width << " does not nest in\n"
            <<"hiercoarse grown by the width of hiercoarse--->supp.\n"
            <<"locally_nest = " << locally_nests
            <<"\nsupp:\n" << supp_mapped_box_level.format("S: ",2)
            <<"\nhiercoarse:\n" << hiercoarse_mapped_box_level.format("H: ",2)
            <<"\nhiercoarse_to_supp:\n"
            << hiercoarse_to_supp.format("HS: ",3));
      }
   }

   TBOX_ASSERT(&supp_to_hiercoarse.getBase() == &supp_mapped_box_level);
   TBOX_ASSERT(&hiercoarse_to_supp.getHead() == &supp_mapped_box_level);
   TBOX_ASSERT(
      supp_to_hiercoarse.getConnectorWidth() >=
      hier::IntVector::ceiling(d_max_stencil_width,
                               d_dst_to_supp.getRatio()));
   TBOX_ASSERT(
      hiercoarse_to_supp.getConnectorWidth() >=
      hier::IntVector::ceiling(d_max_stencil_width,
                               d_dst_to_supp.getRatio()));


   if (d_num_periodic_directions > 0) {
      /*
       * supp_mapped_box_level does not have any periodic images yet.
       * Remove periodic relationships that might have been added when
       * bridging for supp<==>hiercoarse.
       *
       * This makes supp<==>hiercoarse proper transposes and prevent
       * elusive bugs caused by extraneous periodic relationships.
       */
      supp_to_hiercoarse.removePeriodicRelationships();
      hiercoarse_to_supp.removePeriodicRelationships();
   }


   if (s_extra_debug) {
      sanityCheckSupplementalAndHiercoarseLevels(
         supp_to_hiercoarse,
         hiercoarse_to_supp,
         hierarchy,
         next_coarser_ln);
   }

   return;
}




/*
 **************************************************************************
 **************************************************************************
 */
void RefineSchedule::sanityCheckSupplementalAndHiercoarseLevels(
   const hier::Connector &supp_to_hiercoarse,
   const hier::Connector &hiercoarse_to_supp,
   const tbox::Pointer<hier::PatchHierarchy> &hierarchy,
   const int next_coarser_ln)
{

   /*
    * supp_to_hiercoarse and hiercoarse_to_supp should be proper
    * transposes.
    */
   size_t err1 = supp_to_hiercoarse.checkTransposeCorrectness(
      hiercoarse_to_supp);
   if (err1) tbox::perr
      << "supp_to_hiercoarse failed transpose correctness."
      << std::endl;
   size_t err2 = hiercoarse_to_supp.checkTransposeCorrectness(
      supp_to_hiercoarse);
   if (err2) tbox::perr
      << "hiercoarse_to_supp failed transpose correctness."
      << std::endl;

   /*
    * Compute the Connector widths required to properly
    * set up overlap Connectors.
    */

   std::vector<hier::IntVector> self_connector_widths;
   std::vector<hier::IntVector> fine_connector_widths;

   RefineScheduleConnectorWidthRequestor rscwri;
   rscwri.computeRequiredConnectorWidths(self_connector_widths,
                                         fine_connector_widths,
                                         *hierarchy);

   /*
    * To work properly, we must ensure that
    * supplemental^d_max_stencil_width nests in hiercoarse^fine_connector_width.
    *
    * The nesting guarantees that the Connectors sees enough of its
    * surroundings to generate all the necessary relationships in
    * further RefineSchedule recursions.  We know that
    * supp_to_hiercoarse and hiercoarse_to_supp are not complete, but
    * this nesting guarantees we still have enough relationships to
    * avoid miss any relationships when we use these Connectors in
    * bridge operations.
    */
   Connector complete_supp_to_hiercoarse(
      supp_to_hiercoarse.getBase(),
      supp_to_hiercoarse.getHead(),
      supp_to_hiercoarse.getConnectorWidth());
   hier::OverlapConnectorAlgorithm oca;
   oca.findOverlaps(complete_supp_to_hiercoarse);
   MappedBoxLevel external(hiercoarse_to_supp.getBase().getDim());
   Connector supp_to_external;
   hier::MappedBoxLevelConnectorUtils edge_utils;
   edge_utils.computeExternalParts(
      external,
      supp_to_external,
      complete_supp_to_hiercoarse,
      fine_connector_widths[next_coarser_ln]-d_max_stencil_width,
      hierarchy->getPeriodicDomainSearchTree(hier::BlockId(0)));
   supp_to_external.eraseEmptyNeighborSets();
   int err3 = supp_to_external.getGlobalNumberOfRelationships();
   if (err3) {
      tbox::perr << "Some parts of supp lies outside of where we\n"
                 << "guarantee support for recursive RefineSchedule.\n"
                 << supp_to_external.format("SE: ", 2);
   }

   if (err1 || err2 || err3) {
      TBOX_ERROR(
         "supp<==>hiercoarse have problems as reported above.\n"
         << "supp:\n" << supp_to_hiercoarse.getBase().format("SUPP->", 2)
         << "hiercoarse:\n" << supp_to_hiercoarse.getHead().format("HCRS->", 2)
         << "dst_to_supp:\n" << d_dst_to_supp.format("DS->", 2)
         << "supp_to_hiercoarse:\n" << supp_to_hiercoarse.format("SH->", 2)
         << "hiercoarse_to_supp:\n" << hiercoarse_to_supp.format("HS->", 2));
   }
   return;
}



/*
 **************************************************************************
 *
 * Execute the stored communication schedule that copies data into the
 * destination component of the destination level.  The algorithm is as
 * follows:
 *
 * (1) Allocate scratch space on the destination level if it does
 *    not exist.
 * (2) Call the recursive fill function that will fill the scratch
 *    space of the destination level.
 * (3) If the scratch and destination spaces are not the same,
 *    then copy locally from the scratch into the destination.
 * (4) Deallocate any previously allocated data.
 *
 **************************************************************************
 */

void RefineSchedule::fillData(
   double fill_time,
   bool do_physical_boundary_fill) const
{
   t_fill_data->start();

   /*
    * Set the refine items and time for all transactions.  These items will
    * be shared by all transaction objects in the communication schedule.
    */

   d_transaction_factory->setTransactionTime(fill_time);
   d_transaction_factory->setRefineItems(d_refine_items, d_number_refine_items);

   /*
    * Check whether scratch data needs to be allocated on the destination
    * level.  Keep track of those allocated components so that they may be
    * deallocated later.
    */

   hier::ComponentSelector allocate_vector;
   allocateScratchSpace(allocate_vector, d_dst_level, fill_time);

   const int num_blocks = d_dst_level->getGridGeometry()->getNumberBlocks();
   hier::ComponentSelector encon_allocate_vector;
   if (num_blocks > 1) {
      allocateScratchSpace(encon_allocate_vector, d_encon_level, fill_time);
   }

   /*
    * Begin the recursive algorithm that fills from coarser, fills from
    * same, and then fills physical boundaries.
    */

   t_recursive_fill->start();
   recursiveFill(fill_time, do_physical_boundary_fill);
   t_recursive_fill->stop();

   /*
    * Copy the scratch space of the destination level to the destination
    * space.
    */

   copyScratchToDestination();

   /*
    * Deallocate any allocated scratch space on the destination level.
    */

   d_dst_level->deallocatePatchData(allocate_vector);

   if (num_blocks > 1) {
      d_encon_level->deallocatePatchData(encon_allocate_vector);
   }

   /*
    * Unset the refine items for all transactions.  These items are
    * shared by all transaction objects in the communication schedule.
    */

   d_transaction_factory->unsetRefineItems();

   t_fill_data->stop();
}

/*
 **************************************************************************
 *
 * Recursively fill the required fill boxes on the destination level.
 * The scratch component of the level will be filled.  The algorithm
 * is as follows:
 *
 * (1) If we need to get data from a coarser level, then:
 *   (a) allocate scratch data on the coarser level
 *   (b) recursively call this routine to get the data
 *   (c) refine data from the coarse level into this level
 *   (d) deallocate the scratch data on the coarser level
 * (2) Copy data from the same level of refinement
 * (3) Copy data from the physical boundary conditions
 *
 **************************************************************************
 */

void RefineSchedule::recursiveFill(
   double fill_time,
   bool do_physical_boundary_fill) const
{
   /*
    * Copy data from the source interiors of the source level into the ghost
    * cells and interiors of the scratch space on the destination level
    * for data where coarse data takes priority on level boundaries.
    */
   d_coarse_priority_level_schedule->communicate();

   /*
    * If there is a coarser schedule stored in this object, then we will
    * need to get data from a coarser grid level.
    */

   const int num_blocks = d_dst_level->getGridGeometry()->getNumberBlocks();
   if (!d_supp_schedule.isNull()) {

      /*
       * Allocate data on the coarser level and keep track of the allocated
       * components so that they may be deallocated later.
       */

      hier::ComponentSelector allocate_vector;
      allocateScratchSpace(allocate_vector, d_supp_level, fill_time);

      hier::ComponentSelector encon_allocate_vector;
      if (num_blocks > 1) {
         allocateScratchSpace(encon_allocate_vector,
                              d_supp_schedule->d_encon_level,
                              fill_time);
      }

      /*
       * Recursively call the fill routine to fill the required coarse fill
       * boxes on the coarser level.
       */

      d_supp_schedule->recursiveFill(fill_time, do_physical_boundary_fill);

      /*
       * d_supp_level should now be filled.  Now interpolate
       * data from the coarse grid into the fine grid.
       */

      refineScratchData(d_dst_level,
                        d_supp_level,
                        d_supp_to_dst,
                        d_supp_to_unfilled,
                        d_refine_overlaps);

      /*
       * Deallocate the scratch data from the coarse grid.
       */

      d_supp_level->deallocatePatchData(allocate_vector);

      if (num_blocks > 1) {
         d_supp_schedule->d_encon_level->deallocatePatchData(
            encon_allocate_vector);
      }

   }

   if (!d_supp_encon_schedule.isNull()) {

      /*
       * Allocate data on the coarser level and keep track of the allocated
       * components so that they may be deallocated later.
       */

      hier::ComponentSelector allocate_vector;
      allocateScratchSpace(allocate_vector, d_supp_encon_level, fill_time);

      hier::ComponentSelector encon_allocate_vector;
      if (num_blocks > 1) {
         allocateScratchSpace(encon_allocate_vector,
                              d_supp_encon_schedule->d_encon_level,
                              fill_time);
      }


      /*
       * Recursively call the fill routine to fill the required coarse fill
       * boxes on the coarser level.
       */

      d_supp_encon_schedule->recursiveFill(fill_time, do_physical_boundary_fill);

      /*
       * d_supp_encon_level should now be filled.  Now interpolate
       * data from the coarse grid into the fine grid.
       */

      refineScratchData(d_encon_level,
                        d_supp_encon_level,
                        d_supp_encon_to_encon,
                        d_supp_encon_to_unfilled_encon,
                        d_encon_refine_overlaps);

      /*
       * Deallocate the scratch data from the coarse grid.
       */

      d_supp_encon_level->deallocatePatchData(allocate_vector);

      if (num_blocks > 1) {
         d_supp_encon_schedule->d_encon_level->deallocatePatchData(
            encon_allocate_vector);
      }

   }

   /*
    * Copy data from the source interiors of the source level into the ghost
    * cells and interiors of the scratch space on the destination level
    * for data where fine data takes priority on level boundaries.
    */
   d_fine_priority_level_schedule->communicate();

   /*
    * Fill the physical boundaries of the scratch space on the destination
    * level.
    */

   if (do_physical_boundary_fill || d_force_boundary_fill) {
      fillPhysicalBoundaries(fill_time);
   }

   fillSingularityBoundaries(fill_time);
}

/*
 **************************************************************************
 *
 * Fill the physical boundaries of the specified level with data.
 *
 **************************************************************************
 */

void RefineSchedule::fillPhysicalBoundaries(
   double fill_time) const
{
   TBOX_ASSERT(!d_dst_level.isNull());

   d_dst_level->setBoundaryBoxes();

   if (d_refine_patch_strategy) {
      for (hier::PatchLevel::Iterator p(d_dst_level); p; p++) {
         tbox::Pointer<hier::Patch> patch(*p);
         if (patch->getPatchGeometry()->intersectsPhysicalBoundary()) {
            d_refine_patch_strategy->
            setPhysicalBoundaryConditions(*patch,
               fill_time,
               d_boundary_fill_ghost_width);
         }
      }
   }
}

/*
 ********************************************************************
 *
 * Call patch strategy method for filling ghost regions around
 * enhanced connectivity block boundaries.
 *
 ********************************************************************
 */
void RefineSchedule::fillSingularityBoundaries(
   double fill_time) const
{
   TBOX_ASSERT(!d_dst_level.isNull());

   tbox::Pointer<hier::GridGeometry> grid_geometry(
      d_dst_level->getGridGeometry());

   if (grid_geometry->getNumberBlocks() > 1) {

      d_dst_level->setBoundaryBoxes();

      const tbox::Dimension& dim(d_dst_level->getDim());

      hier::IntVector ratio(d_dst_level->getRatioToLevelZero());

      if (d_refine_patch_strategy) {

         for (int bn = 0; bn < grid_geometry->getNumberBlocks(); bn++) {

            hier::BlockId block_id(bn); 

            for (hier::BoxList::Iterator sb(
                    grid_geometry->getSingularityBoxList(bn)); sb; sb++) {

               hier::Box singularity(sb());
               singularity.refine(ratio);

               hier::MappedBoxSetSingleBlockIterator dst_local_iter(
                  d_dst_level->getMappedBoxLevel()->getMappedBoxes(), block_id);

               for ( ; dst_local_iter.isValid(); dst_local_iter++) {

                  const hier::BoxId& mapped_box_id =
                     dst_local_iter->getId();

                  tbox::Pointer<hier::Patch> patch(
                     d_dst_level->getPatch(mapped_box_id));
                  tbox::Pointer<hier::PatchGeometry> pgeom(
                     patch->getPatchGeometry());

                  tbox::Array<hier::BoundaryBox> nboxes =
                     pgeom->getNodeBoundaries();

                  if (nboxes.getSize()) {
                     for (int bb = 0; bb < nboxes.getSize(); bb++) {
                        hier::Box intersection((nboxes[bb].getBox())
                           * singularity);
                        if (!(intersection.empty())) {
                           hier::Box fill_box(
                              pgeom->getBoundaryFillBox(nboxes[bb],
                                 patch->getBox(),
                                 d_boundary_fill_ghost_width));

                           if (!(fill_box.empty())) {

                              d_refine_patch_strategy->
                                 fillSingularityBoundaryConditions(
                                    *patch, *d_encon_level,
                                    d_dst_to_encon,
                                    fill_time, fill_box, nboxes[bb]);

                           }
                        }
                     }
                  }

                  if (dim == tbox::Dimension(3)) {
                     tbox::Array<hier::BoundaryBox> eboxes =
                        pgeom->getEdgeBoundaries();

                     if (eboxes.getSize()) {
                        for (int bb = 0; bb < eboxes.getSize(); bb++) {
                           hier::Box intersection(
                              (eboxes[bb].getBox()) * singularity);
                           if (!(intersection.empty())) {
                              hier::Box fill_box(
                                 pgeom->getBoundaryFillBox(eboxes[bb],
                                    patch->getBox(),
                                    d_boundary_fill_ghost_width));

                              if (!(fill_box.empty())) {
   
                                 d_refine_patch_strategy->
                                    fillSingularityBoundaryConditions(
                                    *patch, *d_encon_level,
                                    d_dst_to_encon,
                                    fill_time, fill_box, eboxes[bb]);
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
}


/*
 **************************************************************************
 *
 * Check whether the scratch data needs to be allocated on the specified
 * level.  Keep track of those allocated components so that they may be
 * deallocated later.
 *
 **************************************************************************
 */

void RefineSchedule::allocateScratchSpace(
   hier::ComponentSelector& allocate_vector,
   tbox::Pointer<hier::PatchLevel> level,
   double fill_time) const
{
   TBOX_ASSERT(!level.isNull());

   allocate_vector.clrAllFlags();

   hier::ComponentSelector preprocess_vector;

   for (int iri = 0; iri < d_number_refine_items; iri++) {
      const int scratch_id = d_refine_items[iri]->d_scratch;
      if (!level->checkAllocated(scratch_id)) {
         allocate_vector.setFlag(scratch_id);
      }
      preprocess_vector.setFlag(scratch_id);
   }

   level->allocatePatchData(allocate_vector, fill_time);

   if (!d_transaction_factory.isNull()) {
      d_transaction_factory->preprocessScratchSpace(level,
         fill_time,
         preprocess_vector);
   }
}

/*
 **************************************************************************
 *
 * Fill the component selector argument with the components needed to
 * allocate source data.
 *
 **************************************************************************
 */

void RefineSchedule::initializeSourceVector(
   hier::ComponentSelector& allocate_vector) const
{
   allocate_vector.clrAllFlags();

   for (int iri = 0; iri < d_number_refine_items; iri++) {
      allocate_vector.setFlag(d_refine_items[iri]->d_src);
   }
}

/*
 * ************************************************************************
 *                                                                        *
 * Fill the component selector argument with the components needed to     *
 * allocate destination data.                                             *
 *                                                                        *
 * ************************************************************************
 */

void RefineSchedule::initializeDestinationVector(
   hier::ComponentSelector& allocate_vector) const
{
   allocate_vector.clrAllFlags();

   for (int iri = 0; iri < d_number_refine_items; iri++) {
      allocate_vector.setFlag(d_refine_items[iri]->d_dst);
   }
}

/*
 **************************************************************************
 *
 * Allocate all destination data and store the destination components
 * in a component selector.
 *
 **************************************************************************
 */

void RefineSchedule::allocateDestinationSpace(
   hier::ComponentSelector& allocate_vector,
   double fill_time) const
{
   TBOX_ASSERT(!d_dst_level.isNull());

   allocate_vector.clrAllFlags();

   for (int iri = 0; iri < d_number_refine_items; iri++) {
      const int dst_id = d_refine_items[iri]->d_dst;
      if (!d_dst_level->checkAllocated(dst_id)) {
         allocate_vector.setFlag(dst_id);
      }
   }

   d_dst_level->allocatePatchData(allocate_vector, fill_time);
}

/*
 **************************************************************************
 *
 * Copy data from the scratch space of the specified level into the
 * destination space.  If the two spaces are the same, then no copy
 * is performed.
 *
 **************************************************************************
 */

void RefineSchedule::copyScratchToDestination() const
{
   TBOX_ASSERT(!d_dst_level.isNull());

   for (hier::PatchLevel::Iterator p(d_dst_level); p; p++) {
      tbox::Pointer<hier::Patch> patch(*p);

      for (int iri = 0; iri < d_number_refine_items; iri++) {
         const int src_id = d_refine_items[iri]->d_scratch;
         const int dst_id = d_refine_items[iri]->d_dst;
         if (src_id != dst_id) {
            TBOX_ASSERT(tbox::MathUtilities<double>::equalEps(patch->
                  getPatchData(dst_id)->getTime(),
                  patch->getPatchData(src_id)->getTime()));
            patch->getPatchData(dst_id)->copy(*patch->getPatchData(src_id));
         }
      }

   }

}

/*
 **************************************************************************
 *
 * Refine data from the coarse level into the fine level on the provided
 * fill box regions.  All operations are performed on the scratch space.
 *
 **************************************************************************
 */

void RefineSchedule::refineScratchData(
   const tbox::Pointer<hier::PatchLevel>& fine_level,
   const tbox::Pointer<hier::PatchLevel>& coarse_level,
   const Connector& coarse_to_fine,
   const Connector& coarse_to_unfilled,
   const tbox::List<tbox::Array<tbox::Pointer<hier::BoxOverlap> > >&
      overlaps) const
{
   t_refine_scratch_data->start();

   const hier::IntVector ratio(fine_level->getRatioToLevelZero()
                               / coarse_level->getRatioToLevelZero());

   tbox::ListIterator<tbox::Array<tbox::Pointer<hier::BoxOverlap> > >
      overlap_iter(overlaps);

   /*
    * Loop over all the coarse patches and find the corresponding
    * destination patch and destination fill boxes.
    */

   const hier::MappedBoxSet& crse_mapped_boxes =
      coarse_level->getMappedBoxLevel()->getMappedBoxes();
   for (hier::MappedBoxSet::const_iterator ni = crse_mapped_boxes.begin();
        ni != crse_mapped_boxes.end(); ++ni) {

      const hier::Box& crse_mapped_box = *ni;
      const hier::MappedBoxSet& dst_nabrs =
         coarse_to_fine.getNeighborSet(crse_mapped_box.getId());
      const hier::Box& dst_mapped_box = *dst_nabrs.begin();
#ifdef DEBUG_CHECK_ASSERTIONS
      /*
       * Each crse_mapped_box can point back to just one dst_mapped_box.
       * All other mapped_boxes in dst_nabrs must be a periodic image of
       * the same dst_mapped_box.
       */
      for (hier::MappedBoxSet::const_iterator na = dst_nabrs.begin();
           na != dst_nabrs.end();
           ++na) {
         TBOX_ASSERT(na->isPeriodicImage() || na == dst_nabrs.begin());
         TBOX_ASSERT(na->getGlobalId() == dst_mapped_box.getGlobalId());
      }
#endif
      tbox::Pointer<hier::Patch> fine_patch(fine_level->getPatch(
            dst_mapped_box.getGlobalId(), crse_mapped_box.getBlockId()));
      tbox::Pointer<hier::Patch> crse_patch(coarse_level->getPatch(
            crse_mapped_box.getGlobalId(), crse_mapped_box.getBlockId()));

      const NeighborSet& unfilled_nabrs = coarse_to_unfilled.getNeighborSet(
            crse_mapped_box.getId());
      TBOX_ASSERT(unfilled_nabrs.size() == 1);
      const hier::Box& unfilled_nabr = *unfilled_nabrs.begin();
      hier::BoxList fill_boxes(unfilled_nabr);

      if (d_refine_patch_strategy) {
         d_refine_patch_strategy->preprocessRefineBoxes(*fine_patch,
            *crse_patch,
            fill_boxes,
            ratio);
      }

      for (int iri = 0; iri < d_number_refine_items; iri++) {
         const xfer::RefineClasses::Data * const ref_item = d_refine_items[iri];
         if (!(ref_item->d_oprefine.isNull())) {

            tbox::Pointer<hier::BoxOverlap> refine_overlap(
               overlap_iter()[ref_item->d_class_index]);

            const int scratch_id = ref_item->d_scratch;

            ref_item->d_oprefine->refine(*fine_patch, *crse_patch,
               scratch_id, scratch_id,
               *refine_overlap, ratio);

         }
      }

      overlap_iter++;

      if (d_refine_patch_strategy) {
         d_refine_patch_strategy->postprocessRefineBoxes(*fine_patch,
            *crse_patch,
            fill_boxes,
            ratio);
      }

   }

   t_refine_scratch_data->stop();
}

/*
 **************************************************************************
 *
 * Compute the overlaps defining where the data will be refined into the
 * fine level.  All operations are performed on the scratch space.
 *
 **************************************************************************
 */
void RefineSchedule::computeRefineOverlaps(
   tbox::List<tbox::Array<tbox::Pointer<hier::BoxOverlap> > >& overlaps,
   const tbox::Pointer<hier::PatchLevel>& fine_level,
   const tbox::Pointer<hier::PatchLevel>& coarse_level,
   const Connector& coarse_to_fine,
   const Connector& coarse_to_unfilled)
{

   tbox::Pointer<hier::PatchDescriptor> fine_patch_descriptor =
      fine_level->getPatchDescriptor();

   const int num_equiv_classes =
      d_refine_classes->getNumberOfEquivalenceClasses();

   /*
    * Loop over all the coarse patches and find the corresponding
    * destination patch and destination fill boxes.
    */

   const hier::MappedBoxSet& coarse_mapped_boxes =
      coarse_level->getMappedBoxLevel()->getMappedBoxes();
   for (hier::MappedBoxSet::const_iterator ni = coarse_mapped_boxes.begin();
        ni != coarse_mapped_boxes.end(); ++ni) {

      const hier::Box& coarse_mapped_box = *ni;
      const hier::MappedBoxSet& fine_nabrs =
         coarse_to_fine.getNeighborSet(coarse_mapped_box.getId());
      const hier::Box& fine_mapped_box = *fine_nabrs.begin();
#ifdef DEBUG_CHECK_ASSERTIONS
      /*
       * Each coarse_mapped_box can point back to just one fine_mapped_box.
       * All other mapped_boxes in fine_nabrs must be a periodic image of
       * the same fine_mapped_box.
       */
      for (hier::MappedBoxSet::const_iterator na = fine_nabrs.begin();
           na != fine_nabrs.end();
           ++na) {
         TBOX_ASSERT(na->isPeriodicImage() || na == fine_nabrs.begin());
         TBOX_ASSERT(na->getGlobalId() == fine_mapped_box.getGlobalId());
      }
#endif
      tbox::Pointer<hier::Patch> fine_patch(fine_level->getPatch(
         fine_mapped_box.getId()));

      const NeighborSet& unfilled_nabrs = coarse_to_unfilled.getNeighborSet(
            coarse_mapped_box.getId());
      TBOX_ASSERT(unfilled_nabrs.size() == 1);
      const hier::Box& unfilled_nabr = *unfilled_nabrs.begin();
      hier::BoxList fill_boxes(unfilled_nabr);

      /*
       * The refine overlap will cover only  the fine fill box regions.
       * of index space.  Note that we restrict the interpolation range
       * to the intersection of the fine fill box and the ghost box of
       * the scratch data component (i.e., the destination of the
       * interpolation).  This is needed for the case where data
       * components treated by the schedule have different ghost
       * cell widths since the fill boxes are generated using the
       * maximum ghost cell width.
       */
      tbox::Array<tbox::Pointer<hier::BoxOverlap> >
         refine_overlaps(num_equiv_classes); 
      for (int ne = 0; ne < num_equiv_classes; ne++) {

         const xfer::RefineClasses::Data& rep_item =
            d_refine_classes->getClassRepresentative(ne);

         if (!(rep_item.d_oprefine.isNull())) {

            const int scratch_id = rep_item.d_scratch;
            tbox::Pointer<hier::PatchDataFactory> fine_pdf(
               fine_patch_descriptor->getPatchDataFactory(scratch_id));
            const hier::IntVector& fine_scratch_gcw =
               fine_pdf->getGhostCellWidth();
            hier::Box scratch_space(fine_patch->getBox());
            scratch_space.grow(fine_scratch_gcw);

            refine_overlaps[ne] =
               rep_item.d_var_fill_pattern->computeFillBoxesOverlap(
                  fill_boxes,
                  fine_patch->getBox(),
                  scratch_space,
                  *fine_pdf);

         }
      }
      overlaps.appendItem(refine_overlaps);
   }
}


/*
 *************************************************************************
 *
 * Generate communication schedule routine creates transactions to move
 * data from interiors of the source space on the source level into the
 * specified fill box regions of the destination level.
 *
 * The resulting transactions will only fill the regions of intersection
 * between the fill_boxes and destination level boxes.  The remaining
 * box regions are returned in unfilled_boxes.
 *
 *********************************************************************
 */

void RefineSchedule::generateCommunicationSchedule(
   tbox::Pointer<MappedBoxLevel>& unfilled_mapped_box_level,
   tbox::Pointer<Connector>& dst_to_unfilled,
   tbox::Pointer<MappedBoxLevel>& unfilled_encon_box_level,
   tbox::Pointer<Connector>& encon_to_unfilled_encon,
   const Connector& dst_to_src,
   const Connector& src_to_dst,
   const Connector& dst_to_fill,
   const FillSet& dst_to_fill_on_src_proc,
   const bool use_time_interpolation,
   const bool create_transactions)
{
   t_gen_comm_sched->start();

   tbox::Pointer<hier::GridGeometry> grid_geometry(
      d_dst_level->getGridGeometry());

   if (s_extra_debug) {
      if (dst_to_src.isInitialized()) {
         dst_to_src.assertTransposeCorrectness(src_to_dst);
         src_to_dst.assertTransposeCorrectness(dst_to_src);
      }
   }

   const MappedBoxLevel& dst_mapped_box_level = dst_to_fill.getBase();

   if (create_transactions) {
      /*
       * We invert the edges so that the src transactions on the local
       * process are created in the same order the dst transactions
       * are created on remote processes.  The ordering is how the
       * source and destination transactions are matched up on the
       * receiving process.  We chose to order the transactions
       * dst-major (all the transactions for one dst Box are
       * grouped together, then all transactions for the next
       * dst Box, and so on).
       */
      const hier::NeighborhoodSet& dst_to_src_edges =
         dst_to_src.getNeighborhoodSets();

      /*
       * Reorder the src_to_dst edge data to arrange neighbors by the
       * dst mapped_boxes, as required to match the transaction ordering
       * on the receiving processors.
       */
      FullNeighborhoodSet src_to_dst_edges_bydst;
      t_invert_edges->start();
      reorderNeighborhoodSetsByDstNodes(src_to_dst_edges_bydst, src_to_dst);
      t_invert_edges->stop();


      t_construct_send_trans->start();

      /*
       * Construct transactions with local source and remote destination.
       */
      for (FullNeighborhoodSet::const_iterator
           ei = src_to_dst_edges_bydst.begin();
           ei != src_to_dst_edges_bydst.end(); ++ei) {

         /*
          * dst_mapped_box can be remote (by definition of FullNeighborhoodSet).
          * local_src_mapped_boxes are the local source mapped_boxes that
          * contribute data to dst_mapped_box.
          */
         const hier::Box& dst_mapped_box = ei->first;
         const hier::MappedBoxSet& local_src_mapped_boxes = ei->second;
         TBOX_ASSERT(!dst_mapped_box.isPeriodicImage());

         FillSet::const_iterator dst_fill_iter =
            dst_to_fill_on_src_proc.find(dst_mapped_box.getId());
         if (dst_fill_iter == dst_to_fill_on_src_proc.end()) {
            /*
             * Missing fill boxes should indicate that the dst mapped_box
             * has no fill box.  One way this is possible is for
             * d_dst_level_fill_pattern to be of type PatchLevelBorderFillPattern
             * and for dst_mapped_box to be away from level borders.
             */
            continue;
         }

         const hier::MappedBoxSet& dst_fill_boxes = dst_fill_iter->second;

         for (hier::MappedBoxSet::const_iterator
              ni = local_src_mapped_boxes.begin();
              ni != local_src_mapped_boxes.end(); ++ni) {
            const hier::Box& src_mapped_box = *ni;

            if (src_mapped_box.getOwnerRank() !=
                dst_mapped_box.getOwnerRank()) {

               constructScheduleTransactions(
                  dst_fill_boxes,
                  dst_mapped_box,
                  src_mapped_box,
                  use_time_interpolation);

            }
         }

      } // end send transactions loop

      t_construct_send_trans->stop();

   } // if create_transactions

   const hier::NeighborhoodSet& dst_to_fill_edges =
      dst_to_fill.getNeighborhoodSets();

   hier::MappedBoxSet level_unfilled_boxes;
   hier::MappedBoxSet level_encon_unfilled_boxes;
   hier::NeighborhoodSet dst_to_unfilled_nbrhood_set;
   hier::NeighborhoodSet encon_to_unfilled_encon_nbrhood_set;

   hier::LocalId last_unfilled_local_id(-1);

   t_construct_recv_trans->start();
   for (hier::NeighborhoodSet::const_iterator cf = dst_to_fill_edges.begin();
        cf != dst_to_fill_edges.end(); ++cf) {

      const hier::BoxId& dst_mapped_box_id(cf->first);
      const hier::Box& dst_mapped_box = *dst_mapped_box_level.getMappedBox(
            dst_mapped_box_id);
      const hier::BlockId& dst_block_id = dst_mapped_box_id.getBlockId();

      const tbox::Dimension& dim = dst_mapped_box.getDim();

      const hier::MappedBoxSet& fill_nabrs = cf->second;

      hier::BoxList fill_boxes_list(dim);
      for (hier::MappedBoxSet::iterator bi = fill_nabrs.begin();
           bi != fill_nabrs.end(); ++bi) {
         fill_boxes_list.appendItem(*bi);
      }

      /*
       * If the destination box touches enhanced connectivity,
       * encon_fill_boxes will hold the portion of fill_boxes_list that lies
       * in the enhanced connectivity ghost region of the destination.
       * Otherwise, encon_fill_boxes will be empty.
       */
      hier::BoxList encon_fill_boxes(dim);
      if (grid_geometry->hasEnhancedConnectivity()) {
         findEnconFillBoxes(encon_fill_boxes,
                            fill_boxes_list,
                            dst_block_id);
      }

      /*
       * unfilled_boxes_for_dst starts out containing the fill boxes
       * for the current dst_mapped_box.  As transactions are created,
       * the space covered by those transactions will be removed from this
       * list, and whatever is left cannot be filled from the source level.
       *
       * The boxes in encon_fill_boxes, if any, are removed from
       * unfilled_boxes_for_dst, because unfilled boxes at enhanced
       * connectivity are handled separately at a later step.
       */
      hier::BoxList unfilled_boxes_for_dst(fill_boxes_list);
      unfilled_boxes_for_dst.removeIntersections(encon_fill_boxes);

      if (create_transactions) {

         const hier::NeighborhoodSet& dst_to_src_edges =
            dst_to_src.getNeighborhoodSets();

         hier::NeighborhoodSet::const_iterator dst_to_src_iter =
            dst_to_src_edges.find(dst_mapped_box_id);

         if (dst_to_src_iter != dst_to_src_edges.end()) {
            const hier::MappedBoxSet& src_mapped_boxes =
               dst_to_src_iter->second; 

            for (hier::MappedBoxSet::const_iterator
                 na = src_mapped_boxes.begin();
                 na != src_mapped_boxes.end(); ++na) {

               const hier::Box& src_mapped_box = *na;
               const hier::BlockId& src_block_id = src_mapped_box.getBlockId();

               /*
                * Remove the source box from the unfilled boxes list.  If
                * necessary, the source box is transformed to the destination
                * coordinate system.
                *
                * The removal of the source box is skipped if the source
                * and destination are enhanced connectivity (singularity)
                * neighbors, since the handling of unfilled boxes in this case
                * is a separate step.
                */
               if (src_block_id != dst_block_id) {

                  if (!grid_geometry->areSingularityNeighbors(dst_block_id,
                                                              src_block_id)) {

                     hier::Box transformed_src_box(src_mapped_box);

                     grid_geometry->transformBox(
                        transformed_src_box,
                        d_dst_level->getRatioToLevelZero(),
                        dst_block_id,
                        src_block_id);

                     unfilled_boxes_for_dst.removeIntersections(
                        transformed_src_box);
 
                  }
               } else {
                  unfilled_boxes_for_dst.removeIntersections(
                     src_mapped_box);
               }
               constructScheduleTransactions(
                  fill_nabrs,
                  dst_mapped_box,
                  src_mapped_box,
                  use_time_interpolation);

            }
         }
      }

      /*
       * Add mapping information to unfilled boxes and add them to
       * containers for the level.
       */
      if (unfilled_boxes_for_dst.size() > 0) {
         NeighborSet& unfilled_nabrs =
            dst_to_unfilled_nbrhood_set[dst_mapped_box_id];

         for (hier::BoxList::Iterator bi(unfilled_boxes_for_dst);
              bi; bi++) {

            hier::Box unfilled_mapped_box(bi(),
                                          ++last_unfilled_local_id,
                                          dst_mapped_box.getOwnerRank(),
                                          dst_block_id);

            level_unfilled_boxes.insert(level_unfilled_boxes.end(),
                                        unfilled_mapped_box);
            unfilled_nabrs.insert(unfilled_nabrs.end(),
                                  unfilled_mapped_box);

         }
      }

      /*
       * Call private method to handle unfilled boxes where the destination
       * box touches enhanced connectivity.
       */
      if (encon_fill_boxes.size() > 0) {
         findEnconUnfilledBoxes(
            level_encon_unfilled_boxes,
            encon_to_unfilled_encon_nbrhood_set,
            last_unfilled_local_id,
            dst_mapped_box,
            dst_to_src,
            encon_fill_boxes);
      } 

   } // End receive/copy transactions loop
   t_construct_recv_trans->stop();

   unfilled_mapped_box_level = new MappedBoxLevel(
      level_unfilled_boxes,
      dst_mapped_box_level.getRefinementRatio(),
      dst_mapped_box_level.getGridGeometry(),
      dst_mapped_box_level.getMPI());

   dst_to_unfilled = new Connector(
      dst_mapped_box_level,
      *unfilled_mapped_box_level,
      dst_to_fill.getConnectorWidth(),
      dst_to_unfilled_nbrhood_set,
      MappedBoxLevel::DISTRIBUTED);

   unfilled_encon_box_level = new MappedBoxLevel(
      level_encon_unfilled_boxes,
      dst_mapped_box_level.getRefinementRatio(),
      dst_mapped_box_level.getGridGeometry(),
      dst_mapped_box_level.getMPI());

   encon_to_unfilled_encon = new Connector(
      *(d_encon_level->getMappedBoxLevel()),
      *unfilled_encon_box_level,
      dst_to_fill.getConnectorWidth(),
      encon_to_unfilled_encon_nbrhood_set,
      MappedBoxLevel::DISTRIBUTED);

   t_gen_comm_sched->stop();
}

/*
 ***********************************************************************
 * Given a BlockId and a list of fill boxes, populate encon_fill_boxes
 * with the portion of the fill boxes that lie across an enhanced
 * connectivity boundary from the destination block.
 ***********************************************************************
 */
void RefineSchedule::findEnconFillBoxes(
   hier::BoxList& encon_fill_boxes,
   const hier::BoxList& fill_boxes_list,
   const hier::BlockId& dst_block_id)
{
   TBOX_ASSERT(encon_fill_boxes.size() == 0);

   tbox::Pointer<hier::GridGeometry> grid_geometry(
      d_dst_level->getGridGeometry());

   const int dst_blk = dst_block_id.getBlockValue();
 
   const tbox::List<hier::GridGeometry::Neighbor>& neighbors =
      grid_geometry->getNeighbors(dst_blk);

   hier::BoxList encon_neighbor_list(encon_fill_boxes.getDim());
   for (tbox::List<hier::GridGeometry::Neighbor>::Iterator
        ni(neighbors); ni; ni++) {

      if(ni().isSingularity()) {
         encon_fill_boxes.unionBoxes(ni().getTransformedDomain());
      }

   }

   encon_fill_boxes.refine(d_dst_level->getRatioToLevelZero());

   encon_fill_boxes.intersectBoxes(fill_boxes_list);

   encon_fill_boxes.coalesceBoxes();
}

/*
 ***********************************************************************
 * Given a list of fill boxes for a destination box at enhanced
 * connectivity, determine which fill boxes cannot be filled from the
 * source level and add those unfilled boxes to the appropriate containers.
 ***********************************************************************
 */
void RefineSchedule::findEnconUnfilledBoxes(
   hier::MappedBoxSet& level_encon_unfilled_boxes,
   hier::NeighborhoodSet& encon_to_unfilled_encon_nbrhood_set,
   hier::LocalId& last_unfilled_local_id,
   const hier::Box& dst_mapped_box,
   const Connector& dst_to_src,
   const hier::BoxList& encon_fill_boxes)
{
   tbox::Pointer<hier::GridGeometry> grid_geometry(
      d_dst_level->getGridGeometry());

   const hier::BoxId& dst_mapped_box_id = dst_mapped_box.getId();
   const hier::BlockId& dst_block_id = dst_mapped_box_id.getBlockId();
   const int dst_blk = dst_block_id.getBlockValue();

   /*
    * map container will hold unfilled boxes for each block that is
    * an enhanced connectivity neighbor of the destination box.
    *
    * To start, each entry in the map will have the intersection of
    * encon_fill_boxes and the neighboring block's domain, represented in
    * the destination coordinate system.
    */
   std::map<hier::BlockId, hier::BoxList> unfilled_encon_nbr_boxes;

   const tbox::List<hier::GridGeometry::Neighbor>& neighbors =
      grid_geometry->getNeighbors(dst_blk);

   for (tbox::List<hier::GridGeometry::Neighbor>::Iterator
        ni(neighbors); ni; ni++) {

      if(ni().isSingularity()) {
         const hier::BlockId nbr_block_id(ni().getBlockNumber());

         unfilled_encon_nbr_boxes[nbr_block_id].unionBoxes(
            ni().getTransformedDomain());
         unfilled_encon_nbr_boxes[nbr_block_id].refine(
            d_dst_level->getRatioToLevelZero());
         unfilled_encon_nbr_boxes[nbr_block_id].intersectBoxes(
            encon_fill_boxes);
      }
   }

   const hier::NeighborhoodSet& dst_to_src_edges =
      dst_to_src.getNeighborhoodSets();

   /*
    * If there are overlapping source boxes found in dst_to_src_edges,
    * and if those source boxes lie across enhanced connectivity from
    * the destination box, then we remove the source box from the
    * source block's entry in the unfilled_encon_nbr_boxes map container.
    */
   hier::NeighborhoodSet::const_iterator dst_to_src_iter =
      dst_to_src_edges.find(dst_mapped_box_id);

   if (dst_to_src_iter != dst_to_src_edges.end()) {
      const hier::MappedBoxSet& src_mapped_boxes = dst_to_src_iter->second;

      /*
       * If at enhanced connectivity, remove source box from container of
       * unfilled boxes
       */
      for (hier::MappedBoxSet::const_iterator na = src_mapped_boxes.begin();
           na != src_mapped_boxes.end(); ++na) {

         const hier::Box& src_mapped_box = *na;
         const hier::BlockId& src_block_id = src_mapped_box.getBlockId();

         if (src_block_id != dst_block_id) {

            if (grid_geometry->areSingularityNeighbors(dst_block_id,
                                                       src_block_id)) {

               hier::Box transformed_src_box(src_mapped_box);
               grid_geometry->transformBox(transformed_src_box,
                                           d_dst_level->getRatioToLevelZero(),
                                           dst_block_id,
                                           src_block_id); 

               unfilled_encon_nbr_boxes[src_block_id].removeIntersections(
                  transformed_src_box);

            }
         }
      }
   }

   /*
    * If any unfilled boxes remain for any block in unfilled_encon_nbr_boxes,
    * they need to be added to the output containers.
    */
   if (unfilled_encon_nbr_boxes.size() > 0) {

      const hier::NeighborhoodSet& dst_to_encon_nbrhood_set =
         d_dst_to_encon.getNeighborhoodSets();

      hier::NeighborhoodSet::const_iterator find_encon_nabrs =
         dst_to_encon_nbrhood_set.find(dst_mapped_box_id);

      if (find_encon_nabrs != dst_to_encon_nbrhood_set.end()) {

         const NeighborSet& dst_encon_nabrs = find_encon_nabrs->second;

         for (NeighborSet::const_iterator de_iter = dst_encon_nabrs.begin();
              de_iter != dst_encon_nabrs.end(); ++de_iter) {

            const hier::BoxId& encon_mapped_box_id =
               de_iter->getId();
            const hier::BlockId& nbr_block_id =
               encon_mapped_box_id.getBlockId();

            NeighborSet& unfilled_nabrs =
               encon_to_unfilled_encon_nbrhood_set[encon_mapped_box_id];

            const hier::BoxList& unfilled_boxes =
               unfilled_encon_nbr_boxes[nbr_block_id];

            if (unfilled_boxes.size() > 0) {

               /*
                * The unfilled boxes are, at this point, represented in
                * the destination coordinates.  Here they are transformed
                * into the neighboring block's coordinates before being
                * added to the output containers. 
                */

               for (hier::BoxList::Iterator bi(unfilled_boxes);
                    bi; bi++) {

                  hier::Box unfilled_box(bi());
                  grid_geometry->transformBox(
                     unfilled_box,
                     d_dst_level->getRatioToLevelZero(),
                     nbr_block_id,
                     dst_block_id); 
 
                  hier::Box unfilled_encon_box(unfilled_box,
                                               ++last_unfilled_local_id,
                                               dst_mapped_box.getOwnerRank(),
                                               nbr_block_id);

                  level_encon_unfilled_boxes.insert(
                     level_encon_unfilled_boxes.end(),
                     unfilled_encon_box);

                  unfilled_nabrs.insert(unfilled_nabrs.end(),
                                        unfilled_encon_box);

               }
            }
         }
      }
   }
}

/*
 ***********************************************************************
 * This method does 2 important things to the src_to_dst edges:
 *
 * 1. It puts the edge data in dst-major order so the src owners can
 * easily loop through the dst-src edges in the same order that dst
 * owners see them.  Transactions must have the same order on the
 * sending and receiving processors.
 *
 * 2. It shifts periodic image dst Boxes back to the zero-shift position,
 * and applies a similar shift to src Boxes so that the overlap is
 * unchanged.  The constructScheduleTransactions method requires all
 * shifts to be absorbed in the src Box.
 ***********************************************************************
 */
void RefineSchedule::reorderNeighborhoodSetsByDstNodes(
   FullNeighborhoodSet& full_inverted_edges,
   const Connector& src_to_dst) const
{
   const tbox::Dimension& dim(d_dst_level->getDim());

   const hier::PeriodicShiftCatalog* shift_catalog =
      hier::PeriodicShiftCatalog::getCatalog(dim);
   const hier::NeighborhoodSet& edges = src_to_dst.getNeighborhoodSets();
   const MappedBoxLevel& src_mapped_box_level = src_to_dst.getBase();
   const hier::IntVector& src_ratio = src_to_dst.getBase().getRefinementRatio();
   const hier::IntVector& dst_ratio = src_to_dst.getHead().getRefinementRatio();

   /*
    * These are the counterparts to shifted dst mapped_boxes and unshifted src
    * mapped_boxes.
    */
   hier::Box shifted_mapped_box(dim), unshifted_nabr(dim);
   full_inverted_edges.clear();
   for (hier::NeighborhoodSet::const_iterator ci = edges.begin();
        ci != edges.end();
        ++ci) {
      const hier::Box& mapped_box =
         *src_mapped_box_level.getMappedBoxStrict(ci->first);
      const NeighborSet& nabrs = ci->second;
      for (NeighborSet::const_iterator na = nabrs.begin();
           na != nabrs.end(); ++na) {
         const hier::Box& nabr = *na;
         if (nabr.isPeriodicImage()) {
            shifted_mapped_box.initialize(
               mapped_box,
               shift_catalog->getOppositeShiftNumber(nabr.getPeriodicId()),
               src_ratio);
            unshifted_nabr.initialize(
               nabr,
               shift_catalog->getZeroShiftNumber(),
               dst_ratio);
            full_inverted_edges[unshifted_nabr].insert(shifted_mapped_box);
         } else {
            full_inverted_edges[nabr].insert(mapped_box);
         }
      }
   }
}

/*
 ****************************************************************
 * Set the fill boxes and the dst--->fill Connector.
 * Initialize this same information on the src processes as well.
 *
 * Compute the fill boxes of the dst level according to the
 * d_dst_level_fill_pattern.  Initialize the fill MappedBoxLevel from
 * the those fill boxes.
 *
 * Set the edges between dst and fill mapped_box_levels (trivial,
 * since fill MappedBoxLevel is related to dst MappedBoxLevel).  Set
 * edges between src and fill mapped_box_levels (by copying the edges
 * between src and dst.
 *
 * If d_dst_level_fill_pattern can compute the relationship between
 * dst boxes and fill boxes on the src processes, have it do that.
 * Otherwise, perform communication to get that information onto the
 * src processes.
 *
 * The Connectors dst_to_dst and dst_to_src (in the arguments)
 * are required only for fill_patterns PatchLevelBorderFillPattern and
 * and PatchLevelBorderAndInteriorFillPattern
 ****************************************************************
 */

void RefineSchedule::setDefaultFillMappedBoxLevel(
   MappedBoxLevel& fill_mapped_box_level,
   Connector& dst_to_fill,
   FillSet& dst_to_fill_on_src_proc,
   const MappedBoxLevel& dst_mapped_box_level,
   const hier::Connector* dst_to_src,
   const hier::Connector* src_to_dst,
   const hier::IntVector& fill_ghost_width)
{
   /*
    * If the destination variable of any registered communication
    * requires that coarse data take precedence over fine data at
    * coarse-fine boundaries, make sure that the ghost cell width is
    * at least one in each direction.  This is necessary so that
    * generateCommunicationSchedule() determines that boundary data for
    * the level needs to be transferred from the next coarser level
    * (i.e. the interiors of the fill boxes overlap the next coarser
    * level in a nontrivial way).
    */

   const tbox::Dimension& dim(d_dst_level->getDim());

   TBOX_DIM_ASSERT_CHECK_DIM_ARGS2(dim, dst_mapped_box_level, fill_ghost_width);

   const hier::IntVector& constant_one_intvector(hier::IntVector::getOne(dim));

   bool need_nontrivial_ghosts = false;
   for (int iri = 0; iri < d_number_refine_items; iri++) {
      if (!(d_refine_items[iri]->d_fine_bdry_reps_var)) {
         need_nontrivial_ghosts = true;
      }
   }

   hier::IntVector fill_gcw = fill_ghost_width;
   if (need_nontrivial_ghosts) {
      fill_gcw.max(constant_one_intvector);
   }

   const Connector dummy_connector;
   const Connector* dst_to_dst =
      dst_mapped_box_level.getPersistentOverlapConnectors().hasConnector(
         dst_mapped_box_level,
         fill_gcw) ?
      &dst_mapped_box_level.getPersistentOverlapConnectors().findConnector(
         dst_mapped_box_level,
         fill_gcw) : &dummy_connector;

   // New data computed here:
   hier::MappedBoxSet fill_mapped_boxes;
   hier::NeighborhoodSet dst_eto_fill;

   /*
    * d_max_fill_boxes is the max number of fill boxes
    * for any dst mapped_box.  d_max_fill_boxes is nominally 1.
    * If more is required (based on fill_pattern),
    * it is incremented below.
    */
   d_max_fill_boxes = 1;

   /*
    * Generate fill mapped_boxes (grown dst mapped_boxes) and
    * get edges between dst and fill mapped_box_levels.
    * These are all local edges by definition.
    * Note that since fill mapped_boxes are simply grown
    * versions of dst_mapped_box (at most), edges between
    * fill and src mapped_box_levels would have gcw of zero.
    *
    * Technically speaking, edges between dst and
    * fill are not complete.  Once the dst box is
    * grown to make the fill box, the fill box may
    * intersect other boxes in the dst mapped_box_level.  We
    * ignore all these overlaps because they are not
    * needed by the algorithm.
    */

   d_dst_level_fill_pattern->computeFillMappedBoxesAndNeighborhoodSets(
      fill_mapped_boxes,
      dst_eto_fill,
      dst_mapped_box_level,
      *dst_to_dst,
      dst_to_src == NULL ? dummy_connector : *dst_to_src,
      src_to_dst == NULL ? dummy_connector : *src_to_dst,
      fill_gcw);

   d_max_fill_boxes = tbox::MathUtilities<int>::Max(
         d_max_fill_boxes,
         d_dst_level_fill_pattern->getMaxFillBoxes());

   fill_mapped_box_level.swapInitialize(
      fill_mapped_boxes,
      dst_mapped_box_level.getRefinementRatio(),
      dst_mapped_box_level.getGridGeometry(),
      dst_mapped_box_level.getMPI());

   /*
    * Note that dst_to_fill is NOT complete.
    * We care only about the intersection between a dst_mapped_box
    * and the fill_mapped_box it created.  In fact, a dst_mapped_box
    * may also intersect the fill_mapped_box of other nearby
    * destination mapped_boxes.
    */
   dst_to_fill.swapInitialize(
      dst_mapped_box_level,
      fill_mapped_box_level,
      fill_gcw,
      dst_eto_fill,
      MappedBoxLevel::DISTRIBUTED);

   if (!d_src_level.isNull()) {
      if (d_dst_level_fill_pattern->needsToCommunicateDestinationFillBoxes()) {

         /*
          * This part assumes src-dst ratio is one.
          * Should be modified if the assumption does not hold.
          */
         TBOX_ASSERT(dst_to_src->getRatio() == hier::IntVector::getOne(dim));

         /*
          * For these fill_pattern, the src owner could not compute fill boxes
          * for all its dst neighbors using local data, so the dst owners
          * must send this information.
          */
         communicateFillBoxes(dst_to_fill_on_src_proc,
            dst_to_fill,
            *dst_to_src,
            *src_to_dst);

      } else {
         d_dst_level_fill_pattern->computeDestinationFillBoxesOnSourceProc(
            dst_to_fill_on_src_proc,
            dst_mapped_box_level,
            *src_to_dst,
            fill_gcw);
      }
   }

   createEnconLevel(fill_gcw);
}

/*
 ****************************************************************
 * Setup the level representing ghost regions of destination patches
 * that touch enhanced connectivity.
 *****************************************************************
 */
void RefineSchedule::createEnconLevel(const hier::IntVector& fill_gcw)
{
   const tbox::Dimension& dim = fill_gcw.getDim();

   tbox::Pointer<hier::GridGeometry> grid_geometry(
      d_dst_level->getGridGeometry()); 
   const int num_blocks = grid_geometry->getNumberBlocks();

   hier::MappedBoxSet encon_mapped_boxes;
   hier::NeighborhoodSet dst_to_encon_nabrs;

   if (num_blocks > 1) {
      hier::IntVector encon_gcw(
         hier::IntVector::max(fill_gcw, hier::IntVector::getOne(dim)));

      hier::LocalId encon_local_id(0);

      /*
       * Loop over blocks to find destination patches on each block
       * that touch enhanced connectivity.
       */
      for (int bn = 0; bn < num_blocks; bn++) {

         hier::BlockId block_id(bn);

         /*
          * Test to see if there are any local destination boxes on this
          * block.  Move on to next block if not.
          */ 
         hier::MappedBoxSetSingleBlockIterator dst_test_iter(
                     d_dst_level->getMappedBoxLevel()->getMappedBoxes(), 
                     block_id);

         if (!dst_test_iter.isValid()) {
            continue;
         } 

         /*
          * An empty singularity box list would mean the block touches no
          * enhanced connectivy boundaries, so no need to go further with
          * this block.
          */ 
         const hier::BoxList& sing_boxes =
            grid_geometry->getSingularityBoxList(bn);

         if (sing_boxes.size() > 0) {
            const tbox::List<hier::GridGeometry::Neighbor>& neighbors =
               grid_geometry->getNeighbors(bn);

            /*
             * Loop over neighboring blocks and find the ones that are
             * singularity neighbors.
             */
            for (tbox::List<hier::GridGeometry::Neighbor>::
                 Iterator ni(neighbors); ni; ni++) {

               if (ni().isSingularity()) {

                  const int nbr_num = ni().getBlockNumber();
                  hier::BlockId nbr_id(nbr_num);

                  /*
                   * Get the transformation from neighbor block to dst
                   * block, and get a representation of the neighbor block
                   * domain in coordinate system of dst block.
                   */
                  hier::Transformation::RotationIdentifier rotation =
                     ni().getRotationIdentifier();
                  hier::IntVector offset(ni().getShift());
                  offset *= (d_dst_level->getRatioToLevelZero());

                  hier::Transformation transformation(rotation, offset);

                  hier::BoxList trans_neighbor_list(dim);
                  grid_geometry->getTransformedBlock(trans_neighbor_list,
                     bn,
                     nbr_num);
                  trans_neighbor_list.refine(
                     d_dst_level->getRatioToLevelZero());

                  /*
                   * Loop over dst mapped boxes for this block and
                   * determine if they touch the current neighbor block
                   * at enhanced connectivity.
                   */
                  hier::MappedBoxSetSingleBlockIterator dst_local_iter(
                     d_dst_level->getMappedBoxLevel()->getMappedBoxes(),
                     block_id);

                  for ( ; dst_local_iter.isValid(); dst_local_iter++) {

                     const hier::BoxId& mapped_box_id =
                        dst_local_iter->getId();
   
                     tbox::Pointer<hier::Patch> patch(
                        d_dst_level->getPatch(mapped_box_id));
                     tbox::Pointer<hier::PatchGeometry> pgeom(
                        patch->getPatchGeometry());

                     if (pgeom->getTouchesRegularBoundary()) {

                        /*
                         * Grow the patch box and intersect with the
                         * representation of the neighboring block.
                         * A non-empty intersection indicates there is
                         * a ghost region across an enhanced connectivity
                         * boundary.
                         */
                        hier::BoxList encon_test_list(patch->getBox());
                        encon_test_list.grow(encon_gcw);
                        encon_test_list.intersectBoxes(trans_neighbor_list);

                        if (encon_test_list.size() > 0) {

                           encon_test_list.coalesceBoxes();
                           TBOX_ASSERT(encon_test_list.size() == 1);

                           /*
                            * Transform the boxes representing the ghost
                            * region back to the neighbor block's
                            * coordinate system, and create a Box
                            * to be added to encon_mapped_boxes. 
                            */ 

                           for (hier::BoxList::Iterator bi(encon_test_list);
                                bi; bi++) {
                              hier::Box encon_box(bi());

                              transformation.inverseTransform(encon_box);

                              /*
                               * If a Box at this location on the
                               * same neighbor block and on the same processor
                               * already exists in encon_mapped_boxes,
                               * do not create another. 
                               */ 
                              hier::Box encon_mapped_box(dim);
                              bool box_exists = false;
                              for (hier::MappedBoxSetSingleBlockIterator
                                   encon_iter(encon_mapped_boxes, nbr_id);
                                   encon_iter.isValid(); ++encon_iter) {

                                 if (encon_box.isSpatiallyEqual(*encon_iter)) {
                                    box_exists = true;
                                    encon_mapped_box = *encon_iter;
                                    break;
                                 }
                              }

                              if (!box_exists) { 
                                 encon_mapped_box =
                                    *encon_mapped_boxes.insert(
                                       encon_mapped_boxes.end(),
                                       hier::Box(encon_box,
                                          ++encon_local_id,
                                          mapped_box_id.getOwnerRank(),
                                          nbr_id));
                              }

                              /*
                               * Add to the neighborhood set for the
                               * d_dst_to_encon connector.
                               */ 
                              dst_to_encon_nabrs[mapped_box_id].insert(
                                 encon_mapped_box);
                           
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }

   /*
    * Create d_encon_level and associated Connectors.
    *
    * Where destination patches have ghost regions across enhanced connectivity
    * boundaries, data communicated by this schedule will not be written
    * directly into those ghost regions, but rather into patches on
    * d_encon_level.  This level, once filled with data, will be provided
    * to RefinePatchStrategy's fillSingularityBoundaryConditions. 
    */
   hier::MappedBoxLevel encon_box_level(encon_mapped_boxes,
                                       d_dst_level->getRatioToLevelZero(),
                                       grid_geometry);

   d_encon_level = new hier::PatchLevel(encon_box_level,
                                        grid_geometry,
                                        d_dst_level->getPatchDescriptor(),
                                        tbox::Pointer<hier::PatchFactory>(NULL),
                                        true);

   d_encon_level->setLevelNumber(d_dst_level->getLevelNumber());

   d_dst_to_encon.initialize(*(d_dst_level->getMappedBoxLevel()),
                             *(d_encon_level->getMappedBoxLevel()),
                             hier::IntVector::getOne(dim),
                             dst_to_encon_nabrs);

   if (!d_src_level.isNull()) {
      d_src_to_encon.initialize(*(d_src_level->getMappedBoxLevel()),
                                *(d_encon_level->getMappedBoxLevel()),
                                hier::IntVector::getZero(dim));
      d_encon_to_src.initialize(*(d_encon_level->getMappedBoxLevel()),
                                *(d_src_level->getMappedBoxLevel()),
                                hier::IntVector::getZero(dim));

      hier::OverlapConnectorAlgorithm oca;

      /*
       * TODO:  Replace with a bridge.
       */
      oca.findOverlaps(d_src_to_encon);
      oca.findOverlaps(d_encon_to_src);
   }
}


/*
 **************************************************************************
 * Calculate the maximum ghost cell width of all destination patch data
 * components.
 *
 * It is possible for dst_to_fill_on_src_proc to omit a visible dst mapped_box
 * if the dst mapped_box has no fill boxes.
 **************************************************************************
 */
void RefineSchedule::communicateFillBoxes(
   FillSet& dst_to_fill_on_src_proc,
   const Connector& dst_to_fill,
   const Connector& dst_to_src,
   const Connector& src_to_dst)
{
   const tbox::Dimension& dim(d_dst_level->getDim());

   const hier::NeighborhoodSet& dst_eto_src = dst_to_src.getNeighborhoodSets();

   std::set<int> src_owners;
   dst_to_src.getNeighborhoodSets().getOwners(src_owners);
   src_owners.erase(dst_to_fill.getBase().getRank());

   std::set<int> dst_owners;
   src_to_dst.getNeighborhoodSets().getOwners(dst_owners);
   dst_owners.erase(dst_to_fill.getBase().getRank());

   std::map<int, std::vector<int> > send_mesgs;
   std::map<int, std::vector<int> > recv_mesgs;

   tbox::AsyncCommStage comm_stage;
   tbox::AsyncCommStage::MemberVec completed;
   tbox::AsyncCommPeer<int>* comms =
      new tbox::AsyncCommPeer<int>[src_owners.size() + dst_owners.size()];

   int mesg_number = 0;

   /*
    * Post receives for fill boxes from dst owners.
    */

   // TODO: should use cyclic ordering for efficiency.
   for (std::set<int>::const_iterator di = dst_owners.begin();
        di != dst_owners.end(); ++di) {
      comms[mesg_number].initialize(&comm_stage);
      comms[mesg_number].setPeerRank(*di);
      // Reuse communicator.  Assuming it has no outstanding messages!
      comms[mesg_number].setMPI(src_to_dst.getMPI());
      comms[mesg_number].setMPITag(0, 1);
      comms[mesg_number].beginRecv();
      if (comms[mesg_number].isDone()) {
         completed.insert(completed.end(), &comms[mesg_number]);
      }
      ++mesg_number;
   }

   /*
    * Pack fill boxes and send messages to src owners.
    * TODO: should use cyclic ordering for efficiency.
    */
   // Pack messages.
   std::vector<int> tmp_mesg;
   hier::MappedBoxSet tmp_fill_boxes;
   const hier::NeighborhoodSet& dst_eto_fill = dst_to_fill.getNeighborhoodSets();
   for (hier::NeighborhoodSet::const_iterator ei = dst_eto_fill.begin();
        ei != dst_eto_fill.end(); ++ei) {
      const hier::BoxId& dst_mapped_box_id = ei->first;
      const NeighborSet& fill_nabrs = ei->second;
      /*
       * Pack dst_mapped_box_id's fill box info into tmp_mesg.
       * - dst_mapped_box_id's LocalId
       * - number of fill neighbors
       * - fill neighbors (could just send box and save 2 ints)
       * Also, create BoxVector object for local use.
       */
      tmp_mesg.clear();
      tmp_mesg.reserve(3 + fill_nabrs.size() * hier::Box::commBufferSize(dim));
      tmp_mesg.insert(tmp_mesg.end(), 3, 0);
      tmp_mesg[0] = dst_mapped_box_id.getLocalId().getValue();
      tmp_mesg[1] = dst_mapped_box_id.getBlockId().getBlockValue();
      tmp_mesg[2] = static_cast<int>(fill_nabrs.size());
      tmp_fill_boxes.clear();
      for (NeighborSet::const_iterator na = fill_nabrs.begin();
           na != fill_nabrs.end(); ++na) {
         tmp_mesg.insert(tmp_mesg.end(), hier::Box::commBufferSize(dim), 0);
         na->putToIntBuffer(&tmp_mesg[tmp_mesg.size()
                                      - hier::Box::commBufferSize(dim)]);
         tmp_fill_boxes.insert(*na);
      }
      // Append tmp_mesg to buffers for sending to src owners.
      hier::NeighborhoodSet::const_iterator di = dst_eto_src.find(dst_mapped_box_id);
      if (di != dst_eto_src.end()) {
         const NeighborSet& src_nabrs = di->second;
         std::set<int> tmp_owners;
         src_nabrs.getOwners(tmp_owners);
         for (std::set<int>::const_iterator so = tmp_owners.begin();
              so != tmp_owners.end(); ++so) {
            const int& src_owner = *so;
            if (src_owner == dst_mapped_box_id.getOwnerRank()) {
               dst_to_fill_on_src_proc[dst_mapped_box_id] = tmp_fill_boxes;
            } else {
               std::vector<int>& send_mesg = send_mesgs[src_owner];
               send_mesg.insert(send_mesg.end(),
                  tmp_mesg.begin(),
                  tmp_mesg.end());
            }
         }
      }
   }
   // Send messages.
   for (std::set<int>::const_iterator si = src_owners.begin();
        si != src_owners.end(); ++si) {
      comms[mesg_number].initialize(&comm_stage);
      comms[mesg_number].setPeerRank(*si);
      // Reuse communicator.  Assuming it has no outstanding messages!
      comms[mesg_number].setMPI(src_to_dst.getMPI());
      comms[mesg_number].setMPITag(0, 1);
      comms[mesg_number].beginSend(&send_mesgs[*si][0],
                                   static_cast<int>(send_mesgs[*si].size()));
      if (comms[mesg_number].isDone()) {
         completed.insert(completed.end(), &comms[mesg_number]);
      }
      ++mesg_number;
   }

   /*
    * Complete communication and unpack messages.
    */
   do {
      for (size_t i = 0; i < completed.size(); ++i) {
         tbox::AsyncCommPeer<int>* peer =
            dynamic_cast<tbox::AsyncCommPeer<int> *>(completed[i]);
         TBOX_ASSERT(completed[i] != NULL);
         TBOX_ASSERT(peer != NULL);
         if (peer < comms + src_owners.size()) {
            // This is a receive.  Unpack it.  (Otherwise, ignore send completion.)
            const int* ptr = peer->getRecvData();
            while (ptr != peer->getRecvData() + peer->getRecvSize()) {
               const hier::BoxId distributed_id(hier::LocalId(ptr[0]),
                                                peer->getPeerRank(),
                                                hier::BlockId(ptr[1]));
               const unsigned int num_fill_mapped_boxes = ptr[2];
               ptr += 3;
               d_max_fill_boxes = tbox::MathUtilities<int>::Max(
                     d_max_fill_boxes,
                     num_fill_mapped_boxes);
               hier::MappedBoxSet& fill_boxes = dst_to_fill_on_src_proc[distributed_id];
               for (size_t ii = 0; ii < num_fill_mapped_boxes; ++ii) {
                  hier::Box tmp_dst_mapped_box(dim);
                  tmp_dst_mapped_box.getFromIntBuffer(ptr);
                  ptr += hier::Box::commBufferSize(dim);
                  fill_boxes.insert(fill_boxes.end(), tmp_dst_mapped_box);
               }
            }
            TBOX_ASSERT(ptr == peer->getRecvData() + peer->getRecvSize());
         }
      }
      completed.clear();
      comm_stage.advanceSome(completed);
   } while (!completed.empty());

   delete[] comms;
}

/*
 **************************************************************************
 * Calculate the maximum ghost cell width of all destination patch data
 * components.
 **************************************************************************
 */

hier::IntVector RefineSchedule::getMaxDestinationGhosts() const
{
   const tbox::Dimension& dim(d_dst_level->getDim());

   hier::IntVector gcw(dim, 0);
   tbox::Pointer<hier::PatchDescriptor> pd(d_dst_level->getPatchDescriptor());

   for (int iri = 0; iri < d_number_refine_items; iri++) {
      const int dst_id = d_refine_items[iri]->d_dst;
      gcw.max(pd->getPatchDataFactory(dst_id)->getGhostCellWidth());
   }

   return gcw;
}

/*
 **************************************************************************
 *
 * Calculate the maximum ghost cell width of all scratch patch data
 * components.
 *
 **************************************************************************
 */

hier::IntVector RefineSchedule::getMaxScratchGhosts() const
{
   const tbox::Dimension& dim(d_dst_level->getDim());

   hier::IntVector gcw(dim, 0);
   tbox::Pointer<hier::PatchDescriptor> pd(d_dst_level->getPatchDescriptor());

   for (int iri = 0; iri < d_number_refine_items; iri++) {
      const int scratch_id = d_refine_items[iri]->d_scratch;
      gcw.max(pd->getPatchDataFactory(scratch_id)->getGhostCellWidth());
   }

   return gcw;
}

/*
 **************************************************************************
 *
 * Calculate the maximum ghost cell width required for all stencils.
 *
 **************************************************************************
 */

hier::IntVector RefineSchedule::getMaxStencilGhosts() const
{
   const tbox::Dimension& dim(d_dst_level->getDim());

   hier::IntVector gcw(dim, 0);
   if (d_refine_patch_strategy) {
      gcw = d_refine_patch_strategy->getRefineOpStencilWidth();
   }

   for (int iri = 0; iri < d_number_refine_items; iri++) {
      if (!(d_refine_items[iri]->d_oprefine.isNull())) {
         gcw.max(d_refine_items[iri]->d_oprefine->getStencilWidth());
      }
   }

   return gcw;
}

/*
 *************************************************************************
 *
 * Private utility function that constructs schedule transactions that
 * move data from source patch on source level to destination patch
 * on destination level on regions defined by list of fil boxes.
 *
 *************************************************************************
 */

void RefineSchedule::constructScheduleTransactions(
   const hier::MappedBoxSet& fill_boxes,
   const hier::Box& dst_mapped_box,
   const hier::Box& src_mapped_box,
   bool use_time_interpolation)
{
   TBOX_ASSERT(!d_dst_level.isNull());
   TBOX_ASSERT(!d_src_level.isNull());
   TBOX_ASSERT(!dst_mapped_box.isPeriodicImage()); // src absorbs the shift, if any.

   const tbox::Dimension& dim(d_dst_level->getDim());
   const int my_rank = d_dst_level->getMappedBoxLevel()->getRank();

   const bool dst_is_local = (dst_mapped_box.getOwnerRank() == my_rank);
   const bool src_is_local = (src_mapped_box.getOwnerRank() == my_rank);

   tbox::Pointer<hier::Patch> dst_patch;
   tbox::Pointer<hier::Patch> src_patch;

   if (dst_is_local) {
      dst_patch = d_dst_level->getPatch(dst_mapped_box.getGlobalId(),
                                        dst_mapped_box.getBlockId());
   }
   if (src_is_local) {
      src_patch = d_src_level->getPatch(src_mapped_box.getGlobalId(),
                                        src_mapped_box.getBlockId());
   }

   const hier::IntVector& constant_zero_intvector(hier::IntVector::getZero(dim));
   const hier::IntVector& constant_one_intvector(hier::IntVector::getOne(dim));

   if (s_extra_debug) {
      tbox::plog << "constructScheduleTransactions: " << use_time_interpolation
                 << std::endl;
      tbox::plog << "  src: L" << d_src_level->getLevelNumber()
                 << "R" << d_src_level->getRatioToLevelZero()
                 << " / " << src_mapped_box << std::endl;
      tbox::plog << "  dst: L" << d_dst_level->getLevelNumber()
                 << "R" << d_dst_level->getRatioToLevelZero()
                 << " / " << dst_mapped_box << std::endl;
      tbox::plog << "  fill_boxes (" << fill_boxes.size() << ")";
      for (hier::MappedBoxSet::const_iterator bi = fill_boxes.begin();
           bi != fill_boxes.end(); ++bi) {
         tbox::plog << " " << *bi;
      }
      tbox::plog << std::endl;
   }

   /*
    * d_src_masks and d_overlaps exist only in this method, but are
    * class members instead of temporaries so that they don't
    * have to be reallocated every time this method is used.
    */
   int max_overlap_array_size = d_max_fill_boxes;

   if (d_src_masks.getNumberOfBoxes() < max_overlap_array_size) {
      for (int i = d_src_masks.getNumberOfBoxes();
           i < max_overlap_array_size; ++i) {
         d_src_masks.appendItem(hier::Box(dim));
      }
   }

   if (d_overlaps.getSize() < max_overlap_array_size) {
      d_overlaps.setNull();
      d_overlaps.resizeArray(max_overlap_array_size);
   }

   tbox::Pointer<hier::PatchDescriptor> dst_patch_descriptor(
      d_dst_level->getPatchDescriptor());
   tbox::Pointer<hier::PatchDescriptor> src_patch_descriptor(
      d_src_level->getPatchDescriptor());

   const hier::Box& dst_box = dst_mapped_box;
#ifdef DEBUG_CHECK_ASSERTIONS
   const hier::Box& src_box = src_mapped_box;
   // src_box used only in debug mode.  Ifdef out in opt mode to avoid warning.
#endif

   const bool same_patch =
      (d_dst_level == d_src_level &&
       dst_mapped_box.getGlobalId() == src_mapped_box.getGlobalId());
   const bool same_patch_no_shift =
      (same_patch && !src_mapped_box.isPeriodicImage());

   const int num_equiv_classes =
      d_refine_classes->getNumberOfEquivalenceClasses();

   const hier::PeriodicShiftCatalog* shift_catalog =
      hier::PeriodicShiftCatalog::getCatalog(dim);

   /*
    * Calculate the shift and the shifted source box.
    */
   hier::IntVector src_shift(dim, 0);
   hier::IntVector dst_shift(dim, 0);
   hier::Box unshifted_src_box = src_mapped_box;
   hier::Box unshifted_dst_box = dst_mapped_box;
   if (src_mapped_box.isPeriodicImage()) {
      TBOX_ASSERT(!dst_mapped_box.isPeriodicImage());
      src_shift = shift_catalog->shiftNumberToShiftDistance(
            src_mapped_box.getPeriodicId());
      src_shift = (d_src_level->getRatioToLevelZero() >
                   constant_zero_intvector) ?
         (src_shift * d_src_level->getRatioToLevelZero()) :
         hier::IntVector::ceiling(src_shift,
            -d_src_level->getRatioToLevelZero());
      unshifted_src_box.shift(-src_shift);
   }
   if (dst_mapped_box.isPeriodicImage()) {
      TBOX_ASSERT(!src_mapped_box.isPeriodicImage());
      dst_shift = shift_catalog->shiftNumberToShiftDistance(
            dst_mapped_box.getPeriodicId());
      dst_shift = (d_dst_level->getRatioToLevelZero() >
                   constant_zero_intvector) ?
         (dst_shift * d_dst_level->getRatioToLevelZero()) :
         hier::IntVector::ceiling(dst_shift,
            -d_dst_level->getRatioToLevelZero());
      unshifted_dst_box.shift(-dst_shift);
   }
   if (s_extra_debug) {
      tbox::plog << "  src_shift: " << src_shift
                 << " unshifted_src_box " << unshifted_src_box << std::endl;
      tbox::plog << "  dst_shift: " << dst_shift
                 << " unshifted_dst_box " << unshifted_dst_box << std::endl;
   }

   /*
    * Transformation initialized to src_shift with no rotation.
    * It will never be modified in single-block runs, nor in multiblock runs
    * when src_mapped_box and dst_mapped_box are on the same block.
    */
   hier::Transformation transformation(src_shift);

   /*
    * When src_mapped_box and dst_mapped_box are on different blocks
    * transformed_src_box is a representation of the source box in the
    * destination coordinate system.
    *
    * For all other cases, transformed_src_box is simply a copy of the
    * box from src_mapped_box.
    */
   hier::Box transformed_src_box(src_mapped_box);

   /*
    * When needed, transform the source box and determine if src and 
    * dst touch at an enhance connectivity singularity.
    */
   bool is_singularity = false;
   if (src_mapped_box.getBlockId() != dst_mapped_box.getBlockId()) {
      const hier::BlockId& dst_block_id = dst_mapped_box.getBlockId();
      const hier::BlockId& src_block_id = src_mapped_box.getBlockId();

      tbox::Pointer<hier::GridGeometry> grid_geometry(
         d_dst_level->getGridGeometry());

      hier::Transformation::RotationIdentifier rotation =
         grid_geometry->getRotationIdentifier(dst_block_id,
                                              src_block_id);
      hier::IntVector offset(
         grid_geometry->getOffset(dst_block_id, src_block_id));

      offset *= d_dst_level->getRatioToLevelZero();

      is_singularity = grid_geometry->areSingularityNeighbors(dst_block_id,
                                                              src_block_id);

      transformation = hier::Transformation(rotation, offset);
      transformation.transform(transformed_src_box);
   }

   /*
    * For any case except when src and dst touch at enhanced connectivity,
    * the transactions use d_dst_level as the destination level and
    * dst_mapped_box as the destination box.  For the enhanced connectivity
    * case, the destination level becomes d_encon_level and the destination
    * box is a member of d_encon_level.
    */
   tbox::Pointer<hier::PatchLevel> transaction_dst_level;
   hier::Box transaction_dst_mapped_box(dim);
   if (is_singularity) {

      /*
       * Determination of transaction_dst_mapped_box is done differently
       * depending of if dst_mapped_box is local.  When it is local
       * (regardless of whether source is also local), the appropriate
       * Box to assign to transaction_dst_mapped_box can be
       * found from the d_dst_to_encon connector.
       *
       * If the destination is not local, then the source is, 
       * and d_src_to_encon is searched to fined the right member of
       * d_encon_level to assign to transaction_dst_mapped_box.
       */
      if (dst_mapped_box.getOwnerRank() == my_rank) {

         const hier::NeighborhoodSet& encon_nbrhood_set =
            d_dst_to_encon.getNeighborhoodSets();

         hier::NeighborhoodSet::const_iterator ei =
            encon_nbrhood_set.find(dst_mapped_box.getId());

         const NeighborSet& encon_nbrs = ei->second;

         const hier::BlockId& src_block_id = src_mapped_box.getBlockId();

         for (NeighborSet::const_iterator en = encon_nbrs.begin();
              en != encon_nbrs.end(); ++en) {

            if (src_block_id == en->getBlockId()) {
               TBOX_ASSERT(transaction_dst_mapped_box.empty()); 
               transaction_dst_mapped_box = *en;
            }
         }

      } else {

         TBOX_ASSERT(src_mapped_box.getOwnerRank() == my_rank);

         hier::IntVector test_gcw(
            hier::IntVector::max(d_boundary_fill_ghost_width,
                                 d_max_stencil_width));
         test_gcw.max(hier::IntVector::getOne(dim));

         hier::Box test_dst_box(dst_mapped_box);
         test_dst_box.grow(test_gcw);

         const hier::NeighborhoodSet& encon_nbrhood_set =
            d_src_to_encon.getNeighborhoodSets();

         hier::NeighborhoodSet::const_iterator ei =
            encon_nbrhood_set.find(src_mapped_box.getId());

         const NeighborSet& encon_nbrs = ei->second;

         for (NeighborSet::const_iterator ni = encon_nbrs.begin();
              ni != encon_nbrs.end(); ++ni) {
            if (ni->getOwnerRank() == dst_mapped_box.getOwnerRank()) {
               hier::Box encon_box(*ni);
               transformation.transform(encon_box);

               if (test_dst_box.contains(encon_box)) {
                  TBOX_ASSERT(transaction_dst_mapped_box.empty()); 
                  transaction_dst_mapped_box = *ni; 
               }
            }
         }
      }

      transaction_dst_level = d_encon_level;

   } else {
      /*
       * All cases except for handling enhance connectivity neighbors
       * go here to do a simple assignment.
       */
      transaction_dst_level = d_dst_level;
      transaction_dst_mapped_box = dst_mapped_box;
   }

   for (int nc = 0; nc < num_equiv_classes; nc++) {

      const xfer::RefineClasses::Data& rep_item =
         d_refine_classes->getClassRepresentative(nc);

      const int rep_item_dst_id = rep_item.d_scratch;
      const int rep_item_src_id = rep_item.d_src;

      tbox::Pointer<hier::PatchDataFactory> src_pdf(
         src_patch_descriptor->getPatchDataFactory(rep_item_src_id));
      tbox::Pointer<hier::PatchDataFactory> dst_pdf(
         dst_patch_descriptor->getPatchDataFactory(rep_item_dst_id));

      const hier::IntVector& dst_gcw = dst_pdf->getGhostCellWidth();

      hier::BoxList::Iterator box_itr(d_src_masks);
      int box_num = 0;
      for (hier::MappedBoxSet::const_iterator bi = fill_boxes.begin();
           bi != fill_boxes.end(); ++bi) {

         const hier::Box& fill_box = *bi;

         /*
          * Get the patch data factories and calculate the overlap.
          * Note that we restrict the domain of the time interpolation
          * to the intersection of the fill box and the ghost box of
          * the destination patch data component.  This is needed for
          * the case where the schedule treats data components with
          * different ghost cell widths since the fill boxes are
          * generated using the largest ghost width.
          */

         hier::Box dst_fill_box(dst_box);
         dst_fill_box.grow(dst_gcw);
         dst_fill_box = dst_fill_box * fill_box; 

         tbox::Pointer<hier::BoxOverlap> overlap;
         hier::Box src_mask(dim);
         if (!is_singularity) {

            /*
             * Create overlap for normal cases (all but enhanced connectivity).
             */

            hier::Box test_mask(dst_fill_box * transformed_src_box);
            if (test_mask.empty() && dst_pdf->dataLivesOnPatchBorder()) {
               if ((dst_gcw == constant_zero_intvector) ||
                   (dst_box.isSpatiallyEqual(fill_box))) {

                  test_mask = dst_fill_box;
                  test_mask.grow(constant_one_intvector);
                  test_mask = test_mask * transformed_src_box;

               }
            }

            src_mask = test_mask;
            transformation.inverseTransform(src_mask);

            overlap =
               rep_item.d_var_fill_pattern->calculateOverlap(
                  *dst_pdf->getBoxGeometry(unshifted_dst_box),
                  *src_pdf->getBoxGeometry(unshifted_src_box),
                  dst_mapped_box,
                  src_mask,
                  fill_box,
                  true, transformation);
         } else {

            /*
             * Create overlap for enhanced connectivity.  This overlap
             * will be used in transaction from d_src_level to d_encon_level.
             */

            hier::Box test_mask(dst_fill_box);
            transformation.inverseTransform(test_mask);
            test_mask = test_mask * src_mapped_box;
            if (test_mask.empty() && dst_pdf->dataLivesOnPatchBorder()) {
               if ((dst_gcw == constant_zero_intvector) ||
                   (dst_box.isSpatiallyEqual(fill_box))) {

                  test_mask = dst_fill_box;
                  test_mask.grow(constant_one_intvector);
                  transformation.inverseTransform(test_mask);
                  test_mask = test_mask * src_mapped_box;

               }
            }

            src_mask = test_mask;
            hier::Box transformed_fill_box(fill_box);
            hier::Box transformed_dst_box(dst_mapped_box);
            transformation.inverseTransform(transformed_fill_box);
            transformation.inverseTransform(transformed_dst_box);

            overlap =
               rep_item.d_var_fill_pattern->calculateOverlap(
                  *dst_pdf->getBoxGeometry(transaction_dst_mapped_box),
                  *src_pdf->getBoxGeometry(unshifted_src_box),
                  transformed_dst_box,
                  src_mask,
                  transformed_fill_box,
                  true, hier::Transformation(hier::IntVector::getZero(dim)));
         }


#ifdef DEBUG_CHECK_ASSERTIONS
         if (overlap.isNull()) {
            TBOX_ERROR("Internal RefineSchedule error..."
               << "\n Overlap is NULL for "
               << "\n src box = " << src_box
               << "\n dst box = " << dst_box
               << "\n src mask = " << src_mask << std::endl);
         }
#endif

         if (s_extra_debug) {
            tbox::plog << "  overlap: ";
            overlap->print(tbox::plog);
         }
         *box_itr = src_mask;
         d_overlaps[box_num] = overlap;
         box_num++;
         box_itr++;

      }

      /*
       * Iterate over components in refine description list
       */
      for (tbox::List<int>::Iterator
           l(d_refine_classes->getIterator(nc)); l; l++) {
         const RefineClasses::Data& item =
            d_refine_classes->getRefineItem(l());
         TBOX_ASSERT(item.d_class_index == nc);

         const int dst_id = item.d_scratch;
         const int src_id = item.d_src;
         const int ritem_count = item.d_tag;

         /*
          * If the src and dst patches, levels, and components are the
          * same, and there is no shift, the data exchange is unnecessary.
          */
         if (!same_patch_no_shift || (dst_id != src_id)) {

            /*
             * Iterate over the fill boxes and create transactions
             * for each box that has a non-empty overlap.
             */
            hier::BoxList::Iterator itr(d_src_masks);
            for (int i = 0; i < box_num; i++, itr++) {

               /*
                * If overlap is not empty, then add the transaction
                * to the appropriate communication schedule.
                * There are two schedules depending on whether
                * coarse or fine data takes precedence at
                * coarse-fine boundaries for communications
                * where the destination variable quantity
                * has data residing on the boundary.
                * There are two types of transactions depending on
                * whether we use time interpolation.
                */

               /*
                * If we only perform the following block only for
                * non-null d_overlap[i], why not just loop through
                * box_num instead of num_fill_boxes?
                */
               if (!d_overlaps[i]->isOverlapEmpty()) {

                  tbox::Pointer<tbox::Transaction> transaction;

                  if (!d_transaction_factory.isNull()) {

                     transaction =
                        d_transaction_factory->allocate(transaction_dst_level,
                           d_src_level,
                           d_overlaps[i],
                           transaction_dst_mapped_box,
                           src_mapped_box,
                           ritem_count);
                  } else if (use_time_interpolation &&
                             item.d_time_interpolate) {

                     transaction = new xfer::RefineTimeTransaction(
                           transaction_dst_level, d_src_level,
                           d_overlaps[i],
                           transaction_dst_mapped_box, src_mapped_box,
                           *itr,
                           ritem_count);

                  } else {  // no time interpolation

                     transaction = new xfer::RefineCopyTransaction(
                           transaction_dst_level, d_src_level,
                           d_overlaps[i],
                           transaction_dst_mapped_box, src_mapped_box,
                           ritem_count);

                  }  // time interpolation conditional

                  if (item.d_fine_bdry_reps_var) {
                     if (same_patch) {
                        d_fine_priority_level_schedule->addTransaction(
                           transaction);
                     } else {
                        d_fine_priority_level_schedule->appendTransaction(
                           transaction);
                     }
                  } else {
                     if (same_patch) {
                        d_coarse_priority_level_schedule->addTransaction(
                           transaction);
                     } else {
                        d_coarse_priority_level_schedule->appendTransaction(
                           transaction);
                     }
                  }

               }  // if overlap not empty

            }  // iterate over fill_boxes

         }

      }  // iterate over refine components in equivalence class

   }  // iterate over refine equivalence classes
}

/*
 *************************************************************************
 *
 * Private member function to initialize data members for hierarchy info.
 *
 *************************************************************************
 */

void
RefineSchedule::initializeDomainAndGhostInformation(
   bool recursive_schedule)
{

   const tbox::Dimension& dim(d_dst_level->getDim());

   d_max_scratch_gcw = getMaxScratchGhosts();
   d_max_stencil_width = getMaxStencilGhosts();

   if (recursive_schedule) {
      d_boundary_fill_ghost_width = d_max_stencil_width;
      d_force_boundary_fill = (d_boundary_fill_ghost_width.max() > 0);
   } else {
      d_boundary_fill_ghost_width = getMaxDestinationGhosts();
      d_force_boundary_fill = false;
   }

   tbox::Pointer<hier::GridGeometry> grid_geom(d_dst_level->getGridGeometry());
   const hier::IntVector& ratio_to_level_zero =
      d_dst_level->getRatioToLevelZero();

   for (int b = 0; b < grid_geom->getNumberBlocks(); b++) {
      d_domain_is_one_box[b] = grid_geom->getDomainIsSingleBox(b);
   }

   d_periodic_shift = grid_geom->getPeriodicShift(ratio_to_level_zero);

   d_num_periodic_directions = 0;
   for (int d = 0; d < dim.getValue(); d++) {
      if (d_periodic_shift(d)) {
         d_num_periodic_directions++;
      }
   }

}

/*
 *************************************************************************
 *
 * Private utility function to set up local array of refine items.
 *
 *************************************************************************
 */

void RefineSchedule::setRefineItems(
   const tbox::Pointer<xfer::RefineClasses> refine_classes)
{

   clearRefineItems();

   d_refine_classes = refine_classes;

   d_number_refine_items = d_refine_classes->getNumberOfRefineItems();

   d_refine_items =
      new const xfer::RefineClasses::Data *[d_number_refine_items];

   int ircount = 0;
   for (int nd = 0; nd < d_number_refine_items; nd++) {
      d_refine_classes->getRefineItem(nd).d_tag = ircount;
      d_refine_items[ircount] = &(d_refine_classes->getRefineItem(nd));
      ircount++;
   }

}

/*
 *************************************************************************
 *
 * Private utility function used during initial schedule set up to
 * check whether patch data entries have proper number of ghost cells.
 * In particular, each scratch data entry must have at least as many
 * ghost cells as the user-defined refine operator stencil width.
 * Other checks are performed in the
 * xfer::RefineClasses::itemIsValid() routine.
 *
 *************************************************************************
 */

void RefineSchedule::initialCheckRefineClassItems() const
{
   const tbox::Dimension& dim(d_dst_level->getDim());

   const hier::IntVector& constant_zero_intvector(hier::IntVector::getZero(dim));

   tbox::Pointer<hier::PatchDescriptor> pd(d_dst_level->getPatchDescriptor());

   hier::IntVector user_gcw(constant_zero_intvector);
   if (d_refine_patch_strategy) {
      user_gcw = d_refine_patch_strategy->getRefineOpStencilWidth();
   }

   if (user_gcw > constant_zero_intvector) {

      for (int iri = 0; iri < d_number_refine_items; iri++) {

         const xfer::RefineClasses::Data * const ref_item = d_refine_items[iri];

#ifdef DEBUG_CHECK_ASSERTIONS
         if (d_refine_classes->itemIsValid(*ref_item, pd)) {
#endif

         const int scratch = ref_item->d_scratch;
         const hier::IntVector& scratch_gcw(pd->getPatchDataFactory(scratch)->
                                            getGhostCellWidth());

         if (user_gcw > scratch_gcw) {
            TBOX_ERROR("Bad data given to RefineSchedule...\n"
               << "User supplied interpolation stencil width = "
               << user_gcw
               << "\nis larger than ghost cell width of `Scratch'\n"
               << "patch data " << pd->mapIndexToName(scratch)
               << " , which is " << scratch_gcw << std::endl);
         }

#ifdef DEBUG_CHECK_ASSERTIONS
      }
#endif

      }

   }

}

/*
 **************************************************************************
 *
 * Private utility function to clear array of refine items.
 *
 **************************************************************************
 */

void RefineSchedule::clearRefineItems()
{
   if (d_refine_items) {
      for (int iri = 0; iri < d_number_refine_items; iri++) {
         d_refine_items[iri] = (xfer::RefineClasses::Data *)NULL;
      }
      delete[] d_refine_items;
      d_refine_items = (const xfer::RefineClasses::Data **)NULL;
      d_number_refine_items = 0;
   }
}

/*
 **************************************************************************
 *
 * Print class data to the specified output stream.
 *
 **************************************************************************
 */

void RefineSchedule::printClassData(
   std::ostream& stream) const
{
   stream << "RefineSchedule::printClassData()\n";
   stream << "--------------------------------------\n";

   d_refine_classes->printClassData(stream);

   stream << "Printing coarse priority refine schedule...\n";
   d_coarse_priority_level_schedule->printClassData(stream);

   stream << "Printing fine priority refine schedule...\n";
   d_fine_priority_level_schedule->printClassData(stream);

   if (!d_supp_schedule.isNull()) {
      stream << "Printing supplemental refine schedule...\n";
      d_supp_schedule->printClassData(stream);
   }
}
/*
 **************************************************************************
 **************************************************************************
 */

void RefineSchedule::initializeCallback()
{
   tbox::Pointer<tbox::Database> idb(tbox::InputManager::getInputDatabase());
   if (idb && idb->isDatabase("RefineSchedule")) {
      tbox::Pointer<tbox::Database> rsdb(idb->getDatabase(
            "RefineSchedule"));
      s_extra_debug =
         rsdb->getBoolWithDefault("extra_debug", s_extra_debug);
      s_barrier_and_time =
         rsdb->getBoolWithDefault("barrier_and_time", s_barrier_and_time);
   }

   t_fill_data = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::fillData()");
   t_recursive_fill = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::recursive_fill");
   t_refine_scratch_data = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::refineScratchData()");
   t_finish_sched_const = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::finishScheduleConstruction()");
   t_finish_sched_const_recurse = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::finishScheduleConstruction()_recurse");
   t_gen_comm_sched = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::generateCommunicationSchedule()");
   t_bridge_connector = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::bridge_connector");
   t_modify_connector = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::modify_connector");
   t_make_seq_map = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::make_seq_map");
   t_shear = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::finish...()_shear");
   t_misc1 = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::finish...()_misc1");
   t_barrier_and_time = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::barrier_and_time");
   t_get_global_mapped_box_count = tbox::TimerManager::getManager()->
      getTimer(
         "xfer::RefineSchedule::finish...()_get_global_mapped_box_count");
   t_coarse_shear = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::finish...()_coarse_shear");
   t_setup_supp_mapped_box_level = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::setupSupplementalMappedBoxLevel()");
   t_misc2 = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::finish...()_misc2");
   t_bridge_supp_hiercoarse = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::finish...()_bridge_supp_hiercoarse");
   t_bridge_dst_hiercoarse = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::finish...()_bridge_dst_hiercoarse");
   t_make_supp_level = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::finish...()_make_supp_level");
   t_make_supp_to_unfilled = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::finish...()_make_supp_to_unfilled");
   t_invert_edges = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::generate...()_invert_edges");
   t_construct_send_trans = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::generate...()_construct_send_trans");
   t_construct_recv_trans = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::generate...()_construct_recv_trans");

}

/*
 ***************************************************************************
 * Release static timers.  To be called by shutdown registry to make sure
 * memory for timers does not leak.
 ***************************************************************************
 */
void RefineSchedule::finalizeCallback()
{
   t_fill_data.setNull();
   t_recursive_fill.setNull();
   t_refine_scratch_data.setNull();
   t_finish_sched_const.setNull();
   t_finish_sched_const_recurse.setNull();
   t_gen_comm_sched.setNull();
   t_bridge_connector.setNull();
   t_modify_connector.setNull();
   t_make_seq_map.setNull();
   t_shear.setNull();
   t_misc1.setNull();
   t_barrier_and_time.setNull();
   t_get_global_mapped_box_count.setNull();
   t_coarse_shear.setNull();
   t_setup_supp_mapped_box_level.setNull();
   t_misc2.setNull();
   t_bridge_supp_hiercoarse.setNull();
   t_bridge_dst_hiercoarse.setNull();
   t_make_supp_level.setNull();
   t_make_supp_to_unfilled.setNull();
   t_invert_edges.setNull();
   t_construct_send_trans.setNull();
   t_construct_recv_trans.setNull();
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif

#endif
