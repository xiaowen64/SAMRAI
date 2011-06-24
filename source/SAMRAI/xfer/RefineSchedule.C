/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
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
#include "SAMRAI/hier/MappedBoxContainerUtils.h"
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
tbox::Pointer<tbox::Timer> RefineSchedule::t_build_supp_mapped_box_level;
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
 * list.  Ony data on the intersection of the two levels will be
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
   const hier::BlockId& block_id,
   const tbox::Pointer<xfer::RefineClasses> refine_classes,
   tbox::Pointer<xfer::RefineTransactionFactory> transaction_factory,
   xfer::RefinePatchStrategy* patch_strategy,
   bool use_time_refinement):
   d_max_stencil_width(dst_level->getDim()),
   d_max_scratch_gcw(dst_level->getDim()),
   d_boundary_fill_ghost_width(dst_level->getDim()),
   d_domain_box(dst_level->getDim()),
   d_periodic_shift(dst_level->getDim()),
   d_unfilled_mapped_box_level(dst_level->getDim()),
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

   d_domain_is_one_box = false;
   d_num_periodic_directions = 0;

   d_coarse_priority_level_schedule = new tbox::Schedule();
   d_fine_priority_level_schedule = new tbox::Schedule();
   d_coarse_priority_level_schedule->setTimerPrefix("xfer::RefineSchedule");
   d_fine_priority_level_schedule->setTimerPrefix("xfer::RefineSchedule");

   d_supp_schedule.setNull();
   d_supp_level.setNull();

   d_max_fill_boxes = 0;
   d_block_id = block_id;

   /*
    * Initialize destination level, ghost cell widths,
    * and domain information data members.
    */

   bool recursive_schedule = false;
   initializeDomainAndGhostInformation(recursive_schedule);

   hier::IntVector min_connector_width = d_max_scratch_gcw;
   min_connector_width.max(d_boundary_fill_ghost_width);

   const Connector dst_to_src =
      dst_level->getMappedBoxLevel()->getPersistentOverlapConnectors().findConnector(
         *src_level->getMappedBoxLevel(),
         min_connector_width);

   const Connector src_to_dst =
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
    * Create the fill box and unfilled box arrays and then the
    * communication schedule for data transfers between source and
    * destination levels.  Note that the fill boxes are initialized here,
    * while the unfilled boxes are set in the generateCommunicationSchedule()
    * and contain the regions for each patch in the destination level
    * that will cannot be filled by the schedule.
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

   MappedBoxLevel unused_unfilled_mapped_box_level(dim);
   Connector unused_dst_to_unfilled;
   generateCommunicationSchedule(
      unused_unfilled_mapped_box_level,
      unused_dst_to_unfilled,
      dst_to_src,
      src_to_dst,
      dst_to_fill,
      dst_to_fill_on_src_proc,
      use_time_refinement);

   if (!d_supp_level.isNull()) {
      computeRefineOverlaps();
   }
}

/*
 **************************************************************************
 *
 * Create a refine schedule that copies data from the source level into
 * the destination level on the components represented by the refine
 * list.  If portions of the destination level remain unfilled, then
 * the algorithm   recursively fills those unfilled portions from coarser
 * levels in the   AMR hierarchy.  It is assumed that the index spaces of
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
   const hier::BlockId& block_id,
   tbox::Pointer<hier::PatchHierarchy> hierarchy,
   const tbox::Pointer<xfer::RefineClasses> refine_classes,
   tbox::Pointer<xfer::RefineTransactionFactory> transaction_factory,
   xfer::RefinePatchStrategy* patch_strategy,
   bool use_time_refinement):
   d_max_stencil_width(dst_level->getDim()),
   d_max_scratch_gcw(dst_level->getDim()),
   d_boundary_fill_ghost_width(dst_level->getDim()),
   d_domain_box(dst_level->getDim()),
   d_periodic_shift(dst_level->getDim()),
   d_unfilled_mapped_box_level(dst_level->getDim()),
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

   d_domain_is_one_box = false;
   d_num_periodic_directions = 0;

   d_coarse_priority_level_schedule.setNull();
   d_fine_priority_level_schedule.setNull();

   d_supp_schedule.setNull();
   d_supp_level.setNull();

   d_max_fill_boxes = 0;

   d_block_id = block_id;

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
    * Create the fill box arrays and then the communication schedule(s)
    * needed to move data from the patch hierarchy to the destination level.
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

   if (!d_supp_schedule.isNull()) {
      computeRefineOverlaps();
   }
}

/*
 **************************************************************************
 *
 * This private constructor creates a refine schedule that copies data
 * into the destination only on the specified fill boxes.
 *
 **************************************************************************
 */

RefineSchedule::RefineSchedule(
   tbox::Pointer<hier::PatchLevel> dst_level,
   tbox::Pointer<hier::PatchLevel> src_level,
   int next_coarser_ln,
   const hier::BlockId& block_id,
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
   d_domain_box(dst_level->getDim()),
   d_periodic_shift(dst_level->getDim()),
   d_unfilled_mapped_box_level(dst_level->getDim()),
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

   d_domain_is_one_box = false;
   d_num_periodic_directions = 0;

   d_coarse_priority_level_schedule.setNull();
   d_fine_priority_level_schedule.setNull();

   d_supp_schedule.setNull();
   d_supp_level.setNull();

   d_max_fill_boxes = 0;
   d_block_id = block_id;

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
    * Finish construction of the communication schedule using the
    * remaining fill boxes.
    *
    * The fill boxes are the dst boxes, grown by the max stencil width
    * (not the ghost width of data).  The reason is that in this
    * private constructor, the dst is always the supplemental level
    * constructed from another RefineSchedule's unfilled boxes.  So,
    * dst is nominally the unfilled boxes.  To interpolate into them,
    * we need data covered by the interpolation stencil width.
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

   if (!d_supp_schedule.isNull()) {
      computeRefineOverlaps();
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
 * Construct schedule with possibility of recursion.
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
   const hier::IntVector zero_vector(hier::IntVector::getZero(dim));
   const hier::IntVector neg1_vector(dim, -1);
   const hier::IntVector huge_vector(dim, SAMRAI::tbox::MathUtilities<int>::getMax());

   hier::MappedBoxLevelConnectorUtils edge_utils;
   hier::OverlapConnectorAlgorithm oca;

   RefineScheduleConnectorWidthRequestor rscwri;
   std::vector<hier::IntVector> self_connector_widths;
   std::vector<hier::IntVector> fine_connector_widths;
   rscwri.computeRequiredConnectorWidths(self_connector_widths,
                                         fine_connector_widths,
                                         *hierarchy);

   const Connector dummy_cnect;
   const MappedBoxLevel& dst_mapped_box_level = dst_to_fill.getBase();
   if (!d_src_level.isNull()) {
      TBOX_ASSERT(dst_to_src.isInitialized());
   }

   /*
    * hiercoarse is the coarse level on the hierarchy.  It is to be
    * differentiated from the supplemental (supp) level, which is at
    * the same refimenent ratio but is not on the hierarchy.
    */
   tbox::Pointer<hier::PatchLevel> hiercoarse_level;
   if (next_coarser_ln >= 0)
      hiercoarse_level = hierarchy->getPatchLevel(next_coarser_ln);

   /*
    * Width of supp<==>hiercoarse must be big enough to
    * (1) generate transactions between supp and hiercoarse, and
    * (2) bridge supp<==>hiercoarser in the next iteration.
    *    (hiercoarser is next coarser level after hiercoarse.)
    *
    * (1) requires width computed as required by
    * generateCommunicationSchedule() for building
    * transaction between supp and hiercoarser:
    * hier::IntVector::max(d_max_scratch_gcw, getMaxDestinationGhosts());
    *
    * (2) requires that hiercoarse^width(hiercoarse->supp) is big
    * enough to nest supp and its ghost data so that we can bridge
    * for supp<==>hiercoarser and guarantee that supp<==>hiercoarser
    * does not miss any critical overlaps.
    */
   hier::IntVector transaction_gcw =
      hier::IntVector::max(d_max_scratch_gcw, getMaxDestinationGhosts());
   transaction_gcw =
      hier::IntVector::max(transaction_gcw, hier::IntVector(dim, 1));


   d_coarse_priority_level_schedule = new tbox::Schedule();
   d_fine_priority_level_schedule = new tbox::Schedule();

   /*
    * If the source level is not null, then generate a communication
    * schedule for moving data from source level to destination level.
    * The schedule generation routine determines the boxes that remain
    * to be filled from coarser levels; i.e., they cannot be filled from
    * the source level.  If the source level is null, then copy all of the
    * fill boxes into the set of boxes to be filled from coarser levels
    * in the hierarchy.
    */

   Connector dst_to_unfilled;

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
   if (d_src_level.isNull() || skip_generate_schedule) {
      d_unfilled_mapped_box_level = fill_mapped_box_level;
      dst_to_unfilled.initialize(
         dst_mapped_box_level,
         d_unfilled_mapped_box_level,
         dst_to_fill.getConnectorWidth(),
         dst_to_fill.getNeighborhoodSets(),
         MappedBoxLevel::DISTRIBUTED);
      dst_to_unfilled.setConnectorType(hier::Connector::BASE_GENERATED);
   } else {
      generateCommunicationSchedule(
         d_unfilled_mapped_box_level,
         dst_to_unfilled,
         dst_to_src,
         src_to_dst,
         dst_to_fill,
         dst_to_fill_on_src_proc,
         use_time_interpolation);
   }
   if (s_extra_debug) {
      tbox::plog << "finishScheduleConstruction in recursion_level="
                 << recursion_level << " next_coarser_ln=" << next_coarser_ln
                 << " after computing unfilled boxes."
                 << "\nd_unfilled_mapped_box_level:\n"
                 << d_unfilled_mapped_box_level.format("UF->", 2)
                 << "\ndst_to_unfilled:\n"
                 << dst_to_unfilled.format("DUF->", 2)
                 << std::endl;
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   const int nblocks = fill_mapped_box_level.getGridGeometry()->getNumberBlocks();
   for ( int bn=0; bn<nblocks; ++bn ) {
      TBOX_ASSERT(fill_mapped_box_level.getLocalBoundingBox(bn).contains(
                     d_unfilled_mapped_box_level.getLocalBoundingBox(bn)));
   }
#endif

   const MappedBoxLevel& periodic_domain_mapped_box_level =
      hierarchy->getDomainMappedBoxLevel();
   TBOX_ASSERT(
      periodic_domain_mapped_box_level.getParallelState() ==
      MappedBoxLevel::GLOBALIZED);

   /*
    * Remove pieces of the boxes to be filled from coarser levels that
    * live outside the physical domain in the non-periodic directions.
    * Then, check whether there still exist boxes to be filled from
    * coarser levels.  Any resulting boxes will be used to generate a
    * supplemental coarse patch level for acquiring data to refine.
    */

   hier::IntVector big_grow_vector(dim, 0);
   if (d_num_periodic_directions > 0) {
      for (int dir = 0; dir < dim.getValue(); dir++) {
         if (d_periodic_shift(dir)) {
            big_grow_vector(dir) = BIG_GHOST_CELL_WIDTH;
         }
      }
   }

   TBOX_ASSERT(d_num_periodic_directions >= 0);
   const bool fully_periodic = d_num_periodic_directions == dim.getValue();

   if (!fully_periodic) {
      /*
       * Shearing is the removal of parts of unfilled_boxes that lie
       * along physical boundaries.  We should not fill these because
       * they would be filled by boundary conditions.
       *
       * We bypass shearing for fully_periodic domains, where it would
       * be a no-op anyway.
       */

      t_shear->start();

      // Shearing for the mapped_box_level.
      Connector unfilled_to_periodic_domain(
         d_unfilled_mapped_box_level,
         periodic_domain_mapped_box_level,
         dst_to_fill.getConnectorWidth());
      oca.findOverlaps(unfilled_to_periodic_domain);
      Connector unfilled_to_sheared;
      MappedBoxLevel sheared_mapped_box_level(dim);
      edge_utils.computeInternalParts(
         sheared_mapped_box_level,
         unfilled_to_sheared,
         unfilled_to_periodic_domain,
         zero_vector);
      t_modify_connector->start();
      hier::MappingConnectorAlgorithm mca;
      mca.modify(dst_to_unfilled,
         unfilled_to_sheared,
         &d_unfilled_mapped_box_level);
      t_modify_connector->stop();
      dst_to_unfilled.eraseEmptyNeighborSets();

      t_shear->stop();

   } // !fully_periodic

   t_get_global_mapped_box_count->barrierAndStart();
   const bool need_to_fill = d_unfilled_mapped_box_level.getGlobalNumberOfBoxes();
   t_get_global_mapped_box_count->stop();

   /*
    * If there remain boxes to be filled from coarser levels, then
    * generate data to do this and recurse to the next coarser level.
    *
    * This big if-block is the last step in this finishScheduleConstruction.
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
       * If there are no coarser levels in the hierarchy or the hierarchy
       * is null, then throw an error.  Something is messed up someplace and
       * code execution cannot proceed.
       */

      if (next_coarser_ln < 0) {
         TBOX_ERROR(
            "Internal error in RefineSchedule::finishScheduleConstruction..."
            << "\n In finishScheduleConstruction() -- "
            << "\n No coarser levels...will not fill from coarser."
            << "\n dst_mapped_box_level:\n" << dst_mapped_box_level.format("DEST->", 2)
            << "\n dst_to_fill:\n" << dst_to_fill.format("DF->", 2)
            << "\n d_unfilled_mapped_box_level:\n" << d_unfilled_mapped_box_level.format("UF->", 2)
            << "\n dst_to_unfilled:\n" << dst_to_unfilled.format("DU->", 2)
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

      const hier::MappedBoxLevel &hiercoarse_mapped_box_level(
         *hiercoarse_level->getMappedBoxLevel());

      /*
       * Calculate the ratio to the next coarser level in the hierarchy
       * and cache the maximum stencil ghost cell width.
       */

      const hier::IntVector dst_hiercoarse_ratio =
         d_dst_level->getRatioToLevelZero()
         / hiercoarse_level->getRatioToLevelZero();

      const hier::BoxList coarser_physical_domain(
         hiercoarse_level->getPhysicalDomain(d_block_id));

      const bool do_coarse_shearing = (!fully_periodic && !d_domain_is_one_box);

      hier::BoxList coarser_shear_domain(coarser_physical_domain);

      if (do_coarse_shearing) {
         t_coarse_shear->start();

         if (d_num_periodic_directions > 0) {
            coarser_shear_domain.grow(big_grow_vector);
         }

         coarser_shear_domain.simplifyBoxes();

         t_coarse_shear->stop();
      }

      /*
       * Convert the unfilled parts of dst into the supplemental MappedBoxLevel
       * supp_mapped_box_level by coarsening it (and shearing if needed).
       * The supp_mapped_box_level is to be filled by the next coarser
       * level in the hierarchy.  The data to be filled on the supp level
       * will be the max stencil width.
       *
       * Build connector d_dst_to_supp using dst_to_unfilled, since
       * supp is came from unfilled.  This Connector is incomplete
       * because each dst MappedBox only has edges to supp MappedBoxes
       * it generated.  Nevertheless, we set its gcw big enough so each
       * dst MappedBox nests its potential supp MappedBoxes so that we
       * bridge to the supp MappedBoxLevel.
       *
       * Build up supp_eto_unfilled, so we can find the unfilled box
       * that the suppplemental box is supposed to fill.
       * d_supp_to_unfilled is used when filling data.
       */

      hier::MappedBoxLevel supp_mapped_box_level(
         hiercoarse_level->getRatioToLevelZero(),
         hiercoarse_level->getGridGeometry(),
         d_unfilled_mapped_box_level.getMPI());

      hier::NeighborhoodSet dst_eto_supp, supp_eto_unfilled;

      const hier::NeighborhoodSet& dst_eto_unfilled = dst_to_unfilled.getNeighborhoodSets();

      t_build_supp_mapped_box_level->start();
      /*
       * This loop builds up supp_mapped_box_level, dst_eto_supp and supp_eto_unfilled.
       */
      for (hier::NeighborhoodSet::const_iterator ei = dst_eto_unfilled.begin();
           ei != dst_eto_unfilled.end(); ++ei) {

         const hier::MappedBoxId& dst_mapped_box_mbid = ei->first;
         const NeighborSet& dst_unfilled_parts = ei->second;

         hier::Box coarser_fill_bounding_box(dim);
         for (hier::MappedBoxSet::const_iterator ni = dst_unfilled_parts.begin();
              ni != dst_unfilled_parts.end(); ++ni) {

            const MappedBox& unfilled_mapped_box = *ni;
            hier::Box supp_box = unfilled_mapped_box.getBox();
            supp_box.coarsen(dst_hiercoarse_ratio);

            if (do_coarse_shearing &&
                (d_dst_level->patchTouchesRegularBoundary(
                    dst_mapped_box_mbid))) {

               hier::BoxList sheared_supp_boxes(supp_box);
               sheared_supp_boxes.intersectBoxes(coarser_shear_domain);
               sheared_supp_boxes.simplifyBoxes();

               (void)hier::BoxUtilities::extendBoxesToDomainBoundary(
                  sheared_supp_boxes,
                  coarser_physical_domain,
                  d_max_stencil_width);
               /*
                * Connector widths must be big enough to make sure
                * we have complete sets after extending mapped_boxes to boundary!
                */
               if (sheared_supp_boxes.size() > 0) {

                  NeighborSet& supp_nabrs = dst_eto_supp[dst_mapped_box_mbid];

                  for (hier::BoxList::Iterator b(sheared_supp_boxes); b; b++) {
                     const MappedBox& supp_mapped_box =
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
                  coarser_physical_domain,
                  d_max_stencil_width);
               if (supp_box != save_supp_box) {
                  TBOX_WARNING("Supplemental box " << save_supp_box
                                                   << " was extended to " << supp_box
                                                   << " at a physical boundary.  This is"
                                                   << " probably ok but a rigorous proof"
                                                   << " that won't cause problem is currently"
                                                   << " lacking.  Expect some sanity checks"
                                                   << " to fail and a slim chance that the"
                                                   << " schedule generate will be bad.");
               }

               const MappedBox& supp_mapped_box = *supp_mapped_box_level.addBox(
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
       *   - width of dst-->fill, but rounded up so it extends equal
       *     to the growth of supp caused by coarsening unfilled.
       *   - extended by the stencil width, where supp has its ghost data.
       */
      const hier::IntVector dst_to_supp_width =
         (hier::IntVector::ceiling(dst_to_fill.getConnectorWidth(),
             dst_hiercoarse_ratio) + d_max_stencil_width)
         * dst_hiercoarse_ratio;

      d_dst_to_supp.swapInitialize(
         dst_mapped_box_level,
         supp_mapped_box_level,
         dst_to_supp_width,
         dst_eto_supp);
      d_dst_to_supp.setConnectorType(hier::Connector::BASE_GENERATED);

      d_supp_to_unfilled.initialize(
         supp_mapped_box_level,
         d_unfilled_mapped_box_level,
         hier::IntVector(dim, 0),
         supp_eto_unfilled);
      d_supp_to_unfilled.setConnectorType(hier::Connector::BASE_GENERATED);

      t_build_supp_mapped_box_level->stop();

      /*
       * Connect the supplemental mapped_box_level (the next recursion's dst)
       * to the hiercoarse mapped_box_level (the next recursion's src).
       */

      /*
       * Get the transpose of d_dst_to_supp, which is simple to compute
       * because we know the edges are all local.
       */
      t_misc2->start();
      d_supp_to_dst.initializeToLocalTranspose(d_dst_to_supp);
      t_misc2->stop();
      if (s_extra_debug) {

         /*
          * Supp should nest in dst^dst_to_supp_width to ensure dst
          * sees all of supp and also supp's ghosts.  Note that supp's
          * relevant ghost data width is d_max_stencil_width.
          *
          * The nesting assures that when bridging across dst<==>supp
          * for supp<==>hiercoarse, we get a complete overlap
          * Connectors.
          */
         bool locally_nests;
         edge_utils.setSanityCheckMethodPreconditions(false);
         edge_utils.setSanityCheckMethodPostconditions(false);
         tbox::plog << "\nsupp_mapped_box_level:\n" << supp_mapped_box_level.format("S-> ", 2)
                    << "\ndst_mapped_box_level:\n" << dst_mapped_box_level.format("D-> ", 2)
                    << "\nd_dst_to_supp:\n" << d_dst_to_supp.format("DS-> ", 2)
                    << "\nd_supp_to_dst:\n" << d_supp_to_dst.format("SD-> ", 2)
                    << std::endl;
         if ( ! edge_utils.baseNestsInHead(
                 &locally_nests,
                 supp_mapped_box_level,
                 dst_mapped_box_level,
                 zero_vector,
                 dst_to_supp_width,
                 zero_vector,
                 NULL) ) {
            TBOX_ERROR("RefineSchedule::finishScheduleConstruction: supp does\n"
                       <<"to nest in dst.\n");
         }

         TBOX_ASSERT(d_supp_to_dst.checkTransposeCorrectness(d_dst_to_supp) == 0);
         TBOX_ASSERT(d_dst_to_supp.checkTransposeCorrectness(d_supp_to_dst) == 0);
      }

      Connector supp_to_hiercoarse;
      Connector hiercoarse_to_supp;
      {
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
          * complete.  However, supp nests in dst+gcw, so they need not be
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
          * bridge for them.
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
                * Connectors dst<==>hiercoarse is not provided.
                * We have to bridge through src for it.
                * (This requires src<==>hiercoarse.)
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

               {
                  /*
                   * hiercoarse has periodic images but dst does not.
                   * Consequently, the bridge for dst<==>hiercoarse
                   * may not catch all periodic neighbors.  Because we
                   * don't need periodic neighbors for
                   * dst<==>hiercoarse, we will remove them to make
                   * dst<==>hiercoarse proper transposes.
                   */
                  hier::NeighborhoodSet tmp_edges;
                  bridged_dst_to_hiercoarse.getNeighborhoodSets().removePeriodicNeighbors(
                     tmp_edges);
                  bridged_dst_to_hiercoarse.swapInitialize(
                     bridged_dst_to_hiercoarse.getBase(),
                     bridged_dst_to_hiercoarse.getHead(),
                     bridged_dst_to_hiercoarse.getConnectorWidth(),
                     tmp_edges,
                     MappedBoxLevel::DISTRIBUTED);
                  tmp_edges.clear();
                  bridged_hiercoarse_to_dst.getNeighborhoodSets().removePeriodicNeighbors(
                     tmp_edges);
                  bridged_hiercoarse_to_dst.swapInitialize(
                     bridged_hiercoarse_to_dst.getBase(),
                     bridged_hiercoarse_to_dst.getHead(),
                     bridged_hiercoarse_to_dst.getConnectorWidth(),
                     tmp_edges,
                     MappedBoxLevel::DISTRIBUTED);
               }

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
                        << next_coarser_ln << ") have problems as reported above.\n"
                        <<"recursion_level=" << recursion_level);
                  }
               }
               dst_to_hiercoarse = &bridged_dst_to_hiercoarse;
               hiercoarse_to_dst = &bridged_hiercoarse_to_dst;
            } // End block bridging for dst<==>hiercoarse.
         } else { /* !dst_is_supplemental_level */
                  /*
                   * dst may be the level next_coarser_ln+1 on the
                   * hierarchy, or it could be a coarsened version.
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
          */

         if (s_barrier_and_time) {
            t_bridge_supp_hiercoarse->barrierAndStart();
         }
         /*
          * Don't use the strict bridge theorem here because it
          * cannot guarantee sufficient width.  We know from how
          * dst nests in hiercoarse what output Connector width
          * can guarantee that all dst MappedBoxes are seen by a
          * hiercoarse MappedBox.
          */
         oca.bridge(
            supp_to_hiercoarse,
            hiercoarse_to_supp,
            d_supp_to_dst,
            *dst_to_hiercoarse,
            *hiercoarse_to_dst,
            d_dst_to_supp,
            fine_connector_widths[next_coarser_ln]);
         if (s_extra_debug) {
            bool locally_nests = false;
            if ( !edge_utils.baseNestsInHead(
                    &locally_nests,
                    supp_to_hiercoarse.getBase(),
                    supp_to_hiercoarse.getHead(),
                    d_max_stencil_width,
                    supp_to_hiercoarse.getConnectorWidth(),
                    zero_vector,
                    &hierarchy->getDomainSearchTree(d_block_id)) ) {
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

         if (s_barrier_and_time) {
            t_bridge_supp_hiercoarse->stop();
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
             * Removing periodic edges should make supp<==>hiercoarse
             * proper transposes.
             */
            hier::NeighborhoodSet tmp_edges;
            supp_to_hiercoarse.getNeighborhoodSets().removePeriodicNeighbors(
               tmp_edges);
            supp_to_hiercoarse.swapInitialize(
               supp_to_hiercoarse.getBase(),
               supp_to_hiercoarse.getHead(),
               supp_to_hiercoarse.getConnectorWidth(),
               tmp_edges,
               MappedBoxLevel::DISTRIBUTED);
            supp_to_hiercoarse.setConnectorType(hier::Connector::INCOMPLETE_OVERLAP);
            tmp_edges.clear();
            hiercoarse_to_supp.getNeighborhoodSets().removePeriodicNeighbors(
               tmp_edges);
            hiercoarse_to_supp.swapInitialize(
               hiercoarse_to_supp.getBase(),
               hiercoarse_to_supp.getHead(),
               hiercoarse_to_supp.getConnectorWidth(),
               tmp_edges,
               MappedBoxLevel::DISTRIBUTED);
            hiercoarse_to_supp.setConnectorType(hier::Connector::INCOMPLETE_OVERLAP);
         }
         if (s_extra_debug) {
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
            bool locally_nests;
            bool err3 = !edge_utils.baseNestsInHead(
               &locally_nests,
               supp_mapped_box_level,
               *hiercoarse_level->getMappedBoxLevel(),
               d_max_stencil_width,
               fine_connector_widths[next_coarser_ln],
               zero_vector,
               &hierarchy->getDomainSearchTree(d_block_id));
            if (err3) tbox::perr
               << "\nsupp does not sufficiently nest in hiercoarse."
               << std::endl;

            Connector complete_supp_to_hiercoarse(
               supp_to_hiercoarse.getBase(),
               supp_to_hiercoarse.getHead(),
               supp_to_hiercoarse.getConnectorWidth());
            oca.findOverlaps(complete_supp_to_hiercoarse);
            MappedBoxLevel external(dim);
            Connector supp_to_external;
            edge_utils.computeExternalParts(
               external,
               supp_to_external,
               complete_supp_to_hiercoarse,
               fine_connector_widths[next_coarser_ln]-d_max_stencil_width,
               hierarchy->getPeriodicDomainSearchTree(d_block_id));
            supp_to_external.eraseEmptyNeighborSets();
            int err4 = supp_to_external.getGlobalNumberOfRelationships();
            if (err4) {
               tbox::perr << "Some parts of supp lies outside of where we\n"
                          << "guarantee support for recursive RefineSchedule.\n"
                          << supp_to_external.format("SE: ", 2);
            }

            if (err1 || err2 || err3) {
               TBOX_ERROR(
                  "supp<==>hiercoarse have problems as reported above."
                  << "dst:\n" << dst_mapped_box_level.format("DEST->", 2)
                  << "supp:\n" << supp_mapped_box_level.format("SUPP->", 2)
                  << "hiercoarse:\n" << hierarchy->getMappedBoxLevel(next_coarser_ln)->format("HCRS->", 2)
                  << "dst_to_hiercoarse:\n" << dst_to_hiercoarse->format("DH->", 2)
                  << "dst_to_supp:\n" << d_dst_to_supp.format("DS->", 2)
                  << "supp_to_hiercoarse:\n" << supp_to_hiercoarse.format("SH->", 2)
                  << "hiercoarse_to_supp:\n" << hiercoarse_to_supp.format("HS->", 2));
            }
         }

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
               hiercoarse_level->getMappedBoxLevel()->getPersistentOverlapConnectors().findConnector(*hiercoarse_level->getMappedBoxLevel(),
               hiercoarse_to_supp.getConnectorWidth());
            edge_utils.setSanityCheckMethodPreconditions(false);
            edge_utils.setSanityCheckMethodPostconditions(false);
            edge_utils.addPeriodicImagesAndRelationships(
               supp_mapped_box_level,
               supp_to_hiercoarse,
               hiercoarse_to_supp,
               hierarchy->getDomainSearchTree(d_block_id),
               hiercoarse_to_hiercoarse);
            edge_utils.setSanityCheckMethodPreconditions(false);
            edge_utils.setSanityCheckMethodPostconditions(false);
         }
         if (s_extra_debug) {
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
            bool locally_nests;
            bool err3 = !edge_utils.baseNestsInHead(
               &locally_nests,
               supp_mapped_box_level,
               *hiercoarse_level->getMappedBoxLevel(),
               d_max_stencil_width,
               fine_connector_widths[next_coarser_ln],
               zero_vector,
               &hierarchy->getDomainSearchTree(d_block_id));
            if (err3) tbox::perr
               << "\nsupp does not sufficiently nest in hiercoarse."
               << std::endl;

            Connector complete_supp_to_hiercoarse(
               supp_to_hiercoarse.getBase(),
               supp_to_hiercoarse.getHead(),
               supp_to_hiercoarse.getConnectorWidth());
            oca.findOverlaps(complete_supp_to_hiercoarse);
            MappedBoxLevel external(dim);
            Connector supp_to_external;
            edge_utils.computeExternalParts(
               external,
               supp_to_external,
               complete_supp_to_hiercoarse,
               fine_connector_widths[next_coarser_ln]-d_max_stencil_width,
               hierarchy->getPeriodicDomainSearchTree(d_block_id));
            supp_to_external.eraseEmptyNeighborSets();
            int err4 = supp_to_external.getGlobalNumberOfRelationships();
            if (err4) {
               tbox::perr << "Some parts of supp lies outside of where we\n"
                          << "guarantee support for recursive RefineSchedule.\n"
                          << supp_to_external.format("SE: ", 2);
            }

            if (err1 || err2 || err3 || err4) {
               TBOX_ERROR(
                  "supp<==>hiercoarse have problems as reported above.\n"
                  << "dst:\n" << dst_mapped_box_level.format("DEST->", 2)
                  << "supp:\n" << supp_mapped_box_level.format("SUPP->", 2)
                  << "hiercoarse:\n" << hierarchy->getMappedBoxLevel(next_coarser_ln)->format("HCRS->", 2)
                  << "dst_to_hiercoarse:\n" << dst_to_hiercoarse->format("DH->", 2)
                  << "dst_to_supp:\n" << d_dst_to_supp.format("DS->", 2)
                  << "supp_to_hiercoarse:\n" << supp_to_hiercoarse.format("SH->", 2)
                  << "hiercoarse_to_supp:\n" << hiercoarse_to_supp.format("HS->", 2));
            }
         }

      } // Block to compute supp<==>hiercoarse.

      t_make_supp_level->start();
      d_supp_level = new hier::PatchLevel(
            supp_mapped_box_level,
            hiercoarse_level->getGridGeometry(),
            hiercoarse_level->getPatchDescriptor());
      t_make_supp_level->stop();
      d_supp_level->setLevelNumber(next_coarser_ln);

      /*
       * Reset supp<==>hiercoarse connectors to use the PatchLevel's
       * MappedBoxLevel.
       */
      supp_to_hiercoarse.initialize( *d_supp_level->getMappedBoxLevel(),
                                     supp_to_hiercoarse.getHead(),
                                     supp_to_hiercoarse.getConnectorWidth(),
                                     supp_to_hiercoarse.getNeighborhoodSets() );
      hiercoarse_to_supp.initialize( hiercoarse_to_supp.getBase(),
                                     *d_supp_level->getMappedBoxLevel(),
                                     hiercoarse_to_supp.getConnectorWidth(),
                                     hiercoarse_to_supp.getNeighborhoodSets()  );


      /*
       * Compute how much hiercoarse has to grow to nest supp.  If dst
       * is a supplemental level (generated by RefineSchedule), we
       * have the info to compute the growth.  If not, we make some
       * assumptions about where dst came from in order to determine
       * how it nests in hiercoarse.
       */
      hier::IntVector hiercoarse_growth_to_nest_supp(dim, -1);
      if (dst_is_supplemental_level) {
         /*
          * Assume that src barely nests in hiercoarse.  (In most
          * places, it nests by a margin equal to the nesting buffer,
          * but we don't count on that because the nesting buffer is
          * relaxed at physical boundaries.
          *
          * FIXME: We may in fact be able to count on the nesting buffer because
          * extending boxes to physical boundaries do not create any
          * extra relationships.  However, we don't currently have
          * access to the size of the nesting buffer.)
          */
         hiercoarse_growth_to_nest_supp =
            src_growth_to_nest_dst + dst_to_fill.getConnectorWidth();
         hiercoarse_growth_to_nest_supp.ceiling(dst_hiercoarse_ratio);
      } else {
         /*
          * dst may be:
          * 1. the hierarchy level just finer than level number next_coarser_ln.
          * 2. a level that nests in level number next_coarser_ln.
          *    a. a new level generated by GriddingAlgorithm.
          *    b. the hierarchy level just finer than level number next_coarser_ln,
          *       coarsened for Richardson extrapolation.
          * In any case, dst should nest in hiercoarse.  Furthermore, it does
          * not grow when coarsened into the hiercoarse index space.
          */
         // hiercoarse_growth_to_nest_supp = hier::IntVector::getZero(dim);
         hiercoarse_growth_to_nest_supp = dst_to_fill.getConnectorWidth();
         hiercoarse_growth_to_nest_supp.ceiling(dst_hiercoarse_ratio);
      }

      /*
       * Recursively fill the coarse schedule using the private
       * refine schedule constructor.
       */

      /*
       * We need to make sure that the coarse schedule uses
       * BoxGeometryVariableFillPattern, so that it fills all needed parts of
       * d_supp_level
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

      t_finish_sched_const_recurse->stop();
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
            d_block_id,
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

   if (s_extra_debug) {
      tbox::plog << "finishScheduleConstruction exiting recursion_level="
                 << recursion_level << " next_coarser_ln=" << next_coarser_ln
                 << std::endl;
      --recursion_level;
   }
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

   if (!d_supp_schedule.isNull()) {

      /*
       * Allocate data on the coarser level and keep track of the allocated
       * components so that they may be deallocated later.
       */

      hier::ComponentSelector allocate_vector;
      allocateScratchSpace(allocate_vector, d_supp_level, fill_time);

      /*
       * Recursively call the fill routine to fill the required coarse fill
       * boxes on the coarser level.
       */

      d_supp_schedule->recursiveFill(fill_time, do_physical_boundary_fill);

      /*
       * The coarse fill boxes should now be filled.  Now interpolate
       * data from the coarse grid into the fine grid.
       */

      refineScratchData();

      /*
       * Deallocate the scratch data from the coarse grid.
       */

      d_supp_level->deallocatePatchData(allocate_vector);
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
         tbox::Pointer<hier::Patch> patch = *p;
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
      tbox::Pointer<hier::Patch> patch = *p;

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

void RefineSchedule::refineScratchData() const
{
   t_refine_scratch_data->start();

   const hier::IntVector ratio(d_dst_level->getRatioToLevelZero()
                               / d_supp_level->getRatioToLevelZero());

   tbox::ListIterator<tbox::Array<tbox::Pointer<hier::BoxOverlap> > >
      overlap_iter(d_refine_overlaps);

   /*
    * Loop over all the supplemental patches and find the corresponding
    * destination patch and destination fill boxes.
    */

   const hier::MappedBoxSet& supp_mapped_boxes =
      d_supp_level->getMappedBoxLevel()->getMappedBoxes();
   for (hier::MappedBoxSet::const_iterator ni = supp_mapped_boxes.begin();
        ni != supp_mapped_boxes.end(); ++ni) {

      const MappedBox& supp_mapped_box = *ni;
      const hier::MappedBoxSet& dst_nabrs =
         d_supp_to_dst.getNeighborSet(supp_mapped_box.getId());
      const hier::MappedBox& dst_mapped_box = *dst_nabrs.begin();
#ifdef DEBUG_CHECK_ASSERTIONS
      /*
       * Each supp_mapped_box can point back to just one dst_mapped_box.
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
      tbox::Pointer<hier::Patch> fine_patch = d_dst_level->getPatch(
            dst_mapped_box.getGlobalId(), supp_mapped_box.getBlockId());
      tbox::Pointer<hier::Patch> supp_patch = d_supp_level->getPatch(
            supp_mapped_box.getGlobalId(), supp_mapped_box.getBlockId());

      const NeighborSet& unfilled_nabrs = d_supp_to_unfilled.getNeighborSet(
            supp_mapped_box.getId());
      TBOX_ASSERT(unfilled_nabrs.size() == 1);
      const MappedBox& unfilled_nabr = *unfilled_nabrs.begin();
      hier::BoxList fill_boxes(unfilled_nabr.getBox());

      if (d_refine_patch_strategy) {
         d_refine_patch_strategy->preprocessRefineBoxes(*fine_patch,
            *supp_patch,
            fill_boxes,
            ratio);
      }

      for (int iri = 0; iri < d_number_refine_items; iri++) {
         const xfer::RefineClasses::Data * const ref_item = d_refine_items[iri];         if (!(ref_item->d_oprefine.isNull())) {

            tbox::Pointer<hier::BoxOverlap> refine_overlap(
               overlap_iter()[ref_item->d_class_id]);

            const int scratch_id = ref_item->d_scratch;

            ref_item->d_oprefine->refine(*fine_patch, *supp_patch,
               scratch_id, scratch_id,
               *refine_overlap, ratio);

         }
      }

      overlap_iter++;

      if (d_refine_patch_strategy) {
         d_refine_patch_strategy->postprocessRefineBoxes(*fine_patch,
            *supp_patch,
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
void RefineSchedule::computeRefineOverlaps()
{

   tbox::Pointer<hier::PatchDescriptor> fine_patch_descriptor =
      d_dst_level->getPatchDescriptor();

   const int num_equiv_classes =
      d_refine_classes->getNumberOfEquivalenceClasses();

   /*
    * Loop over all the supplemental patches and find the corresponding
    * destination patch and destination fill boxes.
    */

   const hier::MappedBoxSet& supp_mapped_boxes =
      d_supp_level->getMappedBoxLevel()->getMappedBoxes();
   for (hier::MappedBoxSet::const_iterator ni = supp_mapped_boxes.begin();
        ni != supp_mapped_boxes.end(); ++ni) {

      const MappedBox& supp_mapped_box = *ni;
      const hier::MappedBoxSet& dst_nabrs =
         d_supp_to_dst.getNeighborSet(supp_mapped_box.getId());
      const hier::MappedBox& dst_mapped_box = *dst_nabrs.begin();
#ifdef DEBUG_CHECK_ASSERTIONS
      /*
       * Each supp_mapped_box can point back to just one dst_mapped_box.
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
      tbox::Pointer<hier::Patch> fine_patch = d_dst_level->getPatch(
            dst_mapped_box.getId());

      const NeighborSet& unfilled_nabrs = d_supp_to_unfilled.getNeighborSet(
            supp_mapped_box.getId());
      TBOX_ASSERT(unfilled_nabrs.size() == 1);
      const MappedBox& unfilled_nabr = *unfilled_nabrs.begin();
      hier::BoxList fill_boxes(unfilled_nabr.getBox());

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
            tbox::Pointer<hier::PatchDataFactory> fine_pdf =
               fine_patch_descriptor->getPatchDataFactory(scratch_id);
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
      d_refine_overlaps.appendItem(refine_overlaps);
   }
}

/*
 *************************************************************************
 *
 * Generate communication schedule routine creates transactions to move
 * data from interiors of the source space on the source level into the
 * specified fill box regions of the destination level.  Each fill box
 * will typically be the interior plus max ghost cells over all data
 * components on the destination region.  If the source and the scratch
 * are the same and the source and destination levels  are the same,
 * then there is no need to copy on the interiors for the same patch.
 *
 * The resulting transactions will only fill the regions of intersection
 * between the fill_boxes and destination level boxes.  The remaining
 * box regions are returned in unfilled_boxes.
 *
 * This main schedule generation routine which passes control to one
 * of the algorithmic variations below based on boolean parameters set
 * to default settings in the constructors and possibly changed via
 * an input file.
 *
 *********************************************************************
 */

void RefineSchedule::generateCommunicationSchedule(
   MappedBoxLevel& unfilled_mapped_box_level,
   Connector& dst_to_unfilled,
   const Connector& dst_to_src,
   const Connector& src_to_dst,
   const Connector& dst_to_fill,
   const FillSet& dst_to_fill_on_src_proc,
   const bool use_time_interpolation)
{
   t_gen_comm_sched->start();

   const tbox::Dimension& dim(d_dst_level->getDim());

   if (s_extra_debug) {
      if (dst_to_src.isInitialized()) {
         dst_to_src.assertTransposeCorrectness(src_to_dst);
         src_to_dst.assertTransposeCorrectness(dst_to_src);
      }
   }
   if (s_extra_debug) {
      tbox::plog
         << "generateCommunicationSchedule: dst_to_src:" << dst_to_src.format("\t", 2)
         << "generateCommunicationSchedule: src_to_dst:" << src_to_dst.format("\t", 2)
         << "generateCommunicationSchedule: dst_to_fill:" << dst_to_fill.format("\t", 2);
   }
   /*
    * We invert the edges so that the src transactions on the local
    * process is created in the same order the dst transactions
    * are created on remote processes.  The ordering is how the
    * source and destination transactions are matched up on the
    * receiving process.  We chose to order the transactions
    * dst-major (all the transactions for one dst MappedBox are
    * grouped together, then all transactions for the next
    * dst MappedBox, and so on).
    */
   const hier::NeighborhoodSet& dst_eto_src = dst_to_src.getNeighborhoodSets();

   const MappedBoxLevel& dst_mapped_box_level = dst_to_src.getBase();

   const bool dst_finer = dst_to_src.getHeadCoarserFlag();
   const bool src_finer = !dst_to_src.getHeadCoarserFlag() &&
      dst_to_src.getRatio() != hier::IntVector::getOne(dim);


   /*
    * Construct sending transactions for local src mapped_boxes (except
    * for transactions with local dst mapped_boxes).
    */
   if (s_extra_debug) {
      tbox::plog << "Constructing send transactions" << std::endl;
   }
   t_construct_send_trans->start();

   /*
    * Restructure the src_to_dst edge data to arange neighbors by the
    * dst mapped_boxes, as required to match the transaction ordering
    * on the receiving processors.  At the same time, shift src-dst
    * pairs to make the dst shifts zero.
    */
   FullNeighborhoodSet src_eto_dst_bydst;
   t_invert_edges->start();
   restructureNeighborhoodSetsByDstNodes(src_eto_dst_bydst, src_to_dst);
   t_invert_edges->stop();

   for (FullNeighborhoodSet::const_iterator ei = src_eto_dst_bydst.begin();
        ei != src_eto_dst_bydst.end(); ++ei) {

      /*
       * dst_mapped_box can be remote (by definition of FullNeighborhoodSet).
       * local_src_mapped_boxes are the local source mapped_boxes that
       * contribute data to dst_mapped_box.
       */
      const MappedBox& dst_mapped_box = ei->first;
      const hier::MappedBoxSet& local_src_mapped_boxes = ei->second;
      TBOX_ASSERT(!dst_mapped_box.isPeriodicImage());

      FillSet::const_iterator dst_fill_boxes_ =
         dst_to_fill_on_src_proc.find(dst_mapped_box.getGlobalId());
      if (dst_fill_boxes_ == dst_to_fill_on_src_proc.end()) {
         /*
          * Missing fill boxes should indicate that the dst mapped_box
          * has no fill box.  One way this is possible is for
          * d_dst_level_fill_pattern to be of type PatchLevelBorderFillPattern
          * and for dst_mapped_box to be away from level borders.
          */
         continue;
      }
      const MappedBoxVector& dst_fill_boxes = dst_fill_boxes_->second;

      /*
       * Construct transactions for data going from local source mapped_boxes
       * to remote dst mapped_boxes.
       */
      for (hier::MappedBoxSet::const_iterator ni = local_src_mapped_boxes.begin();
           ni != local_src_mapped_boxes.end(); ++ni) {
         const MappedBox& src_mapped_box = *ni;
         if (src_mapped_box.getOwnerRank() == dst_mapped_box.getOwnerRank()) {
            /*
             * Disregard local dst_mapped_box to avoid duplicating same
             * transactions created by the second loop below.
             */
            continue;
         }

         if (src_mapped_box.getBlockId() == dst_mapped_box.getBlockId()) {
            constructScheduleTransactions(
               dst_fill_boxes,
               dst_mapped_box,
               src_mapped_box,
               use_time_interpolation);
         }
      }

   }
   t_construct_send_trans->stop();

   /*
    * Construct receiving transactions for local dst mapped_boxes.
    * Set unfilled_mapped_box_level, and dst_eto_unfilled.
    */
   if (s_extra_debug) {
      tbox::plog << "Constructing receive and copy transactions" << std::endl;
   }
   const hier::NeighborhoodSet& dst_eto_fill = dst_to_fill.getNeighborhoodSets();
   hier::MappedBoxSet unfilled_mapped_boxes;
   hier::NeighborhoodSet dst_eto_unfilled;
   hier::LocalId last_unfilled_index(-1);
   t_construct_recv_trans->start();
   for (hier::NeighborhoodSet::const_iterator cf = dst_eto_fill.begin();
        cf != dst_eto_fill.end(); ++cf) {

      const hier::MappedBoxId& dst_mapped_box_id(cf->first);
      const MappedBox& dst_mapped_box = *dst_mapped_box_level.getMappedBox(
            dst_mapped_box_id);

      const NeighborSet& fill_nabrs = cf->second;

      hier::NeighborhoodSet::const_iterator cs = dst_eto_src.find(dst_mapped_box_id);

      if (cs != dst_eto_src.end()) {
         /*
          * dst_mapped_box has overlaps with src mapped_box_level.
          * Fill what we can using src and add the remaining to unfilled.
          */

         const hier::BlockId& dst_block_id = dst_mapped_box.getBlockId();
         MappedBoxVector fill_boxes;
         for (hier::MappedBoxSet::iterator bi = fill_nabrs.begin();
              bi != fill_nabrs.end(); ++bi) {
            if (bi->getBlockId() == dst_block_id) {
               fill_boxes.push_back(*bi);
            }
         }
         hier::BoxList unfilled_box_for_dst;
         if (d_block_id == hier::BlockId::invalidId() ||
             dst_block_id == d_block_id) {
            for (MappedBoxVector::iterator vi = fill_boxes.begin();
                 vi != fill_boxes.end(); ++vi) {
               unfilled_box_for_dst.appendItem(vi->getBox());
            }
         }

         const NeighborSet& src_nabrs = cs->second;
         for (NeighborSet::const_iterator na = src_nabrs.begin();
              na != src_nabrs.end(); ++na) {

            const MappedBox& src_mapped_box = *na;

            if (src_mapped_box.getBlockId() == dst_mapped_box.getBlockId()) {

               constructScheduleTransactions(
                  fill_boxes,
                  dst_mapped_box,
                  src_mapped_box,
                  use_time_interpolation);
               hier::Box src_box_in_dst_space = src_mapped_box.getBox();
               if (dst_finer) {
                  src_box_in_dst_space.refine(src_to_dst.getRatio());
               } else if (src_finer) {
                  src_box_in_dst_space.coarsen(src_to_dst.getRatio());
               }
               unfilled_box_for_dst.removeIntersections(src_box_in_dst_space);
            }
         }

         if (!unfilled_box_for_dst.isEmpty()) {
            NeighborSet& unfilled_nabrs = dst_eto_unfilled[dst_mapped_box_id];
            for (hier::BoxList::Iterator bi(unfilled_box_for_dst);
                 bi; bi++) {
               const MappedBox unfilled_mapped_box =
                  *unfilled_mapped_boxes.insert(unfilled_mapped_boxes.end(),
                     MappedBox(bi(),
                        ++last_unfilled_index,
                        dst_mapped_box.getOwnerRank(),
                        dst_mapped_box.getBlockId()));
               unfilled_nabrs.insert(unfilled_nabrs.end(), unfilled_mapped_box);
            }
         }

      } else {
         /*
          * dst_mapped_box does not overlap with any mapped_boxes in src.
          * All of it goes into unfilled.
          */
         NeighborSet& unfilled_nabrs = dst_eto_unfilled[dst_mapped_box_id];
         for (NeighborSet::const_iterator na = fill_nabrs.begin();
              na != fill_nabrs.end(); ++na) {
            const MappedBox& fill_nabr = *na;
            const MappedBox& unfilled_mapped_box =
               *unfilled_mapped_boxes.insert(unfilled_mapped_boxes.end(),
                  MappedBox(fill_nabr.getBox(),
                     ++last_unfilled_index,
                     fill_nabr.getOwnerRank(),
                     fill_nabr.getBlockId()));
            unfilled_nabrs.insert(unfilled_nabrs.end(), unfilled_mapped_box);
         }
      }

   } // End loop through dst mapped_boxes.
   t_construct_recv_trans->stop();

   unfilled_mapped_box_level.initialize(
      unfilled_mapped_boxes,
      dst_mapped_box_level.getRefinementRatio(),
      dst_mapped_box_level.getGridGeometry(),
      dst_mapped_box_level.getMPI());

   dst_to_unfilled.initialize(
      dst_mapped_box_level,
      unfilled_mapped_box_level,
      dst_to_fill.getConnectorWidth(),
      dst_eto_unfilled,
      MappedBoxLevel::DISTRIBUTED);
   dst_to_unfilled.setConnectorType(hier::Connector::BASE_GENERATED);

   const int nblocks = unfilled_mapped_box_level.getGridGeometry()->getNumberBlocks();
   const hier::MappedBoxLevel &fill_mapped_box_level(dst_to_fill.getHead());
   for ( int bn=0; bn<nblocks; ++bn ) {
      if ( !fill_mapped_box_level.getLocalBoundingBox(bn).isEmpty() &&
           !fill_mapped_box_level.getLocalBoundingBox(bn).contains(
             unfilled_mapped_box_level.getLocalBoundingBox(bn))) {
         TBOX_ERROR(
            "RefineSchedule::generateCommunicationSchedule: library error\n"
            << "Unfilled bounding box "
            << unfilled_mapped_box_level.getLocalBoundingBox(bn)
            << "does not nest inside fill bounding box "
            << dst_to_fill.getHead().getLocalBoundingBox(bn)
            << " in block " << bn << ".\n"
            << "fill:\n" << dst_to_fill.getHead().format("ERR:F-> ", 3)
            << "unfilled:\n" << unfilled_mapped_box_level.format("ERR:U-> ", 3)
            << "dst-->fill:\n" << dst_to_fill.format("ERR:DF->", 3)
            << "dst-->unfilled:\n" << dst_to_fill.format("ERR:DU-> ", 3));
      }
   }

   t_gen_comm_sched->stop();
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
 * 2. It shifts periodic image dst MappedBoxes back to the zero-shift position,
 * and applies a similar shift to src MappedBoxes so that the overlap is
 * unchanged.  The constructScheduleTransactions method requires all
 * shifts to be absorbed in the src MappedBox.
 ***********************************************************************
 */
void RefineSchedule::restructureNeighborhoodSetsByDstNodes(
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
   MappedBox shifted_mapped_box(dim), unshifted_nabr(dim);
   full_inverted_edges.clear();
   for (hier::NeighborhoodSet::const_iterator ci = edges.begin();
        ci != edges.end();
        ++ci) {
      const MappedBox& mapped_box =
         *src_mapped_box_level.getMappedBoxStrict(ci->first);
      const NeighborSet& nabrs = ci->second;
      for (NeighborSet::const_iterator na = nabrs.begin();
           na != nabrs.end(); ++na) {
         const hier::MappedBox& nabr = *na;
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
   dst_to_fill.setConnectorType(hier::Connector::BASE_GENERATED);

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
   MappedBoxVector tmp_fill_boxes;
   const hier::NeighborhoodSet& dst_eto_fill = dst_to_fill.getNeighborhoodSets();
   for (hier::NeighborhoodSet::const_iterator ei = dst_eto_fill.begin();
        ei != dst_eto_fill.end(); ++ei) {
      const hier::MappedBoxId& dst_mapped_box_id = ei->first;
      const NeighborSet& fill_nabrs = ei->second;
      /*
       * Pack dst_mapped_box_id's fill box info into tmp_mesg.
       * - dst_mapped_box_id's LocalId
       * - number of fill neighbors
       * - fill neighbors (could just send box and save 2 ints)
       * Also, create BoxVector object for local use.
       */
      tmp_mesg.clear();
      tmp_mesg.reserve(2 + fill_nabrs.size() * MappedBox::commBufferSize(dim));
      tmp_mesg.insert(tmp_mesg.end(), 2, 0);
      tmp_mesg[0] = dst_mapped_box_id.getLocalId().getValue();
      tmp_mesg[1] = static_cast<int>(fill_nabrs.size());
      tmp_fill_boxes.clear();
      tmp_fill_boxes.reserve(fill_nabrs.size());
      for (NeighborSet::const_iterator na = fill_nabrs.begin();
           na != fill_nabrs.end(); ++na) {
         tmp_mesg.insert(tmp_mesg.end(), MappedBox::commBufferSize(dim), 0);
         na->putToIntBuffer(&tmp_mesg[tmp_mesg.size()
                                      - MappedBox::commBufferSize(dim)]);
         tmp_fill_boxes.insert(tmp_fill_boxes.end(), *na);
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
               dst_to_fill_on_src_proc[dst_mapped_box_id.getGlobalId()] = tmp_fill_boxes;
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
               const hier::GlobalId distributed_id(hier::LocalId(ptr[0]), peer->getPeerRank());
               const unsigned int num_fill_mapped_boxes = ptr[1];
               ptr += 2;
               d_max_fill_boxes = tbox::MathUtilities<int>::Max(
                     d_max_fill_boxes,
                     num_fill_mapped_boxes);
               MappedBoxVector& fill_boxes = dst_to_fill_on_src_proc[distributed_id];
               for (size_t ii = 0; ii < num_fill_mapped_boxes; ++ii) {
                  MappedBox tmp_dst_mapped_box(dim);
                  tmp_dst_mapped_box.getFromIntBuffer(ptr);
                  ptr += MappedBox::commBufferSize(dim);
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
   tbox::Pointer<hier::PatchDescriptor> pd = d_dst_level->getPatchDescriptor();

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
   tbox::Pointer<hier::PatchDescriptor> pd = d_dst_level->getPatchDescriptor();

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
   const MappedBoxVector& fill_boxes,
   const MappedBox& dst_mapped_box,
   const MappedBox& src_mapped_box,
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
      for (MappedBoxVector::const_iterator bi = fill_boxes.begin();
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
      for (int i = d_src_masks.getNumberOfBoxes(); i < max_overlap_array_size; ++i) {
         d_src_masks.appendItem(hier::Box(dim));
      }
   }
   if (d_overlaps.getSize() < max_overlap_array_size) {
      d_overlaps.setNull();
      d_overlaps.resizeArray(max_overlap_array_size);
   }

   tbox::Pointer<hier::PatchDescriptor> dst_patch_descriptor =
      d_dst_level->getPatchDescriptor();
   tbox::Pointer<hier::PatchDescriptor> src_patch_descriptor =
      d_src_level->getPatchDescriptor();

   const hier::Box& dst_box = dst_mapped_box.getBox();
#ifdef DEBUG_CHECK_ASSERTIONS
   const hier::Box& src_box = src_mapped_box.getBox();
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
   hier::Box unshifted_src_box = src_mapped_box.getBox();
   hier::Box unshifted_dst_box = dst_mapped_box.getBox();
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

   for (int nc = 0; nc < num_equiv_classes; nc++) {

      const xfer::RefineClasses::Data& rep_item =
         d_refine_classes->getClassRepresentative(nc);

      const int rep_item_dst_id = rep_item.d_scratch;
      const int rep_item_src_id = rep_item.d_src;

      tbox::Pointer<hier::PatchDataFactory> src_pdf =
         src_patch_descriptor->getPatchDataFactory(rep_item_src_id);
      tbox::Pointer<hier::PatchDataFactory> dst_pdf =
         dst_patch_descriptor->getPatchDataFactory(rep_item_dst_id);

      const hier::IntVector& dst_gcw = dst_pdf->getGhostCellWidth();

      hier::BoxList::Iterator box_itr(d_src_masks);
      int box_num = 0;
      for (MappedBoxVector::const_iterator bi = fill_boxes.begin();
           bi != fill_boxes.end(); ++bi) {

         if (d_block_id == hier::BlockId::invalidId()) {
            if (bi->getBlockId() != dst_mapped_box.getBlockId()) {
               continue;
            }
         } else if (bi->getBlockId() != d_block_id) {
            continue;
         }

         const hier::Box& fill_box = bi->getBox();

         /*
          * Get the patch data factories and calculate the overlap.
          * Note that we restrict the domain of the time interpolation
          * to the intersection of the fill box and the ghost box of
          * the destination patch data component.  This is needed for
          * the case where the schedule treats data components with
          * different ghost cell widths since the fill boxes are
          * generated using the largest ghost width.
          */

         const hier::Box dst_fill_box =
            hier::Box::grow(dst_box, dst_gcw) * fill_box;

         hier::Box test_mask = dst_fill_box * src_mapped_box.getBox();
         if (test_mask.empty() && dst_pdf->dataLivesOnPatchBorder()) {
            if ((dst_gcw == constant_zero_intvector) ||
                (dst_box == fill_box)) { 
               hier::Box tmp_dst_fill_box(
                  hier::Box::grow(dst_fill_box, constant_one_intvector));
               test_mask = tmp_dst_fill_box * src_mapped_box.getBox();
            }
         }
         const hier::Box src_mask = hier::Box::shift(test_mask, -src_shift);

         tbox::Pointer<hier::BoxOverlap> overlap =
            rep_item.d_var_fill_pattern->calculateOverlap(
               *dst_pdf->getBoxGeometry(unshifted_dst_box),
               *src_pdf->getBoxGeometry(unshifted_src_box),
               dst_mapped_box.getBox(),
               src_mask,
               fill_box,
               true, hier::Transformation(src_shift));

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
         TBOX_ASSERT(item.d_class_id == nc);

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
                        d_transaction_factory->allocate(d_dst_level,
                           d_src_level,
                           d_overlaps[i],
                           dst_mapped_box,
                           src_mapped_box,
                           ritem_count);
                  } else if (use_time_interpolation &&
                             item.d_time_interpolate) {

                     transaction = new xfer::RefineTimeTransaction(
                           d_dst_level, d_src_level,
                           d_overlaps[i],
                           dst_mapped_box, src_mapped_box,
                           *itr,
                           ritem_count);

                  } else {  // no time interpolation

                     transaction = new xfer::RefineCopyTransaction(
                           d_dst_level, d_src_level,
                           d_overlaps[i],
                           dst_mapped_box, src_mapped_box,
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

   tbox::Pointer<hier::GridGeometry> grid_geom = d_dst_level->getGridGeometry();
   const hier::IntVector& ratio_to_level_zero =
      d_dst_level->getRatioToLevelZero();

   d_domain_is_one_box = grid_geom->getDomainIsSingleBox(0);

   if (d_domain_is_one_box) {
      hier::BoxList domain(dim);
      grid_geom->computePhysicalDomain(domain,
         ratio_to_level_zero, hier::BlockId::zero());
      hier::BoxList::Iterator itr(domain);
      d_domain_box = *itr;
      itr++;
      for (; itr; itr++) {
         d_domain_box += *itr;
      }
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
 * xfer::RefineClasses::checkRefineItem() routine.
 *
 *************************************************************************
 */

void RefineSchedule::initialCheckRefineClassItems() const
{
   const tbox::Dimension& dim(d_dst_level->getDim());

   const hier::IntVector& constant_zero_intvector(hier::IntVector::getZero(dim));

   tbox::Pointer<hier::PatchDescriptor> pd = d_dst_level->getPatchDescriptor();

   hier::IntVector user_gcw(constant_zero_intvector);
   if (d_refine_patch_strategy) {
      user_gcw = d_refine_patch_strategy->getRefineOpStencilWidth();
   }

   if (user_gcw > constant_zero_intvector) {

      for (int iri = 0; iri < d_number_refine_items; iri++) {

         const xfer::RefineClasses::Data * const ref_item = d_refine_items[iri];

#ifdef DEBUG_CHECK_ASSERTIONS
         if (d_refine_classes->checkRefineItem(*ref_item, pd)) {
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
   tbox::Pointer<tbox::Database> idb = tbox::InputManager::getInputDatabase();
   if (idb && idb->isDatabase("RefineSchedule")) {
      tbox::Pointer<tbox::Database> rsdb = idb->getDatabase(
            "RefineSchedule");
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
   t_build_supp_mapped_box_level = tbox::TimerManager::getManager()->
      getTimer(
         "xfer::RefineSchedule::finish...()_build_supp_mapped_box_level");
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
   t_build_supp_mapped_box_level.setNull();
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
