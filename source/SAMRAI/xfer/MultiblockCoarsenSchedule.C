/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Coarsening schedule for data transfer between AMR levels 
 *
 ************************************************************************/
#ifndef included_xfer_MultiblockCoarsenSchedule_C
#define included_xfer_MultiblockCoarsenSchedule_C

#include "SAMRAI/xfer/MultiblockCoarsenSchedule.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxGeometry.h"
#include "SAMRAI/hier/BoxOverlap.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/MappedBoxSet.h"
#include "SAMRAI/hier/OverlapConnectorAlgorithm.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchDataFactory.h"
#include "SAMRAI/hier/PatchDescriptor.h"
#include "SAMRAI/hier/PeriodicShiftCatalog.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/xfer/PatchInteriorVariableFillPattern.h"
#include "SAMRAI/xfer/PatchLevelFullFillPattern.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace xfer {

/*
 *************************************************************************
 *                                                                       *
 * Initialization for static data members.                               *
 *                                                                       *
 *************************************************************************
 */

std::string
MultiblockCoarsenSchedule::s_schedule_generation_method = "DLBG";

/*
 * ************************************************************************
 *                                                                        *
 * Static function to set box intersection algorithm for schedules.       *
 *                                                                        *
 * ************************************************************************
 */

void MultiblockCoarsenSchedule::setScheduleGenerationMethod(
   const std::string& method)
{
   if (!((method == "ORIG_NSQUARED") ||
         (method == "DLBG"))) {
      TBOX_ERROR("CoarsenSchedule::setScheduleGenerationMethod\n"
         << "Given method std::string "
         << method << " is invalid.\n Options are\n"
         << "'ORIG_NSQUARED' and 'DLBG'."
         << std::endl);
   }

   s_schedule_generation_method = method;
}

/*
 * ************************************************************************
 *                                                                        *
 * Create a coarsening schedule that transfers data from the source       *
 * patch data components of the fine level into the destination patch     *
 * data components of the coarse level.  If the coarsening operators      *
 * require data in ghost cells on the source level, then those ghost	  *
 * cells must be filled before this call.				  *
 *                                                                        *
 * ************************************************************************
 */
MultiblockCoarsenSchedule::MultiblockCoarsenSchedule(
   tbox::Pointer<hier::PatchLevel> crse_level,
   tbox::Pointer<hier::PatchLevel> fine_level,
   const tbox::Pointer<xfer::CoarsenClasses> coarsen_classes,
   tbox::Pointer<hier::PatchHierarchy> hierarchy,
   tbox::Pointer<xfer::CoarsenTransactionFactory> transaction_factory,
   MultiblockCoarsenPatchStrategy* patch_strategy,
   RefinePatchStrategy* refine_strategy,
   bool fill_coarse_data):
   d_ratio_between_levels(crse_level->getDim())
{
   TBOX_ASSERT(!crse_level.isNull());
   TBOX_ASSERT(!fine_level.isNull());
   TBOX_ASSERT(!coarsen_classes.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS2(*crse_level, *fine_level);

   const tbox::Dimension& dim(crse_level->getDim());

   t_coarsen_data = tbox::TimerManager::getManager()->
      getTimer("xfer::MultiblockCoarsenSchedule::coarsenData()");
   t_gen_sched_n_squared = tbox::TimerManager::getManager()->
      getTimer("xfer::MultiblockCoarsenSchedule::generateScheduleNSquared()");
   t_gen_sched_dlbg = tbox::TimerManager::getManager()->
      getTimer("xfer::MultiblockCoarsenSchedule::generateScheduleDLBG()");
   t_modify_connector = tbox::TimerManager::getManager()->
      getTimer("xfer::MultiblockCoarsenSchedule::modify_connector");
   t_invert_edges = tbox::TimerManager::getManager()->
      getTimer("xfer::MultiblockCoarsenSchedule::generate...()_invert_edges");

   /*
    * Initial values; some may change in setup operations.
    */

   d_mblk_crse_level = crse_level;
   d_mblk_fine_level = fine_level;
   d_transaction_factory = transaction_factory;
   d_mblk_temp_crse_level.setNull();

   d_mblk_coarsen_patch_strategy = patch_strategy;
   d_mblk_refine_strategy = refine_strategy;

   d_mblk_hierarchy = hierarchy;

   d_fill_coarse_data = fill_coarse_data;

   d_schedule.setNull();

   d_mblk_fill_coarse_data_alg.setNull();
   d_mblk_fill_coarse_data_sched.setNull();

   d_number_coarsen_items = 0;
   d_coarsen_items = (const xfer::CoarsenClasses::Data **)NULL;

   /*
    * Compute ratio between fine and coarse levels and then check for
    * correctness.
    */

   hier::IntVector fine(d_mblk_fine_level->getRatioToLevelZero());
   hier::IntVector crse(d_mblk_crse_level->getRatioToLevelZero());
   int i;
   for (i = 0; i < dim.getValue(); i++) {
      if (fine(i) > 1) {
         d_ratio_between_levels(i) = fine(i) / crse(i);
      } else {
         d_ratio_between_levels(i) = tbox::MathUtilities<int>::Abs(crse(
                  i) / fine(i));
      }
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   for (i = 0; i < dim.getValue(); i++) {
      TBOX_ASSERT(d_ratio_between_levels(i) != 0);
   }
   if (dim > tbox::Dimension(1))
      for (i = 0; i < dim.getValue(); i++) {
         if (d_ratio_between_levels(i)
             * d_ratio_between_levels((i + 1) % dim.getValue()) < 0) {
            TBOX_ASSERT((d_ratio_between_levels(i) == 1) ||
               (d_ratio_between_levels((i + 1) % dim.getValue()) == 1));
         }
      }
#endif

   setCoarsenItems(coarsen_classes);
   initialCheckCoarsenClassItems();

   const int nblocks = hierarchy->getGridGeometry()->getNumberBlocks();
   d_coarse_to_fine = new const hier::Connector *[nblocks];
   d_fine_to_coarse = new const hier::Connector *[nblocks];
   d_temp_to_coarse.resizeArray(nblocks, hier::Connector());
   d_coarse_to_temp.resizeArray(nblocks, hier::Connector());

   hier::IntVector min_gcw = getMaxGhostsToGrow();
   for (int b = 0; b < nblocks; b++) {
      d_coarse_to_fine[b] = getOverlapConnector(
            *d_mblk_crse_level->getMappedBoxLevel(),
            *d_mblk_fine_level->getMappedBoxLevel(),
            hier::Connector::convertHeadWidthToBase(d_mblk_crse_level->
               getMappedBoxLevel()->getRefinementRatio(),
               d_mblk_fine_level->getMappedBoxLevel()->getRefinementRatio(),
               min_gcw));
      d_fine_to_coarse[b] = getOverlapConnector(
            *d_mblk_fine_level->getMappedBoxLevel(),
            *d_mblk_crse_level->getMappedBoxLevel(),
            min_gcw);
      TBOX_ASSERT(d_coarse_to_fine[b] != NULL);
      TBOX_ASSERT(d_fine_to_coarse[b] != NULL);
   }

   /*
    * Set up refine schedules to transfer coarsened data and to fill temporary
    * coarse level data before coarsening operations, if needed.  Then,
    * generate communication schedules to transfer data.
    */

   setupRefineAlgorithm();

   generateSchedule();

}

/*
 * ************************************************************************
 *                                                                        *
 * The destructor for the coarsen schedule class implicitly deallocates   *
 * all of the data associated with the communication schedule.            *
 *                                                                        *
 * ************************************************************************
 */

MultiblockCoarsenSchedule::~MultiblockCoarsenSchedule()
{
   clearCoarsenItems();

   d_transaction_factory.setNull();
   d_mblk_crse_level.setNull();
   d_mblk_fine_level.setNull();
   d_mblk_temp_crse_level.setNull();

   d_schedule.setNull();

   d_mblk_fill_coarse_data_alg.setNull();
   d_mblk_fill_coarse_data_sched.setNull();

   t_coarsen_data.setNull();
   t_gen_sched_n_squared.setNull();
   t_gen_sched_dlbg.setNull();

   delete[] d_fine_to_coarse;
   delete[] d_coarse_to_fine;
}

/*
 * ************************************************************************
 *                                                                        *
 * Return const pointer to equivalence classes used in schedule.          *
 *                                                                        *
 * ************************************************************************
 */

const tbox::Pointer<xfer::CoarsenClasses>&
MultiblockCoarsenSchedule::getEquivalenceClasses() const
{
   return d_coarsen_classes;
}

/*
 * ************************************************************************
 *                                                                        *
 * Execute the stored communication schedule that copies data into the	  *
 * the destination patch data components of the destination level from	  *
 * the source patch data components of the source level.  The steps	  *
 * to the algorithm are as follows:					  *
 *									  *
 *	(1) Allocate the source space on the temporary patch level.	  *
 *	(2) Coarsen the data from the fine patch level to the temporary	  *
 *	    patch level (local operation).				  *
 *	(3) Copy data from the source space of the temporary patch	  *
 *	    level into the destination space of the destination patch	  *
 *	    level (requires interprocessor communication).		  *
 *	(4) Deallocate the source space on the temporary patch level.	  *
 *                                                                        *
 * ************************************************************************
 */

void MultiblockCoarsenSchedule::coarsenData() const
{
   t_coarsen_data->start();

   /*
    * Set the coarsen items for all transactions.  These items are
    * shared by all transaction objects in the communication schedule.
    */
   d_transaction_factory->setCoarsenItems(d_coarsen_items,
      d_number_coarsen_items);

   /*
    * Allocate the source data space on the temporary patch level.
    * We do not know the current time, so set it to zero.  It should
    * not matter, since the copy routines do not require that
    * the time markers match.
    */

   d_mblk_temp_crse_level->allocatePatchData(d_sources, 0.0);

   if (d_fill_coarse_data) {
      d_mblk_fill_coarse_data_sched->fillData(0.0);
   }

   /*
    * Coarsen the data from the sources on the fine data level into the
    * sources on the temporary data level
    */

   coarsenSourceData(d_mblk_coarsen_patch_strategy);

   /*
    * Copy data from the source interiors of the temporary patch level
    * into the destination interiors of the destination patch level.
    */

   d_mblk_fill_dst_sched->fillData(0.0, false);

   /*
    * Deallocate the source data in the temporary patch level.
    */

   d_mblk_temp_crse_level->deallocatePatchData(d_sources);

   /*
    * Unset the coarsen items for the copy transactions.  These items
    * are shared by all such transaction objects in the communication
    * schedule.
    */
   d_transaction_factory->unsetCoarsenItems();

   t_coarsen_data->stop();
}

/*
 * ************************************************************************
 *                                                                        *
 * Generate the temporary coarse level by coarsening the fine patch	  *
 * level boxes.  Note that no patch data components are allocated until   *
 * they are needed during the coarsening operation.			  *
 *                                                                        *
 * ************************************************************************
 */

void MultiblockCoarsenSchedule::generateTemporaryLevel()
{
   const tbox::Dimension& dim(d_mblk_crse_level->getDim());

   hier::OverlapConnectorAlgorithm oca;

   const int nblocks = d_mblk_fine_level->getGridGeometry()->getNumberBlocks();

   d_mblk_temp_crse_level = new hier::PatchLevel(dim);
   d_mblk_temp_crse_level->setCoarsenedPatchLevel(d_mblk_fine_level,
      d_ratio_between_levels);
   d_mblk_temp_crse_level->setLevelNumber(d_mblk_crse_level->getLevelNumber());
   d_mblk_temp_crse_level->setNextCoarserHierarchyLevelNumber(
      d_mblk_crse_level->getLevelNumber());

   if (d_fill_coarse_data) {
      d_mblk_hierarchy->getGridGeometry()->adjustMultiblockPatchLevelBoundaries(
         *d_mblk_temp_crse_level);
   }

   for (int nb = 0; nb < nblocks; nb++) {
      if (d_fine_to_coarse[nb] != NULL && d_coarse_to_fine[nb] != NULL) {
         /*
          * Generate temporary mapped_box_level and connectors.
          * Use Connector::modify to modify connectors between fine and coarse
          * mapped_box_levels into connector between temp and coarse mapped_box_levels.
          */

#if 0
         /*
          * Assertion fails for aresamr.
          * I think this assertion it too stringent
          * because getMaxGhostsToGrow returns a bigger
          * width than we need.
          * BTNG.
          */
         /*
          * Make sure the provided connectors have enough gcw
          * to cover the growth of the dst level.
          */
         hier::IntVector dst_growth(getMaxGhostsToGrow());
         TBOX_ASSERT((d_coarse_to_fine[nb]->getConnectorWidth() >= dst_growth));
         hier::IntVector fine_gcw(d_fine_to_coarse[nb]->getConnectorWidth());
         TBOX_ASSERT((fine_gcw >= dst_growth * d_ratio_between_levels));
#endif

         /*
          * Compute d_coarse_to_temp and d_temp_to_coarse.
          *
          * We use the fact that d_mblk_temp_crse_level patches are numbered just
          * like the fine level patches.  The Connectors between coarse and
          * temp are very similar to those between coarse and fine.
          */
         hier::NeighborhoodSet coarse_eto_temp;
         d_coarse_to_fine[nb]->getNeighborhoodSets().coarsenNeighbors(
            coarse_eto_temp,
            d_ratio_between_levels);
         d_coarse_to_temp[nb].initialize(
            d_coarse_to_fine[nb]->getBase(),
            *d_mblk_temp_crse_level->getMappedBoxLevel(),
            d_coarse_to_fine[nb]->getConnectorWidth(),
            coarse_eto_temp,
            hier::MappedBoxLevel::DISTRIBUTED);
         coarse_eto_temp.clear();
         d_temp_to_coarse[nb].initialize(
            *d_mblk_temp_crse_level->getMappedBoxLevel(),
            d_coarse_to_fine[nb]->getBase(),
            d_coarse_to_fine[nb]->getConnectorWidth(),
            d_fine_to_coarse[nb]->getNeighborhoodSets(),
            hier::MappedBoxLevel::DISTRIBUTED);
         const hier::IntVector one_vector(dim, 1);
         oca.shrinkConnectorWidth(d_coarse_to_temp[nb], one_vector);
         oca.shrinkConnectorWidth(d_temp_to_coarse[nb], one_vector);

      }
   }
}

/*
 * ************************************************************************
 *                                                                        *
 * Set up refine algorithms to transfer coarsened data and to fill        *
 * temporary coarse level before performing coarsening operations,        *
 * if necessary.                                                          *
 *                                                                        *
 * ************************************************************************
 */

void MultiblockCoarsenSchedule::setupRefineAlgorithm()
{
   const tbox::Dimension& dim(d_mblk_hierarchy->getDim());

   if (d_fill_coarse_data) {

      d_mblk_fill_coarse_data_alg = new RefineAlgorithm(dim);

      for (int ici = 0; ici < d_number_coarsen_items; ici++) {
         const int src_id = d_coarsen_items[ici]->d_src;
         d_mblk_fill_coarse_data_alg->registerRefine(src_id,
            src_id,
            src_id,
            SAMRAI::tbox::Pointer<SAMRAI::hier::RefineOperator>(NULL));
      }
   }

   tbox::Pointer<PatchInteriorVariableFillPattern> interior_fill(
      new PatchInteriorVariableFillPattern(dim));

   d_mblk_fill_dst_alg = new RefineAlgorithm(dim);

   for (int ici = 0; ici < d_number_coarsen_items; ici++) {
      const int dst_id = d_coarsen_items[ici]->d_dst;
      const int src_id = d_coarsen_items[ici]->d_src;
      d_mblk_fill_dst_alg->registerRefine(dst_id,
         src_id,
         dst_id,
         SAMRAI::tbox::Pointer<SAMRAI::hier::RefineOperator>(NULL),
         interior_fill);
   }

//FIXME
//   d_mblk_fill_dst_alg->disableSingularityPatches();
}

/*
 * ************************************************************************
 *                                                                        *
 * Generate communication schedule that copies source patch data          *
 * from the temporary level into the destination patch data of the        *
 * destination (coarse) level.  The source and destination	          *
 * spaces may be the same.						  *
 *                                                                        *
 * ************************************************************************
 */

void MultiblockCoarsenSchedule::generateSchedule()
{

   /*
    * Set up coarsened version of fine level for temporary data storage.
    * Next, create refine algorithm if needed to fill temporary coarse
    * level before coarsen operations occur.  Then, create empty schedule
    * that will hold transactions for moving data.  Finally, generate
    * schedule based on chosen generation method.
    */
   generateTemporaryLevel();

   if (d_fill_coarse_data) {
      d_mblk_fill_coarse_data_sched =
         d_mblk_fill_coarse_data_alg->createSchedule(d_mblk_temp_crse_level,
            d_mblk_crse_level,
            d_mblk_refine_strategy);
   }

   d_schedule = new tbox::Schedule();
   d_schedule->setTimerPrefix("xfer::MultiblockCoarsenSchedule");

   if (s_schedule_generation_method == "ORIG_NSQUARED") {

      generateScheduleNSquared();

   } else if (s_schedule_generation_method == "DLBG") {

      generateScheduleDLBG();

   } else {

      TBOX_ERROR("Internal MultiblockCoarsenSchedule error..."
         << "\n unrecognized schedule generation option: "
         << s_schedule_generation_method << std::endl);

   }

}

/*
 *************************************************************************
 *                                                                       *
 * This version of the schedule generation procedure uses the original   *
 * SAMRAI N^2 algorithms to construct communication schedules.  Here,    *
 * we loop over all of the patches on the source and destination levels. *
 * check to see whether source or destination is local to this processor.*
 * If not, then skip over schedule construction operations.              *
 *                                                                       *
 *************************************************************************
 */

void MultiblockCoarsenSchedule::generateScheduleNSquared()
{

   t_gen_sched_n_squared->start();

   const int nblocks = d_mblk_crse_level->getGridGeometry()->getNumberBlocks();

   TBOX_ASSERT(nblocks == 1);
   for (int nb = 0; nb < nblocks; nb++) {

      if (d_mblk_coarsen_patch_strategy) {
         d_mblk_coarsen_patch_strategy->setCoarsenBlockNumber(nb);
      }

      const int dst_npatches =
         d_mblk_crse_level->getGlobalNumberOfPatches();
      const int src_npatches =
         d_mblk_temp_crse_level->getGlobalNumberOfPatches();

      const hier::ProcessorMapping& dst_mapping =
         d_mblk_crse_level->getProcessorMapping();
      const hier::ProcessorMapping& src_mapping =
         d_mblk_temp_crse_level->getProcessorMapping();

      hier::BoxList::Iterator crse_itr_dp(d_mblk_crse_level->getBoxes());
      for (int dp = 0; dp < dst_npatches; dp++, crse_itr_dp++) {

         const hier::Box dst_mapped_box(
            *crse_itr_dp,
            hier::LocalId(dp),
            dst_mapping.
            getProcessorAssignment(dp),
            hier::BlockId(nb));

         hier::BoxList::Iterator crse_itr_sp(d_mblk_temp_crse_level->getBoxes());
         for (int sp = 0; sp < src_npatches; sp++, crse_itr_sp++) {

            const hier::Box src_mapped_box(
               *crse_itr_sp,
               hier::LocalId(sp),
               src_mapping.
               getProcessorAssignment(sp),
               hier::BlockId(nb));

            if (dst_mapping.isMappingLocal(dp)
                || src_mapping.isMappingLocal(sp)) {

               constructScheduleTransactions(d_mblk_crse_level,
                  dst_mapped_box,
                  d_mblk_temp_crse_level,
                  src_mapped_box);

            }  // if either source or destination patch is local

         } // loop over source patches

      } // loop over destination patches


      if (d_mblk_coarsen_patch_strategy) {
         d_mblk_coarsen_patch_strategy->clearCoarsenBlockNumber();
      }

   } // loop over blocks in the multiblock domain

   t_gen_sched_n_squared->stop();

}

/*
 *************************************************************************
 *************************************************************************
 */
void MultiblockCoarsenSchedule::generateScheduleDLBG()
{

   t_gen_sched_dlbg->start();

   constructDestinationLevelFillSchedule(d_mblk_crse_level,
                                         d_mblk_temp_crse_level);

   t_gen_sched_dlbg->stop();

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
 * 2. It shifts periodic image dst mapped_boxes back to the zero-shift position,
 * and applies a similar shift to src mapped_boxes so that the overlap is
 * unchanged.  The constructScheduleTransactions method requires all
 * shifts to be absorbed in the src mapped_box.
 ***********************************************************************
 */
void MultiblockCoarsenSchedule::restructureNeighborhoodSetsByDstNodes(
   const hier::Connector& src_to_dst,
   FullNeighborhoodSet& full_inverted_edges) const
{
   const tbox::Dimension& dim(d_mblk_crse_level->getDim());

   const hier::PeriodicShiftCatalog* shift_catalog =
      hier::PeriodicShiftCatalog::getCatalog(dim);
   const hier::NeighborhoodSet& edges = src_to_dst.getNeighborhoodSets();
   const hier::MappedBoxLevel& src_mapped_box_level = src_to_dst.getBase();
   const hier::IntVector& src_ratio(src_to_dst.getBase().getRefinementRatio());
   const hier::IntVector& dst_ratio(src_to_dst.getHead().getRefinementRatio());

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
 **************************************************************************
 * Calculate the max ghost cell width to grow boxes to check for
 * overlaps.  Given in source (fine) level's index space.
 **************************************************************************
 */

hier::IntVector MultiblockCoarsenSchedule::getMaxGhostsToGrow() const
{
   const tbox::Dimension& dim(d_mblk_crse_level->getDim());

   /*
    * Box, face and side elements of adjacent cells overlap even though
    * the cells do not overlap.  Therefore, we always grow at least one
    * cell catch overlaps of mapped_box, face and side elements.
    */
   hier::IntVector gcw(dim, 1);

   for (int ici = 0; ici < d_number_coarsen_items; ici++) {

      /*
       * I don't know why we need to grow by ghost width of src_id.
       * Rich seems to recall needing this but can't remember the
       * reason.  I am disabling it, but may have to re-enable it
       * if it causes problems.
       */
      // const int src_id = d_coarsen_items[ici]->d_src;
      // gcw.max(pd->getPatchDataFactory(src_id)->getDefaultGhostCellWidth());

      hier::IntVector gcw1 = d_coarsen_items[ici]->d_gcw_to_coarsen;
      if (!d_coarsen_items[ici]->d_opcoarsen.isNull()) {
         gcw1 += d_coarsen_items[ici]->d_opcoarsen->getStencilWidth();
      }
      gcw.max(gcw1);
   }

   return gcw;
}

/*
 *************************************************************************
 *                                                                       *
 * Private utility function that constructs schedule transactions that   *
 * move data from source patch on source level to destination patch      *
 * on destination level.                                                 *
 *                                                                       *
 *************************************************************************
 */

void MultiblockCoarsenSchedule::constructScheduleTransactions(
   tbox::Pointer<hier::PatchLevel> dst_level,
   const hier::Box& dst_mapped_box,
   tbox::Pointer<hier::PatchLevel> src_level,
   const hier::Box& src_mapped_box)
{
   TBOX_ASSERT(!dst_level.isNull());
   TBOX_ASSERT(!src_level.isNull());

   const tbox::Dimension& dim(d_mblk_crse_level->getDim());

   TBOX_DIM_ASSERT_CHECK_DIM_ARGS4(dim,
      *dst_level,
      *src_level,
      dst_mapped_box,
      src_mapped_box);

   const hier::IntVector& constant_zero_intvector(hier::IntVector::getZero(dim));
   const hier::IntVector& constant_one_intvector(hier::IntVector::getOne(dim));

   tbox::Pointer<hier::PatchDescriptor> dst_patch_descriptor =
      dst_level->getPatchDescriptor();
   tbox::Pointer<hier::PatchDescriptor> src_patch_descriptor =
      src_level->getPatchDescriptor();

   const hier::Box& dst_box = dst_mapped_box;
   const hier::Box& src_box = src_mapped_box;

   const int num_equiv_classes =
      d_coarsen_classes->getNumberOfEquivalenceClasses();

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
      src_shift *= src_level->getRatioToLevelZero();
      unshifted_src_box.shift(-src_shift);
   }
   if (dst_mapped_box.isPeriodicImage()) {
      TBOX_ASSERT(!src_mapped_box.isPeriodicImage());
      dst_shift = shift_catalog->shiftNumberToShiftDistance(
            dst_mapped_box.getPeriodicId());
      dst_shift *= dst_level->getRatioToLevelZero();
      unshifted_dst_box.shift(-dst_shift);
   }

   const int num_coarsen_items = d_coarsen_classes->getNumberOfCoarsenItems();
   tbox::Array<tbox::Pointer<tbox::Transaction> >
   transactions(num_coarsen_items);

   for (int nc = 0; nc < num_equiv_classes; nc++) {

      const xfer::CoarsenClasses::Data& rep_item =
         d_coarsen_classes->getClassRepresentative(nc);

      const int rep_item_dst_id = rep_item.d_dst;
      const int rep_item_src_id = rep_item.d_src;

      tbox::Pointer<hier::PatchDataFactory> src_pdf =
         src_patch_descriptor->getPatchDataFactory(rep_item_src_id);
      tbox::Pointer<hier::PatchDataFactory> dst_pdf =
         dst_patch_descriptor->getPatchDataFactory(rep_item_dst_id);

      const hier::IntVector& dst_gcw(dst_pdf->getGhostCellWidth());

      const hier::Box dst_fill_box(hier::Box::grow(unshifted_dst_box, dst_gcw));

      hier::Box test_mask = dst_fill_box * src_mapped_box;
      if (test_mask.empty() &&
          (dst_gcw == constant_zero_intvector) &&
          dst_pdf->dataLivesOnPatchBorder()) {
         hier::Box tmp_dst_fill_box(
            hier::Box::grow(dst_fill_box, constant_one_intvector));
         test_mask = tmp_dst_fill_box * src_mapped_box;
      }
      hier::Box src_mask = hier::Box::shift(test_mask, -src_shift);

      test_mask = hier::Box::grow(unshifted_src_box,
            hier::IntVector::min(
               rep_item.d_gcw_to_coarsen,
               src_pdf->getGhostCellWidth()));

      src_mask += test_mask;

      hier::Transformation transformation(src_shift);
      tbox::Pointer<hier::BoxOverlap> overlap =
         rep_item.d_var_fill_pattern->calculateOverlap(
            *dst_pdf->getBoxGeometry(unshifted_dst_box),
            *src_pdf->getBoxGeometry(unshifted_src_box),
            dst_mapped_box,
            src_mask,
            dst_fill_box,
            true, transformation);

      if (overlap.isNull()) {
         TBOX_ERROR("Internal MultiblockCoarsenSchedule error..."
            << "\n Overlap is NULL for "
            << "\n src box = " << src_box
            << "\n dst box = " << dst_box
            << "\n src mask = " << src_mask << std::endl);
      }

      if (!overlap->isOverlapEmpty()) {
         for (tbox::List<int>::Iterator l(d_coarsen_classes->getIterator(nc));
              l; l++) {
            const CoarsenClasses::Data& item =
               d_coarsen_classes->getCoarsenItem(l());
            TBOX_ASSERT(item.d_class_id == nc);

            const int citem_count = item.d_tag;
            transactions[citem_count] =
               d_transaction_factory->allocate(dst_level,
                  src_level,
                  overlap,
                  dst_mapped_box,
                  src_mapped_box,
                  citem_count);
         }
      }

   }  // iterate over all coarsen equivalence classes

   for (int i = 0; i < num_coarsen_items; i++) {
      if (!(transactions[i].isNull())) {
         d_schedule->appendTransaction(transactions[i]);
      }
   }
}


/*
 *************************************************************************
 *                                                                       *
 * Private utility function that constructs a schedule to move data      *
 * move data from the temporary coarse level to the destination level.   *
 *                                                                       *
 *************************************************************************
 */
void MultiblockCoarsenSchedule::constructDestinationLevelFillSchedule(
   tbox::Pointer<hier::PatchLevel> dst_level,
   tbox::Pointer<hier::PatchLevel> src_level)
{
   TBOX_ASSERT(!dst_level.isNull());
   TBOX_ASSERT(!src_level.isNull());

#ifdef DEBUG_CHECK_ASSERTIONS
   const tbox::Dimension& dim(d_mblk_crse_level->getDim());
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS2(dim,
      *dst_level,
      *src_level);
#endif

   tbox::Pointer<PatchLevelFullFillPattern> level_fill(
      new PatchLevelFullFillPattern());

   d_mblk_fill_dst_sched = d_mblk_fill_dst_alg->createSchedule(
                              level_fill,
                              dst_level,
                              src_level,
                              d_mblk_refine_strategy);
}


/*
 * ************************************************************************
 *                                                                        *
 * Coarsen data from the source space on the fine patch level into the	  *
 * source space on the coarse temporary patch level.			  *
 *                                                                        *
 * ************************************************************************
 */
void MultiblockCoarsenSchedule::coarsenSourceData(
   xfer::MultiblockCoarsenPatchStrategy* patch_strategy) const
{
   const hier::MappedBoxSet& fine_level_boxes =
      d_mblk_fine_level->getMappedBoxLevel()->getMappedBoxes();

   const int nblocks = d_mblk_fine_level->getGridGeometry()->getNumberBlocks();
   for (int nb = 0; nb < nblocks; nb++) {

      if (d_mblk_coarsen_patch_strategy) {
         d_mblk_coarsen_patch_strategy->setCoarsenBlockNumber(nb);
      }

      if (fine_level_boxes.size()) {

         hier::BlockId block_id(nb);
         hier::MappedBoxSetSingleBlockIterator crse_iter(
            d_mblk_temp_crse_level->getMappedBoxLevel()->getMappedBoxes(),
            block_id);

         if (crse_iter.isValid()) {

            hier::MappedBoxSetSingleBlockIterator fine_iter(
               fine_level_boxes, block_id);

            /*
             * Loop over all local patches (fine and temp have the same mapping)
             */
               for( ; fine_iter.isValid(); fine_iter++) {
                  const hier::BoxId& mapped_box_id = fine_iter->getId();
               tbox::Pointer<hier::Patch> fine_patch(
                  d_mblk_fine_level->getPatch(mapped_box_id));
               tbox::Pointer<hier::Patch> temp_patch(
                  d_mblk_temp_crse_level->getPatch(mapped_box_id));

               const hier::Box& box = temp_patch->getBox();

               /*
                * Coarsen the fine space onto the temporary coarse space
                */

               if (patch_strategy) {
                  patch_strategy->preprocessCoarsen(*temp_patch,
                     *fine_patch, box, d_ratio_between_levels);
               }

               for (int ici = 0; ici < d_number_coarsen_items; ici++) {
                  const xfer::CoarsenClasses::Data * const crs_item =
                     d_coarsen_items[ici];
                  if (!(crs_item->d_opcoarsen.isNull())) {
                     const int source_id = crs_item->d_src;
                     crs_item->d_opcoarsen->coarsen(*temp_patch, *fine_patch,
                        source_id, source_id,
                        box, d_ratio_between_levels);
                  }
               }

               if (patch_strategy) {
                  patch_strategy->postprocessCoarsen(*temp_patch,
                     *fine_patch,
                     box,
                     d_ratio_between_levels);
               }
            }
         }
      }

      if (d_mblk_coarsen_patch_strategy) {
         d_mblk_coarsen_patch_strategy->clearCoarsenBlockNumber();
      }

   }
}

/*
 * ***********************************************************************
 *                                                                       *
 * Private utility function to set up local array of coarsen items.      *
 *                                                                       *
 * ***********************************************************************
 */

void MultiblockCoarsenSchedule::setCoarsenItems(
   const tbox::Pointer<xfer::CoarsenClasses> coarsen_classes)
{

   clearCoarsenItems();

   d_coarsen_classes = coarsen_classes;

   d_number_coarsen_items = d_coarsen_classes->getNumberOfCoarsenItems();

   /*
    * Determine total number of coarsen items and set state of
    * component selector used to manage storage on temporary level.
    */
   d_sources.clrAllFlags();

   int nc;
   for (nc = 0; nc < d_number_coarsen_items; nc++) {
      const CoarsenClasses::Data& item = d_coarsen_classes->getCoarsenItem(nc);
      d_sources.setFlag(item.d_src);
   }

   /*
    * Allocate and initialize array of coarsen items.
    */

   d_coarsen_items =
      new const xfer::CoarsenClasses::Data *[d_number_coarsen_items];

   int ircount = 0;
   for (nc = 0; nc < d_number_coarsen_items; nc++) {
      d_coarsen_classes->getCoarsenItem(nc).d_tag = ircount;
      d_coarsen_items[ircount] = &(d_coarsen_classes->getCoarsenItem(nc));
      ircount++;
   }

}

/*
 * ***********************************************************************
 *                                                                       *
 * Private utility function to check coarsen items in initial setup to   *
 * see whether source and destination patch data components have         *
 * sufficient ghost cell widths to satisfy the "ghost width to coarsen"  *
 * functionality described in the CoarsenAlgorithm class header.   *
 * Specifically, the destination data must have a ghost cell width at    *
 * least as large as the ghost cell width to coarsen.  The source data   *
 * must have a ghost cell width at least as large as the ghost cell      *
 * width to coarsen refined to the source (finer) level index space.     *
 * Other checks are also performed here by calling the                   *
 * CoarsenClasses::checkCoarsenItem() routine.                     *
 *                                                                       *
 * ***********************************************************************
 */

void MultiblockCoarsenSchedule::initialCheckCoarsenClassItems() const
{
   const tbox::Dimension& dim(d_mblk_crse_level->getDim());

   tbox::Pointer<hier::PatchDescriptor> pd =
      hier::VariableDatabase::getDatabase()->getPatchDescriptor();

   hier::IntVector user_gcw(dim, 0);
   if (d_mblk_coarsen_patch_strategy) {
      user_gcw = d_mblk_coarsen_patch_strategy->getCoarsenOpStencilWidth();
   }

   for (int ici = 0; ici < d_number_coarsen_items; ici++) {

      const xfer::CoarsenClasses::Data * const crs_item = d_coarsen_items[ici];

#ifdef DEBUG_CHECK_ASSERTIONS
      if (d_coarsen_classes->checkCoarsenItem(*crs_item, pd)) {
#endif

      const int dst_id = crs_item->d_dst;
      const int src_id = crs_item->d_src;

      tbox::Pointer<hier::PatchDataFactory> dfact =
         pd->getPatchDataFactory(dst_id);
      tbox::Pointer<hier::PatchDataFactory> sfact =
         pd->getPatchDataFactory(src_id);

      const hier::IntVector& dst_gcw(dfact->getGhostCellWidth());
      const hier::IntVector& src_gcw(sfact->getGhostCellWidth());

      if (crs_item->d_gcw_to_coarsen > dst_gcw) {
         TBOX_ERROR("Bad data given to MultiblockCoarsenSchedule...\n"
            << "`ghost cell width to coarsen' specified in\n"
            << "registration of `Destination' patch data "
            << pd->mapIndexToName(dst_id)
            << " with CoarsenAlgorithm\n"
            << " is larger than ghost cell width of data \n"
            << "d_gcw_to_coarsen = " << crs_item->d_gcw_to_coarsen
            << "\n data ghost cell width = " << dst_gcw << std::endl);
      }

      if ((crs_item->d_gcw_to_coarsen * d_ratio_between_levels) > src_gcw) {
         TBOX_ERROR("Bad data given to MultiblockCoarsenSchedule...\n"
            << "`Source' patch data " << pd->mapIndexToName(src_id)
            << " has ghost cell width too small to support the\n"
            << "`ghost cell width to coarsen' specified in"
            << " registration with CoarsenAlgorithm\n"
            << "data ghost cell width = " << src_gcw
            << "d_gcw_to_coarsen = " << crs_item->d_gcw_to_coarsen
            << "\nratio between levels = " << d_ratio_between_levels
            << "\n Thus, data ghost width must be >= "
            << (crs_item->d_gcw_to_coarsen * d_ratio_between_levels)
            << std::endl);
      }

      if (user_gcw > src_gcw) {
         TBOX_ERROR("Bad data given to MultiblockCoarsenSchedule...\n"
            << "User supplied coarsen stencil width = "
            << user_gcw
            << "\nis larger than ghost cell width of `Source'\n"
            << "patch data " << pd->mapIndexToName(src_id)
            << " , which is " << src_gcw << std::endl);
      }

#ifdef DEBUG_CHECK_ASSERTIONS
   }
#endif

   }

}

/*
 * ************************************************************************
 *                                                                        *
 * Private utility function to clear array of coarsen items.              *
 *                                                                        *
 * ************************************************************************
 */

void MultiblockCoarsenSchedule::clearCoarsenItems()
{
   if (d_coarsen_items) {
      for (int ici = 0; ici < d_number_coarsen_items; ici++) {
         d_coarsen_items[ici] = (xfer::CoarsenClasses::Data *)NULL;
      }
      delete[] d_coarsen_items;
      d_coarsen_items = (const xfer::CoarsenClasses::Data **)NULL;
      d_number_coarsen_items = 0;
   }
}

/*
 **************************************************************************
 * Private utility function to access overlap connectors between two
 * MappedBoxLevels.  The Connector must be appropriately registered and
 * be unique.
 **************************************************************************
 */

const hier::Connector *MultiblockCoarsenSchedule::getOverlapConnector(
   const hier::MappedBoxLevel& base,
   const hier::MappedBoxLevel& head,
   const hier::IntVector& min_gcw) const
{
   const hier::Connector* found =
      &base.getPersistentOverlapConnectors().findOrCreateConnector(head,
         min_gcw);
   if (found != NULL) {
      TBOX_ASSERT(found->getConnectorWidth() >= min_gcw);
   }
   return found;
}

/*
 **************************************************************************
 * Private utility function to access overlap connectors between two
 * MappedBoxLevels.  The Connector must be appropriately registered and
 * be unique.
 **************************************************************************
 */

const hier::Connector *MultiblockCoarsenSchedule::getOverlapConnector_strict(
   const hier::MappedBoxLevel& base,
   const hier::MappedBoxLevel& head,
   const hier::IntVector& min_gcw) const
{
   const hier::Connector* found =
      &base.getPersistentOverlapConnectors().findConnector(head, min_gcw);
   TBOX_ASSERT(found != NULL);
   TBOX_ASSERT(found->getConnectorWidth() >= min_gcw);
   return found;
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
