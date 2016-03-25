//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-1/source/transfer/datamovers/locally_active/LocallyActiveDataCoarsenSchedule.C $
// Package:	SAMRAI data transfer
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1846 $
// Modified:	$LastChangedDate: 2008-01-11 09:51:05 -0800 (Fri, 11 Jan 2008) $
// Description:	Coarsening schedule for locally-active data transfer between AMR levels
//

#ifndef included_xfer_LocallyActiveDataCoarsenSchedule_C
#define included_xfer_LocallyActiveDataCoarsenSchedule_C 

#include "LocallyActiveDataCoarsenSchedule.h"

#include "Box.h"
#include "BoxGeometry.h"
#include "BoxGraph.h"
#include "BoxOverlap.h"
#include "BoxTree.h"
#include "LocallyActiveVariableDatabase.h"
#include "Patch.h"
#include "PatchDataFactory.h"
#include "PatchDescriptor.h"
#include "ProcessorMapping.h"
#include "tbox/ArenaManager.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"


#ifndef NULL
#define NULL (0)
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

template<int DIM> const hier::IntVector<DIM>
   LocallyActiveDataCoarsenSchedule<DIM>::s_constant_zero_intvector = 
      hier::IntVector<DIM>(0);
template<int DIM> const hier::IntVector<DIM>
   LocallyActiveDataCoarsenSchedule<DIM>::s_constant_one_intvector = 
      hier::IntVector<DIM>(1);
template<int DIM> std::string
   LocallyActiveDataCoarsenSchedule<DIM>::s_schedule_generation_method = 
      "BOX_TREE";

/*
 * ************************************************************************
 *                                                                        *
 * Static function to set box intersection algorithm for schedules.       *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM>
void LocallyActiveDataCoarsenSchedule<DIM>::setScheduleGenerationMethod(
   const std::string& method)
{
   if ( !((method == "ORIG_NSQUARED") ||
          (method == "BOX_GRAPH") ||
          (method == "BOX_TREE")) ) {
      TBOX_ERROR("xfer::LocallyActiveDataCoarsenSchedule<DIM>::setScheduleGenerationMethod\n"
                 << "Given method string "
                 << method << " is invalid.\n Options are\n"
                 << "'ORIG_NSQUARED', 'BOX_GRAPH', and 'BOX_TREE'."
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
 * require data in ghost cells on the source level, then those ghost      *
 * cells must be filled before this call.                                 *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM>
LocallyActiveDataCoarsenSchedule<DIM>::LocallyActiveDataCoarsenSchedule(
   tbox::Pointer< hier::PatchLevel<DIM> > crse_level,
   tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > crse_level_mgr,
   tbox::Pointer< hier::PatchLevel<DIM> > fine_level,
   tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > fine_level_mgr,
   const tbox::Pointer< xfer::CoarsenClasses<DIM> > coarsen_classes,
   tbox::Pointer< xfer::LocallyActiveDataCoarsenTransactionFactory<DIM> > 
      transaction_factory, 
   xfer::LocallyActiveDataCoarsenPatchStrategy<DIM>* patch_strategy,
   bool fill_coarse_data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!crse_level.isNull());
   TBOX_ASSERT(!crse_level_mgr.isNull());
   TBOX_ASSERT(!fine_level.isNull());
   TBOX_ASSERT(!fine_level_mgr.isNull());
   TBOX_ASSERT(!coarsen_classes.isNull());
   TBOX_ASSERT(!transaction_factory.isNull());
#endif

   t_coarsen_data = tbox::TimerManager::getManager() -> 
      getTimer("xfer::LocallyActiveCoarsenSchedule::coarsenData()");
   t_gen_sched_n_squared = tbox::TimerManager::getManager()->
      getTimer("xfer::LocallyActiveCoarsenSchedule::generateScheduleNSquared()");
   t_gen_sched_box_graph = tbox::TimerManager::getManager()->
      getTimer("xfer::LocallyActiveCoarsenSchedule::generateScheduleBoxGraph()");
   t_gen_sched_box_tree = tbox::TimerManager::getManager()->
      getTimer("xfer::LocallyActiveCoarsenSchedule::generateScheduleBoxTree()");

   /*
    * Initial values; some may change in setup operations.
    */

   d_transaction_factory = transaction_factory;

   d_crse_level             = crse_level;
   d_crse_level_mgr         = crse_level_mgr;
   d_fine_level             = fine_level;
   d_fine_level_mgr         = fine_level_mgr;
   d_temp_crse_level.setNull();
   d_temp_level_active_data_mgr.setNull();

   d_coarsen_patch_strategy = patch_strategy;

   d_fill_coarse_data       = fill_coarse_data;

   d_schedule.setNull();

#if 0
// These are here as a reminder for consistency with regular CoarsenSchedule.
// The details of how refine operations should be used here have note been worked out yet.
   d_precoarsen_refine_algorithm.setNull();
   d_precoarsen_refine_schedule.setNull();
#endif

   d_number_coarsen_items = 0;
   d_coarsen_items = (const typename xfer::CoarsenClasses<DIM>::Data**)NULL;

   /*
    * Compute ratio between fine and coarse levels and then check for
    * correctness.
    */

   hier::IntVector<DIM> fine = d_fine_level->getRatio();
   hier::IntVector<DIM> crse = d_crse_level->getRatio();
   int i;
   for (i = 0; i < DIM; i++) {
      if (fine(i) > 1) {
         d_ratio_between_levels(i) = fine(i) / crse(i);
      } else {
         d_ratio_between_levels(i) = 
            tbox::MathUtilities<int>::Abs( crse(i) / fine(i) );
      }
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   for (i = 0; i < DIM; i++) {
      TBOX_ASSERT( d_ratio_between_levels(i) != 0 );
   }
   if (DIM > 1) { 
      for (i = 0; i < DIM; i++) {
         if ( d_ratio_between_levels(i)*d_ratio_between_levels((i+1)%DIM) < 0 ) {
            TBOX_ASSERT( (d_ratio_between_levels(i) == 1) ||
                     (d_ratio_between_levels((i+1)%DIM) == 1) );
         }
      }
   }
#endif

   bool called_from_reset = false;
   setCoarsenItems(coarsen_classes, called_from_reset);
   initialCheckCoarsenClassItems();

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

template<int DIM>
LocallyActiveDataCoarsenSchedule<DIM>::~LocallyActiveDataCoarsenSchedule()
{
   clearCoarsenItems();

   d_transaction_factory.setNull();
   d_temp_crse_level.setNull();
   d_temp_level_active_data_mgr.setNull();
   d_crse_level.setNull();
   d_fine_level.setNull();
   d_crse_level_mgr.setNull();
   d_fine_level_mgr.setNull();

   d_schedule.setNull();

#if 0
// These are here as a reminder for consistency with regular CoarsenSchedule.
// The details of how refine operations should be used here have note been worked out yet.
   d_precoarsen_refine_algorithm.setNull();
   d_precoarsen_refine_schedule.setNull();
#endif

   t_coarsen_data.setNull();
   t_gen_sched_n_squared.setNull();
   t_gen_sched_box_graph.setNull();
   t_gen_sched_box_tree.setNull();
}

/*
 * ************************************************************************
 *                                                                        *
 * Execute the stored communication schedule that copies data into the    *
 * the destination patch data components of the destination level from    *
 * the source patch data components of the source level.  The steps       *
 * to the algorithm are as follows:                                       *
 *                                                                        *
 *      (1) Allocate the source space on the temporary patch level.       *
 *      (2) Coarsen the data from the fine patch level to the temporary   *
 *          patch level (local operation).                                *
 *      (3) Copy data from the source space of the temporary patch        *
 *          level into the destination space of the destination patch     *
 *          level (requires interprocessor communication).                *
 *      (4) Deallocate the source space on the temporary patch level.     *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM>
void LocallyActiveDataCoarsenSchedule<DIM>::coarsenData() const
{
   t_coarsen_data->start();

   /*
    * Set the coarsen items for all transactions.  These items are
    * shared by all transaction objects in the communication schedule.
    */

   d_transaction_factory->setCoarsenItems(d_coarsen_items, d_number_coarsen_items);

   /*
    * Allocate the source data space on the temporary patch level.
    * We do not know the current time, so set it to zero.  It should
    * not matter, since the copy routines do not require that
    * the time markers match.
    */

   d_temp_level_active_data_mgr->
      allocateAllPatchData(0.0, 
                           tbox::ArenaManager::getManager()->getScratchAllocator());   

   if (d_fill_coarse_data) {
#if 0
// This is here as a reminder for consistency with regular CoarsenSchedule.
// The details of how refine operations should be used here have note been worked out yet.
      d_precoarsen_refine_schedule->fillData(0.0);
#endif
   }

   /*
    * Coarsen the data from the sources on the fine data level into the
    * sources on the temporary data level
    */   

   coarsenSourceData(d_coarsen_patch_strategy);

   /*
    * Copy data from the source interiors of the temporary patch level
    * into the destination interiors of the destination patch level.
    */

   d_schedule->communicate();

   /*
    * Deallocate the source data in the temporary patch level.
    */

   d_temp_level_active_data_mgr->deallocateAllPatchData();

   /*
    * Unset the coarsen items for the transactions.  These items are 
    * shared by all such transaction objects in the communication
    * schedule.
    */

   d_transaction_factory->unsetCoarsenItems();

   t_coarsen_data->stop();
}

/*
 * ************************************************************************
 *                                                                        *
 * Generate the temporary coarse level by coarsening the fine patch       *
 * level boxes.  Note that no patch data components are allocated until   *
 * they are needed during the coarsening operation.  However, the         *
 * required source data on the temporary level is recorded so that it     *
 * may be allocated when needed.                                          *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM>
void LocallyActiveDataCoarsenSchedule<DIM>::generateTemporaryLevel()
{
   d_temp_crse_level = new hier::PatchLevel<DIM>();
   d_temp_crse_level->setCoarsenedPatchLevel(d_fine_level,
                                             d_ratio_between_levels);
   d_temp_crse_level->setLevelNumber(d_crse_level->getLevelNumber());
   d_temp_crse_level->setNextCoarserHierarchyLevelNumber(
                      d_crse_level->getLevelNumber());

   resetTempCoarseLevelDataManager();
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

template<int DIM>
void LocallyActiveDataCoarsenSchedule<DIM>::setupRefineAlgorithm()
{

   if (d_fill_coarse_data) {
#if 0
// This is here as a reminder for consistency with regular CoarsenSchedule.
// The details of how refine operations should be used here have note been worked out yet.
      d_precoarsen_refine_algorithm.setNull();
      d_precoarsen_refine_algorithm = 
         new xfer::LocallyActiveDataRefineAlgorithm<DIM>();

      for (int ici = 0; ici < d_number_coarsen_items; ici++) {
         const int src_id = d_coarsen_items[ici]->d_src;
         d_precoarsen_refine_algorithm->registerRefine(src_id,
                                                       src_id,
                                                       src_id,
                                                       NULL);
      }
#endif
   }

}

/*
 * ************************************************************************
 *                                                                        *
 * Generate communication schedule that copies source patch data          *
 * from the temporary level into the destination patch data of the        *
 * destination (coarse) level.  The source and destination                *
 * spaces may be the same.                                                *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM>
void LocallyActiveDataCoarsenSchedule<DIM>::generateSchedule()
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
#if 0
// This is here as a reminder for consistency with regular CoarsenSchedule.
// The details of how refine operations should be used here have note been worked out yet.
      d_precoarsen_refine_schedule =
         d_precoarsen_refine_algorithm->createSchedule(d_temp_crse_level,
                                                       d_crse_level);
#endif
   }

   d_schedule = new tbox::Schedule();

   if (s_schedule_generation_method == "ORIG_NSQUARED") {

       generateScheduleNSquared();

   } else if (s_schedule_generation_method == "BOX_GRAPH") {

       generateScheduleBoxGraph();

   } else if (s_schedule_generation_method == "BOX_TREE") {

       generateScheduleBoxTree();

   } else {

      TBOX_ERROR("Internal LocallyActiveDataCoarsenSchedule<DIM> error..."
                 << "\n unrecognized schedule generation option: "
                 << s_schedule_generation_method << std::endl);

   }

}

/*
*************************************************************************
*                                                                       *
* This version of the schedule generation procedure uses the original   *
* SAMRAI N^2 algorithms to construct communication schedules.   Here,   *
* we loop over all of the patches on the source and destination levels. *
* check to see whether source or destination is local to this processor.*
* If not, then skip over schedule construction operations.              *
*                                                                       *
*************************************************************************
*/

template<int DIM>
void LocallyActiveDataCoarsenSchedule<DIM>::generateScheduleNSquared()
{
   t_gen_sched_n_squared->start();

   const int dst_npatches = d_crse_level->getNumberOfPatches();
   const int src_npatches = d_temp_crse_level->getNumberOfPatches();

   const hier::ProcessorMapping& dst_mapping =
      d_crse_level->getProcessorMapping();
   const hier::ProcessorMapping& src_mapping =
      d_temp_crse_level->getProcessorMapping();

   for (int dp = 0; dp < dst_npatches; dp++) {

      for (int sp = 0; sp < src_npatches; sp++) {

         if (dst_mapping.isMappingLocal(dp)
             || src_mapping.isMappingLocal(sp)) {

            constructScheduleTransactions(
               d_crse_level, d_crse_level_mgr, dp,
               d_temp_crse_level, d_temp_level_active_data_mgr, sp);

         }  // if either source or destination patch is local

      } // loop over source patches

   } // loop over destination patches

   t_gen_sched_n_squared->stop();
}

/*
*************************************************************************
*                                                                       *
* This version of the schedule generation procedure uses a bipartite    *
* graph algorithm to determine which source patches contribute data to  *
* each destination patch.   Once we have generated the graph, we        *
* only perform schedule construction operations between source and      *
* destination patches where one is local to processor and their         *
* overlap is non-empty.                                                 *
*                                                                       *
*************************************************************************
*/

template<int DIM>
void LocallyActiveDataCoarsenSchedule<DIM>::generateScheduleBoxGraph()
{
   t_gen_sched_box_graph->start();

   hier::IntVector<DIM> dst_growth = getMaxGhostsToGrow();

   tbox::Pointer< hier::BoxGraph<DIM> > box_graph =
      new hier::BoxGraph<DIM>(d_temp_crse_level->getBoxes(),
                              d_temp_crse_level->getShiftsForLevel(),
                              d_temp_crse_level->getProcessorMapping(),
                              d_crse_level->getBoxes(),
                              dst_growth);

   const int dst_npatches = d_crse_level->getNumberOfPatches();
   const hier::ProcessorMapping& dst_mapping = d_crse_level->getProcessorMapping();

   for (int dp = 0; dp < dst_npatches; dp++) {

      tbox::Array<int> src_nabor_indices;
      if (dst_mapping.isMappingLocal(dp)) {
         src_nabor_indices = box_graph->getSrcOverlapIndices(dp);
      } else {
         src_nabor_indices = box_graph->getLocalSrcOverlapIndices(dp);
      }

      int src_len = src_nabor_indices.getSize();
      for (int spp = 0; spp < src_len; spp++) {

         int sp = src_nabor_indices[spp];

         constructScheduleTransactions(
            d_crse_level, d_crse_level_mgr, dp,
            d_temp_crse_level, d_temp_level_active_data_mgr, sp);

      } // loop over source patches

   } // loop over destination patches

   t_gen_sched_box_graph->stop();
}

/*
*************************************************************************
*                                                                       *
* This version of the schedule generation procedure uses a recursive    *
* binary box tree algorithm to determine which source patches           *
* contribute data to each destination patch.   Once we have generated   *
* the graph, we only perform schedule construction operations between   *
* source and destination patches where one is local to processor and    *
* their overlap is non-empty.                                           *
*                                                                       *
*************************************************************************
*/

template<int DIM>
void LocallyActiveDataCoarsenSchedule<DIM>::generateScheduleBoxTree()
{
   t_gen_sched_box_tree->start();

   tbox::Pointer< hier::BoxTree<DIM> > box_tree = d_temp_crse_level->getBoxTree();

   hier::IntVector<DIM> dst_growth = getMaxGhostsToGrow();

   const int dst_npatches = d_crse_level->getNumberOfPatches();
   const hier::BoxArray<DIM>& dst_boxes = d_crse_level->getBoxes();
   const hier::ProcessorMapping& dst_mapping = d_crse_level->getProcessorMapping();

   for (int dp = 0; dp < dst_npatches; dp++) {

      const hier::Box<DIM>& dst_box = dst_boxes(dp);
      hier::Box<DIM> dst_box_plus_ghosts = dst_box;
      dst_box_plus_ghosts.grow(dst_growth);

      tbox::Array<int> src_nabor_indices;
      if (dst_mapping.isMappingLocal(dp)) {
         box_tree->findOverlapIndices(
            src_nabor_indices, dst_box_plus_ghosts);
      } else {
         box_tree->findLocalOverlapIndices(
            src_nabor_indices, dst_box_plus_ghosts);
      }

      int src_len = src_nabor_indices.getSize();
      for (int spp = 0; spp < src_len; spp++) {

         int sp = src_nabor_indices[spp];

         constructScheduleTransactions(
            d_crse_level, d_crse_level_mgr, dp,
            d_temp_crse_level, d_temp_level_active_data_mgr, sp);

      } // loop over source patches

   } // loop over destination patches

   t_gen_sched_box_tree->stop();
}

/*
**************************************************************************
*                                                                        *
* Calculate the max ghost cell width to grow boxes to check for overlaps.*
*                                                                        *
**************************************************************************
*/

template<int DIM>
hier::IntVector<DIM> 
LocallyActiveDataCoarsenSchedule<DIM>::getMaxGhostsToGrow() const
{
   hier::IntVector<DIM> gcw(0);
   tbox::Pointer< hier::PatchDescriptor<DIM> > pd = 
      d_temp_crse_level->getPatchDescriptor();

   for (int ici = 0; ici < d_number_coarsen_items; ici++) {
      const int src_id = d_coarsen_items[ici]->d_src;
      gcw.max(pd->getPatchDataFactory(src_id)->getGhostCellWidth());
      gcw.max(d_coarsen_items[ici]->d_gcw_to_coarsen);
   }

   hier::IntVector<DIM> growth(1);
   growth.max(gcw);
   return(growth);
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

template<int DIM>
void LocallyActiveDataCoarsenSchedule<DIM>::constructScheduleTransactions(
   tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
   tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > dst_level_mgr,
   int dst_patch_id,
   tbox::Pointer< hier::PatchLevel<DIM> > src_level,
   tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > src_level_mgr,
   int src_patch_id)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst_level.isNull());
   TBOX_ASSERT(!src_level.isNull());
#endif

   const int num_equiv_classes =
      d_coarsen_classes->getNumberEquivalenceClasses();

   int nc;

   /*
    * Determine which equivalence classes are active; i.e., which
    * contain data items to be copied from given source patch to
    * given destination patch.  For an equivalence class to be active,
    * both the source data must be active on the source patch and the
    * destination data must be active on the destination patch.
    */

   bool any_equivalence_classes_active = false;

   tbox::Array<bool> active_equivalence_class(num_equiv_classes);
   for (nc = 0; nc < num_equiv_classes; nc++) {
      active_equivalence_class[nc] = false;
      for (typename tbox::List<typename xfer::CoarsenClasses<DIM>::Data>::Iterator 
           l(d_coarsen_classes->getIterator(nc)); l; l++) {
         if ( src_level_mgr->getPatchDataActive( hier::PatchDataId(l().d_src), 
                                                 hier::PatchNumber(src_patch_id) ) &&
              dst_level_mgr->getPatchDataActive( hier::PatchDataId(l().d_dst), 
                                                 hier::PatchNumber(dst_patch_id) ) ) {
            active_equivalence_class[nc] = true;
         }
      }
      any_equivalence_classes_active |= active_equivalence_class[nc];
   }

   if (any_equivalence_classes_active) {

      tbox::Pointer< hier::PatchDescriptor<DIM> > dst_patch_descriptor =
         dst_level->getPatchDescriptor();
      tbox::Pointer< hier::PatchDescriptor<DIM> > src_patch_descriptor =
         src_level->getPatchDescriptor();

      const hier::Box<DIM>& dst_box = dst_level->getBoxes().getBox(dst_patch_id);
      const hier::Box<DIM>& src_box = src_level->getBoxes().getBox(src_patch_id);

      /*
       * Test all potential intersections between source box and fill
       * boxes including boxes shifted via periodic boundary conditions.
       * If there are periodic shifts, the source is always shifted
       * relative to the destination.
       */

      typename tbox::List< hier::IntVector<DIM> >::Iterator
         sh(src_level->getShiftsForPatch(src_patch_id));

      bool zero_shift = true;

      while (sh || zero_shift) {

         /*
          * Calculate the shift and the shifted source box.
          */

         hier::IntVector<DIM> shift(0);
         if (!zero_shift) {
            shift = sh();
         }

         const hier::Box<DIM> shifted(hier::Box<DIM>::shift(src_box, shift));

         for (nc = 0; nc < num_equiv_classes; nc++) {

            if (active_equivalence_class[nc]) {

               const typename xfer::CoarsenClasses<DIM>::Data& rep_item =
                  d_coarsen_classes->getClassRepresentative(nc);

               const int rep_item_dst_id = rep_item.d_dst;
               const int rep_item_src_id = rep_item.d_src;

               tbox::Pointer< hier::PatchDataFactory<DIM> > src_pdf =
                  src_patch_descriptor->getPatchDataFactory(rep_item_src_id);
               tbox::Pointer< hier::PatchDataFactory<DIM> > dst_pdf =
                  dst_patch_descriptor->getPatchDataFactory(rep_item_dst_id);

               const hier::IntVector<DIM>& dst_gcw = dst_pdf->getGhostCellWidth();

               hier::Box<DIM> dst_fill_box(hier::Box<DIM>::grow(dst_box, dst_gcw));

               hier::Box<DIM> test_mask(dst_fill_box*shifted);
               if ( test_mask.empty() &&
                    (dst_gcw == s_constant_zero_intvector) &&
                    dst_pdf->dataLivesOnPatchBorder() ) {
                  hier::Box<DIM> tmp_dst_fill_box(
                            hier::Box<DIM>::grow(dst_fill_box,
                                            s_constant_one_intvector));
                  test_mask = tmp_dst_fill_box * shifted;
               }
               hier::Box<DIM> src_mask( hier::Box<DIM>::shift(test_mask,-shift) );

               test_mask = hier::Box<DIM>::grow(src_box,
                                           hier::IntVector<DIM>::min(
                                              rep_item.d_gcw_to_coarsen,
                                              src_pdf->getGhostCellWidth()) );

               src_mask += test_mask;

               tbox::Pointer< hier::BoxOverlap<DIM> > overlap =
                  dst_pdf->getBoxGeometry(dst_box)
                              ->calculateOverlap(
                                *src_pdf->getBoxGeometry(src_box),
                                src_mask,
                                true, shift);

               if (overlap.isNull()) {
                  TBOX_ERROR("Internal LocallyActiveDataCoarsenSchedule<DIM> error..."
                             << "\n Overlap is NULL for "
                             << "\n src box = " << src_box
                             << "\n dst box = " << dst_box
                             << "\n src mask = " << src_mask << std::endl);
               }

               if (!overlap->isOverlapEmpty()) {
                  for (typename 
                       tbox::List<typename xfer::CoarsenClasses<DIM>::Data>::Iterator 
                       l(d_coarsen_classes->getIterator(nc)); l; l++) {

                     if (src_level_mgr->
                            getPatchDataActive( hier::PatchDataId(l().d_src), 
                                                hier::PatchNumber(src_patch_id) ) &&
                         dst_level_mgr->
                            getPatchDataActive( hier::PatchDataId(l().d_dst),
                                                hier::PatchNumber(dst_patch_id) ) ) {

                        d_schedule->addTransaction(
                           d_transaction_factory->allocate(dst_level,
                                                           src_level,
                                                           overlap,
                                                           dst_patch_id,
                                                           src_patch_id,
                                                           l().d_tag) );
                     }

                  } // iterate over coarsen components in equivalence class
               }

            }  // if equivalence class is active

         }  // iterate over all coarsen equivalence classes

         if (!zero_shift) {
            sh++;
         } else {
            zero_shift = false;
         }

      }  // iterate over valid shifts of source patch

   }  // if any equivalence class is active

}

/*
 * ************************************************************************
 *                                                                        *
 * Coarsen data from the source space on the fine patch level into the    *
 * source space on the coarse temporary patch level.                      *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM>
void LocallyActiveDataCoarsenSchedule<DIM>::coarsenSourceData(
   xfer::LocallyActiveDataCoarsenPatchStrategy<DIM>* patch_strategy) const
{
   NULL_USE(patch_strategy);

   const tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > 
      fine_level_mgr =
         hier::LocallyActiveVariableDatabase<DIM>::getDatabase()->
         getLocallyActiveDataPatchLevelManager(d_fine_level);

   /*
    * Loop over all local patches (fine and temp have the same mapping)
    */

   for (typename hier::PatchLevel<DIM>::Iterator p(d_fine_level); p; p++) {
      tbox::Pointer< hier::Patch<DIM> > fine_patch = d_fine_level->getPatch(p());
      tbox::Pointer< hier::Patch<DIM> > temp_patch = d_temp_crse_level->getPatch(p());

      tbox::List<const typename xfer::CoarsenClasses<DIM>::Data*> active_src_data;
      tbox::List<int> active_src_data_ids;
      for (int ici = 0; ici < d_number_coarsen_items; ici++) {
         const int src_data_id = d_coarsen_items[ici]->d_src;
         if ( fine_level_mgr->getPatchDataActive( hier::PatchDataId(src_data_id), 
                                                  hier::PatchNumber(p()) ) ) {
            /*
             * IMPORTANT! Proper generation of refine data lists requires
             *            appendItem(); addItem() is incorrect!
             */
            active_src_data.appendItem(d_coarsen_items[ici]);
            active_src_data_ids.appendItem(src_data_id);
         }
      }

      const hier::Box<DIM>& box = temp_patch->getBox();

      /*
       * Coarsen the fine space onto the temporary coarse space
       */

      if (d_coarsen_patch_strategy) {
         d_coarsen_patch_strategy->preprocessCoarsen(*temp_patch,
                                                     *fine_patch,
                                                     active_src_data_ids,
                                                     box,
                                                     d_ratio_between_levels);
      }

      for (typename tbox::List<const typename xfer::CoarsenClasses<DIM>::Data*>::Iterator 
           vi(active_src_data); vi; vi++) {
         int src_id = vi()->d_src;

         if ( !vi()->d_opcoarsen.isNull() ) {
            vi()->d_opcoarsen->coarsen(*temp_patch,
                                       *fine_patch,
                                       src_id,
                                       src_id,
                                       box,
                                       d_ratio_between_levels);
         }

      } // loop over active variables 

      if (d_coarsen_patch_strategy) {
         d_coarsen_patch_strategy->postprocessCoarsen(*temp_patch,
                                                     *fine_patch,
                                                     active_src_data_ids,
                                                     box,
                                                     d_ratio_between_levels);
      }

   } // loop over patches

}

/*
 * ***********************************************************************
 *                                                                       *
 * Private utility function to set up local array of coarsen items.      *
 *                                                                       *
 * ***********************************************************************
 */

template<int DIM>
void LocallyActiveDataCoarsenSchedule<DIM>::setCoarsenItems(
   const tbox::Pointer< xfer::CoarsenClasses<DIM> > coarsen_classes,
   bool called_from_reset)
{

   clearCoarsenItems();

   d_coarsen_classes        = coarsen_classes;

   const int num_coarsen_classes =
      d_coarsen_classes->getNumberEquivalenceClasses();

   d_number_coarsen_items = 0;

   /*
    * Determine total number of coarsen items and set state of
    * level manager used to manage storage on temporary level.
    */

   int nc;
   for (nc = 0; nc < num_coarsen_classes; nc++) {
      for (typename tbox::List<typename xfer::CoarsenClasses<DIM>::Data>::Iterator 
           l(d_coarsen_classes->getIterator(nc)); l; l++) {
         d_number_coarsen_items++;
      }
   }

   /*
    * Allocate and initialize array of coarsen items.
    */

   d_coarsen_items =
      new const typename xfer::CoarsenClasses<DIM>::Data*[d_number_coarsen_items];

   int ircount = 0;
   for (nc = 0; nc < num_coarsen_classes; nc++) {
      for (typename tbox::List<typename xfer::CoarsenClasses<DIM>::Data>::Iterator 
           l(d_coarsen_classes->getIterator(nc)); l; l++) {
         l().d_tag = ircount;
         d_coarsen_items[ircount] = &(l());
         ircount++;
      }
   }

   if (called_from_reset) {
      resetTempCoarseLevelDataManager();
   }

}

/*
 * ************************************************************************
 *                                                                        *
 * Private utility function to clear array of coarsen items.              *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM>
void LocallyActiveDataCoarsenSchedule<DIM>::clearCoarsenItems()
{
   if (d_coarsen_items) {
      for (int ici = 0; ici < d_number_coarsen_items; ici++) {
         d_coarsen_items[ici] = (typename xfer::CoarsenClasses<DIM>::Data*)NULL;
      }
      delete [] d_coarsen_items;
      d_coarsen_items = (const typename xfer::CoarsenClasses<DIM>::Data**)NULL;
      d_number_coarsen_items   = 0;
   }
}


/*
 * ***********************************************************************
 *                                                                       *
 * Private utility function to set up locally-active data manager for    *
 * temporary coarse level.                                               *
 *                                                                       *
 * ***********************************************************************
 */

template<int DIM>
void LocallyActiveDataCoarsenSchedule<DIM>::resetTempCoarseLevelDataManager()
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_temp_crse_level.isNull());
#endif

   d_temp_level_active_data_mgr.setNull();
   d_temp_level_active_data_mgr = 
      new hier::LocallyActiveDataPatchLevelManager<DIM>(d_temp_crse_level);

   const tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > 
      fine_level_mgr = hier::LocallyActiveVariableDatabase<DIM>::getDatabase()->
                       getLocallyActiveDataPatchLevelManager(d_fine_level);

   const int num_temp_patches = d_temp_crse_level->getNumberOfPatches();
   for (int tp = 0; tp < num_temp_patches; tp++) {
      for (int ici = 0; ici < d_number_coarsen_items; ici++) {
         int data_id = d_coarsen_items[ici]->d_src;
         if ( fine_level_mgr->getPatchDataActive( hier::PatchDataId(data_id), 
                                                  hier::PatchNumber(tp) ) ) {
            d_temp_level_active_data_mgr->
               setPatchDataActive( hier::PatchDataId(data_id), 
                                   hier::PatchNumber(tp) );
         }
      }
   }

}

/*
 * ***********************************************************************
 *                                                                       *
 * Private utility function to check coarsen items in initial setup to   *
 * see whether source and destination patch data components have         *
 * sufficient ghost cell widths to satisfy the "ghost width to coarsen"  *
 * functionality described in the CoarsenAlgorithm<DIM> class header.    *
 * Specifically, the destination data must have a ghost cell width at    *
 * least as large as the ghost cell width to coarsen.  The source data   *
 * must have a ghost cell width at least as large as the ghost cell      *
 * width to coarsen refined to the source (finer) level index space.     *
 * Other checks are also performed here by calling the                   *
 * xfer::CoarsenClasses<DIM>::checkCoarsenItem() routine.                *
 *                                                                       *
 * ***********************************************************************
 */

template<int DIM>
void LocallyActiveDataCoarsenSchedule<DIM>::initialCheckCoarsenClassItems() const
{
   tbox::Pointer< hier::PatchDescriptor<DIM> > pd = d_crse_level->getPatchDescriptor();

   hier::IntVector<DIM> user_gcw(0);
   if (d_coarsen_patch_strategy) {
      user_gcw = d_coarsen_patch_strategy->getCoarsenOpStencilWidth();
   }

   for (int ici = 0; ici < d_number_coarsen_items; ici++) {

      const typename xfer::CoarsenClasses<DIM>::Data* const crs_item = 
         d_coarsen_items[ici];

#ifdef DEBUG_CHECK_ASSERTIONS
      if (d_coarsen_classes->checkCoarsenItem(*crs_item, pd)) {
#endif

         const int dst_id = crs_item->d_dst;
         const int src_id = crs_item->d_src;

         tbox::Pointer< hier::PatchDataFactory<DIM> > dfact =
            pd->getPatchDataFactory(dst_id);
         tbox::Pointer< hier::PatchDataFactory<DIM> > sfact =
            pd->getPatchDataFactory(src_id);

         const hier::IntVector<DIM>& dst_gcw = dfact->getGhostCellWidth();
         const hier::IntVector<DIM>& src_gcw = sfact->getGhostCellWidth();

         if (crs_item->d_gcw_to_coarsen > dst_gcw) {
            TBOX_ERROR("Bad data given to LocallyActiveDataCoarsenSchedule<DIM>...\n"
                       << "`ghost cell width to coarsen' specified in\n"
                       << "registration of `Destination' patch data "
                       << pd->mapIndexToName(dst_id)
                       << " with CoarsenAlgorithm<DIM>\n"
                       << " is larger than ghost cell width of data \n"
                       << "d_gcw_to_coarsen = " << crs_item->d_gcw_to_coarsen
                       << "\n data ghost cell width = " << dst_gcw << std::endl);
         }

         if ( (crs_item->d_gcw_to_coarsen * d_ratio_between_levels) > src_gcw ) {
            TBOX_ERROR("Bad data given to LocallyActiveDataCoarsenSchedule<DIM>...\n"
                       << "`Source' patch data " << pd->mapIndexToName(src_id)
                       << " has ghost cell width too small to support the\n"
                       << "`ghost cell width to coarsen' specified in"
                       << " registration with CoarsenAlgorithm<DIM>\n"
                       << "data ghost cell width = " << src_gcw
                       << "d_gcw_to_coarsen = " << crs_item->d_gcw_to_coarsen
                       << "\nratio between levels = " << d_ratio_between_levels
                       << "\n Thus, data ghost width must be >= "
                       << (crs_item->d_gcw_to_coarsen * d_ratio_between_levels)
                       << std::endl);
         }
         if ( user_gcw > src_gcw) {
            TBOX_ERROR("Bad data given to LocallyActiveDataCoarsenSchedule<DIM>...\n"
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
 * Print coarsen schedule data to the specified output stream.            *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM>
void LocallyActiveDataCoarsenSchedule<DIM>::printClassData(std::ostream& stream) const
{
   stream << "xfer::LocallyActiveDataCoarsenSchedule<DIM>::printClassData()" << std::endl;
   stream << "---------------------------------------" << std::endl;
   stream << "s_schedule_generation_method = "
          << s_schedule_generation_method << std::endl;
   stream << "d_fill_coarse_data = " << d_fill_coarse_data << std::endl;

   d_coarsen_classes->printClassData(stream);

   d_schedule->printClassData(stream);

   if (d_fill_coarse_data) {
#if 0
// This is here as a reminder for consistency with regular CoarsenSchedule.
// The details of how refine operations should be used here have note been worked out yet.
      stream <<
      "Printing pre-coarsen refine algorithm that fills data before coarsening...\n";
      d_precoarsen_refine_algorithm->printClassData(stream);
      stream <<
      "Printing pre-coarsen refine schedule that fills data before coarsening...\n";
      d_precoarsen_refine_schedule->printClassData(stream);
#endif
   }
}

}
}

#endif
