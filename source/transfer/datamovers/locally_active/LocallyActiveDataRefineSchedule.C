//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-1/source/transfer/datamovers/locally_active/LocallyActiveDataRefineSchedule.C $
// Package:     SAMRAI data transfer
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1846 $
// Modified:	$LastChangedDate: 2008-01-11 09:51:05 -0800 (Fri, 11 Jan 2008) $
// Description:	Refine schedule for locally-active data transfer between AMR levels
//

#ifndef included_xfer_RefineSchedule_C
#define included_xfer_RefineSchedule_C

#include "LocallyActiveDataRefineSchedule.h"

#include "BoxArray.h"
#include "BoxGeometry.h"
#include "BoxOverlap.h"
#include "BoxTree.h"
#include "BoxUtilities.h"
#include "LocallyActiveVariableDatabase.h"
#include "Patch.h"
#include "PatchDataFactory.h"
#include "PatchDescriptor.h"
#include "ProcessorMapping.h"
#include "tbox/ArenaManager.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"
#include "LocallyActiveDataFillBox.h"
#include "RefineCopyTransaction.h"

#include "tbox/PIO.h"


namespace SAMRAI {
   namespace xfer {

#ifndef NULL
#define NULL (0)
#endif

#define BIG_GHOST_CELL_WIDTH  (10)

/*
*************************************************************************
*                                                                       *
* Initialization for static data members.                               *
*                                                                       *
*************************************************************************
*/
template<int DIM> const hier::IntVector<DIM>
   LocallyActiveDataRefineSchedule<DIM>::s_constant_zero_intvector = 
      hier::IntVector<DIM>(0);
template<int DIM> const hier::IntVector<DIM>
   LocallyActiveDataRefineSchedule<DIM>::s_constant_one_intvector = 
      hier::IntVector<DIM>(1);
template<int DIM> std::string
   LocallyActiveDataRefineSchedule<DIM>::s_schedule_generation_method = 
      "BOX_TREE";

/*
 * ************************************************************************
 *                                                                        *
 * Static function to set box intersection algorithm for schedules.       *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM>
void LocallyActiveDataRefineSchedule<DIM>::setScheduleGenerationMethod(
   const std::string& method)
{
   if ( !((method == "ORIG_NSQUARED") ||
          (method == "BOX_GRAPH") ||
          (method == "BOX_TREE")) ) {
      TBOX_ERROR("LocallyActiveDataRefineSchedule<DIM>::setScheduleGenerationMethod\n"
                 << "Given method string "
                 << method << " is invalid.\n Options are\n"
                 << "'ORIG_NSQUARED', 'BOX_GRAPH', and 'BOX_TREE'."
                 << std::endl);
   }

   s_schedule_generation_method = method;
}

/*
**************************************************************************
*
* Create a refine schedule that moves data from the source level into
* the destination level on the components represented by the refine
* list.  Ony data on the intersection of the two levels will be
* copied. It is assumed that the index spaces of the source and
* destination levels are "consistent"; i.e., they represent the same
* grid resolution.  The levels do not have to be part of the same
* AMR patch hierarchy, however.
*
**************************************************************************
*/

template<int DIM>
LocallyActiveDataRefineSchedule<DIM>::LocallyActiveDataRefineSchedule(
   tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
   tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > dst_level_mgr,
   tbox::Pointer< hier::PatchLevel<DIM> > src_level,
   tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > src_level_mgr,
   const tbox::Pointer< xfer::RefineClasses<DIM> > refine_classes,
   tbox::Pointer< xfer::LocallyActiveDataRefineTransactionFactory<DIM> > transaction_factory,
   xfer::LocallyActiveDataRefinePatchStrategy<DIM>* patch_strategy,
   bool use_time_interpolation)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst_level.isNull());
   TBOX_ASSERT(!dst_level_mgr.isNull());
   TBOX_ASSERT(dst_level_mgr->checkLevel(dst_level));
   TBOX_ASSERT(!src_level.isNull());
   TBOX_ASSERT(!src_level_mgr.isNull());
   TBOX_ASSERT(src_level_mgr->checkLevel(src_level));
   TBOX_ASSERT(!refine_classes.isNull());
   TBOX_ASSERT(!transaction_factory.isNull());
#endif

   t_fill_data = tbox::TimerManager::getManager() ->
      getTimer("xfer::LocallyActiveDataRefineSchedule::fillData()");
   t_gen_comm_sched = tbox::TimerManager::getManager()->
      getTimer("xfer::LocallyActiveDataRefineSchedule::generateCommunicationSchedule()");
   t_finish_sched_const = tbox::TimerManager::getManager()->
      getTimer("xfer::LocallyActiveDataRefineSchedule::finishScheduleConstruction()");

   t_gen_comm_sched_unfilled = tbox::TimerManager::getManager()->
      getTimer("xfer::LocallyActiveDataRefineSchedule::gen_comm_sched_unfilled");
   t_gen_comm_sched_trans = tbox::TimerManager::getManager()->
      getTimer("xfer::LocallyActiveDataRefineSchedule::gen_comm_sched_trans");

   /*
    * Initial values; some will change in setup operations.
    */

   d_transaction_factory = transaction_factory;

   d_dst_level = dst_level;
   d_dst_level_mgr = dst_level_mgr;

   d_refine_patch_strategy = patch_strategy;

   d_number_refine_items = 0;
   d_refine_items = (const typename xfer::RefineClasses<DIM>::Data**)NULL;

   setRefineItems(refine_classes);
   initialCheckRefineClassItems();

   d_force_boundary_fill       = false;

   d_domain_is_one_box         = false;
   d_num_periodic_directions   = 0;

   d_coarse_priority_level_schedule = new tbox::Schedule();
   d_fine_priority_level_schedule = new tbox::Schedule();

   d_coarse_schedule.setNull();
   d_coarse_level.setNull();
   d_coarse_level_mgr.setNull();

   d_max_fill_boxes = 0;

   /*
    * Initialize destination level, ghost cell widths,
    * and domain information data members.
    */

   bool recursive_schedule = false;
   initializeDomainAndGhostInformation(recursive_schedule);

   /*
    * Create the fill box and unfilled box arrays and then the
    * communication schedule for data transfers between source and
    * destination levels.  Note that the fill boxes are initialized here,
    * while the unfilled boxes are set in the generateCommunicationSchedule()
    * and contain the regions for each patch in the destination level
    * that will cannot be filled by the schedule.
    */

   tbox::Array< xfer::LocallyActiveDataFillBoxSet<DIM> > la_fill_boxes;
   allocateDefaultFillBoxes(la_fill_boxes,
                            d_dst_level,
                            d_dst_level_mgr,
                            d_boundary_fill_ghost_width);

   tbox::Array< xfer::LocallyActiveDataFillBoxSet<DIM> >
      la_unfilled_boxes(dst_level->getNumberOfPatches());

   tbox::Pointer<tbox::Timer> t_gen_comm_sched =
      tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::generate_comm_schedule");

   t_gen_comm_sched->start();
   generateCommunicationSchedule(d_coarse_priority_level_schedule,
                                 d_fine_priority_level_schedule,
                                 d_dst_level,
                                 d_dst_level_mgr, 
                                 src_level,
                                 src_level_mgr,
                                 la_fill_boxes,
                                 la_unfilled_boxes,
                                 use_time_interpolation);
   t_gen_comm_sched->stop();
}

/*
*************************************************************************
*
* Create a refine schedule that moves data from the source level into
* the destination level for the components represented by the refine
* list and which are defined on both destination and source levels.
* If portions of the destination level remain unfilled, then the algorithm
* recursively fills those unfilled portions from coarser levels in the AMR
* hierarchy, again where data is defined on both destination and source levels.
* It is assumed that the index spaces of the source and destination levels are
* "consistent"; i.e., they represent the same grid resolution.  Also, the next
* coarser level integer argument must be the number of level in the specified
* hierarchy representing the next coarser level of mesh resolution to
* the destination level.
*
* IMPORTANT NOTES: The source level may be NULL, in which case the
* destination level will be filled only using data interpolated from
* coarser levels in the AMR hierarchy.  The hierarchy may be NULL only
* if the next coarser level is < 0 (that is, there is no coarser level).
*
*************************************************************************
*/

template<int DIM>
LocallyActiveDataRefineSchedule<DIM>::LocallyActiveDataRefineSchedule(
   tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
   tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > dst_level_mgr,
   tbox::Pointer< hier::PatchLevel<DIM> > src_level,
   tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > src_level_mgr,
   int next_coarser_level,
   tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
   const tbox::Pointer< xfer::RefineClasses<DIM> > refine_classes,
   tbox::Pointer< xfer::LocallyActiveDataRefineTransactionFactory<DIM> > transaction_factory,
   xfer::LocallyActiveDataRefinePatchStrategy<DIM>* patch_strategy,
   bool use_time_interpolation)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst_level.isNull());
   TBOX_ASSERT(!dst_level_mgr.isNull());
   TBOX_ASSERT(dst_level_mgr->checkLevel(dst_level));
   if (!src_level.isNull()) {
      TBOX_ASSERT(!src_level_mgr.isNull());
      TBOX_ASSERT(src_level_mgr->checkLevel(src_level));
   }
   TBOX_ASSERT((next_coarser_level < 0) || !hierarchy.isNull());
   TBOX_ASSERT(!refine_classes.isNull());
   TBOX_ASSERT(!transaction_factory.isNull());
#endif

   t_fill_data = tbox::TimerManager::getManager() ->
      getTimer("xfer::LocallyActiveDataRefineSchedule::fillData()");
   t_gen_comm_sched = tbox::TimerManager::getManager()->
      getTimer("xfer::LocallyActiveDataRefineSchedule::generateCommunicationSchedule()");
   t_finish_sched_const = tbox::TimerManager::getManager()->
      getTimer("xfer::LocallyActiveDataRefineSchedule::finishScheduleConstruction()");

   t_gen_comm_sched_unfilled = tbox::TimerManager::getManager()->
      getTimer("xfer::LocallyActiveDataRefineSchedule::gen_comm_sched_unfilled");
   t_gen_comm_sched_trans = tbox::TimerManager::getManager()->
      getTimer("xfer::LocallyActiveDataRefineSchedule::gen_comm_sched_trans");

   /*
    * Initial values; some will change in setup operations.
    */

   d_transaction_factory = transaction_factory;

   d_dst_level = dst_level;
   d_dst_level_mgr = dst_level_mgr;

   d_refine_patch_strategy = patch_strategy;

   d_number_refine_items = 0;
   d_refine_items = (const typename xfer::RefineClasses<DIM>::Data**)NULL;

   setRefineItems(refine_classes);
   initialCheckRefineClassItems();

   d_force_boundary_fill       = false;

   d_domain_is_one_box         = false;
   d_num_periodic_directions   = 0;

   d_coarse_priority_level_schedule.setNull();
   d_fine_priority_level_schedule.setNull();

   d_coarse_schedule.setNull();
   d_coarse_level.setNull();
   d_coarse_level_mgr.setNull();

   d_max_fill_boxes = 0;

   /*
    * Initialize destination level, ghost cell widths,
    * and domain information data members.
    */

   bool recursive_schedule = false;
   initializeDomainAndGhostInformation(recursive_schedule);

   /*
    * Create the fill box arrays and then the communication schedule(s)
    * needed to move data from the patch hierarchy to the destination level.
    */

   tbox::Array< xfer::LocallyActiveDataFillBoxSet<DIM> > la_fill_boxes;
   allocateDefaultFillBoxes(la_fill_boxes,
                            d_dst_level,
                            d_dst_level_mgr,
                            d_boundary_fill_ghost_width);

   t_finish_sched_const->start();
   finishScheduleConstruction(src_level,
                              src_level_mgr,
                              next_coarser_level,
                              hierarchy,
                              la_fill_boxes,
                              use_time_interpolation);
   t_finish_sched_const->stop();

}

/*
*************************************************************************
*									*
* This private constructor creates a refine schedule that moves data    *
* into the destination only on the specified fill boxes.                *
*									*
*************************************************************************
*/

template<int DIM>
LocallyActiveDataRefineSchedule<DIM>::LocallyActiveDataRefineSchedule(
   tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
   tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > dst_level_mgr,
   tbox::Pointer< hier::PatchLevel<DIM> > src_level,
   tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > src_level_mgr,
   int next_coarser_level,
   tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
   const tbox::Pointer< xfer::RefineClasses<DIM> > refine_classes,
   tbox::Pointer< xfer::LocallyActiveDataRefineTransactionFactory<DIM> > transaction_factory, 
   const tbox::Array< xfer::LocallyActiveDataFillBoxSet<DIM> >& la_fill_boxes,
   xfer::LocallyActiveDataRefinePatchStrategy<DIM>* patch_strategy)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst_level.isNull());
   TBOX_ASSERT(!dst_level_mgr.isNull());
   TBOX_ASSERT(dst_level_mgr->checkLevel(dst_level));
   if (!src_level.isNull()) {
      TBOX_ASSERT(!src_level_mgr.isNull());
      TBOX_ASSERT(src_level_mgr->checkLevel(src_level));
   }
   TBOX_ASSERT((next_coarser_level < 0) || !hierarchy.isNull());
   TBOX_ASSERT(!refine_classes.isNull());
   TBOX_ASSERT(!transaction_factory.isNull());
   TBOX_ASSERT(la_fill_boxes.getSize() == dst_level->getNumberOfPatches());
#endif

   t_fill_data = tbox::TimerManager::getManager() ->
      getTimer("xfer::LocallyActiveDataRefineSchedule::fillData()");
   t_gen_comm_sched = tbox::TimerManager::getManager()->
      getTimer("xfer::LocallyActiveDataRefineSchedule::generateCommunicationSchedule()");
   t_finish_sched_const = tbox::TimerManager::getManager()->
      getTimer("xfer::LocallyActiveDataRefineSchedule::finishScheduleConstruction()");

   t_gen_comm_sched_unfilled = tbox::TimerManager::getManager()->
      getTimer("xfer::LocallyActiveDataRefineSchedule::gen_comm_sched_unfilled");
   t_gen_comm_sched_trans = tbox::TimerManager::getManager()->
      getTimer("xfer::LocallyActiveDataRefineSchedule::gen_comm_sched_trans");

   /*
    * Initial values; some will change in setup operations.
    * Note that we do not check refine items here, since this
    * constructor is private and called recursively (i.e., the
    * items have been checked already).
    */

   d_transaction_factory = transaction_factory;

   d_dst_level = dst_level;

   d_dst_level_mgr = dst_level_mgr;

   d_refine_patch_strategy = patch_strategy;

   d_number_refine_items = 0;
   d_refine_items = (const typename xfer::RefineClasses<DIM>::Data**)NULL;

   setRefineItems(refine_classes);

   d_force_boundary_fill       = false;

   d_domain_is_one_box         = false;
   d_num_periodic_directions   = 0;

   d_coarse_priority_level_schedule.setNull();
   d_fine_priority_level_schedule.setNull();

   d_coarse_schedule.setNull();
   d_coarse_level.setNull();
   d_coarse_level_mgr.setNull();

   d_max_fill_boxes = 0;   

   /*
    * Initialize destination level, ghost cell widths,
    * and domain information data members.
    */

   bool recursive_schedule = true;
   initializeDomainAndGhostInformation(recursive_schedule);

   /*
    * Finish construction of the communication schedule using the remaining
    * fill boxes.
    */

   const int nfill_boxes = la_fill_boxes.getSize();
   for (int p = 0; p < nfill_boxes; p++) {
      d_max_fill_boxes =
         tbox::MathUtilities<int>::Max(d_max_fill_boxes,
                                       la_fill_boxes[p].getNumberOfBoxes());
   }
 
   bool use_time_interpolation = true;

   finishScheduleConstruction(src_level,
                              src_level_mgr,
                              next_coarser_level,
                              hierarchy,
                              la_fill_boxes,
                              use_time_interpolation);

}

/*
*************************************************************************
*									*
* The destructor for the refine schedule class implicitly deallocates	*
* all of the data associated with the communication schedule.		*
*									*
*************************************************************************
*/

template<int DIM>
LocallyActiveDataRefineSchedule<DIM>::~LocallyActiveDataRefineSchedule()
{
   clearRefineItems();

   d_dst_level.setNull();
   d_dst_level_mgr.setNull();

   d_refine_patch_strategy = (xfer::LocallyActiveDataRefinePatchStrategy<DIM>*)NULL;

   d_transaction_factory.setNull();

   d_coarse_priority_level_schedule.setNull();
   d_fine_priority_level_schedule.setNull();

   d_coarse_schedule.setNull();
   d_coarse_level.setNull();
   d_coarse_level_mgr.setNull();

   t_fill_data.setNull();
   t_gen_comm_sched.setNull();
   t_finish_sched_const.setNull();
   t_gen_comm_sched_unfilled.setNull();
   t_gen_comm_sched_trans.setNull();
}

/*
*************************************************************************
*									*
* Finish construction of the communication schedule when filling from	*
* coarser levels.  This routine guarantees that the fill boxes will be	*
* filled from the same level or coarser levels.  If the fill boxes can	*
* not be filled, then the routine aborts with an internal error.	*
*									*
*************************************************************************
*/

template<int DIM>
void LocallyActiveDataRefineSchedule<DIM>::finishScheduleConstruction(
   tbox::Pointer< hier::PatchLevel<DIM> > src_level,
   tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > src_level_mgr,
   int next_coarser_level,
   tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
   const tbox::Array< xfer::LocallyActiveDataFillBoxSet<DIM> >& la_fill_boxes,
   bool use_time_interpolation)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (!src_level.isNull()) {
      TBOX_ASSERT(!src_level_mgr.isNull());
      TBOX_ASSERT(src_level_mgr->checkLevel(src_level));
   }
   TBOX_ASSERT((next_coarser_level < 0) || !hierarchy.isNull());
   TBOX_ASSERT(la_fill_boxes.getSize() == d_dst_level->getNumberOfPatches());
#endif

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
   const int dst_npatches = d_dst_level->getNumberOfPatches();
   tbox::Array< xfer::LocallyActiveDataFillBoxSet<DIM> > la_unfilled_boxes(dst_npatches);

   if (!src_level.isNull()) {

      t_gen_comm_sched->start();
      generateCommunicationSchedule(d_coarse_priority_level_schedule,
                                    d_fine_priority_level_schedule,
                                    d_dst_level,
                                    d_dst_level_mgr,
                                    src_level,
                                    src_level_mgr,
                                    la_fill_boxes,
                                    la_unfilled_boxes,
                                    use_time_interpolation);
      t_gen_comm_sched->stop();

   } else {
      for (int p = 0; p < dst_npatches; p++) {
         la_unfilled_boxes[p].setTo(la_fill_boxes[p]);
      }
   }

   /*
    * Remove pieces of the boxes to be filled from coarser levels that 
    * live outside the physical domain in the non-periodic directions.  
    * Then, check whether there still exist boxes to be filled from 
    * coarser levels.  Any resulting boxes will be used to generate a 
    * temporary coarse patch level for acquiring data to refine.
    */

   hier::IntVector<DIM> big_grow_vector(s_constant_zero_intvector);
   if (d_num_periodic_directions > 0) {
      for (int dir=0; dir < DIM; dir++) {
         if (d_periodic_shift(dir)) {
            big_grow_vector(dir) = BIG_GHOST_CELL_WIDTH;
         }
      }
   }

   const bool do_shearing = ( (d_num_periodic_directions >= 0) &&
                              (d_num_periodic_directions < DIM) );

   bool need_to_fill = false;
   if (do_shearing) {

      hier::BoxList<DIM> shear_domain(d_dst_level->getPhysicalDomain());
      if (d_domain_is_one_box) {
        shear_domain.clearItems();
        shear_domain.addItem(d_domain_box);
      }

      if (d_num_periodic_directions > 0) {
         shear_domain.grow(big_grow_vector);
      }

      shear_domain.simplifyBoxes();

      for (int p = 0; p < dst_npatches; p++) {
         if (la_unfilled_boxes[p].getNumberOfBoxes() > 0) {
            if (d_dst_level->patchTouchesRegularBoundary(p)) {
               la_unfilled_boxes[p].intersectBoxes(shear_domain);
               if (la_unfilled_boxes[p].getNumberOfBoxes() > 0 ) {
                  need_to_fill = true;
               }
            } else {
               need_to_fill = true;
            }
         }
      }

   } else {  // !do_shearing

      for (int p = 0; p < dst_npatches; p++) {
         if (la_unfilled_boxes[p].getNumberOfBoxes() > 0) {
            need_to_fill = true;
         }
      }

   } 

   /*
    * If there remain boxes to be filled from coarser levels, then 
    * generate data to do this and recurse to the next coarser level. 
    */

   if (need_to_fill) {

      /*
       * If there are no coarser levels in the hierarchy or the hierarchy pointer
       * is null, then throw an error.  Something is messed up someplace and
       * code execution cannot proceed.
       */

      if (next_coarser_level < 0) {
         TBOX_ERROR("Internal xfer::LocallyActiveDataRefineSchedule<DIM> error..."
                    << "\n In finishScheduleConstruction() -- "
                    << "\n No coarser levels...could not fill from coarser."
                    << std::endl);
      } else {
         if (hierarchy.isNull()) {
            TBOX_ERROR("Internal xfer::LocallyActiveDataRefineSchedule<DIM> error..."
                       << "\n In finishScheduleConstruction() -- "
                       << "\n Need to fill from coarser hierarchy level and \n"
                       << "hierarchy is unavailable." << std::endl);
         }
      }

      tbox::Pointer< hier::PatchLevel<DIM> > coarser_level_from_hierarchy =
         hierarchy->getPatchLevel(next_coarser_level);
      const hier::IntVector<DIM> ratio =
         d_dst_level->getRatio() / coarser_level_from_hierarchy->getRatio();

      const hier::BoxList<DIM> coarser_physical_domain(
         coarser_level_from_hierarchy->getPhysicalDomain());

      const bool do_coarse_shearing = (do_shearing && !d_domain_is_one_box);

      hier::BoxList<DIM> coarser_shear_domain(coarser_physical_domain);

      if (do_coarse_shearing) {

         if (d_num_periodic_directions > 0) {
            coarser_shear_domain.grow(big_grow_vector);
         }

         coarser_shear_domain.simplifyBoxes();

      }

      /*
       * Create the box description list and fine-to-coarse mapping for
       * the next coarser level in the hierarchy.  Each fine patch will
       * have zero or more corresponding coarse patches needed to obtain
       * the data to fill the unfilled regions of the patch.  The fine patch
       * will have zero coarse patches if its unfilled box list is empty.
       * It may have more than one if the generated coarse level intersects
       * the physical domain in funny ways.
       */

      hier::BoxList<DIM> coarse_box_list;
      tbox::List<int> coarse_to_fine_mapping;

      for (int fp = 0; fp < dst_npatches; fp++) {
         if (la_unfilled_boxes[fp].getNumberOfBoxes() > 0) {

            /*
             * Create boxes that will be used to construct temporary
             * coarse patch level.  We first bound the unfilled boxes
             * with a single box.  Then, we coarsen this bounding box.
             * After intersecting the bounding box with the physical
             * domain (either coarse level or sheared for periodic shifts),
             * we extend any of the resulting boxes to the physical
             * domain to avoid bizarre intersections of the ghost cell
             * regions with the domain boundary.  Finally, we add the
             * coarse box(es) and coarse-to-fine mapping information
             * to the proper lists for making the temporary coarse level.
             */

            hier::Box<DIM> coarser_fill_bounding_box(
               la_unfilled_boxes[fp].getBoundingBox());

            coarser_fill_bounding_box.coarsen(ratio);

            if ( do_coarse_shearing &&
                 (d_dst_level->patchTouchesRegularBoundary(fp)) ) {

               hier::BoxList<DIM> coarse_boxes(coarser_fill_bounding_box); 
               coarse_boxes.intersectBoxes(coarser_shear_domain);
               coarse_boxes.simplifyBoxes();

               (void) hier::BoxUtilities<DIM>::extendBoxesToDomainBoundary(
                                          coarse_boxes, 
                                          coarser_physical_domain, 
                                          d_max_scratch_gcw);

               for (typename hier::BoxList<DIM>::Iterator b(coarse_boxes); b; b++) {
                  coarse_box_list.appendItem(b()); 
                  coarse_to_fine_mapping.appendItem(fp);
               }

            } else {

               (void) hier::BoxUtilities<DIM>::extendBoxToDomainBoundary(
                                          coarser_fill_bounding_box,
                                          coarser_physical_domain,
                                          d_max_scratch_gcw);

               coarse_box_list.appendItem(coarser_fill_bounding_box);
               coarse_to_fine_mapping.appendItem(fp);

            }

         }
      }

      /*
       * Create a temporary patch level from the coarse boxes and
       * in doing so, use the coarse-to-fine mapping array.
       */

      const int ncoarse = coarse_to_fine_mapping.getNumberItems();

      d_coarse_to_fine_mapping.resizeArray(ncoarse);
      int coarse_index = 0;
      for (tbox::List<int>::Iterator ci(coarse_to_fine_mapping); ci; ci++) {
         d_coarse_to_fine_mapping[coarse_index++] = ci();
      }

      hier::ProcessorMapping coarse_processor_mapping(ncoarse);
      for (int c = 0; c < ncoarse; c++) {
         coarse_processor_mapping.setProcessorAssignment(
            c, d_dst_level->getMappingForPatch(d_coarse_to_fine_mapping[c]));
      }

      hier::BoxArray<DIM> coarse_box_array(coarse_box_list);
      
      d_coarse_level = new hier::PatchLevel<DIM>(
         coarse_box_array,
         coarse_processor_mapping,
         coarser_level_from_hierarchy->getRatio(),
         coarser_level_from_hierarchy->getGridGeometry(),
         coarser_level_from_hierarchy->getPatchDescriptor());
      d_coarse_level->setLevelNumber(next_coarser_level);
      d_coarse_level->setNextCoarserHierarchyLevelNumber(next_coarser_level);

      /*
       * Cache the fine fill boxes in this schedule object for moving
       * interpolated data to this level and generate boxes to fill on
       * the temporary coarser level.  The fine fill boxes represent
       * the valid range of the coarse to fine data interpolation.
       * The coarse fill boxes represent the coarse regions where data
       * must be obtained to perform the appropriate interpolation.
       */

      d_la_fine_fill_boxes.resizeArray(ncoarse);
      tbox::Array< xfer::LocallyActiveDataFillBoxSet<DIM> > la_coarse_fill_boxes(ncoarse);
      la_coarse_fill_boxes.resizeArray(ncoarse);

      d_coarse_level_mgr = new hier::LocallyActiveDataPatchLevelManager<DIM>(d_coarse_level);

      for (int cp = 0; cp < ncoarse; cp++) {
         const int fp = d_coarse_to_fine_mapping[cp];

         const hier::Box<DIM> coarse_box = coarse_box_array.getBox(cp);
         const hier::Box<DIM> fine_box = hier::Box<DIM>::refine(coarse_box, ratio);

         d_la_fine_fill_boxes[cp].setTo(la_unfilled_boxes[fp]);

         d_la_fine_fill_boxes[cp].intersectBoxes(fine_box);

#if 0
         hier::BoxList<DIM> coarse_to_fill;
         if (d_coarse_level->patchTouchesRegularBoundary(cp)) {
            coarse_to_fill.addItem(coarse_box);
         } else {
            coarse_to_fill = d_la_fine_fill_boxes[cp].getBoxList();
            coarse_to_fill.coarsen(ratio);
         }
#else
         hier::BoxList<DIM> coarse_to_fill = d_la_fine_fill_boxes[cp].getBoxList();
         coarse_to_fill.coarsen(ratio);
#endif

         coarse_to_fill.grow(d_max_stencil_gcw);

         const tbox::List<const typename xfer::RefineClasses<DIM>::Data*>& var_data = 
            d_la_fine_fill_boxes[cp].getUnionActiveRefineVarData();
         for (typename tbox::List<const typename xfer::RefineClasses<DIM>::Data*>::Iterator 
              vi(var_data); vi; vi++) {
            d_coarse_level_mgr->setPatchDataActive( 
                                hier::PatchDataId(vi()->d_scratch), 
                                hier::PatchNumber(cp) );
         }
 
         for (typename hier::BoxList<DIM>::Iterator cfb(coarse_to_fill); cfb; cfb++) {
            la_coarse_fill_boxes[cp].addLocallyActiveFillBox(cfb(), var_data);
         }

      }

      /*
       * Recursively fill the coarse schedule using the private
       * refine schedule constructor.
       */
      tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> >
         coarser_level_from_hierarchy_mgr =
            hier::LocallyActiveVariableDatabase<DIM>::getDatabase()->
            getLocallyActiveDataPatchLevelManager(coarser_level_from_hierarchy);
      d_coarse_schedule =
         new LocallyActiveDataRefineSchedule<DIM>(d_coarse_level,
                                                  d_coarse_level_mgr,
                                                  coarser_level_from_hierarchy,
                                                  coarser_level_from_hierarchy_mgr,
                                                  next_coarser_level-1,
                                                  hierarchy,
                                                  d_refine_classes,
                                                  d_transaction_factory,
                                                  la_coarse_fill_boxes,
                                                  d_refine_patch_strategy);

   }

}

/*
*************************************************************************
*									*
* Execute the stored communication schedule that moves data into the	*
* destination component of the destination level.  The algorithm is as	*
* follows:								*
*									*
*	(1) Allocate scratch space on the destination level if it does	*
*	    not exist.							*
*	(2) Call the recursive fill function that will fill the scratch	*
*	    space of the destination level.				*
*	(3) If the scratch and destination spaces are not the same,	*
*	    then copy locally from the scratch into the destination.	*
*	(4) Deallocate any previously allocated data.			*
*									*
*************************************************************************
*/

template<int DIM>
void LocallyActiveDataRefineSchedule<DIM>::fillData(
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
 
   hier::LocallyActiveDataPatchLevelManager<DIM> allocate_mgr(d_dst_level);
   allocateScratchSpace(d_dst_level, d_dst_level_mgr, fill_time, allocate_mgr); 
   
   /*
    * Begin the recursive algorithm that fills from coarser, fills from
    * same, and then fills physical boundaries.
    */

   recursiveFill(fill_time, do_physical_boundary_fill);

   /*
    * Copy the scratch space of the destination level to the destination
    * space.
    */

   copyScratchToDestination(d_dst_level, d_dst_level_mgr);

   /*
    * Deallocate any allocated scratch space on the destination level.
    */

   allocate_mgr.deallocateAllPatchData();

   /*
    * Unset the refine items for all transactions.  These items are
    * shared by all transaction objects in the communication schedule.
    */

   d_transaction_factory->unsetRefineItems();

   t_fill_data->stop();
}

/*
**************************************************************************
*                                                                        *
* Recursively fill the required fill boxes on the destination level.     *
* The scratch component of the level will be filled.  The algorithm      *
* is as follows:                                                         *
*                                                                        *
*   (1) If we need to get data from a coarser level, then:               *
*      (a) allocate scratch data on the coarser level                    *
*      (b) recursively call this routine to get the data                 *
*      (c) refine data from the coarse level into this level             *
*      (d) deallocate the scratch data on the coarser level              *
*   (2) Copy data from the same level of refinement                      *
*   (3) Copy data from the physical boundary conditions                  *
*                                                                        *
**************************************************************************
*/

template<int DIM>
void LocallyActiveDataRefineSchedule<DIM>::recursiveFill(
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

   if (!d_coarse_schedule.isNull()) {

      /*
       * Allocate data on the coarser level and keep track of the allocated
       * components so that they may be deallocated later.
       */

      hier::LocallyActiveDataPatchLevelManager<DIM> allocate_mgr(d_coarse_level);
      allocateScratchSpace(d_coarse_level, d_coarse_level_mgr, fill_time, allocate_mgr);
      
      /*
       * Recursively call the fill routine to fill the required coarse fill
       * boxes on the coarser level.
       */

      d_coarse_schedule->recursiveFill(fill_time, do_physical_boundary_fill);

      /*
       * The coarse fill boxes should now be filled.  Now interpolate
       * data from the coarse grid into the fine grid.
       */

      refineScratchData();

      /*
       * Deallocate the scratch data from the coarse grid.
       */

      allocate_mgr.deallocateAllPatchData(); 
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
      fillPhysicalBoundaries(d_dst_level, fill_time);
   }
}

/*
*************************************************************************
*									*
* Fill the physical boundaries of the specified level with data.	*
*									*
*************************************************************************
*/

template<int DIM>
void LocallyActiveDataRefineSchedule<DIM>::fillPhysicalBoundaries(
   tbox::Pointer< hier::PatchLevel<DIM> > level,
   double fill_time) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!level.isNull());
   TBOX_ASSERT(!d_dst_level_mgr.isNull());
   TBOX_ASSERT(d_dst_level_mgr->checkLevel(level));
#endif

   if (d_refine_patch_strategy) {

      for (typename hier::PatchLevel<DIM>::Iterator p(level); p; p++) {
         tbox::Pointer< hier::Patch<DIM> > patch = level->getPatch(p());

         if (patch->getPatchGeometry()->intersectsPhysicalBoundary()) {

            tbox::List<const typename xfer::RefineClasses<DIM>::Data*> active_data;
            for (int iri = 0; iri < d_number_refine_items; iri++) {
               if ( d_dst_level_mgr->getPatchDataActive(
                                     hier::PatchDataId(d_refine_items[iri]->d_scratch), 
                                     hier::PatchNumber(p()) ) ) {
                   /*
                    * IMPORTANT! Proper generation of refine data lists requires
                    *            appendItem(); addItem() is incorrect since
                    *            it results in incorrect order. 
                    */
                   active_data.appendItem(d_refine_items[iri]);
               }
            }

            d_refine_patch_strategy->
               setPhysicalBoundaryConditions(*patch,
                                             active_data,
                                             fill_time,
                                             d_boundary_fill_ghost_width);

         }
      }

   }

}

/*
**************************************************************************
*                                                                        *
* Check whether the scratch data needs to be allocated on the specified  *
* level.  Keep track of those allocated components so that they may be   *
* deallocated later.                                                     *
*                                                                        *
**************************************************************************
*/

template<int DIM>
void LocallyActiveDataRefineSchedule<DIM>::allocateScratchSpace(
   tbox::Pointer< hier::PatchLevel<DIM> > level,
   tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > level_mgr,
   double fill_time,
   hier::LocallyActiveDataPatchLevelManager<DIM>& allocate_mgr) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!level.isNull());
   TBOX_ASSERT(!level_mgr.isNull());
   TBOX_ASSERT(level_mgr->checkLevel(level));
#endif

   hier::LocallyActiveDataPatchLevelManager<DIM> preprocess_mgr(level);

   for (typename hier::PatchLevel<DIM>::Iterator p(level); p; p++) {
      tbox::Pointer< hier::Patch<DIM> > patch = level->getPatch(p());

      for (int iri = 0; iri < d_number_refine_items; iri++) {
         const int data_id = d_refine_items[iri]->d_scratch;
         if ( level_mgr->getPatchDataActive( hier::PatchDataId(data_id), 
                                             hier::PatchNumber(p()) ) ) {
            if (!patch->checkAllocated(data_id)) {
               allocate_mgr.setPatchDataActive( hier::PatchDataId(data_id), 
                                                hier::PatchNumber(p()) );
            }
            preprocess_mgr.setPatchDataActive( hier::PatchDataId(data_id), 
                                               hier::PatchNumber(p()) );
         }
      }
   }

   allocate_mgr.allocateAllPatchData(fill_time,
                tbox::ArenaManager::getManager()->getScratchAllocator());

   d_transaction_factory->preprocessScratchSpace(level,
                                                 fill_time,
                                                 preprocess_mgr);

}

/*
**************************************************************************
*                                                                        *
* Copy data from the scratch space of the specified level into the       *
* destination space.  If the two spaces are the same, then no copy       *
* is performed.                                                          *
*                                                                        *
**************************************************************************
*/

template<int DIM>
void LocallyActiveDataRefineSchedule<DIM>::copyScratchToDestination(
   tbox::Pointer< hier::PatchLevel<DIM> > level,
   tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > level_mgr) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!level.isNull());
   TBOX_ASSERT(!level_mgr.isNull());
   TBOX_ASSERT(level_mgr->checkLevel(level));
#endif

   for (typename hier::PatchLevel<DIM>::Iterator p(level); p; p++) {
      tbox::Pointer< hier::Patch<DIM> > patch = level->getPatch(p());

      for (int iri = 0; iri < d_number_refine_items; iri++) {
         const int src_id = d_refine_items[iri]->d_scratch;
         const int dst_id = d_refine_items[iri]->d_dst;
         if ( level_mgr->getPatchDataActive( hier::PatchDataId(dst_id), 
                                             hier::PatchNumber(p()) ) &&
              (src_id != dst_id) ) {
#ifdef DEBUG_CHECK_ASSERTIONS
            const double dst_time = patch->getPatchData(dst_id)->getTime();
            const double src_time = patch->getPatchData(src_id)->getTime();
            TBOX_ASSERT(dst_time == src_time);
#endif
            patch->getPatchData(dst_id)->copy(*patch->getPatchData(src_id));
         }
      }

   }

}

/*
*************************************************************************
*									*
* Refine data from the coarse level into the fine level on the provided	*
* fill box regions.  All operations are performed on the scratch space.	*
*									*
*************************************************************************
*/

template<int DIM>
void LocallyActiveDataRefineSchedule<DIM>::refineScratchData() const
{

   const hier::IntVector<DIM> ratio =
      d_dst_level->getRatio() / d_coarse_level->getRatio();

   /*
    * Loop over all the coarse patches and find the corresponding destination
    * patch and destination fill boxes.
    */

   for (typename hier::PatchLevel<DIM>::Iterator p(d_coarse_level); p; p++) {
      const int mapping = d_coarse_to_fine_mapping[p()];
      tbox::Pointer< hier::Patch<DIM> > fine_patch = d_dst_level->getPatch(mapping);
      tbox::Pointer< hier::Patch<DIM> > crse_patch = d_coarse_level->getPatch(p());
      const xfer::LocallyActiveDataFillBoxSet<DIM>& fill_boxes = 
         d_la_fine_fill_boxes[p()];

      /*
       * Refine only on the fine fill box regions of index space.
       * Note that we restrict the interpolation range to the
       * intersection of the fine fill box and the ghost box of
       * the scratch data component (i.e., the destination of the
       * interpolation).  This is needed for the case where data
       * components treated by the schedule have different ghost
       * cell widths since the fill boxes are generated using the
       * maximum ghost cell width.
       */

      if (d_refine_patch_strategy) {
         d_refine_patch_strategy->preprocessRefineBoxes(*fine_patch,
                                                        *crse_patch,
                                                        fill_boxes,
                                                        ratio);
      }

      typename tbox::List< xfer::LocallyActiveDataFillBox<DIM> >::Iterator 
         fbi(fill_boxes.getLocallyActiveDataBoxes());
      for ( ; fbi; fbi++) {
         for (typename tbox::List<const typename xfer::RefineClasses<DIM>::Data*>::Iterator
              vi(fbi().getActiveRefineVarData()); vi; vi++) {
            if ( !(vi()->d_oprefine).isNull() ) {
               const int scratch_id = vi()->d_scratch;
               if ( d_coarse_level_mgr->getPatchDataActive( hier::PatchDataId(scratch_id),
                                                            hier::PatchNumber(p()) ) ) {
#ifdef DEBUG_CHECK_ASSERTIONS
                  TBOX_ASSERT(fine_patch->checkAllocated(scratch_id));
                  TBOX_ASSERT(crse_patch->checkAllocated(scratch_id));
#endif
                  const hier::Box<DIM>& scratch_space =
                     fine_patch->getPatchData(scratch_id)->getGhostBox();
                  vi()->d_oprefine->refine(*fine_patch,
                                           *crse_patch,
                                           scratch_id,
                                           scratch_id,
                                           fbi().getBox()*scratch_space,
                                           ratio);
               }
            }
         } // loop over active variables
      } // loop over fill boxes
      
      if (d_refine_patch_strategy) {
         d_refine_patch_strategy->postprocessRefineBoxes(*fine_patch,
                                                         *crse_patch,
                                                         fill_boxes,
                                                         ratio);
      }
   }
}

/*
*************************************************************************
*									*
* Generate a communication schedule that moves data from the interiors	*
* of the source box on the source level into the specified fill box	*
* region of the destination level.  The fill box will typically be the	*
* interior plus ghost cells of the destination region. The schedule     *
* does NOT copy patch interiors from source to destination, it only     *
* fills the ghost region.                                               *
*									*
* This routine manages periodic boundary conditions using the offset	*
* shift capability of the box overlap calculation.  The periodic shift	*
* vector contains a nonzero entry for every dimension that is periodic	*
* and zero entries in the other dimensions.  The nonzero entry gives	*
* the periodic shift in that dimension.					*
*									*
* The source level should always be a level in the MFHierarchy.  It may *
* NOT be a temporary level.                                             *
*									*
*************************************************************************
*/

template<int DIM>
void LocallyActiveDataRefineSchedule<DIM>::generateCommunicationSchedule(
   tbox::Pointer<tbox::Schedule> coarse_priority_schedule,
   tbox::Pointer<tbox::Schedule> fine_priority_schedule,
   tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
   tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > dst_level_mgr,  
   tbox::Pointer< hier::PatchLevel<DIM> > src_level,
   tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > src_level_mgr,  
   const tbox::Array< xfer::LocallyActiveDataFillBoxSet<DIM> >& la_fill_boxes,
   tbox::Array< xfer::LocallyActiveDataFillBoxSet<DIM> >& la_unfilled_boxes,
   bool use_time_interpolation)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!coarse_priority_schedule.isNull());
   TBOX_ASSERT(!fine_priority_schedule.isNull());
   TBOX_ASSERT(!dst_level.isNull());
   TBOX_ASSERT(!dst_level_mgr.isNull());
   TBOX_ASSERT(dst_level_mgr->checkLevel(dst_level));
   TBOX_ASSERT(!src_level.isNull());
   TBOX_ASSERT(src_level_mgr->checkLevel(src_level));
   TBOX_ASSERT(la_fill_boxes.getSize() == dst_level->getNumberOfPatches());
   TBOX_ASSERT(la_unfilled_boxes.getSize() == dst_level->getNumberOfPatches());
#endif

   const int dst_npatches = dst_level->getNumberOfPatches();
   const hier::BoxArray<DIM>& dst_boxes = dst_level->getBoxes();
   const hier::ProcessorMapping& dst_mapping = dst_level->getProcessorMapping();
   tbox::Pointer< hier::PatchDescriptor<DIM> > dst_patch_descriptor =
      dst_level->getPatchDescriptor();

   const hier::BoxArray<DIM>& src_boxes = src_level->getBoxes();
   tbox::Pointer< hier::BoxTree<DIM> > src_box_tree = src_level->getBoxTree();
   tbox::Pointer< hier::PatchDescriptor<DIM> > src_patch_descriptor =
      src_level->getPatchDescriptor();

   hier::IntVector<DIM> growth =
      hier::IntVector<DIM>::max(d_max_scratch_gcw, getMaxDestinationGhosts());
   int max_gcw = tbox::MathUtilities<int>::Max(growth.max(), 1);
   hier::IntVector<DIM> dst_growth(max_gcw);

   /*
    * Determine boxes and data that cannot be filled from source level.
    */

   t_gen_comm_sched_unfilled->start();
   for (int dst_patch_id = 0; dst_patch_id < dst_npatches; dst_patch_id++) {

      hier::Box<DIM> dst_box_plus_ghosts = dst_boxes(dst_patch_id);
      dst_box_plus_ghosts.grow(dst_growth);
 
      tbox::Array<int> src_nabor_indices;
      src_box_tree->findOverlapIndices(src_nabor_indices, dst_box_plus_ghosts);

      const tbox::List<const typename xfer::RefineClasses<DIM>::Data*>& 
          dst_patch_active_data = la_fill_boxes[dst_patch_id].getUnionActiveRefineVarData();

      /*
       * Loop over source patches with possible overlap with current destination patch
       */

      hier::BoxList<DIM> tmp_unfilled_boxes = la_fill_boxes[dst_patch_id].getBoxList();

      const int src_len = src_nabor_indices.getSize();
      for (int spp = 0; spp < src_len; spp++) {
         int src_patch_id = src_nabor_indices[spp];

         const hier::Box<DIM>& src_box = src_boxes(src_patch_id);

         tbox::List<const typename xfer::RefineClasses<DIM>::Data*> missing_src_data;
         bool some_src_data_is_missing = false;
         for (typename tbox::List<const typename xfer::RefineClasses<DIM>::Data*>::Iterator
              vi(dst_patch_active_data); vi; vi++) {
            if (!src_level_mgr->getPatchDataActive( hier::PatchDataId(vi()->d_src), 
                                                    hier::PatchNumber(src_patch_id) ) ) {
               /*
                * IMPORTANT! Proper generation of refine data lists requires
                *            appendItem(); addItem() is incorrect!
                */
               missing_src_data.appendItem(vi());
               some_src_data_is_missing = true;
            }
         }

         typename tbox::List< hier::IntVector<DIM> >::Iterator
            sh(src_level->getShiftsForPatch(src_patch_id));

         bool zero_shift = true;

         while (sh || zero_shift) {

            hier::IntVector<DIM> shift(s_constant_zero_intvector);
            if (!zero_shift) {
               shift = sh();
            } 

            const hier::Box<DIM> shifted(hier::Box<DIM>::shift(src_box, shift));

            if (some_src_data_is_missing) {
               hier::BoxList<DIM> fill_boxes = la_fill_boxes[dst_patch_id].getBoxList();
               fill_boxes.intersectBoxes(shifted);
               for (typename hier::BoxList<DIM>::Iterator b(fill_boxes); b; b++) {
                  la_unfilled_boxes[dst_patch_id].
                     addLocallyActiveFillBox(b(), missing_src_data);
               }
            }

            tmp_unfilled_boxes.removeIntersections(shifted);

            if (!zero_shift) {
               sh++;
            } else {
               zero_shift = false;
            }

         } // while loop (sh || zero_shift)

      } // loop over source patches

      for (typename hier::BoxList<DIM>::Iterator b(tmp_unfilled_boxes); b; b++) {
         la_unfilled_boxes[dst_patch_id].
            addLocallyActiveFillBox(b(), dst_patch_active_data);
      }

   } // loop over destination patches
   t_gen_comm_sched_unfilled->stop();

   t_gen_comm_sched_trans->start();

   if (d_src_masks.getNumberOfBoxes() < d_max_fill_boxes) {
      d_src_masks.resizeBoxArray(0);
      d_src_masks.resizeBoxArray(d_max_fill_boxes);
   }
   if (d_overlaps.getSize() < d_max_fill_boxes) {
      d_overlaps.setNull();
      d_overlaps.resizeArray(d_max_fill_boxes);
   }

   const int num_equiv_classes =
      d_refine_classes->getNumberEquivalenceClasses();

   const bool same_level = (dst_level == src_level);

   for (int dst_patch_id = 0; dst_patch_id < dst_npatches; dst_patch_id++) {

      const hier::Box<DIM>& dst_box = dst_boxes(dst_patch_id); 

      hier::Box<DIM> dst_box_plus_ghosts = dst_box;
      dst_box_plus_ghosts.grow(dst_growth);

      tbox::Array<int> src_nabor_indices;
      if (dst_mapping.isMappingLocal(dst_patch_id)) {
         src_box_tree->findOverlapIndices(
            src_nabor_indices, dst_box_plus_ghosts);
      } else {
         src_box_tree->findLocalOverlapIndices(
            src_nabor_indices, dst_box_plus_ghosts);
      }

      const xfer::LocallyActiveDataFillBoxSet<DIM>& fill_boxes = la_fill_boxes[dst_patch_id];
      const int num_fill_boxes = fill_boxes.getNumberOfBoxes();

      int src_len = src_nabor_indices.getSize();
      for (int spp = 0; spp < src_len; spp++) {

         int src_patch_id = src_nabor_indices[spp];

         /*
          * Determine which equivalence classes are active; i.e., which
          * contain data items to be copied from given source patch to
          * given destination patch.  For an equivalence class to be active,
          * both the source data must be active on the source patch and the
          * destination data must be active on the destination patch.
          */

         bool any_equivalence_classes_active = false;

         tbox::Array<bool> active_equivalence_class(num_equiv_classes);
         tbox::Array<bool> active_ritem(d_number_refine_items);

         for (int nc = 0; nc < num_equiv_classes; nc++) {
            active_equivalence_class[nc] = false;
            for (typename tbox::List< typename xfer::RefineClasses<DIM>::Data >::Iterator
                    l(d_refine_classes->getIterator(nc)); l; l++) {
               int ritem_count = l().d_tag;
               active_ritem[ritem_count] = false;
               if ( src_level_mgr->getPatchDataActive( 
                                   hier::PatchDataId(l().d_src), 
                                   hier::PatchNumber(src_patch_id) ) &&
                    dst_level_mgr->getPatchDataActive( 
                                   hier::PatchDataId(l().d_dst), 
                                   hier::PatchNumber(dst_patch_id) ) ) {
                  active_equivalence_class[nc] = true;
                  active_ritem[ritem_count] = true;
               }
               
            }
            any_equivalence_classes_active |= active_equivalence_class[nc];
         }

         if (any_equivalence_classes_active) {

            const hier::Box<DIM>& src_box = src_boxes(src_patch_id);

            const bool same_patch = ( same_level && (dst_patch_id == src_patch_id) );

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

               bool same_patch_no_shift = (same_patch && zero_shift);            

               /*
                * Calculate the shift and the shifted source box.
                */

               hier::IntVector<DIM> shift(0);
               if (!zero_shift) {
                  shift = sh();
               }
   
               const hier::Box<DIM> shifted(hier::Box<DIM>::shift(src_box, shift));

               for (int nc = 0; nc < num_equiv_classes; nc++) {

                  if (active_equivalence_class[nc]) {

                     const typename xfer::RefineClasses<DIM>::Data& rep_item =
                        d_refine_classes->getClassRepresentative(nc);

                     const int rep_item_dst_id = rep_item.d_scratch;
                     const int rep_item_src_id = rep_item.d_src;

                     tbox::Pointer< hier::PatchDataFactory<DIM> > src_pdf =
                        src_patch_descriptor->getPatchDataFactory(rep_item_src_id);
                     tbox::Pointer< hier::PatchDataFactory<DIM> > dst_pdf =
                        dst_patch_descriptor->getPatchDataFactory(rep_item_dst_id);
   
                     const hier::IntVector<DIM>& dst_gcw = 
                        dst_pdf->getGhostCellWidth();

                     /*
                      * Iterate over boxes in fill box list for destination patch.
                      * Then, calculate the src_mask and overlap for each fill box.
                      * For each equivalence class, this loop is executed once.
                      */

                     int box_num = 0;
                     for (typename hier::BoxList<DIM>::Iterator b(fill_boxes.getBoxList()); 
                             b; b++) {

                        const hier::Box<DIM>& fill_box = b();

                        hier::Box<DIM> dst_fill_box(hier::Box<DIM>::grow(dst_box, dst_gcw));
                        dst_fill_box = dst_fill_box * fill_box;
   
                        hier::Box<DIM> test_mask(dst_fill_box*shifted);
                        if ( test_mask.empty() &&
                             (dst_gcw == s_constant_zero_intvector) &&
                             dst_pdf->dataLivesOnPatchBorder() ) {
                           hier::Box<DIM> tmp_dst_fill_box(
                                          hier::Box<DIM>::grow(dst_fill_box,
                                                               s_constant_one_intvector));
                           test_mask = tmp_dst_fill_box * shifted;
                        }
                        hier::Box<DIM> src_mask( hier::Box<DIM>::shift( test_mask,-shift) );
   
                        tbox::Pointer< hier::BoxOverlap<DIM> > overlap =
                           dst_pdf->getBoxGeometry(dst_box)
                                  ->calculateOverlap(
                                     *src_pdf->getBoxGeometry(src_box),
                                     src_mask,
                                     true, shift);

#ifdef DEBUG_CHECK_ASSERTIONS
                        if (overlap.isNull()) {
                           TBOX_ERROR(
                              "Internal xfer::LocallyActiveDataRefineSchedule<DIM> error..."
                              << "\n Overlap is NULL for "
                              << "\n src box = " << src_box
                              << "\n dst box = " << dst_box
                              << "\n src mask = " << src_mask << std::endl);
                        }
#endif

                        d_src_masks.getBox(box_num) = src_mask;
                        d_overlaps[box_num] = overlap;
                        box_num++;

                     }  // iterate over fill boxes

                     /*
                      * Iterate over components in refine description list
                      */

                     for (typename tbox::List<typename xfer::RefineClasses<DIM>::Data >::Iterator
                             l(d_refine_classes->getIterator(nc)); l; l++) {
                        const int dst_id = l().d_scratch;
                        const int src_id = l().d_src;
                        int ritem_count = l().d_tag;

                        /*
                         * If the src and dst patches, levels, and components are the
                         * same, and there is no shift, the data exchange is unnecessary.
                         */
                        if ( (!same_patch_no_shift || (dst_id != src_id)) &&
                              active_ritem[ritem_count] ) {

                           /*
                            * Iterate over the fill boxes and create transactions
                            * for each box that has a non-empty overlap.
                            */
                           for (int i = 0; i < num_fill_boxes; i++) {

                              if (!d_overlaps[i]->isOverlapEmpty()) {

                                 bool do_time_interpolation = 
                                    (use_time_interpolation &&
                                     l().d_time_interpolate);

                                 tbox::Pointer<tbox::Transaction> transaction =
                                    d_transaction_factory->allocate(dst_level,
                                                                    src_level,
                                                                    d_overlaps[i],
                                                                    dst_patch_id,
                                                                    src_patch_id,
                                                                    ritem_count,
                                                                    d_src_masks.getBox(i),
                                                                    do_time_interpolation);

                                 if (l().d_fine_bdry_reps_var) {
                                    if (same_patch) {
                                       fine_priority_schedule->
                                          addTransaction(transaction);
                                    } else {
                                       fine_priority_schedule->
                                          appendTransaction(transaction);
                                    }
                                 } else {
                                    if (same_patch) {
                                       coarse_priority_schedule->
                                          addTransaction(transaction);
                                    } else {
                                       coarse_priority_schedule-> 
                                          appendTransaction(transaction);
                                    }
                                 }

                              }  // if overlap not empty

                           }  // iterate over fill_boxes

                        }  // if copy transaction is needed 
                           // (i.e., not same patch data object)

                     }  // iterate over refine components in equivalence class

                  }  // if equivalence class is active

               }  // iterate over refine equivalence classes

               if (!zero_shift) {
                  sh++;
               } else {
                  zero_shift = false;
               }

            }  // while (sh || zero_shift)

         } // if any equivalence class active

      } // loop over source patches

   } // loop over destination patches
   t_gen_comm_sched_trans->stop();

}

/*
*************************************************************************
*									*
* Allocate the default fill boxes for the specified destination array.	*
* The fill box for a particular patch involves the box associated with	*
* that patch grown by the maximum ghost cell width of all of the patch	*
* data destination components on that patch.  Also, if any              *
* communications have destination variables quantities that have coarse *
* data taking precedence at course-fine boundaries, make sure that the  *
* box associated with each patch is grown by at least one cell in each  *
* coordinate direction.  Finally, create an array of integer lists      *
* specifying the active data on destination level patches to be used    *
* later to determine which patches on the dest level should be filled.  *
*									*
*************************************************************************
*/

template<int DIM>
void LocallyActiveDataRefineSchedule<DIM>::allocateDefaultFillBoxes(
   tbox::Array< xfer::LocallyActiveDataFillBoxSet<DIM> >& la_fill_boxes,
   tbox::Pointer< hier::PatchLevel<DIM> > level,
   tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > level_mgr,
   const hier::IntVector<DIM>& fill_ghost_width)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(la_fill_boxes.getSize() == 0);
   TBOX_ASSERT(!level.isNull());
   TBOX_ASSERT(!level_mgr.isNull());
   TBOX_ASSERT(!(fill_ghost_width < s_constant_zero_intvector));
#endif

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

   bool need_nontrivial_ghosts = false;
   for (int iri = 0; iri < d_number_refine_items; iri++) {
      if (!(d_refine_items[iri]->d_fine_bdry_reps_var)) {
         need_nontrivial_ghosts = true;
      }
   }

   hier::IntVector<DIM> gcw = fill_ghost_width;
   if (need_nontrivial_ghosts) {
      gcw.max(s_constant_one_intvector);
   }

   d_max_fill_boxes = 0;

   const hier::BoxArray<DIM>& boxes  = level->getBoxes();
   const int nboxes                  = boxes.getNumberOfBoxes();

   la_fill_boxes.resizeArray(nboxes);

   for (int p = 0; p < nboxes; p++) {
      tbox::List<const typename xfer::RefineClasses<DIM>::Data*> var_data;
      for (int iri = 0; iri < d_number_refine_items; iri++) {
         if ( level_mgr->getPatchDataActive( 
                         hier::PatchDataId(d_refine_items[iri]->d_dst), 
                         hier::PatchNumber(p) ) ) {
            /*
             * IMPORTANT! Proper generation of refine data lists requires
             *            appendItem(); addItem() is incorrect!
             */
            var_data.appendItem(d_refine_items[iri]);
         }
      }

      la_fill_boxes[p].resetLocallyActiveFillBoxes(hier::Box<DIM>::grow(boxes(p), gcw),
                                                   var_data);

      d_max_fill_boxes =
         tbox::MathUtilities<int>::Max(d_max_fill_boxes,
                                       la_fill_boxes[p].getNumberOfBoxes());

   }

}

/*
**************************************************************************
*                                                                        *
* Calculate the maximum ghost cell width of all destination patch data   *
* components.                                                            *
*                                                                        *
**************************************************************************
*/

template<int DIM>
hier::IntVector<DIM>
LocallyActiveDataRefineSchedule<DIM>::getMaxDestinationGhosts() const
{
   hier::IntVector<DIM> gcw(s_constant_zero_intvector);
   tbox::Pointer< hier::PatchDescriptor<DIM> > pd = d_dst_level->getPatchDescriptor();

   for (int iri = 0; iri < d_number_refine_items; iri++) {
      const int dst_id = d_refine_items[iri]->d_dst;
      gcw.max(pd->getPatchDataFactory(dst_id)->getGhostCellWidth());
   }

   return(gcw);
}

/*
**************************************************************************
*                                                                        *
* Calculate the maximum ghost cell width of all scratch patch data       *
* components.                                                            *
*                                                                        *
**************************************************************************
*/

template<int DIM>
hier::IntVector<DIM>
LocallyActiveDataRefineSchedule<DIM>::getMaxScratchGhosts() const
{
   hier::IntVector<DIM> gcw(s_constant_zero_intvector);
   tbox::Pointer< hier::PatchDescriptor<DIM> > pd = d_dst_level->getPatchDescriptor();

   for (int iri = 0; iri < d_number_refine_items; iri++) {
      const int scratch_id = d_refine_items[iri]->d_scratch;
      gcw.max(pd->getPatchDataFactory(scratch_id)->getGhostCellWidth());
   }

   return(gcw);
}

/*
**************************************************************************
*                                                                        *
* Calculate the maximum ghost cell width required for all stencils.      *
*                                                                        *
**************************************************************************
*/

template<int DIM>
hier::IntVector<DIM> LocallyActiveDataRefineSchedule<DIM>::getMaxStencilGhosts() const
{
   hier::IntVector<DIM> gcw(s_constant_zero_intvector);
   if (d_refine_patch_strategy) {
      gcw = d_refine_patch_strategy->getRefineOpStencilWidth();
   }

   for (int iri = 0; iri < d_number_refine_items; iri++) {
      if (!(d_refine_items[iri]->d_oprefine.isNull())) {
         gcw.max(d_refine_items[iri]->d_oprefine->getStencilWidth());
      }
   }

   return(gcw);
}

/*
*************************************************************************
*                                                                       *
* Private member function to initialize data members for hierarchy
* information and refine data items.
*                                                                       *
*************************************************************************
*/

template<int DIM>
void
LocallyActiveDataRefineSchedule<DIM>::initializeDomainAndGhostInformation(
   bool recursive_schedule)
{

   d_max_scratch_gcw = getMaxScratchGhosts();
   d_max_stencil_gcw = getMaxStencilGhosts();

   if (recursive_schedule) {
      d_boundary_fill_ghost_width = d_max_stencil_gcw;
      d_force_boundary_fill = (d_boundary_fill_ghost_width.max() > 0);
   } else {
      d_boundary_fill_ghost_width = getMaxDestinationGhosts();
      d_force_boundary_fill = false;
   }

   tbox::Pointer< hier::GridGeometry<DIM> > grid_geom = d_dst_level->getGridGeometry();
   const hier::IntVector<DIM>& ratio_to_level_zero = d_dst_level->getRatio();

   d_domain_is_one_box = grid_geom->getDomainIsSingleBox();

   if (d_domain_is_one_box) {
      hier::BoxArray<DIM> domain;
      grid_geom->computePhysicalDomain(domain, ratio_to_level_zero);
      d_domain_box = domain.getBox(0);
      for (int i = 1; i < domain.getNumberOfBoxes(); i++) {
         d_domain_box += domain.getBox(i);
      }
   }

   d_periodic_shift = grid_geom->getPeriodicShift(ratio_to_level_zero);

   d_num_periodic_directions = 0;
   for (int d = 0; d < DIM; d++) {
      if (d_periodic_shift(d)) {
        d_num_periodic_directions++;
      }
   }

}

/*
*************************************************************************
*                                                                       *
* Private utility function to set up local array of refine items.       *
*                                                                       *
*************************************************************************
*/

template<int DIM>
void LocallyActiveDataRefineSchedule<DIM>::setRefineItems(
   const tbox::Pointer< xfer::RefineClasses<DIM> > refine_classes)
{

   clearRefineItems();

   d_refine_classes = refine_classes;

   const int num_refine_classes =
      d_refine_classes->getNumberEquivalenceClasses();

   int nc;
   for (nc = 0; nc < num_refine_classes; nc++) {
      for (typename tbox::List<typename xfer::RefineClasses<DIM>::Data>::Iterator 
              l(d_refine_classes->getIterator(nc)); l; l++) {
         d_number_refine_items++;
      }
   }

   d_refine_items =
      new const typename xfer::RefineClasses<DIM>::Data*[d_number_refine_items];

   int ircount = 0;
   for (nc = 0; nc < num_refine_classes; nc++) {
      for (typename tbox::List<typename xfer::RefineClasses<DIM>::Data>::Iterator 
              l(d_refine_classes->getIterator(nc)); l; l++) {
         l().d_tag = ircount;
         d_refine_items[ircount] = &(l());
         ircount++;
      }
   }

}

/*
*************************************************************************
*                                                                       *
* Private utility function used during initial schedule set up to       *
* check whether patch data entries have proper number of ghost cells.   *
* In particular, each scratch data entry must have at least as many     *
* ghost cells as the user-defined refine operator stencil width.        *
* Other checks are performed in the                                     *
* xfer::RefineClasses<DIM>::checkRefineItem() routine.                       *
*                                                                       *
*************************************************************************
*/

template<int DIM>
void LocallyActiveDataRefineSchedule<DIM>::initialCheckRefineClassItems() const
{
   tbox::Pointer< hier::PatchDescriptor<DIM> > pd = d_dst_level->getPatchDescriptor();

   hier::IntVector<DIM> user_gcw(s_constant_zero_intvector);
   if (d_refine_patch_strategy) {
      user_gcw = d_refine_patch_strategy->getRefineOpStencilWidth();
   }

   if (user_gcw > s_constant_zero_intvector) {

      for (int iri = 0; iri < d_number_refine_items; iri++) {

         const typename xfer::RefineClasses<DIM>::Data* const ref_item = d_refine_items[iri];

#ifdef DEBUG_CHECK_ASSERTIONS
         if (d_refine_classes->checkRefineItem(*ref_item, pd)) {
#endif

            const int scratch = ref_item->d_scratch;
            const hier::IntVector<DIM>& scratch_gcw(pd->getPatchDataFactory(scratch)->
                                               getGhostCellWidth());

            if (user_gcw > scratch_gcw) {
               TBOX_ERROR("Bad data given to xfer::LocallyActiveDataRefineSchedule<DIM>...\n"
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
*                                                                        *
* Private utility function to clear array of refine items.               *
*                                                                        *
**************************************************************************
*/

template<int DIM>
void LocallyActiveDataRefineSchedule<DIM>::clearRefineItems()
{
   if (d_refine_items) {
      for (int iri = 0; iri < d_number_refine_items; iri++) {
         d_refine_items[iri] = (typename xfer::RefineClasses<DIM>::Data*)NULL;
      }
      delete [] d_refine_items;
      d_refine_items = (const typename xfer::RefineClasses<DIM>::Data**)NULL;
      d_number_refine_items   = 0;
   }
}

/*
*************************************************************************
*									*
* Print class data to the specified output stream.			*
*									*
*************************************************************************
*/

template<int DIM>
void LocallyActiveDataRefineSchedule<DIM>::printClassData(std::ostream& stream) const
{
   stream << "LocallyActiveDataRefineSchedule<DIM>::printClassData()\n";
   stream << "--------------------------------------\n";
   stream << "s_schedule_generation_method = "
          << s_schedule_generation_method << std::endl;

   stream << "d_refine_classes ..." << std::endl;
   d_refine_classes->printClassData(stream);

   stream << "d_dst_level = " << (hier::PatchLevel<DIM>*)d_dst_level << std::endl;
   stream << "d_dst_level_mgr = " 
          << (hier::LocallyActiveDataPatchLevelManager<DIM>*)d_dst_level_mgr << std::endl;

   stream << "d_refine_patch_strategy = " 
          << (xfer::LocallyActiveDataRefinePatchStrategy<DIM>*)d_refine_patch_strategy << std::endl;

   stream << "d_max_stencil_gcw = " << d_max_stencil_gcw << std::endl; 
   stream << "d_max_scratch_gcw = " << d_max_scratch_gcw << std::endl; 
   stream << "d_boundary_fill_ghost_width = " << d_boundary_fill_ghost_width << std::endl; 

   stream << "d_force_boundary_fill = " << d_force_boundary_fill << std::endl;
   stream << "d_domain_is_one_box = " << d_domain_is_one_box << std::endl;
   stream << "d_domain_box = " << d_domain_box << std::endl;
   stream << "d_num_periodic_directions = " << d_num_periodic_directions << std::endl;
   stream << "d_periodic_shift = " << d_periodic_shift << std::endl;

   stream << "d_coarse_level = " << (hier::PatchLevel<DIM>*)d_coarse_level << std::endl;
   stream << "d_coarse_level_mgr = " 
          << (hier::LocallyActiveDataPatchLevelManager<DIM>*)d_coarse_level_mgr << std::endl;

   stream << "d_coarse_to_fine_mapping size = "
          << d_coarse_to_fine_mapping.size() << std::endl;

   stream << "d_la_fine_fill_boxes size = "
          << d_la_fine_fill_boxes.size() << std::endl;

   if (!d_coarse_priority_level_schedule.isNull()) {
      stream << "Printing coarse priority schedule..." << std::endl;
      d_coarse_priority_level_schedule->printClassData(stream);
   } else {
      stream << "coarse priority schedule is null" << std::endl;
   }
   if (!d_fine_priority_level_schedule.isNull()) {
      stream << "Printing fine priority schedule..." << std::endl;
      d_fine_priority_level_schedule->printClassData(stream);
   } else {
      stream << "fine priority schedule is null" << std::endl;
   }
   if (!d_coarse_schedule.isNull()) {
      stream << "Printing coarse refine schedule..." << std::endl;
      d_coarse_schedule->printClassData(stream);
   } else {
      stream << "coarse refine schedule is null" << std::endl;
   }
}

}
}

#endif
