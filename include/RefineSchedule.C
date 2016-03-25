//
// File:         $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/transfer/datamovers/standard/RefineSchedule.C $
// Package:      SAMRAI data transfer
// Copyright:    (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:     $LastChangedRevision: 2141 $
// Modified:     $LastChangedDate: 2008-04-23 08:36:33 -0700 (Wed, 23 Apr 2008) $
// Description:  Refine schedule for data transfer between AMR levels
//

#ifndef included_xfer_RefineSchedule_C
#define included_xfer_RefineSchedule_C

#include "RefineSchedule.h"
#include "BoxGeometry.h"
#include "BoxOverlap.h"
#include "BoxTop.h"
#include "BoxGraph.h"
#include "BoxTree.h"
#include "BoxUtilities.h"
#include "ComponentSelector.h"
#include "Patch.h"
#include "PatchData.h"
#include "PatchGeometry.h"
#include "tbox/ArenaManager.h"
#include "tbox/InputManager.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"
#include "RefinePatchStrategy.h"


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
   RefineSchedule<DIM>::s_constant_zero_intvector = hier::IntVector<DIM>(0);
template<int DIM> const hier::IntVector<DIM>
   RefineSchedule<DIM>::s_constant_one_intvector = hier::IntVector<DIM>(1);
template<int DIM> std::string
   RefineSchedule<DIM>::s_schedule_generation_method = "BOX_TREE";

/*
 * ************************************************************************
 *                                                                        *
 * Static function to set box intersection algorithm for schedules.       *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM> void RefineSchedule<DIM>::setScheduleGenerationMethod(
   const std::string& method)
{
   if ( !((method == "ORIG_NSQUARED") ||
          (method == "BOX_GRAPH") ||
          (method == "BOX_TREE")) ) {
      TBOX_ERROR("RefineSchedule<DIM>::setScheduleGenerationMethod\n"
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

template<int DIM>  RefineSchedule<DIM>::RefineSchedule(
   const std::string& fill_pattern,
   tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
   tbox::Pointer< hier::PatchLevel<DIM> > src_level,
   const tbox::Pointer< xfer::RefineClasses<DIM> > refine_classes,
   tbox::Pointer< xfer::RefineTransactionFactory<DIM> > transaction_factory,
   xfer::RefinePatchStrategy<DIM>* patch_strategy,
   bool use_time_interpolation)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst_level.isNull());
   TBOX_ASSERT(!src_level.isNull());
   TBOX_ASSERT(!refine_classes.isNull());
   TBOX_ASSERT(!transaction_factory.isNull());
#endif

   t_fill_data =tbox::TimerManager::getManager() ->
      getTimer("xfer::RefineSchedule::fillData()");
   t_gen_sched_n_squared = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::generateCommunicationScheduleNSquared()");
   t_gen_sched_box_graph = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::generateCommunicationScheduleBoxGraph()");
   t_gen_sched_box_tree = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::generateCommunicationScheduleBoxTree()");

   /*
    * Initial values; some will change in setup operations.
    */

   d_transaction_factory = transaction_factory;

   d_dst_level                 = dst_level;

   d_refine_patch_strategy     = patch_strategy;

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

   tbox::Array< xfer::FillBoxSet<DIM> > fill_boxes;
   allocateFillBoxes(fill_pattern,
                     fill_boxes,
                     d_dst_level,
                     d_boundary_fill_ghost_width);

   tbox::Array< xfer::FillBoxSet<DIM> >
      unfilled_boxes(dst_level->getNumberOfPatches());

   tbox::Pointer<tbox::Timer> t_gen_comm_sched =
      tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::generate_comm_schedule");

   t_gen_comm_sched->start();
   generateCommunicationSchedule(d_coarse_priority_level_schedule,
				 d_fine_priority_level_schedule,
				 d_dst_level,
				 src_level,
				 fill_boxes,
				 unfilled_boxes,
				 use_time_interpolation);
   t_gen_comm_sched->stop();
}

/*
**************************************************************************
*
* Create a refine schedule that moves data from the source level into
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
* if the next coarser level is < 0 (that is, there is no coarser level).
*
**************************************************************************
*/

template<int DIM>  RefineSchedule<DIM>::RefineSchedule(
   const std::string &fill_pattern,
   tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
   tbox::Pointer< hier::PatchLevel<DIM> > src_level,
   int next_coarser_level,
   tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
   const tbox::Pointer< xfer::RefineClasses<DIM> > refine_classes,
   tbox::Pointer< xfer::RefineTransactionFactory<DIM> > transaction_factory,
   xfer::RefinePatchStrategy<DIM>* patch_strategy,
   bool use_time_interpolation)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst_level.isNull());
   TBOX_ASSERT((next_coarser_level < 0) || !hierarchy.isNull());
   TBOX_ASSERT(!refine_classes.isNull());
   TBOX_ASSERT(!transaction_factory.isNull());
#endif

   t_fill_data =tbox::TimerManager::getManager() ->
      getTimer("xfer::RefineSchedule::fillData()");
   t_gen_sched_n_squared = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::generateCommunicationScheduleNSquared()");
   t_gen_sched_box_graph = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::generateCommunicationScheduleBoxGraph()");
   t_gen_sched_box_tree = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::generateCommunicationScheduleBoxTree()");

   /*
    * Initial values; some will change in setup operations.
    */

   d_transaction_factory = transaction_factory;

   d_dst_level                 = dst_level;

   d_refine_patch_strategy     = patch_strategy;

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

   tbox::Array< xfer::FillBoxSet<DIM> > fill_boxes;
   allocateFillBoxes(fill_pattern,
                     fill_boxes,
                     d_dst_level,
                     d_boundary_fill_ghost_width);

   tbox::Pointer<tbox::Timer> t_finish_sched_const =
      tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::finish_schedule_const");

   bool skip_first_generate_schedule =
      (fill_pattern == "FILL_LEVEL_BORDERS_ONLY");

   t_finish_sched_const->start();
   finishScheduleConstruction(src_level,
			      next_coarser_level,
			      hierarchy,
			      fill_boxes,
                              use_time_interpolation,
                              skip_first_generate_schedule);
   t_finish_sched_const->stop();
}

/*
**************************************************************************
*                                                                        *
* This private constructor creates a refine schedule that moves data     *
* into the destination only on the specified fill boxes.                 *
*                                                                        *
**************************************************************************
*/

template<int DIM>  RefineSchedule<DIM>::RefineSchedule(
   tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
   tbox::Pointer< hier::PatchLevel<DIM> > src_level,
   int next_coarser_level,
   tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
   const tbox::Pointer< xfer::RefineClasses<DIM> > refine_classes,
   tbox::Pointer< xfer::RefineTransactionFactory<DIM> > transaction_factory,
   const tbox::Array< xfer::FillBoxSet<DIM> >& fill_boxes,
   xfer::RefinePatchStrategy<DIM>* patch_strategy)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst_level.isNull());
   TBOX_ASSERT((next_coarser_level < 0) || !hierarchy.isNull());
   TBOX_ASSERT(!refine_classes.isNull());
   TBOX_ASSERT(!transaction_factory.isNull());
   TBOX_ASSERT(fill_boxes.getSize() == dst_level->getNumberOfPatches());
#endif

   t_fill_data =tbox::TimerManager::getManager() ->
      getTimer("xfer::RefineSchedule::fillData()");
   t_gen_sched_n_squared = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::generateCommunicationScheduleNSquared()");
   t_gen_sched_box_graph = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::generateCommunicationScheduleBoxGraph()");
   t_gen_sched_box_tree = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::generateCommunicationScheduleBoxTree()");

   /*
    * Initial values; some will change in setup operations.
    * Note that we do not check refine items here, since this
    * constructor is private and called recursively (i.e., the
    * items have been checked already).
    */

   d_transaction_factory = transaction_factory;

   d_dst_level                 = dst_level;

   d_refine_patch_strategy     = patch_strategy;

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

   const int nfill_boxes = fill_boxes.getSize();
   for (int p = 0; p < nfill_boxes; p++) {
      d_max_fill_boxes =
	    tbox::MathUtilities<int>::Max( d_max_fill_boxes,
				           fill_boxes[p].getNumberOfBoxes() );
   }

   bool use_time_interpolation = true;

   finishScheduleConstruction(src_level,
			      next_coarser_level,
			      hierarchy,
			      fill_boxes,
			      use_time_interpolation,
                              false);

}

/*
**************************************************************************
*                                                                        *
* The destructor for the refine schedule class implicitly deallocates    *
* all of the data associated with the communication schedule.            *
*                                                                        *
**************************************************************************
*/

template<int DIM>  RefineSchedule<DIM>::~RefineSchedule()
{
   clearRefineItems();

   d_dst_level.setNull();

   d_refine_patch_strategy = (xfer::RefinePatchStrategy<DIM>*)NULL;

   d_transaction_factory.setNull();

   d_coarse_priority_level_schedule.setNull();
   d_fine_priority_level_schedule.setNull();

   d_coarse_schedule.setNull();
   d_coarse_level.setNull();

   t_fill_data.setNull();
   t_gen_sched_n_squared.setNull();
   t_gen_sched_box_graph.setNull();
   t_gen_sched_box_tree.setNull();
}

     /*
*************************************************************************
*                                                                       *
* Reset schedule with new set of refine items.                          *
*                                                                      *
************************************************************************
*/

template<int DIM> void RefineSchedule<DIM>::reset(
   const tbox::Pointer< xfer::RefineClasses<DIM> > refine_classes)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!refine_classes.isNull());
#endif
   setRefineItems(refine_classes);
   if (!d_coarse_schedule.isNull()) {
      d_coarse_schedule->reset(refine_classes);
   }
}

/*
**************************************************************************
*                                                                        *
* Return const pointer to equivalence classes used in schedule.          *
*                                                                        *
**************************************************************************
*/

template<int DIM> const tbox::Pointer< xfer::RefineClasses<DIM> >&
RefineSchedule<DIM>::getEquivalenceClasses() const
{
   return(d_refine_classes);
}

template<int DIM>
const hier::IntVector<DIM>&
RefineSchedule<DIM>::getBoundaryFillGhostWidth() const
{
   return(d_boundary_fill_ghost_width);
}


/*
**************************************************************************
*                                                                        *
* Finish construction of the communication schedule when filling from    *
* coarser levels.  This routine guarantees that the fill boxes will be   *
* filled from the same level or coarser levels.  If the fill boxes can   *
* not be filled, then the routine aborts with an internal error.         *
*                                                                        *
**************************************************************************
*/

template<int DIM> void RefineSchedule<DIM>::finishScheduleConstruction(
   tbox::Pointer< hier::PatchLevel<DIM> > src_level,
   int next_coarser_level,
   tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
   const tbox::Array< xfer::FillBoxSet<DIM> >& fill_boxes,
   bool use_time_interpolation,
   bool skip_generate_schedule)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((next_coarser_level < 0) || !hierarchy.isNull());
   TBOX_ASSERT(fill_boxes.getSize() == d_dst_level->getNumberOfPatches());
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
   tbox::Array< xfer::FillBoxSet<DIM> > unfilled_boxes(dst_npatches);

   if ( !src_level.isNull() && !skip_generate_schedule ) {

      tbox::Pointer<tbox::Timer> t_gen_comm_sched =
         tbox::TimerManager::getManager()->
         getTimer("xfer::RefineSchedule::generate_comm_schedule");

      t_gen_comm_sched->start();
      generateCommunicationSchedule(d_coarse_priority_level_schedule,
				    d_fine_priority_level_schedule,
				    d_dst_level,
				    src_level,
				    fill_boxes,
				    unfilled_boxes,
				    use_time_interpolation);
      t_gen_comm_sched->stop();

   } else {
      for (int p = 0; p < dst_npatches; p++) {
	 unfilled_boxes[p] = fill_boxes[p];
      }
   }

   /*
    * Remove pieces of the boxes to be filled from coarser levels that
    * live outside the physical domain in the non-periodic directions.
    * Then, check whether there still exist boxes to be filled from
    * coarser levels.  Any resulting boxes will be used to generate a
    * temporary coarse patch level for acquiring data to refine.
    */

   hier::IntVector<DIM> big_grow_vector(0);
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
         if (unfilled_boxes[p].getNumberOfBoxes() > 0) {
            if (d_dst_level->patchTouchesRegularBoundary(p)) {
               unfilled_boxes[p].intersectBoxes(shear_domain);
               if (!d_domain_is_one_box) {
                  unfilled_boxes[p].simplifyBoxes();
               }
               if (unfilled_boxes[p].getNumberOfBoxes() > 0 ) {
                  need_to_fill = true;
               }
            } else {
               need_to_fill = true;
            }
         }
      }

   } else {  // !do_shearing

      for (int p = 0; p < dst_npatches; p++) {
	 if (unfilled_boxes[p].getNumberOfBoxes() > 0) {
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
       * If there are no coarser levels in the hierarchy or the hierarchy
       * is null, then throw an error.  Something is messed up someplace and
       * code execution cannot proceed.
       */

      if (next_coarser_level < 0) {
	 TBOX_ERROR("Internal RefineSchedule<DIM> error..."
		    << "\n In finishScheduleConstruction() -- "
		    << "\n No coarser levels...could not fill from coarser."
		    << std::endl);
      } else {
         if (hierarchy.isNull()) {
	    TBOX_ERROR("Internal RefineSchedule<DIM> error..."
                       << "\n In finishScheduleConstruction() -- "
                       << "\n Need to fill from coarser hierarchy level and \n"
                       << "hierarchy is unavailable." << std::endl);
         }
      }

      /*
       * Calculate the ratio to the next coarser level in the hierarchy
       * and cache the maximum stencil ghost cell width.
       */

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
	 if (unfilled_boxes[fp].getNumberOfBoxes() > 0) {

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
		      unfilled_boxes[fp].getBoundingBox());

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

      const int ncoarse = coarse_to_fine_mapping.getNumberOfItems();

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

      d_fine_fill_boxes.resizeArray(ncoarse);
      tbox::Array< xfer::FillBoxSet<DIM> > coarse_fill_boxes(ncoarse);

      for (int cp = 0; cp < ncoarse; cp++) {
	 const int fp = d_coarse_to_fine_mapping[cp];

	 const hier::Box<DIM> coarse_box = coarse_box_array[cp];
	 const hier::Box<DIM> fine_box = hier::Box<DIM>::refine(coarse_box, ratio);

	 d_fine_fill_boxes[cp] = unfilled_boxes[fp];

	 d_fine_fill_boxes[cp].intersectBoxes(fine_box);

	 hier::BoxList<DIM> coarse_to_fill;
	 if (d_coarse_level->patchTouchesRegularBoundary(cp)) {
	    coarse_to_fill.addItem(coarse_box);
	 } else {
	    coarse_to_fill = d_fine_fill_boxes[cp].getBoxList();
	    coarse_to_fill.coarsen(ratio);
	 }

	 coarse_to_fill.grow(d_max_stencil_gcw);

	 coarse_fill_boxes[cp].resetFillBoxes(coarse_to_fill);
      }

      /*
       * Recursively fill the coarse schedule using the private
       * refine schedule constructor.
       */

      d_coarse_schedule = new RefineSchedule<DIM>(d_coarse_level,
						   coarser_level_from_hierarchy,
						   next_coarser_level-1,
						   hierarchy,
						   d_refine_classes,
                                                   d_transaction_factory,
						   coarse_fill_boxes,
						   d_refine_patch_strategy);

   }

}

/*
**************************************************************************
*                                                                        *
* Execute the stored communication schedule that moves data into the    *
* destination component of the destination level.  The algorithm is as   *
* follows:                                                               *
*                                                                        *
*   (1) Allocate scratch space on the destination level if it does       *
*       not exist.                                                       *
*   (2) Call the recursive fill function that will fill the scratch      *
*       space of the destination level.                                  *
*   (3) If the scratch and destination spaces are not the same,          *
*       then copy locally from the scratch into the destination.         *
*   (4) Deallocate any previously allocated data.                        *
*                                                                        *
**************************************************************************
*/

template<int DIM> void RefineSchedule<DIM>::fillData(
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
   allocateScratchSpace(d_dst_level, fill_time, allocate_vector);

   /*
    * Begin the recursive algorithm that fills from coarser, fills from
    * same, and then fills physical boundaries.
    */

   recursiveFill(fill_time, do_physical_boundary_fill);

   /*
    * Copy the scratch space of the destination level to the destination
    * space.
    */

   copyScratchToDestination(d_dst_level);

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

template<int DIM> void RefineSchedule<DIM>::recursiveFill(
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

      hier::ComponentSelector allocate_vector;
      allocateScratchSpace(d_coarse_level, fill_time, allocate_vector);

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

      d_coarse_level->deallocatePatchData(allocate_vector);
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
**************************************************************************
*                                                                        *
* Fill the physical boundaries of the specified level with data.         *
*                                                                        *
**************************************************************************
*/

template<int DIM> void RefineSchedule<DIM>::fillPhysicalBoundaries(
   tbox::Pointer< hier::PatchLevel<DIM> > level,
   double fill_time) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!level.isNull());
#endif

   level->setBoundaryBoxes();

   if (d_refine_patch_strategy) {
      for (typename hier::PatchLevel<DIM>::Iterator p(level); p; p++) {
	 tbox::Pointer< hier::Patch<DIM> > patch = level->getPatch(p());
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
*                                                                        *
* Check whether the scratch data needs to be allocated on the specified  *
* level.  Keep track of those allocated components so that they may be   *
* deallocated later.                                                     *
*                                                                        *
**************************************************************************
*/

template<int DIM> void RefineSchedule<DIM>::allocateScratchSpace(
   tbox::Pointer< hier::PatchLevel<DIM> > level,
   double fill_time,
   hier::ComponentSelector& allocate_vector) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!level.isNull());
#endif

   allocate_vector.clrAllFlags();

   hier::ComponentSelector preprocess_vector;

   for (int iri = 0; iri < d_number_refine_items; iri++) {
      const int scratch_id = d_refine_items[iri]->d_scratch;
      if (!level->checkAllocated(scratch_id)) {
	 allocate_vector.setFlag(scratch_id);
      }
      preprocess_vector.setFlag(scratch_id);
   }

   level->allocatePatchData(allocate_vector, fill_time,
      tbox::ArenaManager::getManager()->getScratchAllocator());

   d_transaction_factory->preprocessScratchSpace(level,
                                                 fill_time,
                                                 preprocess_vector);
}

/*
**************************************************************************
*                                                                        *
* Fill the component selector argument with the components needed to     *
* allocate source data.                                                  *
*                                                                        *
**************************************************************************
*/

template<int DIM> void RefineSchedule<DIM>::initializeSourceVector(
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

template<int DIM>
void RefineSchedule<DIM>::initializeDestinationVector(
   hier::ComponentSelector& allocate_vector) const
{
   allocate_vector.clrAllFlags();

   for (int iri = 0; iri < d_number_refine_items; iri++) {
      allocate_vector.setFlag(d_refine_items[iri]->d_dst);
   }
}


/*
**************************************************************************
*                                                                        *
* Allocate all destination data and store the destination components     *
* in a component selector.                                               *
*                                                                        *
**************************************************************************
*/

template<int DIM> void RefineSchedule<DIM>::allocateDestinationSpace(
   double fill_time,
   hier::ComponentSelector& allocate_vector) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_dst_level.isNull());
#endif
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
*                                                                        *
* Copy data from the scratch space of the specified level into the       *
* destination space.  If the two spaces are the same, then no copy       *
* is performed.                                                          *
*                                                                        *
**************************************************************************
*/

template<int DIM> void RefineSchedule<DIM>::copyScratchToDestination(
   tbox::Pointer< hier::PatchLevel<DIM> > level) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!level.isNull());
#endif

   for (typename hier::PatchLevel<DIM>::Iterator p(level); p; p++) {
      tbox::Pointer< hier::Patch<DIM> > patch = level->getPatch(p());

      for (int iri = 0; iri < d_number_refine_items; iri++) {
	 const int src_id = d_refine_items[iri]->d_scratch;
	 const int dst_id = d_refine_items[iri]->d_dst;
	 if (src_id != dst_id) {
#ifdef DEBUG_CHECK_ASSERTIONS
	    const double dst_time = patch->getPatchData(dst_id)->getTime();
	    const double src_time = patch->getPatchData(src_id)->getTime();
	    TBOX_ASSERT(tbox::MathUtilities<double>::equalEps(dst_time, src_time));
#endif
	    patch->getPatchData(dst_id)->copy(*patch->getPatchData(src_id));
	 }
      }

   }

}


/*
**************************************************************************
*                                                                        *
* Refine data from the coarse level into the fine level on the provided  *
* fill box regions.  All operations are performed on the scratch space.  *
*                                                                        *
**************************************************************************
*/

template<int DIM> void RefineSchedule<DIM>::refineScratchData() const
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
      const hier::BoxList<DIM>& fill_boxes = d_fine_fill_boxes[p()].getBoxList();

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

      for (typename hier::BoxList<DIM>::Iterator b(fill_boxes); b; b++) {
	 for (int iri = 0; iri < d_number_refine_items; iri++) {
	    const typename xfer::RefineClasses<DIM>::Data* const ref_item = d_refine_items[iri];
	    if (!(ref_item->d_oprefine.isNull())) {
	       const int scratch_id = ref_item->d_scratch;
	       const hier::Box<DIM>& scratch_space =
		     fine_patch->getPatchData(scratch_id)->getGhostBox();
	       ref_item->d_oprefine->refine(*fine_patch, *crse_patch,
					    scratch_id, scratch_id,
					    b()*scratch_space, ratio);
	    }
	 }
      }

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
*                                                                       *
* Generate communication schedule routine creates transactions to move  *
* data from interiors of the source space on the source level into the  *
* specified fill box regions of the destination level.  Each fill box   *
* will typically be the interior plus max ghost cells over all data     *
* components on the destination region.  If the source and the scratch  *
* are the same and the source and destination levels  are the same,     *
* then there is no need to copy on the interiors for the same patch.    *
*                                                                       *
* The resulting transactions will only fill the regions of intersection *
* between the fill_boxes and destination level boxes.  The remaining    *
* box regions are returned in unfilled_boxes.                           *
*                                                                       *
* This main schedule generation routine which passes control to one     *
* of the algorithmic variations below based on boolean parameters set   *
* to default settings in the constructors and possibly changed via      *
* an input file.                                                        *
*                                                                       *
** ***********************************************************************
*/

template<int DIM> void RefineSchedule<DIM>::generateCommunicationSchedule(
   tbox::Pointer<tbox::Schedule> coarse_priority_schedule,
   tbox::Pointer<tbox::Schedule> fine_priority_schedule,
   tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
   tbox::Pointer< hier::PatchLevel<DIM> > src_level,
   const tbox::Array< xfer::FillBoxSet<DIM> >& fill_boxes,
   tbox::Array< xfer::FillBoxSet<DIM> >& unfilled_boxes,
   bool use_time_interpolation)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!coarse_priority_schedule.isNull());
   TBOX_ASSERT(!fine_priority_schedule.isNull());
   TBOX_ASSERT(!dst_level.isNull());
   TBOX_ASSERT(!src_level.isNull());
   TBOX_ASSERT(fill_boxes.getSize() == dst_level->getNumberOfPatches());
   TBOX_ASSERT(unfilled_boxes.getSize() == dst_level->getNumberOfPatches());
#endif

   if (s_schedule_generation_method == "ORIG_NSQUARED") {

      generateCommunicationScheduleNSquared(coarse_priority_schedule,
					    fine_priority_schedule,
					    dst_level,
					    src_level,
					    fill_boxes,
					    unfilled_boxes,
					    use_time_interpolation);


   } else if (s_schedule_generation_method == "BOX_GRAPH") {

      generateCommunicationScheduleBoxGraph(coarse_priority_schedule,
					    fine_priority_schedule,
					    dst_level,
					    src_level,
					    fill_boxes,
					    unfilled_boxes,
					    use_time_interpolation);

   } else if (s_schedule_generation_method == "BOX_TREE") {

       generateCommunicationScheduleBoxTree(coarse_priority_schedule,
					   fine_priority_schedule,
					   dst_level,
					   src_level,
					   fill_boxes,
					   unfilled_boxes,
					   use_time_interpolation);

   } else {

      TBOX_ERROR("Internal RefineSchedule error..."
		 << "\n unrecognized schedule generation option: "
		 << s_schedule_generation_method << std::endl);

   }

}

/*
*************************************************************************
*                                                                       *
* This version of the schedule generation procedure uses the original   *
* SAMRAI N^2 algorithms to construct communication schedules.  First,   *
* construct the set of "unfilled boxes" for each destination which      *
* indicates the regions that cannot filled from the source level (i.e., *
* they require filling from coarser levels in a patch hierarchy.  Then, *
* loop over all of the patches on the source and destination levels.    *
* check to see whether source or destination is local to this processor.*
* If not, then skip over schedule construction operations.              *
*                                                                       *
*************************************************************************
*/

template<int DIM> void RefineSchedule<DIM>::generateCommunicationScheduleNSquared(
   tbox::Pointer<tbox::Schedule> coarse_priority_schedule,
   tbox::Pointer<tbox::Schedule> fine_priority_schedule,
   tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
   tbox::Pointer< hier::PatchLevel<DIM> > src_level,
   const tbox::Array< xfer::FillBoxSet<DIM> >& fill_boxes,
   tbox::Array< xfer::FillBoxSet<DIM> >& unfilled_boxes,
   bool use_time_interpolation)
{

   t_gen_sched_n_squared->start();

   const int dst_npatches = dst_level->getNumberOfPatches();
   const int src_npatches = src_level->getNumberOfPatches();

   const hier::ProcessorMapping& dst_mapping = dst_level->getProcessorMapping();
   const hier::ProcessorMapping& src_mapping = src_level->getProcessorMapping();

    makeUnfilledBoxesNSquared(unfilled_boxes,
			      fill_boxes,
			      dst_level,
			      src_level);

   for (int dp = 0; dp < dst_npatches; dp++) {

      for (int sp = 0; sp < src_npatches; sp++) {

	 if (   dst_mapping.isMappingLocal(dp)
	     || src_mapping.isMappingLocal(sp) ) {

	    constructScheduleTransactions(fine_priority_schedule,
					  coarse_priority_schedule,
					  fill_boxes[dp].getBoxList(),
					  dst_level, dp,
					  src_level, sp,
					  use_time_interpolation);

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
* each destination patch.  First, construct the set of "unfilled boxes" *
* for each destination which indicates the regions that cannot filled   *
* from the source level (i.e., they require filling from coarser levels *
* in a patch hierarchy.  Then, generated the graph and only perform     *
* schedule construction operations between source and destination       *
* patches where one is local to processor and their overlap is non-     *
* empty.                                                                *
*                                                                       *
*************************************************************************
*/

template<int DIM> void RefineSchedule<DIM>::generateCommunicationScheduleBoxGraph(
   tbox::Pointer<tbox::Schedule> coarse_priority_schedule,
   tbox::Pointer<tbox::Schedule> fine_priority_schedule,
   tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
   tbox::Pointer< hier::PatchLevel<DIM> > src_level,
   const tbox::Array< xfer::FillBoxSet<DIM> >& fill_boxes,
   tbox::Array< xfer::FillBoxSet<DIM> >& unfilled_boxes,
   bool use_time_interpolation)
{

   t_gen_sched_box_graph->start();

   const int dst_npatches = dst_level->getNumberOfPatches();
   const hier::BoxArray<DIM>& src_boxes = src_level->getBoxes();
   const hier::BoxArray<DIM>& dst_boxes = dst_level->getBoxes();
   const hier::ProcessorMapping& dst_mapping = dst_level->getProcessorMapping();

   tbox::Pointer< hier::BoxTop<DIM> > box_top = src_level->getBoxTop();

   for (int dp = 0; dp < dst_npatches; dp++) {
      unfilled_boxes[dp] = fill_boxes[dp];
      box_top->removeIntersections(unfilled_boxes[dp].getBoxListToChange());
   }

   tbox::Pointer< hier::BoxGraph<DIM> > box_graph;

   if (dst_level == src_level) {

      box_graph = dst_level->getBoxGraph();

   } else {

      hier::IntVector<DIM> growth =
	 hier::IntVector<DIM>::max(d_max_scratch_gcw, getMaxDestinationGhosts());
      int max_gcw = tbox::MathUtilities<int>::Max(growth.max(),1);
      hier::IntVector<DIM> dst_growth(max_gcw);

      box_graph = new hier::BoxGraph<DIM>(src_boxes,
				     src_level->getShiftsForLevel(),
				     src_level->getProcessorMapping(),
				     dst_boxes,
				     dst_growth);

   }

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

	 constructScheduleTransactions(fine_priority_schedule,
				       coarse_priority_schedule,
				       fill_boxes[dp].getBoxList(),
				       dst_level, dp,
				       src_level, sp,
				       use_time_interpolation);

      } // loop over source patches

   } // loop over destination patches

   t_gen_sched_box_graph->stop();

}

/*
*************************************************************************
*                                                                       *
* This version of the schedule generation procedure uses a recursive    *
* binary box tree algorithm to determine which source patches           *
* contribute data to each destination patch and to compute              *
* unfilled_boxes.  First, construct the set of "unfilled boxes"         *
* for each destination which indicates the regions that cannot filled   *
* from the source level (i.e., they require filling from coarser levels *
* in a patch hierarchy.  Then, generated the graph and only perform     *
* schedule construction operations between source and destination       *
* patches where one is local to processor and their overlap is non-     *
* empty.                                                                *
*                                                                       *
*************************************************************************
*/

template<int DIM> void RefineSchedule<DIM>::generateCommunicationScheduleBoxTree(
   tbox::Pointer<tbox::Schedule> coarse_priority_schedule,
   tbox::Pointer<tbox::Schedule> fine_priority_schedule,
   tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
   tbox::Pointer< hier::PatchLevel<DIM> > src_level,
   const tbox::Array< xfer::FillBoxSet<DIM> >& fill_boxes,
   tbox::Array< xfer::FillBoxSet<DIM> >& unfilled_boxes,
   bool use_time_interpolation)
{

   t_gen_sched_box_tree->start();

   const int dst_npatches = dst_level->getNumberOfPatches();
   const hier::BoxArray<DIM>& dst_boxes = dst_level->getBoxes();
   const hier::ProcessorMapping& dst_mapping = dst_level->getProcessorMapping();

   tbox::Pointer< hier::BoxTree<DIM> > box_tree = src_level->getBoxTree();

   for (int dp = 0; dp < dst_npatches; dp++) {
      unfilled_boxes[dp] = fill_boxes[dp];
      box_tree->removeIntersections(unfilled_boxes[dp].getBoxListToChange());
   }

   hier::IntVector<DIM> growth =
      hier::IntVector<DIM>::max(d_max_scratch_gcw, getMaxDestinationGhosts());
   int max_gcw = tbox::MathUtilities<int>::Max(growth.max(),1);
   hier::IntVector<DIM> dst_growth(max_gcw);

   for (int dp = 0; dp < dst_npatches; dp++) {

      const hier::Box<DIM>& dst_box = dst_boxes[dp];
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

	 constructScheduleTransactions(fine_priority_schedule,
				       coarse_priority_schedule,
				       fill_boxes[dp].getBoxList(),
				       dst_level, dp,
				       src_level, sp,
				       use_time_interpolation);

      } // loop over source patches

   } // loop over destination patches

   t_gen_sched_box_tree->stop();

}



/*
**************************************************************************
*
* New fill box methods...
*
**************************************************************************
*/

template<int DIM> void RefineSchedule<DIM>::allocateFillBoxes(
   const std::string& fill_pattern,
   tbox::Array< xfer::FillBoxSet<DIM> >& fill_boxes,
   tbox::Pointer< hier::PatchLevel<DIM> > level,
   const hier::IntVector<DIM>& fill_ghost_width)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!level.isNull());
   TBOX_ASSERT(fill_boxes.getSize() == 0);
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

   const hier::BoxArray<DIM>& boxes = level->getBoxes();
   const int nboxes = boxes.getNumberOfBoxes();

   fill_boxes.resizeArray(nboxes);

   if ( fill_pattern == "DEFAULT_FILL" ) {
      /*
       * The default option fills the patch boxes grown by gcw.
       */
      for (int p = 0; p < nboxes; p++) {
         fill_boxes[p].resetFillBoxes(hier::Box<DIM>::grow(boxes[p], gcw));
         d_max_fill_boxes =
            tbox::MathUtilities<int>::Max( d_max_fill_boxes,
                                           fill_boxes[p].getNumberOfBoxes() );
      }
   }

   else if ( fill_pattern == "FILL_INTERIORS_ONLY" ) {
      /*
       * Fill just the interior.  Disregard gcw.
       */
      for (int p = 0; p < nboxes; p++) {
         fill_boxes[p].resetFillBoxes(boxes[p]);
         d_max_fill_boxes =
            tbox::MathUtilities<int>::Max( d_max_fill_boxes,
                                           fill_boxes[p].getNumberOfBoxes() );
      }
   }

   else if ( fill_pattern == "FILL_LEVEL_BORDERS_ONLY" ) {
      /*
       * To get the level border, grow each patch box and remove
       * the level from it.
       */
      const hier::BoxTree<DIM> &boxtree = *level->getBoxTree();
      for (int p = 0; p < nboxes; p++) {
         hier::BoxList<DIM> levelborder( hier::Box<DIM>::grow(boxes[p], gcw) );
         boxtree.removeIntersections( levelborder );
         if ( ! levelborder.isEmpty() ) {
            fill_boxes[p].resetFillBoxes( levelborder );
            d_max_fill_boxes =
               tbox::MathUtilities<int>::Max( d_max_fill_boxes,
                                              fill_boxes[p].getNumberOfBoxes() );
         }
      }
   }

   else if ( fill_pattern == "FILL_LEVEL_BORDERS_AND_INTERIORS" ) {
      /*
       * Grow each patch box and remove the level from it, except the
       * patch box itself.
       */
      const hier::BoxTree<DIM> &boxtree = *level->getBoxTree();
      for (int p = 0; p < nboxes; p++) {
         hier::Box<DIM> ghostbox = hier::Box<DIM>::grow(boxes[p], gcw);
         hier::BoxList<DIM> nofill;
         boxtree.findOverlapBoxes( nofill, ghostbox );
         for ( typename hier::BoxList<DIM>::Iterator li(nofill); li; li++ ) {
            if ( *li == boxes[p] ) {
               // Exclude the patch box itself.
               nofill.removeItem(li);
               break;
            }
         }
         hier::BoxList<DIM> tofill( ghostbox );
         tofill.removeIntersections( nofill );
         fill_boxes[p].resetFillBoxes( tofill );
         d_max_fill_boxes =
            tbox::MathUtilities<int>::Max( d_max_fill_boxes,
                                           fill_boxes[p].getNumberOfBoxes() );
      }
   }

   else {
      TBOX_ERROR("RefineSchedule<DIM>::allocateFillBoxes\n"
                 << "Given communication pattern string "
                 << fill_pattern << " is invalid.\n Valid options are\n"
                 << "'DEFAULT_FILL', 'FILL_LEVEL_BORDERS_ONLY',\n"
                 << "'FILL_INTERIORS_ONLY', 'FILL_LEVEL_BORDERS_AND_INTERIORS'\n"
                 << std::endl);
   }

   return;
}


/*
**************************************************************************
*                                                                        *
* Calculate the maximum ghost cell width of all destination patch data   *
* components.                                                            *
*                                                                        *
**************************************************************************
*/

template<int DIM> hier::IntVector<DIM> RefineSchedule<DIM>::getMaxDestinationGhosts() const
{
   hier::IntVector<DIM> gcw(0);
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

template<int DIM> hier::IntVector<DIM> RefineSchedule<DIM>::getMaxScratchGhosts() const
{
   hier::IntVector<DIM> gcw(0);
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

template<int DIM> hier::IntVector<DIM> RefineSchedule<DIM>::getMaxStencilGhosts() const
{
   hier::IntVector<DIM> gcw(0);
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
*
* Private utility function to construct the unfilled_boxes list
* associated with each destination patch box which specifies those
* regions that cannot be filled from the source level (i.e they
* require filling from coarser levels in a patch hierarchy).
*
* NOTE: This is the original N^2 implementation in SAMRAI.
*
*************************************************************************
*/

template<int DIM> void RefineSchedule<DIM>::makeUnfilledBoxesNSquared(
   tbox::Array< xfer::FillBoxSet<DIM> >& unfilled_boxes,
   const tbox::Array< xfer::FillBoxSet<DIM> >& fill_boxes,
   tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
   tbox::Pointer< hier::PatchLevel<DIM> > src_level) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst_level.isNull());
   TBOX_ASSERT(!src_level.isNull());
   TBOX_ASSERT(fill_boxes.getSize() == dst_level->getNumberOfPatches());
   TBOX_ASSERT(unfilled_boxes.getSize() == dst_level->getNumberOfPatches());
#endif

   const int dst_npatches = dst_level->getNumberOfPatches();
   const int src_npatches = src_level->getNumberOfPatches();

   const hier::BoxArray<DIM>& src_boxes = src_level->getBoxes();

   for (int dp = 0; dp < dst_npatches; dp++) {

      unfilled_boxes[dp] = fill_boxes[dp];

      for (int sp = 0; sp < src_npatches; sp++) {

	 const hier::Box<DIM>& src_box = src_boxes[sp];

	 typename tbox::List< hier::IntVector<DIM> >::Iterator
	    sh(src_level->getShiftsForPatch(sp));

	 bool zero_shift = true;

	 while (sh || zero_shift) {

	    hier::IntVector<DIM> shift(0);
	    if (!zero_shift) {
	       shift = sh();
	    }

	    const hier::Box<DIM> shifted(hier::Box<DIM>::shift(src_box, shift));

	    unfilled_boxes[dp].removeIntersections(shifted);
	    if (!zero_shift) {
	       sh++;
	    } else {
	       zero_shift = false;
	    }

	 }  // loop over valid shifts of source box

      }  // loop over source boxes

   }

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

template<int DIM> void RefineSchedule<DIM>::constructScheduleTransactions(
   tbox::Pointer<tbox::Schedule> fine_priority_schedule,
   tbox::Pointer<tbox::Schedule> coarse_priority_schedule,
   const hier::BoxList<DIM>& fill_boxes,
   tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
   int dst_patch_id,
   tbox::Pointer< hier::PatchLevel<DIM> > src_level,
   int src_patch_id,
   bool use_time_interpolation)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!coarse_priority_schedule.isNull());
   TBOX_ASSERT(!fine_priority_schedule.isNull());
   TBOX_ASSERT(!dst_level.isNull());
   TBOX_ASSERT(!src_level.isNull());
#endif

   if (d_src_masks.getNumberOfBoxes() < d_max_fill_boxes) {
      d_src_masks.resizeBoxArray(0);
      d_src_masks.resizeBoxArray(d_max_fill_boxes);
   }
   if (d_overlaps.getSize() < d_max_fill_boxes) {
      d_overlaps.setNull();
      d_overlaps.resizeArray(d_max_fill_boxes);
   }

   tbox::Pointer< hier::PatchDescriptor<DIM> > dst_patch_descriptor =
      dst_level->getPatchDescriptor();
   tbox::Pointer< hier::PatchDescriptor<DIM> > src_patch_descriptor =
      src_level->getPatchDescriptor();

   const hier::Box<DIM>& dst_box = dst_level->getBoxes()[dst_patch_id];
   const hier::Box<DIM>& src_box = src_level->getBoxes()[src_patch_id];

   const bool same_patch = ( (dst_level == src_level) &&
			     (dst_patch_id == src_patch_id) );

   const int num_fill_boxes = fill_boxes.getNumberOfBoxes();
   const int num_equiv_classes =
      d_refine_classes->getNumberOfEquivalenceClasses();

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

      for (int nc = 0; nc < num_equiv_classes; nc++) {

	 const typename xfer::RefineClasses<DIM>::Data& rep_item =
	    d_refine_classes->getClassRepresentative(nc);

	 const int rep_item_dst_id = rep_item.d_scratch;
	 const int rep_item_src_id = rep_item.d_src;

	 tbox::Pointer< hier::PatchDataFactory<DIM> > src_pdf =
	    src_patch_descriptor->getPatchDataFactory(rep_item_src_id);
	 tbox::Pointer< hier::PatchDataFactory<DIM> > dst_pdf =
	    dst_patch_descriptor->getPatchDataFactory(rep_item_dst_id);

	 const hier::IntVector<DIM>& dst_gcw = dst_pdf->getGhostCellWidth();

         /*
          * Iterate over boxes in fill box list for destination patch.
          * Then, calculate the src_mask and overlap for each fill box.
          * For each equivalence class, this loop is executed once.
          */

	 int box_num = 0;
	 for (typename hier::BoxList<DIM>::Iterator b(fill_boxes); b; b++) {

	    const hier::Box<DIM>& fill_box = b();

            /*
             * Get the patch data factories and calculate the
             * overlap.  Note that we restrict the domain of the
             * time interpolation to the intersection of the
             * fill box and the ghost box of the destination
             * patch data component.  This is needed for the
             * case where the schedule treats data components
             * with different ghost cell widths since the fill
             * boxes are generated using the largest ghost width.
             */

            hier::Box<DIM> dst_fill_box(hier::Box<DIM>::grow(dst_box, dst_gcw));
            dst_fill_box = dst_fill_box * fill_box;

            hier::Box<DIM> test_mask(dst_fill_box*shifted);
            if ( test_mask.empty() &&
                 (dst_gcw == s_constant_zero_intvector) &&
                 dst_pdf->dataLivesOnPatchBorder() ) {
               hier::Box<DIM> tmp_dst_fill_box(hier::Box<DIM>::grow(dst_fill_box,
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
	       TBOX_ERROR("Internal RefineSchedule<DIM> error..."
			  << "\n Overlap is NULL for "
			  << "\n src box = " << src_box
			  << "\n dst box = " << dst_box
			  << "\n src mask = " << src_mask << std::endl);
	    }
#endif

	    d_src_masks[box_num] = src_mask;
	    d_overlaps[box_num] = overlap;
	    box_num++;

	 }

         /*
          * Iterate over components in refine description list
          */
	 bool same_patch_no_shift = (same_patch && zero_shift);
	 for (typename tbox::List<typename xfer::RefineClasses<DIM>::Data>::Iterator
                 l(d_refine_classes->getIterator(nc)); l; l++) {
	    const int dst_id = l().d_scratch;
	    const int src_id = l().d_src;
            const int ritem_count = l().d_tag;

            /*
             * If the src and dst patches, levels, and components are the
             * same, and there is no shift, the data exchange is unnecessary.
             */
	    if ( !same_patch_no_shift || (dst_id != src_id) ) {

               /*
                * Iterate over the fill boxes and create transactions
                * for each box that has a non-empty overlap.
                */
	       for (int i = 0; i < num_fill_boxes; i++) {

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

		  if (!d_overlaps[i]->isOverlapEmpty()) {

		     tbox::Pointer<tbox::Transaction> transaction;

                     bool do_time_interpolation = (use_time_interpolation &&
                                                   l().d_time_interpolate);
                     transaction =
                        d_transaction_factory->allocate(dst_level,
                                                        src_level,
                                                        d_overlaps[i],
                                                        dst_patch_id,
                                                        src_patch_id,
                                                        ritem_count,
                                                        d_src_masks[i],
                                                        do_time_interpolation);

		     if (l().d_fine_bdry_reps_var) {
			if (same_patch) {
			   fine_priority_schedule->addTransaction(transaction);
			} else {
			   fine_priority_schedule->appendTransaction(transaction);
			}
		     } else {
			if (same_patch) {
			   coarse_priority_schedule->addTransaction(transaction);
			} else {
			   coarse_priority_schedule->appendTransaction(transaction);
			}
		     }

		  }  // if overlap not empty

	       }  // iterate over fill_boxes

	    } // if copy transaction is needed
              // (i.e., not same patch data object)

	 }  // iterate over refine components in equivalence class

      }  // iterate over refine equivalence classes

      if (!zero_shift) {
	 sh++;
      } else {
	 zero_shift = false;
      }

   }  // while (sh || zero_shift)

}

/*
*************************************************************************
*
* Private member function to initialize data members for hierarchy info.
*
*************************************************************************
*/

template<int DIM> void
RefineSchedule<DIM>::initializeDomainAndGhostInformation(
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
      d_domain_box = domain[0];
      for (int i = 1; i < domain.getNumberOfBoxes(); i++) {
	 d_domain_box += domain[i];
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

template<int DIM> void RefineSchedule<DIM>::setRefineItems(
   const tbox::Pointer< xfer::RefineClasses<DIM> > refine_classes)
{

   clearRefineItems();

   d_refine_classes            = refine_classes;

   const int num_refine_classes =
      d_refine_classes->getNumberOfEquivalenceClasses();

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

template<int DIM> void RefineSchedule<DIM>::initialCheckRefineClassItems() const
{
   tbox::Pointer< hier::PatchDescriptor<DIM> > pd = d_dst_level->getPatchDescriptor();

   hier::IntVector<DIM> user_gcw(0);
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
	       TBOX_ERROR("Bad data given to RefineSchedule<DIM>...\n"
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

template<int DIM> void RefineSchedule<DIM>::clearRefineItems()
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
**************************************************************************
*
* Print class data to the specified output stream.
*
**************************************************************************
*/

template<int DIM> void RefineSchedule<DIM>::printClassData(std::ostream& stream) const
{
   stream << "RefineSchedule<DIM>::printClassData()\n";
   stream << "--------------------------------------\n";
   stream << "s_schedule_generation_method = "
	  << s_schedule_generation_method << std::endl;

   d_refine_classes->printClassData(stream);

   stream << "d_dst_level = " << (hier::PatchLevel<DIM>*)d_dst_level << std::endl;

   stream << "d_refine_patch_strategy = "
          << (xfer::RefinePatchStrategy<DIM>*)d_refine_patch_strategy << std::endl;

   stream << "d_max_stencil_gcw = " << d_max_stencil_gcw << std::endl;
   stream << "d_max_scratch_gcw = " << d_max_scratch_gcw << std::endl;
   stream << "d_boundary_fill_ghost_width = " << d_boundary_fill_ghost_width << std::endl;

   stream << "d_force_boundary_fill = " << d_force_boundary_fill << std::endl;
   stream << "d_domain_is_one_box = " << d_domain_is_one_box << std::endl;
   stream << "d_domain_box = " << d_domain_box << std::endl;
   stream << "d_num_periodic_directions = " << d_num_periodic_directions << std::endl;
   stream << "d_periodic_shift = " << d_periodic_shift << std::endl;

   stream << "d_coarse_level = " << (hier::PatchLevel<DIM>*)d_coarse_level << std::endl;

   stream << "d_coarse_to_fine_mapping size = "
          << d_coarse_to_fine_mapping.size() << std::endl;

   stream << "d_fine_fill_boxes size = "
          << d_fine_fill_boxes.size() << std::endl;

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
