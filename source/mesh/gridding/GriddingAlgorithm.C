//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-0/source/mesh/gridding/GriddingAlgorithm.C $
// Package:     SAMRAI mesh
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2297 $
// Modified:    $LastChangedDate: 2008-07-14 17:06:37 -0700 (Mon, 14 Jul 2008) $
// Description: AMR hierarchy generation and regridding routines.
//

#ifndef included_mesh_GriddingAlgorithm_C
#define included_mesh_GriddingAlgorithm_C

#include "GriddingAlgorithm.h"

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include "Box.h"
#include "BoxArray.h"
#include "BoxList.h"
#include "BoxUtilities.h"
#include "GridGeometry.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchDataFactory.h"
#include "PatchDescriptor.h"
#include "PatchCellDataBasicOps.h"
#include "VariableDatabase.h"
#include "CellData.h"
#include "CellIterator.h"
#include "CellVariable.h"
#include "tbox/List.h"
#include "tbox/PIO.h"
#include "tbox/RestartManager.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"
#include "RefineAlgorithm.h"
#include "RefineSchedule.h"

#define ALGS_GRIDDING_ALGORITHM_VERSION (2)

#ifndef NULL
#define NULL (0)
#endif

#ifdef DEBUG_NO_INLINE
#include "GriddingAlgorithm.I"
#endif
namespace SAMRAI {
    namespace mesh {

/*
*************************************************************************
*                                                                       *
* Initialize the static data members.                                   *
*                                                                       *
*************************************************************************
*/
template<int DIM> int GriddingAlgorithm<DIM>::s_instance_counter = 0;
template<int DIM> int GriddingAlgorithm<DIM>::s_tag_indx = -1;
template<int DIM> int GriddingAlgorithm<DIM>::s_buf_tag_indx = -1;

/*
*************************************************************************
*                                                                       *
* Constructor and destructor for GriddingAlgorithm<DIM>.                *
*                                                                       *
*************************************************************************
*/
template<int DIM> GriddingAlgorithm<DIM>::GriddingAlgorithm(
   const std::string& object_name,
   tbox::Pointer<tbox::Database> input_db,
   tbox::Pointer< TagAndInitializeStrategy<DIM> > tag_init_strategy,
   tbox::Pointer< BoxGeneratorStrategy<DIM> > generator,
   tbox::Pointer< LoadBalanceStrategy<DIM> > balancer,
   bool register_for_restart)
   :
   d_check_nonrefined_tags('w'),
   d_check_overlapping_patches('i')
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(!input_db.isNull());
    TBOX_ASSERT(!tag_init_strategy.isNull());
    TBOX_ASSERT(!generator.isNull());
    TBOX_ASSERT(!balancer.isNull());
#endif

   d_object_name     = object_name;
   d_registered_for_restart = register_for_restart;

   if (d_registered_for_restart) {
      tbox::RestartManager::getManager()->
         registerRestartItem(d_object_name, this);
   }

   d_tag_init_strategy  = tag_init_strategy;
   d_box_generator   = generator;
   d_load_balancer   = balancer;

   /*
    * Construct integer tag variables and add to variable database.  Note that
    * variables and patch data indices are shared among all GriddingAlgorithm
    * instances.  The VariableDatabase holds the variables, once contructed and
    * registered via the VariableDatabase::registerInternalSAMRAIVariable()
    * function call.  Note that variables are registered and patch data indices
    * are made only for the first time through the constructor.
    */
   hier::VariableDatabase<DIM>* var_db = hier::VariableDatabase<DIM>::getDatabase();

   static std::string tag_interior_variable_name("GriddingAlgorithm__tag-interior");
   static std::string tag_buffer_variable_name("GriddingAlgorithm__tag-buffer");

   d_tag = var_db->getVariable(tag_interior_variable_name);
   if (d_tag.isNull()) {
      d_tag = new pdat::CellVariable<DIM,int>(tag_interior_variable_name, 1);
   }

   d_buf_tag = var_db->getVariable(tag_buffer_variable_name);
   if (d_buf_tag.isNull()) {
      d_buf_tag = new pdat::CellVariable<DIM,int>(tag_buffer_variable_name, 1);
   }

   if (s_tag_indx < 0) {
      s_tag_indx =
         var_db->registerInternalSAMRAIVariable(d_tag,
                                                hier::IntVector<DIM>(0));
   }
   if (s_buf_tag_indx < 0) {
      s_buf_tag_indx =
         var_db->registerInternalSAMRAIVariable(d_buf_tag,
                                                hier::IntVector<DIM>(1));
   }

   d_tag_indx = s_tag_indx;
   d_buf_tag_indx = s_buf_tag_indx;

   /*
    * Tag value for refined cells is one; others are zero.
    */
   d_true_tag = 1;
   d_false_tag = 0;

   /*
    * Initialize communication algorithm for exchanging buffered tag data.
    */
   d_bdry_fill_tags = new xfer::RefineAlgorithm<DIM>();

   d_bdry_fill_tags->
     registerRefine(d_buf_tag_indx,
                    d_buf_tag_indx,
                    d_buf_tag_indx,
                    ((tbox::Pointer<xfer::RefineOperator<DIM> >)NULL));

   /*
    * Set default state for gridding parameters.  These values will be
    * reset to valid values when reading from input or restart.
    */
   d_max_levels = 1;

   d_ratio_to_coarser.resizeArray(d_max_levels);

   d_efficiency_tolerance.resizeArray(d_max_levels);
   d_combine_efficiency.resizeArray(d_max_levels);

   d_smallest_patch_size.resizeArray(d_max_levels);
   d_largest_patch_size.resizeArray(d_max_levels);

   d_proper_nesting_buffer.resizeArray(d_max_levels);
   d_allow_patches_smaller_than_ghostwidth = false;
   d_allow_patches_smaller_than_minimum_size_to_prevent_overlaps = false;

   for (int ln = 0; ln < d_max_levels; ln++) {
     d_ratio_to_coarser[ln] = hier::IntVector<DIM>(1);
     d_efficiency_tolerance[ln] = 0.8e0;
     d_combine_efficiency[ln] = 0.8e0;

     d_smallest_patch_size[ln] = hier::IntVector<DIM>(0);
     d_largest_patch_size[ln] = hier::IntVector<DIM>(0);
     d_proper_nesting_buffer[ln] = 1;
   }

   d_write_dumped_level_boxes = false;
   d_read_dumped_level_boxes = false;

   d_barrier_before_clustering = false;
   d_sort_boxes_after_clustering = false;
   d_coalesce_boxes = true;

   d_extend_tags_to_bdry = false;

   /*
    * Timers:  for gathering performance information about box
    * calculus and other regridding operations.
    */
   t_find_domain_complement = tbox::TimerManager::getManager() ->
     getTimer("mesh::GriddingAlgorithm::findDomainComplement()");
   t_intersect_boxes_find_refinement = tbox::TimerManager::getManager() ->
     getTimer("mesh::GriddingAlgorithm::intersect_boxes_find_refinement");
   t_load_balance = tbox::TimerManager::getManager() ->
     getTimer("mesh::GriddingAlgorithm::load_balance");
   t_bdry_fill_tags_create = tbox::TimerManager::getManager() ->
     getTimer("mesh::GriddingAlgorithm::bdry_fill_tags_create");
   t_make_coarsest = tbox::TimerManager::getManager() ->
     getTimer("mesh::GriddingAlgorithm::makeCoarsestLevel()");
   t_make_finer = tbox::TimerManager::getManager() ->
     getTimer("mesh::GriddingAlgorithm::makeFinerLevel()");
   t_remove_intersections_make_finer = tbox::TimerManager::getManager() ->
     getTimer("mesh::GriddingAlgorithm::remove_intersections_make_finer");
   t_regrid_all_finer = tbox::TimerManager::getManager() ->
     getTimer("mesh::GriddingAlgorithm::regridAllFinerLevels()");
   t_remove_intersections_regrid_all = tbox::TimerManager::getManager()->
     getTimer("mesh::GriddingAlgorithm::remove_intersections_regrid_all");
   t_remove_intersections_find_proper =  tbox::TimerManager::getManager()->
     getTimer("mesh::GriddingAlgorithm::remove_intersections_find_proper");
   t_intersect_boxes_find_proper = tbox::TimerManager::getManager() ->
     getTimer("mesh::GriddingAlgorithm::intersect_boxes_find_proper");
   t_set_tags = tbox::TimerManager::getManager() ->
     getTimer("mesh::GriddingAlgorithm::setTagsOnLevel()");
   t_buffer_tags = tbox::TimerManager::getManager()->
     getTimer("mesh::GriddingAlgorithm::bufferTagsOnLevel()");
   t_bdry_fill_tags_comm = tbox::TimerManager::getManager() ->
     getTimer("mesh::GriddingAlgorithm::bdry_fill_tags_comm");
   t_find_refinement = tbox::TimerManager::getManager() ->
     getTimer("mesh::GriddingAlgorithm::findRefinementBoxes()");
   t_find_boxes_containing_tags = tbox::TimerManager::getManager()->
     getTimer("mesh::GriddingAlgorithm::find_boxes_containing_tags");
   t_find_nesting_restriction = tbox::TimerManager::getManager()->
     getTimer("mesh::GriddingAlgorithm::find_nesting_restriction");
   t_apply_nesting_restriction = tbox::TimerManager::getManager()->
     getTimer("mesh::GriddingAlgorithm::apply_nesting_restriction");
   t_coalesce_boxes = tbox::TimerManager::getManager()->
     getTimer("mesh::GriddingAlgorithm::coalesce_boxes");
   t_before_load_balance = tbox::TimerManager::getManager()->
     getTimer("mesh::GriddingAlgorithm::before_load_balance");
   t_grow_boxes_within_domain = tbox::TimerManager::getManager()->
     getTimer("mesh::GriddingAlgorithm::grow_boxes_within_domain");
   t_reset_hier = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::reset_hierarchy_config");
   t_tag_cells_for_refinement = tbox::TimerManager::getManager() ->
     getTimer("mesh::GriddingAlgorithm::tag_cells_for_refinement");
   t_box_massage = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::box_massage");
   t_enforce_nesting = tbox::TimerManager::getManager()->
     getTimer("mesh::GriddingAlgorithm::enforce_nesting");
   t_extend_to_domain_boundary = tbox::TimerManager::getManager()->
     getTimer("mesh::GriddingAlgorithm::extend_to_domain_boundary");
   t_regrid_finer_create = tbox::TimerManager::getManager() ->
     getTimer("mesh::GriddingAlgorithm::regridFinerLevel()_create");

  /*
   * Initialize object with data read from input and restart databases.
   */
  bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();
  if ( is_from_restart ) {
     getFromRestart();
  }

  getFromInput(input_db, is_from_restart);

  d_proper_nesting_boxes.resizeArray(d_max_levels);
  d_proper_nesting_complement.resizeArray(d_max_levels);
  d_bdry_sched_tags.resizeArray(d_max_levels);

  if (d_read_dumped_level_boxes) {
     d_regrid_box_utility =
        new hier::BoxIOUtility<DIM>(d_regrid_boxes_filename,
                               hier::BoxIOUtility<DIM>::READ);
  }

  if (d_write_dumped_level_boxes) {
     d_regrid_box_utility =
        new hier::BoxIOUtility<DIM>(d_regrid_boxes_filename,
                               hier::BoxIOUtility<DIM>::WRITE);
  }

  if (d_read_dumped_level_boxes || d_write_dumped_level_boxes) {
     d_regrid_box_counter.resizeArray(d_max_levels);
     for (int il = 0; il < d_regrid_box_counter.getSize(); il++) {
        d_regrid_box_counter[il] = 0;
     }
  }

  s_instance_counter++;

}

/*
*************************************************************************
*                                                                       *
* Destructor tells the tbox::RestartManager to remove this object from  *
* the list of restart items.                                            *
*                                                                       *
*************************************************************************
*/
template<int DIM> GriddingAlgorithm<DIM>::~GriddingAlgorithm()
{

   s_instance_counter--;

   /*
    * If we are writing boxes, have the box utility generate the file.
    */
   if (d_write_dumped_level_boxes) {
      d_regrid_box_utility->writeLevelBoxesDatabase();
   }

   if (d_registered_for_restart) {
      tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
   }

   if (s_instance_counter == 0) {
      hier::VariableDatabase<DIM>::getDatabase()->
         removeInternalSAMRAIVariablePatchDataIndex(s_tag_indx);
      s_tag_indx = -1;
      hier::VariableDatabase<DIM>::getDatabase()->
         removeInternalSAMRAIVariablePatchDataIndex(s_buf_tag_indx);
      s_buf_tag_indx = -1;
   }

}

/*
*************************************************************************
*                                                                       *
* Construct coarsest level in the patch hierarchy (i.e., level 0).      *
* This routine operates in two modes:                                   *
*                                                                       *
* (1) If level 0 does not exist in the hierarchy, then a new level      *
*     will be created based on the domain description provided by       *
*     the grid geometry associated with the hierarchy.                  *
* (2) If level 0 exists in the hierarchy, then a new level is made      *
*     by re-applying the load balancing routines to the existing level. *
*     The pre-existing level will be removed from the hierarchy.        *
*                                                                       *
* In either case, the new level is placed in the hierarchy as level 0.  *
* If level 0 can be refined (i.e., it will be subject to tagging cells  *
* for refinement), the domain boxes are checked against the constraints *
* of regridding process.                                                *
*                                                                       *
*************************************************************************
*/

template<int DIM> void GriddingAlgorithm<DIM>::makeCoarsestLevel(
   tbox::Pointer< hier::BasePatchHierarchy<DIM> > hierarchy,
   const double level_time,
   const hier::BoxArray<DIM>& override_boxes,
   const hier::ProcessorMapping& override_mapping)
{
   tbox::Pointer< hier::PatchHierarchy<DIM> > patch_hierarchy = hierarchy;

   bool override_load_balance = false;
   if (override_boxes.size() != 0) {
      override_load_balance = true;
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!(patch_hierarchy.isNull()));
   TBOX_ASSERT(!(patch_hierarchy->getGridGeometry().isNull()));
   TBOX_ASSERT(d_max_levels > 0);
#endif
   /*
    * Start timer for this method.
    */
   t_make_coarsest->start();

   const int level_number = 0;

   bool level_zero_exists = false;
   if ( (patch_hierarchy->getNumberOfLevels() > 0) ) {
      if ( !(patch_hierarchy->getPatchLevel(level_number).isNull()) ) {
         level_zero_exists = true;
      }
   }


   hier::BoxArray<DIM> domain = patch_hierarchy->getGridGeometry()
                                          ->getPhysicalDomain();

#ifdef DEBUG_CHECK_ASSERTIONS
   if (override_load_balance) {
      TBOX_ASSERT(override_boxes.size() == override_mapping.getProcessorMapping().size());

      hier::BoxList<DIM> override_list(override_boxes);
      override_list.removeIntersections(domain);
      TBOX_ASSERT(override_list.size() == 0);
      hier::BoxList<DIM> test_domain_list(domain);
      test_domain_list.removeIntersections(override_boxes);
      TBOX_ASSERT(test_domain_list.size() == 0);
   }
#endif

   /*
    * Read/write coarse level boxes from/to file.
    */
   if (d_read_dumped_level_boxes) {
      d_regrid_box_utility->getLevelBoxes(domain,
                                          level_number,
                                          d_regrid_box_counter[level_number]);
      d_regrid_box_counter[level_number]++;
   } else {
      if (d_write_dumped_level_boxes) {
         d_regrid_box_utility->putLevelBoxes(domain,
                                             level_number,
                                             d_regrid_box_counter[level_number]);
         d_regrid_box_counter[level_number]++;
      }
   }

   hier::BoxList<DIM> domain_list(domain);
   hier::IntVector<DIM> smallest_patch;
   hier::IntVector<DIM> smallest_box_to_refine;
   hier::IntVector<DIM> largest_patch;
   hier::IntVector<DIM> extend_ghosts;
   // "true" argument: for_building_finer level = true
   getGriddingParameters(smallest_patch,
                         smallest_box_to_refine,
                         largest_patch,
                         extend_ghosts,
                         patch_hierarchy,
                         level_number,
                         true);

   /*
    * If there is no level 0 in the patch hierarchy, then check
    * constraints on domain description.
    */

   if ( !level_zero_exists ) {

      for (typename hier::BoxList<DIM>::Iterator b(domain_list); b; b++) {
         hier::Box<DIM> test_box = b();
         for (int dir = 0; dir < DIM; dir++) {
            if (test_box.numberCells(dir) < smallest_patch(dir)) {
              int error_coarsen_ratio =
                  d_tag_init_strategy->getErrorCoarsenRatio();
               if (error_coarsen_ratio > 1) {
                  TBOX_ERROR(d_object_name << ": "
                             << "\nA box from the input file violates"
                             << "the minimum patch size constraints."
                             << "\nVerify that boxes are larger than"
                             << "the maximum ghost width and/or"
                             << "\nthe specified minimum patch size."
                             << "\nNOTE: to assure the constraints are"
                             << "properly enforced during coarsening for"
                             << "\nerror computation, the minimum patch"
                             << "size is the smallest patch size multiplied"
                             << "\nby the error coarsen ratio, which is "
                             << error_coarsen_ratio << " in this case."
                             << std::endl);
               } else {
                  TBOX_ERROR(d_object_name << ": "
                             << "\nA box from the input file violates"
                             << "the minimum patch size constraints."
                             << "\nVerify that boxes are larger than"
                             << "the maximum ghost width and/or"
                             << "\nthe specified minimum patch size."
                             << std::endl);
               }
            }
         }
      }

      if ( domain_list.boxesIntersect() ) {
         TBOX_ERROR(d_object_name << ":  "
                 << "Boxes specified for coarsest level "
                 << "contain intersections with each other!");
      }

      if ( (d_max_levels > 1)
           && (!d_tag_init_strategy->coarsestLevelBoxesOK(domain)) ) {
         TBOX_ERROR(d_object_name << ":  "
                 << "level gridding strategy encountered"
                 << " a problem with the domain boxes!");
      }

   }

   /*
    * Apply load balancing algorithm to boxes describing domain
    * to determine configuration of patches on level 0.
    */

   hier::IntVector<DIM> patch_cut_factor(d_tag_init_strategy->
                                    getErrorCoarsenRatio());

   const hier::IntVector<DIM> ratio_to_level_zero(1);
   hier::BoxArray<DIM> level_boxes;
   hier::ProcessorMapping mapping;

   if (!override_load_balance) {
      tbox::SAMRAI_MPI::barrier();
      t_load_balance->start();
      d_load_balancer->loadBalanceBoxes(level_boxes, mapping,
                                        domain_list,
                                        patch_hierarchy, level_number, domain,
                                        ratio_to_level_zero,
                                        smallest_patch,
                                        largest_patch,
                                        patch_cut_factor, extend_ghosts);
      t_load_balance->stop();
   } else {
      level_boxes = override_boxes;
      mapping.setProcessorMapping(override_mapping.getProcessorMapping());
   }

   if ( !level_zero_exists ) {

      patch_hierarchy->makeNewPatchLevel(level_number,
                                         ratio_to_level_zero,
                                         level_boxes,
                                         mapping);

      // "true" argument: const bool initial_time = true;
      d_tag_init_strategy->initializeLevelData(patch_hierarchy, level_number,
                                               level_time,
                                               levelCanBeRefined(level_number),
                                               true);

   } else {

      tbox::Pointer< hier::PatchLevel<DIM> > old_level =
         patch_hierarchy->getPatchLevel(level_number);

      patch_hierarchy->removePatchLevel(level_number);

      patch_hierarchy->makeNewPatchLevel(level_number,
                                         ratio_to_level_zero,
                                         level_boxes,
                                         mapping);

      // "false" argument: const bool initial_time = false;
      d_tag_init_strategy->initializeLevelData(patch_hierarchy, level_number,
                                               level_time,
                                               levelCanBeRefined(level_number),
                                               false,
                                               old_level);

      old_level.setNull();

   }

   tbox::SAMRAI_MPI::barrier();
   t_reset_hier->start();
   d_tag_init_strategy->resetHierarchyConfiguration(patch_hierarchy,
                                                    level_number,
                                                    level_number);
   tbox::SAMRAI_MPI::barrier();
   t_reset_hier->stop();

   t_make_coarsest->stop();

}

/*
*************************************************************************
*                                                                       *
* Perform error estimation process on the finest hierarchy level to     *
* determine if and where a new finest level is needed.  If it is        *
* determined  that a new finest level should exist, it is created and   *
* its patch data is allocated and initialized.  The algorithm is        *
* summarized as follows:                                                *
*                                                                       *
* (1) Compute boxes whose union constitutes the region within the level *
*     in which the next finer level may reside (i.e., proper nesting).  *
*                                                                       *
* (2) Set tags on the level to indicate which cells should be refined.  *
*                                                                       *
* (3) Buffer the tags.  This prevents disturbances from moving off      *
*     refined regions before the next remesh occurs.                    *
*                                                                       *
* (4) Determine boxes for new patches that will cover the tagged cells. *
*     This step includes load balancing of the patches.                 *
*                                                                       *
* (5) If there exist some regions to refine, construct the next finer   *
*     level and insert it in the hierarchy.  Then, initialize the data  *
*     on that level.                                                    *
*                                                                       *
*************************************************************************
*/

template<int DIM> void GriddingAlgorithm<DIM>::makeFinerLevel(
   tbox::Pointer< hier::BasePatchHierarchy<DIM> > hierarchy,
   const double level_time,
   const bool initial_time,
   const int tag_buffer,
   const double regrid_start_time)
{
   tbox::Pointer< hier::PatchHierarchy<DIM> > patch_hierarchy = hierarchy;
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!(patch_hierarchy.isNull()));
   TBOX_ASSERT(!(patch_hierarchy->getPatchLevel(
                       patch_hierarchy->getFinestLevelNumber()).isNull()));
   TBOX_ASSERT(tag_buffer >= 0);
#endif

   t_make_finer->start();

   const int level_number = patch_hierarchy->getFinestLevelNumber();
   const int fine_level_number = level_number + 1;

   if ( levelCanBeRefined(level_number) ) {

      const tbox::Pointer< hier::PatchLevel<DIM> >
         patch_level = patch_hierarchy->getPatchLevel(level_number);

      hier::BoxArray<DIM> fine_boxes;
      hier::ProcessorMapping mapping;

      /*
       * The boolean "do_tagging" specifies whether or not tagging will
       * be performed. This will be true except in two circumstances:
       *    1) only user supplied refine boxes are used
       *    2) the boxes are read from a previously dumped file.
       *
       * If either of these circumstances is true, tagging operations
       * are NOT necessary so do_tagging will be set to false.
       */
      bool do_tagging = true;
      if (d_read_dumped_level_boxes ||
          d_tag_init_strategy->refineUserBoxInputOnly()) do_tagging = false;

      /*
       * Tag cells, determine refine boxes from tagged regions, and
       * load balance in preparation for constructing new refined level.
       */
      if (do_tagging) {

         t_enforce_nesting->start();
         t_find_nesting_restriction->start();
         /*
          * Determine proper nesting box array for specified level.  If the
          * level is the coarsest, the proper nesting box array covers same
          * region as the physical domain.  Otherwise, the nesting boxes are
          * computed from a proper interior of the level.
          */
         findProperNestingData(hierarchy, level_number);
         t_find_nesting_restriction->stop();
         t_enforce_nesting->stop();

         /*
          * Create communication schedule for buffer tags on this level.
          */
         t_bdry_fill_tags_create->start();
         d_bdry_sched_tags[level_number] =
         d_bdry_fill_tags->createSchedule(patch_level);
         t_bdry_fill_tags_create->stop();

         /*
          * Initialize integer tag arrays on level to false.
          */

         patch_level->allocatePatchData(d_tag_indx);
         setTagsOnLevel(d_false_tag,
                        patch_level,
                        d_tag_indx,
                        patch_level->getBoxes());
         /*
          * Perform pre-processing of error estimation data, if appropriate.
          */
         if (errorEstimationUsesTimeIntegration() ) {
            d_tag_init_strategy->
                preprocessErrorEstimation(patch_hierarchy,
                                          level_number,
                                          level_time,
                                          regrid_start_time,
                                          initial_time);
         }

         /*
          * Determine cells needing refinement on level and set tags to true.
          * Because we are constructing a new level, not regridding the level,
          * the coarsest_sync_level argument is always false.
          */
         bool coarsest_sync_level = false;
         tbox::SAMRAI_MPI::barrier();
         t_tag_cells_for_refinement->start();
         d_tag_init_strategy->
            tagCellsForRefinement(patch_hierarchy,
                                  level_number,
                                  level_time,
                                  d_tag_indx,
                                  initial_time,
                                  coarsest_sync_level,
                                  levelCanBeRefined(level_number),
                                  regrid_start_time);
         tbox::SAMRAI_MPI::barrier();
         t_tag_cells_for_refinement->stop();

         /*
          * Check for user-tagged cells that violate proper nesting.
          * except if user specified that the violating tags be ignored.
          */
         if ( d_check_nonrefined_tags != 'i' ) {
            checkNonrefinedTags(patch_level,
                                *d_proper_nesting_complement[level_number]);
         }

         /*
          * Buffer true tagged cells by specified amount which should be
          * sufficient to keep disturbance on refined region until next regrid
          * of the level occurs.
          */
         patch_level->allocatePatchData(d_buf_tag_indx);
         bufferTagsOnLevel(d_true_tag, patch_level, tag_buffer);

         if ( d_extend_tags_to_bdry ) {
            /*
             * Extend tags to physical boundaries whereever tags extend to
             * within a the ghost cell width of the physical boundaries.  This
             * avoids construction of patches whose with ghost cell regions that
             * are only partially outside the physical boundary.
             */
            extendTagsToBoundary(d_true_tag, patch_level,patch_hierarchy);
         }
         patch_level->deallocatePatchData(d_buf_tag_indx);

         /*
          * Determine box array and processor mapping for new fine level.
          */
         findRefinementBoxes(fine_boxes,
                             mapping,
                             patch_hierarchy,
                             level_number);

         /*
          * Deallocate tag arrays and schedule -- no longer needed.
          */
         patch_level->deallocatePatchData(d_tag_indx);
         d_bdry_sched_tags[level_number].setNull();

         /*
          * Write refine boxes to file if requested.
          */
         if (d_write_dumped_level_boxes) {
            d_regrid_box_utility->
              putLevelBoxes(fine_boxes,
                            fine_level_number,
                            d_regrid_box_counter[fine_level_number]);
            d_regrid_box_counter[fine_level_number]++;
         }

     } else {

         /*
          * If tagging is not necessary (do_tagging = false) we simply
          * need to access the level boxes, either from a dumped file or
          * from user-supplied refine boxes, and load balance them before
          * constructing the finer level.
          */
         bool remove_old_fine_level = false;
         readLevelBoxes(fine_boxes,
                        mapping,
                        hierarchy,
                        level_number,
                        level_time,
                        remove_old_fine_level);
     }


     /*
      * Make new finer level (fine_level_number == level_number+1), if appropriate.
      */
     if ( !(fine_boxes.getNumberOfBoxes() == 0) ) {

         hier::IntVector<DIM> ratio_to_level_zero =
            patch_level->getRatio() * getRatioToCoarserLevel(fine_level_number);

         patch_hierarchy->makeNewPatchLevel(fine_level_number,
                                      ratio_to_level_zero,
                                      fine_boxes, mapping);

         if (d_check_overlapping_patches != 'i') {
            checkOverlappingPatches(
               patch_hierarchy->getPatchLevel(fine_level_number));
         }

         d_tag_init_strategy->initializeLevelData(patch_hierarchy,
                                                  fine_level_number,
                                                  level_time,
                                                  levelCanBeRefined(
                                                  fine_level_number),
                                                  initial_time);



         tbox::SAMRAI_MPI::barrier();
         t_reset_hier->start();
         d_tag_init_strategy->resetHierarchyConfiguration(patch_hierarchy,
                                                          fine_level_number,
                                                          fine_level_number);
         tbox::SAMRAI_MPI::barrier();
         t_reset_hier->stop();
      }

   }  // if level cannot be refined, the routine drops through...

   t_make_finer->stop();

}

/*
*************************************************************************
*                                                                       *
* Regrid each level in the hierarchy which is finer than the specified  *
* level.  First, we recursively compute proper nesting boxes for each   *
* level that will be subject to regridding. If the regridding procedure *
* employs time integration, we perform any pre-processing necessary     *
* to regrid the levels.  Then, each level finer than the specified      *
* level is regridded from fine to coarse.  The recursive regridding     *
* procedure is performed by the function regridFinerLevel().  Finally,  *
* after the new hierarchy configuration is set, the application-        *
* specific operations for resetting hierarchy-dependent infomation is   *
* called.                                                               *
*                                                                       *
*************************************************************************
*/

template<int DIM> void GriddingAlgorithm<DIM>::regridAllFinerLevels(
   tbox::Pointer< hier::BasePatchHierarchy<DIM> > hierarchy,
   const int level_number,
   const double regrid_time,
   const tbox::Array<int>& tag_buffer,
   const tbox::Array<double> regrid_start_time,
   const bool level_is_coarsest_sync_level)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!(hierarchy.isNull()));
   TBOX_ASSERT( (level_number>=0)
           && (level_number <= hierarchy->getFinestLevelNumber()) );
   TBOX_ASSERT(!(hierarchy->getPatchLevel(level_number).isNull()));
   TBOX_ASSERT(tag_buffer.getSize() >= level_number+1);
   for (int i = 0; i < tag_buffer.getSize(); i++) {
      TBOX_ASSERT(tag_buffer[i] >= 0);
   }
#endif

   t_regrid_all_finer->start();

   if ( levelCanBeRefined(level_number) ) {

      const tbox::Pointer< hier::PatchLevel<DIM> > patch_level =
         hierarchy->getPatchLevel(level_number);

      t_find_nesting_restriction->start();
      findProperNestingData(hierarchy, level_number);
      t_find_nesting_restriction->stop();

      /*
       * Perform pre-processing of error estimation data, if
       * appropriate.
       */
      if (errorEstimationUsesTimeIntegration() ) {
         for (int ln = level_number;
                  ln <= hierarchy->getFinestLevelNumber(); ln++) {
            if (levelCanBeRefined(ln) ) {
               bool initial_time = false;
               double level_regrid_start_time = 0.;
               if (regrid_start_time.getSize() < ln+1) {
                  level_regrid_start_time =
                     tbox::MathUtilities<double>::getSignalingNaN();
               } else {
                  level_regrid_start_time = regrid_start_time[ln];
               }

               d_tag_init_strategy->
                  preprocessErrorEstimation(hierarchy,
                                            ln,
                                            regrid_time,
                                            level_regrid_start_time,
                                            initial_time);
            }
         }
      }


      /*
       * Recursively regrid each finer level.
       */
      const int finest_level_not_regridded = level_number;
      regridFinerLevel(hierarchy,
                       level_number,
                       regrid_time,
                       finest_level_not_regridded,
                       level_is_coarsest_sync_level,
                       tag_buffer,
                       regrid_start_time);

      /*
       * Invoke application-specific routines to reset information for those
       * levels which have been modified.
       */

      if ( hierarchy->getFinestLevelNumber() >= (level_number+1) ) {
         tbox::SAMRAI_MPI::barrier();
         t_reset_hier->start();
         d_tag_init_strategy->
            resetHierarchyConfiguration(hierarchy,
                                        level_number+1,
                                        hierarchy->
                                        getFinestLevelNumber());
         tbox::SAMRAI_MPI::barrier();
         t_reset_hier->stop();

      }

   } //  if level cannot be refined, the routine drops through...

   t_regrid_all_finer->stop();

}

/*
*************************************************************************
*                                                                       *
* Recursively, regrid each AMR hierarchy level finer than the specified *
* level (indexed by level_number).  The process is as follows:          *
*                                                                       *
* (1) Initialize tags to false on the level.                            *
*                                                                       *
* (2) If a finer level exists, set tag to true on level for each cell   *
*     that is refined.                                                  *
*                                                                       *
* (3) Tag cells for refinement on level by applying application-        *
*     specific error estimation routines.                               *
*                                                                       *
* (4) If a finer level exists, invoke process recursively (i.e.,        *
*     invoke step 1 on next finer level).                               *
*                                                                       *
* (5) (Note we have popped out of recursion at this point).  Buffer     *
*     true tags on current level to keep disturbances on fine grids     *
*     until regridding occurs next.                                     *
*                                                                       *
* (6) Determine box configuration for new finer level, by calling       *
*     findRefinementBoxes() function.                                   *
*                                                                       *
* (7) If a finer level should exist in the hierarchy, create its        *
*     patches from the box description and initialize its data.  If     *
*     necessary, discard old level.                                     *
*                                                                       *
*************************************************************************
*/

template<int DIM> void GriddingAlgorithm<DIM>::regridFinerLevel(
   tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
   const int level_number,
   const double regrid_time,
   const int finest_level_not_regridded,
   const bool level_is_coarsest_sync_level,
   const tbox::Array<int>& tag_buffer,
   const tbox::Array<double>& regrid_start_time)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!(hierarchy.isNull()));
   TBOX_ASSERT( (level_number >= 0)
          && (level_number <= hierarchy->getFinestLevelNumber()) );
   TBOX_ASSERT(!(hierarchy->getPatchLevel(level_number).isNull()));
   TBOX_ASSERT(finest_level_not_regridded >= 0
          && finest_level_not_regridded <= level_number);
   TBOX_ASSERT(tag_buffer.getSize() >= level_number+1);
   for (int i = 0; i < tag_buffer.getSize(); i++) {
      TBOX_ASSERT(tag_buffer[i] >= 0);
   }
#endif

   if ( levelCanBeRefined(level_number) ) {

      int fine_level_number = level_number+1;

      tbox::Pointer< hier::PatchLevel<DIM> >
         patch_level = hierarchy->getPatchLevel(level_number);

      hier::BoxArray<DIM> fine_boxes;
      hier::ProcessorMapping mapping;

      /*
       * The boolean "do_tagging" specifies whether or not tagging will
       * be performed. This will be true except in two circumstances:
       *    1) only user supplied refine boxes are used
       *    2) the boxes are read from a previously dumped file.
       *
       * If either of these circumstances is true, tagging operations
       * are NOT necessary so do_tagging will be set to false.
       *
       * The old level is generally removed when regridding, but
       * some circumstances may warrant keeping the old level.  For
       * example, if the refine region has not changed, there is no
       * need to regenerate the finer level.  The boolean
       * "remove_old_fine_level" specifies if the old level should
       * be removed.
       */
      bool do_tagging = true;
      if (d_read_dumped_level_boxes ||
          d_tag_init_strategy->refineUserBoxInputOnly()) do_tagging = false;
      bool remove_old_fine_level = true;

      /*
       * Tag cells, determine refine boxes from tagged regions, and
       * load balance in preparation for constructing new refined level.
       */
      if (do_tagging) {

         /*
          * Create communication schedule for buffer tags and set tags to
          * false.
          */

        t_bdry_fill_tags_create->start();
        d_bdry_sched_tags[level_number] =
          d_bdry_fill_tags->createSchedule(patch_level);
        t_bdry_fill_tags_create->stop();

        patch_level->allocatePatchData(d_tag_indx);
        setTagsOnLevel(d_false_tag,
                       patch_level, d_tag_indx, patch_level->getBoxes());

        /*
         * Set tags to true for cells that currently cover next finer level.
         * Note that this is not needed for all regridding strategies.  But
         * knowledge of currently refined cells is generally useful to avoid
         * repeated refining and coarsening of cells near boundaries of
         * refined regions where error estimates may hover around the error
         * tolerance. For example, the regridding scheme may require that
         * the error in a cell that is currently refined fall below a
         * certain tolerance  (generally different than the tolerance to
         * refine the cell in the first place) before the cell will be
         * de-refined.
         */

        if ( hierarchy->finerLevelExists(level_number) ) {
           tbox::Pointer< hier::PatchLevel<DIM> > fine_level =
              hierarchy->getPatchLevel(fine_level_number);
           hier::BoxArray<DIM> fine_box_array(fine_level->getBoxes());
           fine_box_array.coarsen(getRatioToCoarserLevel(fine_level_number));
           setTagsOnLevel(d_true_tag, patch_level, d_tag_indx, fine_box_array);
        }

        /*
         * Determine cells needing refinement according to a specific
         * error estimation procedure and set to true.
         *
         * The "level_is_coarsest_sync_level" is provided as an argument
         * to this method.  Provide the additional check of whether the
         * level is not the coarsest level and that it is not a new level
         * in the hierarchy.  If all three conditions are true, the
         * "coarsest_sync_level" argument passed into the tagCells method
         * will be true.  Otherwise, it will be false.
         */

        bool coarsest_sync_level =
           level_is_coarsest_sync_level &&
           level_number > 0 &&
           level_number <= finest_level_not_regridded;

        bool initial_time = false;
        double level_regrid_start_time = 0.;
        if (regrid_start_time.getSize() < level_number+1) {
           level_regrid_start_time =
              tbox::MathUtilities<double>::getSignalingNaN();
        } else {
           level_regrid_start_time = regrid_start_time[level_number];
        }
        d_tag_init_strategy->
           tagCellsForRefinement(hierarchy,
                                 level_number,
                                 regrid_time,
                                 d_tag_indx,
                                 initial_time,
                                 coarsest_sync_level,
                                 levelCanBeRefined(level_number),
                                 level_regrid_start_time);

         /*
          * Check for user-tagged cells that violate proper nesting.
          * except if user specified that the violating tags be ignored.
          */
         if ( d_check_nonrefined_tags != 'i' ) {
            checkNonrefinedTags(patch_level,
                                *d_proper_nesting_complement[level_number]);
         }

        /*
         * Perform regridding recursively on finer levels, if appropriate.
         */
        if ( hierarchy->finerLevelExists(level_number)
             && levelCanBeRefined(fine_level_number) ) {
           regridFinerLevel(hierarchy,
                            fine_level_number,
                            regrid_time,
                            finest_level_not_regridded,
                            false,
                            tag_buffer,
                            regrid_start_time);
        }

         /*
          * Buffer true tagged cells by specified amount which should be
          * sufficient to keep disturbance on refined region until next
          * regrid of the level occurs.
          */
         patch_level->allocatePatchData(d_buf_tag_indx);
         bufferTagsOnLevel(d_true_tag, patch_level, tag_buffer[level_number]);

        /*
         * Add tags built from two levels above to ensure that the the fine
         * level nests it.
         */
 	 if ( hierarchy->finerLevelExists(fine_level_number) ) {
            tbox::Pointer< hier::PatchLevel<DIM> > second_fine_level =
               hierarchy->getPatchLevel(fine_level_number+1);
            hier::BoxArray<DIM> second_fine_boxes(second_fine_level->getBoxes());
            second_fine_boxes.grow(
               getRatioToCoarserLevel(fine_level_number+1)*
               getProperNestingBuffer(fine_level_number));
            second_fine_boxes.coarsen(getRatioToCoarserLevel(fine_level_number+1)*
                                        getRatioToCoarserLevel(fine_level_number));
            setTagsOnLevel(d_true_tag, patch_level, d_tag_indx, second_fine_boxes);
	 }

         if ( d_extend_tags_to_bdry ) {
            /*
             * Extend tags to physical boundaries whereever tags extend to
             * within a the ghost cell width of the physical boundaries.  This
             * avoids construction of patches whose with ghost cell regions that
             * are only partially outside the physical boundary.
             */
            extendTagsToBoundary(d_true_tag, patch_level, hierarchy);
         }

         patch_level->deallocatePatchData(d_buf_tag_indx);

         /*
          * Determine box array containing cells on level with a true tag
          * value.  The box array must be contained in array of proper
          * nesting boxes.
          */
         findRefinementBoxes(fine_boxes, mapping,
                             hierarchy, level_number);

         /*
          * Deallocate tag arrays and schedule; no longer needed on current
          * level.
          */

         patch_level->deallocatePatchData(d_tag_indx);
         d_bdry_sched_tags[level_number].setNull();


         /*
          * Now that we have determined refine boxes, write those to
          * file.
          */
         if (d_write_dumped_level_boxes) {
            d_regrid_box_utility->putLevelBoxes(fine_boxes,
                                                fine_level_number,
                                                d_regrid_box_counter[fine_level_number]);
            d_regrid_box_counter[fine_level_number]++;
         }


      } else {

 	 if ( hierarchy->finerLevelExists(level_number)
	     && levelCanBeRefined(fine_level_number) ) {
	    regridFinerLevel(hierarchy,
                             fine_level_number,
                             regrid_time,
                             finest_level_not_regridded,
                             false,
                             tag_buffer,
                             regrid_start_time);
         }

         /*
          * If tagging is not necessary (do_tagging = false) we simply
          * need to access the level boxes, either from a dumped file or
          * from user-supplied refine boxes, and load balance them before
          * constructing the finer level.
          */
         readLevelBoxes(fine_boxes,
                        mapping,
                        hierarchy,
                        level_number,
                        regrid_time,
                        remove_old_fine_level);
      }

      /*
       * Make new finer level (fine_level_number) if necessary, or remove
       * next finer level if it is no longer needed.
       */

      if ( !(fine_boxes.getNumberOfBoxes() == 0) ) {

         tbox::SAMRAI_MPI::barrier();
         t_regrid_finer_create->start();

         /*
          * Either remove pre-existing fine level from hierarchy and make
          * a new level, or just make a new fine level for hierarchy.
          */

         tbox::Pointer< hier::PatchLevel<DIM> > old_fine_level;

         hier::IntVector<DIM> ratio(patch_level->getRatio()
                               * getRatioToCoarserLevel(fine_level_number));

         if ( hierarchy->finerLevelExists(level_number) ) {
            old_fine_level = hierarchy->getPatchLevel(fine_level_number);
            hierarchy->removePatchLevel(fine_level_number);
            ratio = old_fine_level->getRatio();
         }

         hierarchy->makeNewPatchLevel(fine_level_number, ratio,
                                      fine_boxes, mapping);

         // "false" argument": const bool initial_time = false;
         d_tag_init_strategy->initializeLevelData(hierarchy,
                                                  fine_level_number,
                                                  regrid_time,
                                                  levelCanBeRefined(
                                                     fine_level_number),
                                                  false,
                                                  old_fine_level);

         /*
          * Destroy old patch level, if such a level existed prior to regrid.
          */
         old_fine_level.setNull();
         tbox::SAMRAI_MPI::barrier();
         t_regrid_finer_create->stop();


         if ( do_tagging ) {

            /*
             * If current level (level_number) is coarser than the finest
             * level subject to regridding, then adjust tags for next coarser
             * level (i.e., level_number-1) so that the next finer level
             * (fine_level_number) will be nested within the new level
             * (level_number) on return from recursion.
             */

            if (level_number > finest_level_not_regridded) {
               fine_boxes.coarsen(getRatioToCoarserLevel(fine_level_number)
                                  * getRatioToCoarserLevel(level_number));
               setTagsOnLevel(d_true_tag,
                              hierarchy->getPatchLevel(level_number-1),
                              d_tag_indx,
                              fine_boxes);
            }

         }

      } else {

         /*
          * If there are no boxes for the new fine level, remove the
          * pre-existing fine level if it existed.
          */

         if ( hierarchy->finerLevelExists(level_number)
              && remove_old_fine_level) {
            hierarchy->removePatchLevel(fine_level_number);
         }

      } // if we are not re-gridding the level

   } //  if level cannot be refined, the routine drops through...

}


/*
*************************************************************************
*************************************************************************
*/
template<int DIM> void GriddingAlgorithm<DIM>::checkNonrefinedTags(
   const tbox::Pointer<hier::PatchLevel<DIM> > &level,
   const hier::BoxTree<DIM> &proper_nesting_complement )
{
   /*
    * Check for user-tagged cells that violate proper nesting.
    */
   math::PatchCellDataBasicOps<DIM,int> dataop;
   int maxval = 0;
   for ( typename hier::PatchLevel<DIM>::Iterator pi(level); pi; pi++ ) {
      tbox::Pointer<hier::Patch<DIM> > patch = level->getPatch(*pi);
      const hier::Box<DIM> &pbox = patch->getBox();
      tbox::Pointer<pdat::CellData<DIM,int> > tag_data =
         patch->getPatchData(d_tag_indx);
      hier::BoxList<DIM> overlap_boxes;
      proper_nesting_complement.findOverlapBoxes(
         overlap_boxes, patch->getBox() );
      for ( typename hier::BoxList<DIM>::Iterator bi(overlap_boxes); bi; bi++ ) {
         hier::Box<DIM> ovlap = bi() * pbox;
         maxval = dataop.max( tag_data, ovlap );
         if ( maxval > 0 ) break;
      }
      if ( maxval > 0 ) break;
   }
   maxval = tbox::SAMRAI_MPI::maxReduction(maxval);

   if ( maxval > 0 ) {
      if ( d_check_nonrefined_tags == 'w' ) {
         TBOX_WARNING("User code has tagged cells in\n"
                      <<"violation of nesting requirements.\n"
                      <<"Violating tags will be discarded.\n"
                      <<"See GriddingAlgorithm::checkNonrefinedTags()\n");
      }
      else if ( d_check_nonrefined_tags == 'e' ) {
         TBOX_ERROR("User code has tagged cells in\n"
                    <<"violation of nesting requirements.\n"
                    <<"See GriddingAlgorithm::checkNonrefinedTags()\n");
      }
   }

   return;
}

/*
*************************************************************************
*************************************************************************
*/
template<int DIM> void GriddingAlgorithm<DIM>::checkOverlappingPatches(
   const tbox::Pointer<hier::PatchLevel<DIM> > &level )
{
   hier::BoxList<DIM> level_boxes(level->getBoxes());

   if (level_boxes.boxesIntersect()) {
      if ( d_check_overlapping_patches == 'w' ) {
         TBOX_WARNING("PatchLevel has patches which overlap in index space\n"
                      <<"See GriddingAlgorithm::checkOverlappingPatches()\n");
      }
      else if ( d_check_nonrefined_tags == 'e' ) {
         TBOX_ERROR("PatchLevel has patches which overlap in index space\n"
                    <<"See GriddingAlgorithm::checkOverlappingPatches()\n");
      }
   }

   return;
}

/*
*************************************************************************
*                                                                       *
*   For cases where tagging is not performed read the new level boxes   *
*   either from user input or from stored level boxes.                  *
*                                                                       *
*************************************************************************
*/
template<int DIM> void GriddingAlgorithm<DIM>::readLevelBoxes(
   hier::BoxArray<DIM>& new_level_boxes,
   hier::ProcessorMapping& mapping,
   const tbox::Pointer<hier::PatchHierarchy<DIM> > hierarchy,
   const int level_number,
   const double regrid_time,
   bool& remove_old_fine_level)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!hierarchy.isNull());
   TBOX_ASSERT( (level_number >= 0)
          && (level_number <= hierarchy->getFinestLevelNumber()) );
#endif
   int fine_level_number = level_number + 1;
   hier::BoxArray<DIM> boxes_to_refine;

   if (d_read_dumped_level_boxes) {

     d_regrid_box_utility->
       getLevelBoxes(boxes_to_refine,
                     fine_level_number,
                     d_regrid_box_counter[fine_level_number]);

     d_regrid_box_counter[fine_level_number]++;

   }

   /*
    * Access the user supplied refine boxes.  The
    * "new_level_has_new_boxes" boolean specifies whether the
    * level boxes have changed from the last time
    * getUserSuppliedRefineBoxes() was called.  If they have changed,
    * it returns true.  If they are unchanged, it returns false.
    */
   bool new_level_has_new_boxes = true;
   if (d_tag_init_strategy->refineUserBoxInputOnly()) {

     new_level_has_new_boxes = d_tag_init_strategy->
       getUserSuppliedRefineBoxes(boxes_to_refine,
                                  level_number,
                                  regrid_time);

     boxes_to_refine.refine(getRatioToCoarserLevel(fine_level_number));

   }

   /*
    * If "new_level_has_new_boxes" is false we wish to keep the
    * existing fine level intact.  Avoid further work by setting
    * the parameter "compute_load_balanced_level_boxes" to false
    * and indicate that we want to avoid removing the old fine level
    * by setting "remove_old_fine_level" to false.
    */
   bool compute_load_balanced_level_boxes = true;
   if (!new_level_has_new_boxes) {
     compute_load_balanced_level_boxes = false;
     remove_old_fine_level = false;
   }

   /*
    * If we are using the nonuniform load balance option, we
    * still need to redo the load balance and construct a new level,
    * even if the level boxes have not changed.
    */

   if (d_load_balancer->getLoadBalanceDependsOnPatchData(fine_level_number)
       && boxes_to_refine.getNumberOfBoxes() > 0) {
     compute_load_balanced_level_boxes = true;
     remove_old_fine_level = true;
   }

   /*
    * If the boxes_to_refine are empty, this implies that no
    * refinement is desired so a new finer level will NOT be
    * constructed.  In this case, avoid load balance steps and
    * specify that we want to remove the old fine level.
    */
   if (boxes_to_refine.getNumberOfBoxes() == 0) {
     compute_load_balanced_level_boxes = false;
     remove_old_fine_level = true;
   }

   if (compute_load_balanced_level_boxes) {

     hier::BoxList<DIM> fine_box_list(boxes_to_refine);
     tbox::Pointer<hier::PatchLevel<DIM> > patch_level =
        hierarchy->getPatchLevel(level_number);

     hier::IntVector<DIM> ratio_to_level_zero =
       patch_level->getRatio() * getRatioToCoarserLevel(fine_level_number);

     hier::BoxArray<DIM> physical_domain;
     hierarchy->getGridGeometry()->
       computePhysicalDomain(physical_domain,
                             ratio_to_level_zero);

     hier::IntVector<DIM>
       patch_cut_factor(getRatioToCoarserLevel(fine_level_number));

     hier::IntVector<DIM> smallest_patch;
     hier::IntVector<DIM> smallest_box_to_refine;
     hier::IntVector<DIM> largest_patch;
     hier::IntVector<DIM> extend_ghosts;
     // "false" argument: for_building_finer level = false
     getGriddingParameters(smallest_patch,
                           smallest_box_to_refine,
                           largest_patch,
                           extend_ghosts,
                           hierarchy,
                           fine_level_number,
                           false);

     tbox::SAMRAI_MPI::barrier();
     t_load_balance->start();
     d_load_balancer->loadBalanceBoxes(new_level_boxes, mapping,
                                       fine_box_list, hierarchy,
                                       fine_level_number, physical_domain,
                                       ratio_to_level_zero,
                                       smallest_patch,
                                       getLargestPatchSize(fine_level_number),
                                       patch_cut_factor, extend_ghosts);
     t_load_balance->stop();
   }
}



/*
*************************************************************************
*                                                                       *
* Flat (non-recursive) version of findProperNestingBoxes.               *
*                                                                       *
* Compute d_proper_nesting_complement for the given level and           *
* finer levels.  The proper nesting for the given level is its          *
* level domain squeezed in by its proper nesting buffer.                *
* The proper nesting for each finer level is the previous one,          *
* squeezed in by the buffer of the finer level.  This method            *
* computes the complements of those proper nestings.                    *
*                                                                       *
*************************************************************************
*/

template<int DIM> void GriddingAlgorithm<DIM>::findProperNestingData(
   const tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
   const int base_ln)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!(hierarchy.isNull()));
#endif

   tbox::Pointer< hier::PatchLevel<DIM> > base_level =
      hierarchy->getPatchLevel(base_ln);
   const hier::BoxArray<DIM> &base_domain = base_level->getPhysicalDomain();

   tbox::Array< hier::BoxList<DIM> >
      proper_nesting_complement( d_proper_nesting_complement.size() );

   int ln = base_ln;

   do {

      tbox::Pointer< hier::PatchLevel<DIM> > level =
         hierarchy->getPatchLevel(ln);

      // Shorthands definitions.
      hier::BoxList<DIM> &complement = proper_nesting_complement[ln];
      hier::BoxList<DIM> &propernest = d_proper_nesting_boxes[ln];

      /*
       * The base level computes the complement from its own boxes
       * while finer levels start from the previous level's complement.
       */
      if ( ln == base_ln ) {
         complement = hier::BoxList<DIM>( base_domain );
         base_level->getBoxTree()->removeIntersections(complement);
      }
      else {
         complement = proper_nesting_complement[ln-1];
         complement.refine(getRatioToCoarserLevel(ln));
      }

      /*
       * Grow boxes in complement and throw away pieces outside of physical
       * domain.  Then, proper nesting boxes are complement of the complement.
       */

      complement.grow(getProperNestingBuffer(ln));
      t_intersect_boxes_find_proper->start();
      complement.intersectBoxes( level->getPhysicalDomain() );
      t_intersect_boxes_find_proper->stop();

      propernest = hier::BoxList<DIM>( level->getPhysicalDomain() );
      t_remove_intersections_find_proper->start();
      propernest.removeIntersections(complement);
      t_remove_intersections_find_proper->stop();

      ++ln;

   } while ( hierarchy->finerLevelExists(ln-1) &&
             levelCanBeRefined(ln) );

   /*
    * Above computation cannot contain domain complement.
    * But it is needed so add it now.
    */

   hier::Box<DIM> universe;
   for ( int i=0; i<base_domain.size(); ++i ) {
      universe += base_domain[i];
   }
   universe.grow(universe.numberCells());

   hier::BoxList<DIM> domain_complement = universe;
   domain_complement.removeIntersections( base_domain );
   proper_nesting_complement[base_ln].copyItems( domain_complement );
   d_proper_nesting_complement[base_ln] =
      new hier::BoxTree<DIM>( hier::BoxArray<DIM>(proper_nesting_complement[base_ln]) );
   const int end_ln = ln;
   for ( ln=base_ln+1; ln<end_ln; ++ln ) {
      tbox::Pointer< hier::PatchLevel<DIM> > level =
         hierarchy->getPatchLevel(ln);
      domain_complement.refine( level->getRatioToCoarserLevel() );
      proper_nesting_complement[ln].copyItems( domain_complement );
      d_proper_nesting_complement[ln] =
         new hier::BoxTree<DIM>( hier::BoxArray<DIM>(proper_nesting_complement[ln]) );
   }

   return;
}




/*
*************************************************************************
*                                                                       *
* Set each integer value in specified tag array to tag_value where      *
* patch level intersects given box array.                               *
*                                                                       *
*************************************************************************
*/

template<int DIM> void GriddingAlgorithm<DIM>::setTagsOnLevel(
   const int tag_value,
   const tbox::Pointer< hier::PatchLevel<DIM> > level,
   const int index,
   const hier::BoxArray<DIM>& boxes,
   const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( (tag_value == d_true_tag) || (tag_value == d_false_tag) );
   TBOX_ASSERT(!(level.isNull()));
   TBOX_ASSERT(index == d_tag_indx || index == d_buf_tag_indx);
#endif
   /*
    * Start timer for this method.
    */
   t_set_tags->start();

   const int n_boxes = boxes.getNumberOfBoxes();
   for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
      tbox::Pointer< hier::Patch<DIM> > patch = level->getPatch(ip());

      tbox::Pointer< pdat::CellData<DIM,int> >
         tag_data = patch->getPatchData(index);
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT( !(tag_data.isNull()) );
#endif

      hier::Box<DIM> set_box( (interior_only ? tag_data->getBox() :
                                          tag_data->getGhostBox()) );

      for (int ib = 0; ib < n_boxes; ib++) {
         hier::Box<DIM> intersection = boxes[ib] * set_box;
         if ( !(intersection.empty()) ) {
            tag_data->fill(tag_value, intersection);
         }
      }
   }
   t_set_tags->stop();
}

/*
*************************************************************************
*                                                                       *
* Buffer each integer tag with given value on the patch level by the    *
* specified buffer size.  Note that the patch data indexed by           *
* d_buf_tag_indx is used temporarily to buffer the tag data. The        *
* communication of ghost cell (i.e., buffer) information forces all     *
* tags on all patch interiors to represent a consistent buffering of    *
* the original configuration of tagged cells.                           *
*                                                                       *
*************************************************************************
*/

template<int DIM> void GriddingAlgorithm<DIM>::bufferTagsOnLevel(
   const int tag_value,
   const tbox::Pointer< hier::PatchLevel<DIM> > level,
   const int buffer_size) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( (tag_value == d_true_tag) || (tag_value == d_false_tag) );
   TBOX_ASSERT(!(level.isNull()));
   TBOX_ASSERT(buffer_size >= 0);
#endif
   /*
    * Start timer for this method.
    */
   t_buffer_tags->start();

   /*
    * Set temporary buffered tags based on buffer width and
    * distance from actual tags.
    */
   const int not_tag = ((tag_value == d_true_tag) ? d_false_tag : d_true_tag);
   for (typename hier::PatchLevel<DIM>::Iterator ip1(level); ip1; ip1++) {
      tbox::Pointer< hier::Patch<DIM> > patch = level->getPatch(ip1());

      tbox::Pointer< pdat::CellData<DIM,int> >
         buf_tag_data = patch->getPatchData(d_buf_tag_indx);
      tbox::Pointer< pdat::CellData<DIM,int> >
         tag_data = patch->getPatchData(d_tag_indx);

      buf_tag_data->fillAll(not_tag);

      hier::Box<DIM> interior = patch->getBox();

      for (int bc = buffer_size; bc >= 0; bc--) {

         int fill_val = buffer_size - bc + 1;

         for (pdat::CellIterator<DIM> ic(interior); ic; ic++) {
            if ( (*tag_data)(ic()) == tag_value ) {
               hier::Box<DIM> buf_box(ic()-bc, ic()+bc);
               buf_tag_data->fill(fill_val, buf_box * interior);
            }
         }

      }

   }

   /*
    * Communicate boundary data for buffered tag array so that tags
    * near patch boundaries will become buffered properly.
    */
   const double dummy_time = 0.0;

   t_bdry_fill_tags_comm->start();
   d_bdry_sched_tags[level->getLevelNumber()]->fillData(dummy_time);
   t_bdry_fill_tags_comm->stop();

   /*
    * Buffer tags on patch interior according to buffered tag data.
    */
   for (typename hier::PatchLevel<DIM>::Iterator ip2(level); ip2; ip2++) {
      tbox::Pointer< hier::Patch<DIM> > patch = level->getPatch(ip2());

      tbox::Pointer< pdat::CellData<DIM,int> > buf_tag_data =
                                 patch->getPatchData(d_buf_tag_indx);
      tbox::Pointer< pdat::CellData<DIM,int> > tag_data =
                                 patch->getPatchData(d_tag_indx);

      hier::Box<DIM> buf_tag_box = buf_tag_data->getGhostBox();
      hier::Box<DIM> tag_box = tag_data->getBox();

      /*
       * Set all interior tags to tag value where buffer tags non-zero.
       */
      for (pdat::CellIterator<DIM> ic(tag_box); ic; ic++) {
         (*tag_data)(ic()) = ( (*buf_tag_data)(ic()) ? tag_value : not_tag );
      }

      /*
       * Set all interior tags in buffers around tags in ghosts.
       */
      for (pdat::CellIterator<DIM> ic2(buf_tag_box); ic2; ic2++) {
         int tval = (*buf_tag_data)(ic2());
         if ( tval > 1 ) {
            int buf_size = tval - 1;
            hier::Box<DIM> buf_box(ic2()-buf_size, ic2()+buf_size);
            tag_data->fill(tag_value, buf_box);
         }
      }
   }

   t_buffer_tags->stop();
}

/*
*************************************************************************
*                                                                       *
*   If any tags are set that are within "extend_ghosts" cells from the  *
*   boundary, fill them to the boundary using the following procedure:  *
*                                                                       *
*   1) Grow physical domain boxes by "extend_ghosts"                    *
*   2) Intersect patch box with the grown physical domain boxes         *
*   3) On this intersected box:                                         *
*   a) for any tag that lies in this box, create a small 1 cell box     *
*   around it                                                           *
*   b) grow the 1 cell box to the domain boundaries                     *
*   c) set tags inside this grown box.                                  *
*                                                                       *
*************************************************************************
*/
template<int DIM> void GriddingAlgorithm<DIM>::extendTagsToBoundary(
   const int tag_value,
   const tbox::Pointer< hier::PatchLevel<DIM> > level,
   const tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( (tag_value == d_true_tag) || (tag_value == d_false_tag) );
   TBOX_ASSERT(!(level.isNull()));
#endif

   int level_number = level->getLevelNumber();
   int fine_level_number = level_number + 1;

   hier::IntVector<DIM> smallest_patch;
   hier::IntVector<DIM> smallest_box_to_refine;
   hier::IntVector<DIM> largest_patch;
   hier::IntVector<DIM> extend_ghosts;
   // "true" argument: for_building_finer level = true
   getGriddingParameters(smallest_patch,
                         smallest_box_to_refine,
                         largest_patch,
                         extend_ghosts,
                         hierarchy,
                         level_number,
                         true);

   /*
    * Make sure coarsened version of level will accomodate the ghost cell
    * needs of the variables.
    */
   for (int i = 0; i < DIM; i++) {
     int extend = extend_ghosts(i)
       / d_ratio_to_coarser[fine_level_number](i);
     if ( !((extend_ghosts(i)
             - extend * d_ratio_to_coarser[fine_level_number](i)) == 0) ) {
       extend++;
     }
     extend_ghosts(i) = extend;
   }

   /*
    * If any tags are set that are within "extend_ghosts" cells from the
    * boundary, fill them to the boundary using the following procedure:
    *
    * 1) Grow physical domain boxes by "extend_ghosts"
    * 2) Intersect patch box with the grown physical domain boxes
    * 3) On this intersected box:
    *    a) for any tag that lies in this box, create a small 1 cell box
    *       around it
    *    b) grow the 1 cell box to the domain boundaries
    *    c) set tags inside this grown box.
    */
   const int not_tag = ((tag_value == d_true_tag) ? d_false_tag : d_true_tag);
   hier::BoxList<DIM> level_domain = level->getPhysicalDomain();
   hier::BoxList<DIM> grown_level_domain(level_domain);
   grown_level_domain.grow(extend_ghosts);

   for (typename hier::PatchLevel<DIM>::Iterator ip1(level); ip1; ip1++) {
     tbox::Pointer<hier::Patch<DIM> > patch = level->getPatch(ip1());

     if (patch->getPatchGeometry()->intersectsPhysicalBoundary()) {

        tbox::Pointer< pdat::CellData<DIM,int> >
           buf_tag_data = patch->getPatchData(d_buf_tag_indx);
        tbox::Pointer< pdat::CellData<DIM,int> >
           tag_data = patch->getPatchData(d_tag_indx);

        buf_tag_data->fillAll(not_tag);

        hier::Box<DIM> pbox = patch->getBox();
        hier::BoxList<DIM> patch_extend_ghost_boxes(grown_level_domain);
        patch_extend_ghost_boxes.intersectBoxes(pbox);

        for (typename hier::BoxList<DIM>::Iterator
                b(patch_extend_ghost_boxes); b; b++) {
           hier::Box<DIM> extend_box = b();

           for (typename pdat::CellIterator<DIM> ic1(extend_box); ic1; ic1++) {
              if ( (*tag_data)(ic1()) != not_tag ) {
                 hier::Box<DIM> tag_box(ic1(),ic1());

                 hier::BoxUtilities<DIM>::extendBoxToDomainBoundary(tag_box,
                                                               level_domain,
                                                               extend_ghosts);
                 buf_tag_data->fill(tag_value,tag_box);
              }
           }
        }

        /*
         * Set all interior tags to tag value where buffer tags non-zero.
         */
        hier::Box<DIM> tag_box = tag_data->getBox();
        for (typename pdat::CellIterator<DIM> ic2(tag_box); ic2; ic2++) {
           (*tag_data)(ic2()) = ( (*buf_tag_data)(ic2()) ? tag_value : not_tag );
        }

     } // if patch touches boundary

   } // loop over patches


}

/*
*************************************************************************
*                                                                       *
* Given a patch level, determine an appropriate array of boxes from     *
* which a new finer level may be constructed.  That is, find an array   *
* of boxes that covers all tags having the specified tag value.  Note   *
* that it is assumed that the integer tag arrays have been set          *
* properly; i.e., cells have been tagged through error estimation and   *
* the tags have been buffered to ensure disturbances remain on fine     *
* level until next regrid occurs.  Note that load balancing is          *
* performed once an appropriate list of boxes containing the tags is    *
* found.  This procedure massages the list of boxes further and then    *
* assigns each to a single processor (i.e., the mapping).               *
*                                                                       *
*************************************************************************
*/

template<int DIM> void GriddingAlgorithm<DIM>::findRefinementBoxes(
   hier::BoxArray<DIM>& fine_boxes,
   hier::ProcessorMapping& mapping,
   const tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
   const int coarse_level_number ) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( (coarse_level_number >= 0)
           && (coarse_level_number <= hierarchy->getFinestLevelNumber()) );
   TBOX_ASSERT(!(hierarchy.isNull()));
   TBOX_ASSERT(!(hierarchy->getPatchLevel(coarse_level_number).isNull()));
#endif
   /*
    * Start timer for this method.
    */
   t_find_refinement->start();

   const int fine_level_number = coarse_level_number + 1;

   /*
    * Construct list of boxes covering the true tags on the level.
    * Note that box list will be contained in the bounding box
    * but will not be contained in the list of proper nesting boxes,
    * in general.  So we intersect the box list against the list of
    * nesting boxes.  Note that this may produce boxes which are too
    * small.  Thus, boxes are regrown later.
    */

   hier::IntVector<DIM> smallest_patch;
   hier::IntVector<DIM> smallest_box_to_refine;
   hier::IntVector<DIM> largest_patch;
   hier::IntVector<DIM> extend_ghosts;
   // "true" argument: for_building_finer level = true
   getGriddingParameters(smallest_patch,
                         smallest_box_to_refine,
                         largest_patch,
                         extend_ghosts,
                         hierarchy,
                         coarse_level_number,
                         true);

   hier::BoxList<DIM> box_list;
   tbox::Pointer< hier::PatchLevel<DIM> > level =
   hierarchy->getPatchLevel(coarse_level_number);

   /*
    * Determine single smallest bounding box for all nesting boxes.
    */
   hier::Box<DIM> bounding_box;
   {
      const hier::BoxArray<DIM> &level_boxes = level->getBoxes();
      for ( int i=0; i<level_boxes.size(); ++i ) {
         bounding_box += level_boxes[i];
      }
   }

   if ( d_barrier_before_clustering ) {
      tbox::SAMRAI_MPI::barrier();
   }

   t_find_boxes_containing_tags->start();
   d_box_generator->findBoxesContainingTags(
         box_list, level, d_tag_indx, d_true_tag, bounding_box,
         smallest_box_to_refine,
         getEfficiencyTolerance(coarse_level_number),
         getCombineEfficiency(coarse_level_number));
   t_find_boxes_containing_tags->stop();

   if (!box_list.isEmpty()) {

      t_box_massage->start();

      if(d_sort_boxes_after_clustering) {
         /*
          * Temporarily stop timer and sort box for
          * performance comparison between ABR and BR.
          * Theoretically, this step should take the
          * same ammount of time on all processors.
          */
         t_find_refinement->stop();
         hier::BoxArray<DIM> box_array(box_list);
         hier::Box<DIM> *box_ptr = &box_array[0];
         qsort( (void*)box_ptr,
                box_array.size(),
                sizeof(hier::Box<DIM>),
                qsortBoxCompare );
         box_list = hier::BoxList<DIM>(box_array);
         t_find_refinement->start();
      }

      t_enforce_nesting->start();
      t_apply_nesting_restriction->start();

      t_intersect_boxes_find_refinement->start();
      box_list.intersectBoxes(d_proper_nesting_boxes[coarse_level_number]);
      // d_proper_nesting_complement[coarse_level_number]->removeIntersections(box_list);
      t_intersect_boxes_find_refinement->stop();
      t_enforce_nesting->start();

      /*
       * If overlaps are to be avoided, at the possible expense of constructing
       * boxes smaller than the minimum patch size, set the
       * "grow_after_nesting" argument to false in input.  By default, this
       * argument is true.  If set false, the boxes will not be grown within
       * the nesting domain and no overlaps will be introduced.
       */

      t_extend_to_domain_boundary->start();
      hier::BoxList<DIM> level_domain(level->getPhysicalDomain());
      for (int i = 0; i < DIM; i++) {
         int extend = extend_ghosts(i)
            / d_ratio_to_coarser[fine_level_number](i);
         if ( !((extend_ghosts(i)
                 - extend * d_ratio_to_coarser[fine_level_number](i)) == 0) ) {
            extend++;
         }
         extend_ghosts(i) = extend;
      }

      (void) hier::BoxUtilities<DIM>::extendBoxesToDomainBoundary(
         box_list,
         level_domain,
         extend_ghosts);
      t_extend_to_domain_boundary->stop();

      t_apply_nesting_restriction->stop();

      if ( d_coalesce_boxes ) {
         t_coalesce_boxes->start();
         box_list.coalesceBoxes();
         t_coalesce_boxes->stop();
      }

      if (!d_allow_patches_smaller_than_minimum_size_to_prevent_overlaps) {
         t_grow_boxes_within_domain->start();
         hier::BoxUtilities<DIM>::growBoxesWithinDomain(
            box_list,
            d_proper_nesting_boxes[coarse_level_number],
            smallest_box_to_refine);
         t_grow_boxes_within_domain->stop();
      } else {
         hier::IntVector<DIM> periodic_dirs(
            hierarchy->getGridGeometry()->getPeriodicShift());

         bool need_to_grow = false;
         hier::IntVector<DIM> min_size(1);
         for (int i = 0; i < DIM; i++) {
            if (periodic_dirs(i)) {
               need_to_grow = true;
               min_size(i) = smallest_box_to_refine(i);
            }
         }

         if (need_to_grow) {
            hier::BoxUtilities<DIM>::growBoxesWithinDomain(
               box_list,
               d_proper_nesting_boxes[coarse_level_number],
               min_size);
         }
      }

      t_box_massage->stop();

      if ( !box_list.isEmpty() ) {
         /*
          * Refine boxes to index space of next finer level and load balance.
          */

         t_before_load_balance->start();

         box_list.refine(getRatioToCoarserLevel(fine_level_number));

         hier::IntVector<DIM> ratio_to_level_zero =
            level->getRatio() * getRatioToCoarserLevel(fine_level_number);

         hier::BoxArray<DIM> physical_domain;
         hierarchy->getGridGeometry()->computePhysicalDomain(physical_domain,
                                                             ratio_to_level_zero);

         hier::IntVector<DIM> patch_cut_factor(getRatioToCoarserLevel(fine_level_number));
         // "false" argument: for_building_finer level = false
         getGriddingParameters(smallest_patch,
                               smallest_box_to_refine,
                               largest_patch,
                               extend_ghosts,
                               hierarchy,
                               fine_level_number,
                               false);
         t_before_load_balance->stop();

         tbox::SAMRAI_MPI::barrier();
         t_load_balance->start();
         d_load_balancer->loadBalanceBoxes(
            fine_boxes, mapping,
            box_list, hierarchy, fine_level_number, physical_domain,
            ratio_to_level_zero,
            smallest_patch,
            getLargestPatchSize(fine_level_number),
            patch_cut_factor, extend_ghosts);
         t_load_balance->stop();
      }

   } else {

      fine_boxes = hier::BoxArray<DIM>(box_list);

   }

   t_find_refinement->stop();
}

/*
*************************************************************************
*                                                                       *
* Set patch size and ghost cell information needed to create new        *
* mesh levels.  The maximum number of ghost cells over all variables    *
* is used to compute the smallest patch size allowed and the extent to  *
* which patches may be extended to touch the physical boundary.  This   *
* avoids problems in setting ghost cell data that may occur when ghost  *
* cell regions intersect the physical boundary in strange ways.         *
*                                                                       *
*************************************************************************
*/

template<int DIM> void GriddingAlgorithm<DIM>::getGriddingParameters(
   hier::IntVector<DIM>& smallest_patch,
   hier::IntVector<DIM>& smallest_box_to_refine,
   hier::IntVector<DIM>& largest_patch,
   hier::IntVector<DIM>& extend_ghosts,
   const tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
   const int level_number,
   const bool for_building_finer) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!(hierarchy.isNull()));
   TBOX_ASSERT( (level_number >= 0) && (level_number < d_max_levels) );
#endif
   /*
    * Determine maximum ghost cell width needed over all variables
    * currently known to the patch descriptor, and set the smallest
    * patch size.  The maximum number of ghosts is multiplied by the
    * error coarsen ratio (which should always be 1 unless regridding
    * uses error estimation).  This assures that when levels are
    * coarsened during error estimation, the coarser level patches
    * will meet the ghost cell constraint.
    */
   smallest_patch = getSmallestPatchSize(level_number);
   hier::IntVector<DIM> max_ghosts =
      hierarchy->getPatchDescriptor()->getMaxGhostWidth();
   max_ghosts = max_ghosts * d_tag_init_strategy->getErrorCoarsenRatio();
   smallest_patch = getSmallestPatchSize(level_number);
   if (!d_allow_patches_smaller_than_ghostwidth) {
      smallest_patch.max(max_ghosts);
   } else {
      hier::IntVector<DIM> periodic_dirs(
         hierarchy->getGridGeometry()->getPeriodicShift());

      for (int i = 0; i < DIM; i++) {
         if (periodic_dirs(i)) {
            smallest_patch(i) =
               tbox::MathUtilities<int>::Max(smallest_patch(i), max_ghosts(i));
         }
      }
   }

   /*
    * Set largest patch size.
    */
   largest_patch = getLargestPatchSize(level_number);
   largest_patch.max(smallest_patch);

   /*
    * Set the smallest box to refine based on the number of cells that
    * coarsened patches must accomodate to meet ghost cell needs of variables.
    * On the finest level, the smallest box to refine is the smallest patch.
    * On coarser levels, it is a function of the error coarsen ratio and
    * the ratio to the next finer level.
    *
    * If we are accessing gridding parameters for a level that is being
    * reconstructed, the smallest box to refine is not applicable so we
    * set it to -1 to indicate an invalid entry in case it is used.
    */
   if (for_building_finer) {

      if (level_number < d_max_levels-1) {

         int fine_level_number = level_number + 1;

         for (int i = 0; i < DIM; i++) {
            int den = d_ratio_to_coarser[fine_level_number](i)
               / d_tag_init_strategy->getErrorCoarsenRatio();
            int sz = max_ghosts(i)/den;
            if ( max_ghosts(i) - sz * den ) sz++;

            smallest_box_to_refine(i) = tbox::MathUtilities<int>::Max( sz,
                smallest_patch(i)/d_ratio_to_coarser[fine_level_number](i) );
         }

      } else {

        smallest_box_to_refine = smallest_patch;

      }

   } else {

      smallest_box_to_refine = hier::IntVector<DIM>(-1);

   }

   /*
    * Determine number of cells box may be extended to physical
    * domain boundary to accomodate ghost cells.
    */
   extend_ghosts = max_ghosts
                   * d_tag_init_strategy->getErrorCoarsenRatio();

}


/*
*************************************************************************
*                                                                       *
* Print out all attributes of class instance for debugging.             *
*                                                                       *
*************************************************************************
*/

template<int DIM> void GriddingAlgorithm<DIM>::printClassData(std::ostream& os) const
{
   os << "\nGriddingAlgorithm<DIM>::printClassData..." << std::endl;
   os << "   static data members:" << std::endl;
   os << "      s_tag_indx = " << s_tag_indx << std::endl;
   os << "      s_buf_tag_indx = " << s_buf_tag_indx << std::endl;
   os << "GriddingAlgorithm<DIM>: this = "
      << (GriddingAlgorithm<DIM>*)this << std::endl;
   os << "d_object_name = " << d_object_name << std::endl;
   os << "d_tag_init_strategy = "
      << (TagAndInitializeStrategy<DIM>*)d_tag_init_strategy << std::endl;
   os << "d_box_generator = "
      << (BoxGeneratorStrategy<DIM>*)d_box_generator << std::endl;
   os << "d_load_balancer = "
      << (LoadBalanceStrategy<DIM>*)d_load_balancer << std::endl;
   os << "d_tag = " << d_tag.getPointer() << std::endl;
   os << "d_tag_indx = " << d_tag_indx << std::endl;
   os << "d_buf_tag_indx = " << d_buf_tag_indx << std::endl;
   os << "d_true_tag = " << d_true_tag << std::endl;
   os << "d_false_tag = " << d_false_tag << std::endl;
   os << "d_max_levels = " << d_max_levels << std::endl;

   int ln;

   os << "d_proper_nesting_buffer..." << std::endl;
   for (ln = 0; ln < d_proper_nesting_buffer.getSize(); ln++) {
      os << "    d_proper_nesting_buffer[" << ln << "] = "
         << d_proper_nesting_buffer[ln] << std::endl;
   }

   os << "d_efficiency_tolerance..." << std::endl;
   for (ln = 0; ln < d_efficiency_tolerance.getSize(); ln++) {
      os << "    d_efficiency_tolerance[" << ln << "] = "
         << d_efficiency_tolerance[ln] << std::endl;
   }
   os << "d_combine_efficiency..." << std::endl;
   for (ln = 0; ln < d_combine_efficiency.getSize(); ln++) {
      os << "    d_combine_efficiency[" << ln << "] = "
         << d_combine_efficiency[ln] << std::endl;
   }

   os << "d_smallest_patch_size..." << std::endl;
   for (ln = 0; ln < d_smallest_patch_size.getSize(); ln++) {
      os << "    d_smallest_patch_size[" << ln << "] = "
         << d_smallest_patch_size[ln] << std::endl;
   }

   os << "d_largest_patch_size..." << std::endl;
   for (ln = 0; ln < d_largest_patch_size.getSize(); ln++) {
      os << "    d_largest_patch_size[" << ln << "] = "
         << d_largest_patch_size[ln] << std::endl;
   }

   os << "d_write_dumped_level_boxes = " << d_write_dumped_level_boxes << std::endl;
   os << "d_read_dumped_level_boxes = " << d_read_dumped_level_boxes << std::endl;
   os << "d_regrid_boxes_filename = " << d_regrid_boxes_filename << std::endl;
   for (int il = 0; il < d_regrid_box_counter.getSize(); il++) {
      os << "d_regrid_box_counter[" << il << "] = "
         << d_regrid_box_counter[il] << std::endl;
   }
}

/*
*************************************************************************
*                                                                       *
* Write out class version number and data members to database.          *
*                                                                       *
*************************************************************************
*/

template<int DIM> void GriddingAlgorithm<DIM>::putToDatabase(tbox::Pointer<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!db.isNull());
#endif

   db->putInteger("ALGS_GRIDDING_ALGORITHM_VERSION",
                   ALGS_GRIDDING_ALGORITHM_VERSION);

   db->putInteger("d_true_tag", d_true_tag);
   db->putInteger("d_false_tag", d_false_tag);
   db->putInteger("d_max_levels", d_max_levels);

   int ln;

   tbox::Pointer<tbox::Database> ratio_to_coarser_db =
      db->putDatabase("d_ratio_to_coarser");
   for (ln = 1; ln < d_max_levels; ln++) {
      int* temp_ratio_to_coarser = d_ratio_to_coarser[ln];
      std::string level_name = "level_" + tbox::Utilities::intToString(ln);
      ratio_to_coarser_db->putIntegerArray(level_name,
                                           temp_ratio_to_coarser, DIM);
   }

   db->putDoubleArray("d_efficiency_tolerance", d_efficiency_tolerance);
   db->putDoubleArray("d_combine_efficiency", d_combine_efficiency);

   tbox::Pointer<tbox::Database> largest_patch_db =
      db->putDatabase("d_largest_patch_size");
   for (ln = 0; ln < d_largest_patch_size.getSize(); ln++) {
      int* tmp_array = d_largest_patch_size[ln];
      std::string level_name = "level_" + tbox::Utilities::intToString(ln);
      largest_patch_db->putIntegerArray(level_name, tmp_array, DIM);
   }


   tbox::Pointer<tbox::Database> smallest_patch_db =
      db->putDatabase("d_smallest_patch_size");
   for (ln = 0; ln < d_smallest_patch_size.getSize(); ln++) {
      int* tmp_array = d_smallest_patch_size[ln];
      std::string level_name = "level_" + tbox::Utilities::intToString(ln);
      smallest_patch_db->putIntegerArray(level_name, tmp_array, DIM);
   }

   db->putIntegerArray("d_proper_nesting_buffer", d_proper_nesting_buffer);
}

/*
*************************************************************************
*                                                                       *
* If simulation is not from restart, read data from input database.     *
* Otherwise, override data members initialized from restart with        *
* values in the input database.                                         *
*                                                                       *
*************************************************************************
*/

template<int DIM> void GriddingAlgorithm<DIM>::getFromInput(
   tbox::Pointer<tbox::Database> db,
   bool is_from_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!db.isNull());
#endif

   int ln;

   /*
    * Read input for maximum number of levels.
    */

   int new_max_levels = d_max_levels;

   if (db->keyExists("max_levels")) {
      new_max_levels = db->getInteger("max_levels");
   } else {
      if (!is_from_restart) {
         TBOX_ERROR(d_object_name << ":  "
                    << "Key data `max_levels' not found in input.");
      }
   }
   if (new_max_levels < 1) {
      TBOX_ERROR(d_object_name << ":  "
                 << "Key data `max_levels' found in input is < 1.");
   }

   int ln_lo = 1;
   int ln_hi = new_max_levels-1;

   if (is_from_restart) {
      if (new_max_levels > d_max_levels) {
         ln_lo = d_max_levels;
      } else {
         ln_hi = ln_lo - 1;
      }
   }

   /*
    * Read input for ratios between levels.  Default is 1 for level zero
    * ratio to coarser.  If starting from restart, coarsen ratio data
    * is only read in for new levels if new max levels is greater than
    * old max levels.
    */

   d_ratio_to_coarser.resizeArray(new_max_levels);
   d_ratio_to_coarser[0] = hier::IntVector<DIM>(1);

   tbox::Pointer<tbox::Database> ratio_to_coarser_db;
   if ( (is_from_restart && (new_max_levels > d_max_levels))
        || (new_max_levels > 1)) {
      ratio_to_coarser_db = db->getDatabase("ratio_to_coarser");
      if (ratio_to_coarser_db.isNull()) {
         TBOX_ERROR(d_object_name << ":  "
                    << "Key data `ratio_to_coarser' not found in input.");
      }
   }

   tbox::Array<std::string> ratio_to_coarser_patch_keys;
   int num_ratio_to_coarser_keys = 0;

   if(!ratio_to_coarser_db.isNull()) {
      ratio_to_coarser_patch_keys = ratio_to_coarser_db -> getAllKeys();
      num_ratio_to_coarser_keys = ratio_to_coarser_patch_keys.getSize();
   }

   d_ratio_to_coarser.resizeArray(new_max_levels);
   d_ratio_to_coarser[0] = hier::IntVector<DIM>(1);

   for (ln = ln_lo; ln <= tbox::MathUtilities<int>::Min(num_ratio_to_coarser_keys, ln_hi); ln++) {

#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(!ratio_to_coarser_db.isNull());
#endif
      std::string level_name = "level_" + tbox::Utilities::intToString(ln);

      if (!ratio_to_coarser_db->keyExists(level_name)) {
          TBOX_ERROR(d_object_name << ":  "
                     <<"Key data `" << level_name
                     << "' not found in ratio_to_coarser input.");
      }

      int* temp_ratio_to_coarser = d_ratio_to_coarser[ln];
      ratio_to_coarser_db->getIntegerArray(level_name,
                                           temp_ratio_to_coarser, DIM);

      hier::IntVector<DIM> test =
         hier::IntVector<DIM>::min(d_ratio_to_coarser[ln], hier::IntVector<DIM>(1));
      if (test.max() < 1) {
         TBOX_ERROR(d_object_name << ":  "
                    <<"Key data `" << level_name
                     << "' in ratio_to_coarser input has entries less than 1.");
      }
   }


   /*
    * If user did not supply explicit values for each level use the last one supplied for higher levels.
    */
   for(ln = tbox::MathUtilities<int>::Min(num_ratio_to_coarser_keys, ln_hi) + 1; ln <= ln_hi; ln++) {
      d_ratio_to_coarser[ln] = d_ratio_to_coarser[num_ratio_to_coarser_keys];
   }

   /*
    * Read input for largest patch sizes.
    */

   tbox::Pointer<tbox::Database> largest_patch_db =
      db->getDatabase("largest_patch_size");

   if (largest_patch_db.isNull()) {
       TBOX_ERROR(d_object_name << ":  "
                  << "Key data `largest_patch_size' not found in input.");
   }

   tbox::Array<std::string> largest_patch_keys = largest_patch_db->getAllKeys();
   int num_larg_patch_keys = largest_patch_keys.getSize();

   if ((num_larg_patch_keys < 1) && !is_from_restart) {
      TBOX_ERROR(d_object_name << ":  "
                 <<"No keys in `largest_patch_size' input section.");
   } else {
      d_largest_patch_size.resizeArray(num_larg_patch_keys);
      for (ln = 0; ln < num_larg_patch_keys; ln++) {
         int* temp_array = d_largest_patch_size[ln];

	 std::string level_name = "level_" + tbox::Utilities::intToString(ln);

	 if (!largest_patch_db->keyExists(level_name)) {
	    TBOX_ERROR(d_object_name << ":  "
		       <<"Key data `" << level_name
		       << "' not found in smallest_patch_size input.");
	 }

         largest_patch_db->getIntegerArray(largest_patch_keys[ln],
                                           temp_array, DIM);
      }
   }

   /*
    * Read input for smallest patch sizes.
    */

   if (db->keyExists("smallest_patch_size")) {
      tbox::Pointer<tbox::Database> smallest_patch_db =
         db->getDatabase("smallest_patch_size");

      tbox::Array<std::string> smallest_patch_keys = smallest_patch_db->getAllKeys();
      int num_smal_patch_keys = smallest_patch_keys.getSize();

      if (num_smal_patch_keys < 1) {
         TBOX_ERROR(d_object_name << ":  "
                    <<"No keys in `smallest_patch_size' input section.");
      } else {
         d_smallest_patch_size.resizeArray(num_smal_patch_keys);
         for (ln = 0; ln < num_smal_patch_keys; ln++) {
            int* temp_array = d_smallest_patch_size[ln];
	    std::string level_name = "level_" + tbox::Utilities::intToString(ln);

	    if (!smallest_patch_db->keyExists(level_name)) {
	       TBOX_ERROR(d_object_name << ":  "
			  <<"Key data `" << level_name
			  << "' not found in smallest_patch_size input.");
	    }

            smallest_patch_db->getIntegerArray(level_name, temp_array, DIM);
         }
      }
   }

   /*
    * Read input for proper nesting buffer and option of whether to grow
    * boxes after nesting is applied.
    */

   if (db->keyExists("proper_nesting_buffer")) {
      d_proper_nesting_buffer = db->getIntegerArray("proper_nesting_buffer");

      for (ln = 0; ln < d_proper_nesting_buffer.getSize(); ln++) {
         if (d_proper_nesting_buffer[ln] < 0) {
            TBOX_ERROR(d_object_name << ":  "
                       << "Key data `proper_nesting_buffer' has values < 0.");
         }
         if (d_proper_nesting_buffer[ln] == 0) {
            TBOX_WARNING(d_object_name << ":  "
                         << "Using zero `proper_nesting_buffer' values.");

         }
      }

   }

   if (db->keyExists("allow_patches_smaller_than_ghostwidth")) {
      d_allow_patches_smaller_than_ghostwidth =
         db->getBoolWithDefault("allow_patches_smaller_than_ghostwidth",
                                d_allow_patches_smaller_than_ghostwidth);

      if (d_allow_patches_smaller_than_ghostwidth) {
            TBOX_WARNING(d_object_name << ":  "
                         << "Allowing patches smaller than the max "
                         << "ghostwidth.  Note:  If periodic boundary "
                         << "conditions are used, this flag is ignored "
                         << "in the periodic directions.");

      }
   }

   if (db->keyExists("allow_patches_smaller_than_minimum_size_to_prevent_overlaps")) {
      d_allow_patches_smaller_than_minimum_size_to_prevent_overlaps =
         db->getBoolWithDefault(
            "allow_patches_smaller_than_minimum_size_to_prevent_overlaps",
            d_allow_patches_smaller_than_minimum_size_to_prevent_overlaps);

      if (d_allow_patches_smaller_than_minimum_size_to_prevent_overlaps) {
            TBOX_WARNING(d_object_name << ":  "
                         << "Allowing patches smaller than the given "
                         << "smallest patch size.  Note:  If periodic "
                         << "boundary conditions are used, this flag is "
                         << "ignored in the periodic directions.");
                                                                                
      }
   }



   /*
    * Read input for efficiency tolerance.
    */

   if (db->keyExists("efficiency_tolerance")) {
      d_efficiency_tolerance = db->getDoubleArray("efficiency_tolerance");

      for (ln = 0; ln < d_efficiency_tolerance.getSize(); ln++) {
         if ( (d_efficiency_tolerance[ln] <= 0.0e0)
             || (d_efficiency_tolerance[ln] >= 1.0e0) ) {
            TBOX_ERROR(d_object_name << ":  "
                       << "Key data `efficiency_tolerance' has values"
                       << " out of range 0.0 < tol < 1.0.");

         }
      }

   }

   /*
    * Read input for combine efficiency.
    */

   if (db->keyExists("combine_efficiency")) {
      d_combine_efficiency = db->getDoubleArray("combine_efficiency");

      for (ln = 0; ln < d_combine_efficiency.getSize(); ln++) {
         if ( (d_combine_efficiency[ln] <= 0.0e0)
             || (d_combine_efficiency[ln] >= 1.0e0) ) {
            TBOX_ERROR(d_object_name << ":  "
                       << "Key data `combine_efficiency' has values"
                       << " out of range 0.0 < tol < 1.0.");

         }
      }

   }

   d_max_levels = new_max_levels;

   if (d_max_levels > 1) {
      d_tag_init_strategy->checkCoarsenRatios(d_ratio_to_coarser);
   }

   d_check_nonrefined_tags =
      db->getCharWithDefault( "check_nonrefined_tags",
                              d_check_nonrefined_tags );
   if ( d_check_nonrefined_tags != 'i' &&
        d_check_nonrefined_tags != 'w' &&
        d_check_nonrefined_tags != 'e' ) {
      TBOX_ERROR("GriddingAlgorithm: input parameter check_nonrefined_tags\n"
                 <<"can only be i (ignore), w (warn) or e (error)");
   }

   d_check_overlapping_patches =
      db->getCharWithDefault( "check_overplapping_patches",
                              d_check_overlapping_patches );
   if ( d_check_overlapping_patches != 'i' &&
        d_check_overlapping_patches != 'w' &&
        d_check_overlapping_patches != 'e' ) {
      TBOX_ERROR("GriddingAlgorithm: input parameter check_overlapping_patches\n"
                 <<"can only be i (ignore), w (warn) or e (error)");
   }

   /*
    * Read input options for writing or reading refinement boxes.
    */
   if (db->keyExists("write_regrid_boxes")) {
      d_write_dumped_level_boxes = db->getBool("write_regrid_boxes");
   }
   if (db->keyExists("read_regrid_boxes")) {
      d_read_dumped_level_boxes = db->getBool("read_regrid_boxes");
   }
   if (d_write_dumped_level_boxes && d_read_dumped_level_boxes) {
      TBOX_ERROR("both 'write_regrid_boxes' and 'read_regrid_boxes' \n"
                 << "are specified in input... can only do one or the \n"
                 << "other.");
   }

   if (db->keyExists("regrid_boxes_filename")) {
      d_regrid_boxes_filename = db->getString("regrid_boxes_filename");
   } else {
      if (d_write_dumped_level_boxes || d_read_dumped_level_boxes) {
         TBOX_ERROR("In order to read or write refine boxes in the"
                    << "\nGriddingAlgorithm, you must specify the"
                    << "\n'regrid_boxes_filename' in the input file.");
      }

   }

   if (db->keyExists("barrier_before_clustering")) {
      d_barrier_before_clustering = db->getBool("barrier_before_clustering");
   }
   if (db->keyExists("sort_boxes_after_clustering")) {
      d_sort_boxes_after_clustering = db->getBool("sort_boxes_after_clustering");
   }
   if (db->keyExists("coalesce_boxes")) {
      d_coalesce_boxes = db->getBool("coalesce_boxes");
   }

   d_extend_tags_to_bdry = db->getBoolWithDefault("extend_tags_to_bdry",
                                                  d_extend_tags_to_bdry);

}

/*
*************************************************************************
*                                                                       *
* Gets the database in the root database that corresponds to the object *
* name.  This method then checks to make sure that the version number   *
* of the class is that same as the version number in the restart file.  *
* If these values are equal, the data members are read in from the      *
* restart database.                                                     *
*                                                                       *
*************************************************************************
*/

template<int DIM> void GriddingAlgorithm<DIM>::getFromRestart()
{
   tbox::Pointer<tbox::Database> root_db =
      tbox::RestartManager::getManager()->getRootDatabase();

   tbox::Pointer<tbox::Database> db;
   if ( root_db->isDatabase(d_object_name) ) {
      db = root_db->getDatabase(d_object_name);
   } else {
      TBOX_ERROR("Restart database corresponding to "
              << d_object_name << " not found in restart file.");
   }

   int ver = db->getInteger("ALGS_GRIDDING_ALGORITHM_VERSION");
   if (ver != ALGS_GRIDDING_ALGORITHM_VERSION) {
      TBOX_ERROR(d_object_name << ":  "
              << "Restart file version different than class version.");
   }

   int ln;

   d_max_levels = db->getInteger("d_max_levels");

   d_ratio_to_coarser.resizeArray(d_max_levels);

   tbox::Pointer<tbox::Database> ratio_to_coarser_db =
      db->getDatabase("d_ratio_to_coarser");

   for (ln = 1; ln < d_max_levels; ln++) {
      std::string level_name = "level_" + tbox::Utilities::intToString(ln);

      int* temp_ratio_to_coarser = d_ratio_to_coarser[ln];
      ratio_to_coarser_db->getIntegerArray(level_name,
                                           temp_ratio_to_coarser, DIM);
   }

   d_efficiency_tolerance = db->getDoubleArray("d_efficiency_tolerance");
   d_combine_efficiency = db->getDoubleArray("d_combine_efficiency");

   tbox::Pointer<tbox::Database> smallest_patch_db =
      db->getDatabase("d_smallest_patch_size");

   tbox::Array<std::string> smallest_patch_keys = smallest_patch_db->getAllKeys();
   d_smallest_patch_size.resizeArray(smallest_patch_keys.getSize());
   for (ln = 0; ln < smallest_patch_keys.getSize(); ln++) {
      int* tmp_array = d_smallest_patch_size[ln];
      smallest_patch_db->getIntegerArray(smallest_patch_keys[ln],
                                         tmp_array, DIM);
   }

   tbox::Pointer<tbox::Database> largest_patch_db =
      db->getDatabase("d_largest_patch_size");

   tbox::Array<std::string> largest_patch_keys = largest_patch_db->getAllKeys();
   d_largest_patch_size.resizeArray(largest_patch_keys.getSize());
   for (ln = 0; ln < largest_patch_keys.getSize(); ln++) {
      int* tmp_array = d_largest_patch_size[ln];
      largest_patch_db->getIntegerArray(largest_patch_keys[ln],
                                        tmp_array, DIM);
   }

   d_proper_nesting_buffer = db->getIntegerArray("d_proper_nesting_buffer");
}



/*
 * ************************************************************************
 *
 *  for use when sorting integers using the C-library qsort
 *
 * ************************************************************************
 */
template<int DIM> int GriddingAlgorithm<DIM>::qsortBoxCompare(const void *v, const void *w)
{
   const hier::IntVector<DIM> &lowv = ((const hier::Box<DIM> *)v)->lower();
   const hier::IntVector<DIM> &loww = ((const hier::Box<DIM> *)w)->lower();
   for ( int i=0; i<DIM; ++i ) {
      if ( lowv[i] > loww[i] ) return  1;
      if ( lowv[i] < loww[i] ) return -1;
   }
   return (0);
}



}
}

#endif
