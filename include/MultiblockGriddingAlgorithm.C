//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-0/source/mesh/multiblock/MultiblockGriddingAlgorithm.C $
// Package:     SAMRAI multiblock package
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2292 $
// Modified:    $LastChangedDate: 2008-07-11 11:31:57 -0700 (Fri, 11 Jul 2008) $
// Description: AMR hierarchy generation and regridding routines.
//

#ifndef included_mesh_MultiblockGriddingAlgorithm_C
#define included_mesh_MultiblockGriddingAlgorithm_C

#include "MultiblockGriddingAlgorithm.h"

#include "BoxUtilities.h"
#include "CellData.h"
#include "tbox/RestartManager.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"

#include <stdio.h>
#include <stdlib.h>
#include <fstream>

#define ALGS_GRIDDING_ALGORITHM_VERSION (2)

#ifndef NULL
#define NULL (0)
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

template<int DIM> int MultiblockGriddingAlgorithm<DIM>::s_instance_counter = 0;
template<int DIM> int MultiblockGriddingAlgorithm<DIM>::s_tag_indx = -1;
template<int DIM> int MultiblockGriddingAlgorithm<DIM>::s_buf_tag_indx = -1;
template<int DIM> int MultiblockGriddingAlgorithm<DIM>::s_buf_tag_src_indx = -1;

/*
*************************************************************************
*                                                                       *
* Constructor and destructor for MultiblockGriddingAlgorithm<DIM>.      *
*                                                                       *
*************************************************************************
*/
template<int DIM>
MultiblockGriddingAlgorithm<DIM>::MultiblockGriddingAlgorithm(
   const std::string& object_name,
   tbox::Pointer<tbox::Database> input_db,
   tbox::Pointer< hier::MultiblockPatchHierarchy<DIM> > multiblock,
   tbox::Pointer< mesh::TagAndInitializeStrategy<DIM> > tag_init_strategy,
   tbox::Pointer< mesh::BoxGeneratorStrategy<DIM> > generator,
   tbox::Pointer< mesh::LoadBalanceStrategy<DIM> > balancer,
   MultiblockGriddingTagger<DIM>* mb_tagger_strategy,
   bool register_for_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(!input_db.isNull());
    TBOX_ASSERT(!multiblock.isNull());
    TBOX_ASSERT(!tag_init_strategy.isNull());
    TBOX_ASSERT(!generator.isNull());
    TBOX_ASSERT(!balancer.isNull());
#endif

   d_object_name = object_name;
   d_registered_for_restart = register_for_restart;

   if (d_registered_for_restart) {
      tbox::RestartManager::getManager()->
         registerRestartItem(d_object_name, this);
   }

   d_tag_init_strategy = tag_init_strategy;
   d_box_generator = generator;
   d_load_balancer = balancer;
   d_mb_tagger_strategy = mb_tagger_strategy;

   if (mb_tagger_strategy) {
      d_mb_tagger_strategy = mb_tagger_strategy;
      d_internal_tagger_strategy = false;
   } else {
      d_mb_tagger_strategy = new MultiblockGriddingTagger<DIM>(); 
      d_internal_tagger_strategy = true;
   } 

   /*
    * Construct integer tag variables and add to variable database.  Note that 
    * variables and patch data indices are shared among all MultiblockGriddingAlgorithm 
    * instances.  The VariableDatabase holds the variables, once contructed and registered
    * via the VariableDatabase::registerInternalSAMRAIVariable()
    * function call.  Note that variables are registered and patch data indices
    * are made only for the first time through the constructor.
    */
   /*
    * Construct integer tag variables and add to database
    */
   hier::VariableDatabase<DIM>* var_db =
      hier::VariableDatabase<DIM>::getDatabase();

   static std::string tag_interior_variable_name("MultiblockGriddingAlgorithm__tag-interior");
   static std::string tag_buffer_variable_name("MultiblockGriddingAlgorithm__tag-buffer");
   static std::string tag_src_variable_name("MultiblockGriddingAlgorithm__tag-source");

   d_tag = var_db->getVariable(tag_interior_variable_name);
   if (d_tag.isNull()) {
      d_tag = new pdat::CellVariable<DIM,int>(tag_interior_variable_name, 1);
   }

   d_buf_tag = var_db->getVariable(tag_buffer_variable_name);
   if (d_buf_tag.isNull()) {
      d_buf_tag = new pdat::CellVariable<DIM,int>(tag_buffer_variable_name, 1);
   }

   d_buf_src_tag = var_db->getVariable(tag_src_variable_name);
   if (d_buf_src_tag.isNull()) {
      d_buf_src_tag = new pdat::CellVariable<DIM,int>(tag_src_variable_name, 1);
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
   if (s_buf_tag_src_indx < 0) {
      s_buf_tag_src_indx =
         var_db->registerInternalSAMRAIVariable(d_buf_src_tag,
                                                hier::IntVector<DIM>(0));
   }

   d_tag_indx = s_tag_indx;
   d_buf_tag_indx = s_buf_tag_indx;
   d_buf_tag_src_indx = s_buf_tag_src_indx;

   d_mb_tagger_strategy->setScratchTagPatchDataIndex(d_buf_tag_indx);

   /*
    * Tag value for refined cells is one; others are zero.
    */
   d_true_tag = 1;
   d_false_tag = 0;

   /*
    * Initialize communication algorithm for exchanging buffered tag data.
    */

   tbox::Pointer< xfer::RefineAlgorithm<DIM> > refine_alg = new xfer::RefineAlgorithm<DIM>();

   d_bdry_fill_tags = new xfer::MultiblockRefineAlgorithm<DIM>(refine_alg,
                                                          multiblock);

   d_bdry_fill_tags->
      registerRefine(d_buf_tag_indx,
                     d_buf_tag_src_indx,
                     d_buf_tag_indx,
                     ((tbox::Pointer< xfer::RefineOperator<DIM> >)NULL));

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

   /*
    * calculus and other regridding operations.
    */
   t_find_proper_nesting = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm<DIM>::findProperNestingBoxes()");
   t_intersect_boxes_find_refinement = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm<DIM>::intersect_boxes_find_refinement");
   t_load_balance_boxes = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm<DIM>::load_balance_boxes");
   t_bdry_fill_tags_create = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm<DIM>::bdry_fill_tags_create");

   t_make_coarsest = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm<DIM>::makeCoarsestLevel()");

   t_regrid_all_finer = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm<DIM>::regridAllFinerLevels()");
   t_remove_intersections_regrid_all =
      tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm<DIM>::remove_intersections_regrid_all");

   t_make_finer = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm<DIM>::makeFinerLevel()");
   t_remove_intersections_make_finer = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm<DIM>::remove_intersections_make_finer");

   t_set_tags = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm<DIM>::setTagsOnLevel()");

   t_buffer_tags = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm<DIM>::bufferTagsOnLevel()");
   t_bdry_fill_tags_comm = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm<DIM>::bdry_fill_tags_comm");

   t_find_refinement = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm<DIM>::findRefinementBoxes()");
   t_find_boxes_containing_tags = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm<DIM>::find_boxes_containing_tags");
   t_remove_intersections_find_proper = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm<DIM>::remove_intersections_find_proper");
   t_intersect_boxes_find_proper = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm<DIM>::intersect_boxes_find_proper");

   /*
    * Initialize object with data read from input and restart databases.
    */
   bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();
   if ( is_from_restart ) {
      getFromRestart();
   }

   getFromInput(input_db, is_from_restart);

   d_bdry_sched_tags.resizeArray(d_max_levels);

   d_proper_nesting_boxes.resizeArray(multiblock->getNumberOfBlocks());
   for (int nb = 0; nb < multiblock->getNumberOfBlocks(); nb++) {
      d_proper_nesting_boxes[nb].resizeArray(d_max_levels);
   }

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
* Destructor tells the RestartManager to remove this object from   *
* the list of restart items.                                            *
*                                                                       *
*************************************************************************
*/
template<int DIM>
MultiblockGriddingAlgorithm<DIM>::~MultiblockGriddingAlgorithm()
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

   if (d_internal_tagger_strategy && d_mb_tagger_strategy) {
      delete d_mb_tagger_strategy;
   }

   if (s_instance_counter == 0) {
      hier::VariableDatabase<DIM>::getDatabase()->
         removeInternalSAMRAIVariablePatchDataIndex(s_tag_indx);
      s_tag_indx = -1;
      hier::VariableDatabase<DIM>::getDatabase()->
         removeInternalSAMRAIVariablePatchDataIndex(s_buf_tag_indx);
      s_buf_tag_indx = -1;
      hier::VariableDatabase<DIM>::getDatabase()->
         removeInternalSAMRAIVariablePatchDataIndex(s_buf_tag_src_indx);
      s_buf_tag_src_indx = -1;
   }

}

template<int DIM>
void MultiblockGriddingAlgorithm<DIM>::makeCoarsestLevel(
   tbox::Pointer< hier::BasePatchHierarchy<DIM> > hierarchy,
   const double level_time,
   const hier::BoxArray<DIM>& override_boxes,
   const hier::ProcessorMapping& override_mapping)
{
   NULL_USE(override_boxes);
   NULL_USE(override_mapping);

   tbox::Pointer< hier::MultiblockPatchHierarchy<DIM> > multiblock = hierarchy;

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!(multiblock.isNull()));
   TBOX_ASSERT(d_max_levels > 0);
#endif

   /*
    * Start timer for this method.
    */
   t_make_coarsest->start();

   const int level_number = 0;
   const int nblocks = multiblock->getNumberOfBlocks();

   bool level_zero_exists = false;
   for (int b = 0; b < nblocks; b++) {
      tbox::Pointer< hier::PatchHierarchy<DIM> > block_hierarchy =
         multiblock->getHierarchy(b);
      if ( (block_hierarchy->getNumberOfLevels() > 0) ) {
         if ( !(block_hierarchy->getPatchLevel(level_number).isNull()) ) {
            level_zero_exists = true;
         }
      }
   }

   tbox::Array< tbox::Pointer< hier::PatchLevel<DIM> > > old_level(nblocks);
 
   for (int i = 0; i < nblocks; i++) {
 
      tbox::Pointer< hier::PatchHierarchy<DIM> > block_hierarchy =
         multiblock->getHierarchy(i);

      hier::BoxArray<DIM> domain = block_hierarchy->getGridGeometry()->
                                                  getPhysicalDomain();
      /*
       * Read/write coarse level boxes from/to file.
       */
      if (d_read_dumped_level_boxes) {
         d_regrid_box_utility->getLevelBoxes(
                                  domain,
                                  level_number,
                                  d_regrid_box_counter[i][level_number]);
         d_regrid_box_counter[i][level_number]++;
      } else {
         if (d_write_dumped_level_boxes) {
            d_regrid_box_utility->putLevelBoxes(
                                     domain,
                                     level_number,
                                     d_regrid_box_counter[i][level_number]);
            d_regrid_box_counter[i][level_number]++;
         }
      }

      hier::BoxList<DIM> domain_list(domain);
      hier::IntVector<DIM> smallest_patch;
      hier::IntVector<DIM> smallest_box_to_refine;
      hier::IntVector<DIM> largest_patch;
      hier::IntVector<DIM> extend_ghosts;
      getGriddingParameters(smallest_patch,
                            smallest_box_to_refine,
                            largest_patch,
                            extend_ghosts,
                            block_hierarchy, level_number);

      /*
       * If there is no level 0 in the patch hierarchy, then check
       * constraints on domain description.
       */

      if ( !level_zero_exists ) {

         for (typename hier::BoxList<DIM>::Iterator b(domain_list); b; b++) {
            hier::Box<DIM> test_box = b();
            for (int dir = 0; dir < DIM; dir++) {
               if (test_box.numberCells(dir) < smallest_patch(dir)) {
                  TBOX_ERROR(d_object_name << ":  "
                          << "A Box from input file is smaller than max ghost "
                          << "width or minimum patch size");
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

      const hier::IntVector<DIM> ratio(1);
      hier::BoxArray<DIM> level_boxes;
      hier::ProcessorMapping mapping;

      t_load_balance_boxes->start();
      d_load_balancer->loadBalanceBoxes(level_boxes, mapping,
                                        domain_list,
                                        block_hierarchy, level_number, domain,
                                        hier::IntVector<DIM>(1), 
                                        smallest_patch,
                                        largest_patch,
                                        patch_cut_factor, extend_ghosts);
      t_load_balance_boxes->stop();

      if ( !level_zero_exists ) {

         block_hierarchy->makeNewPatchLevel(level_number,
                                      ratio,
                                      level_boxes,
                                      mapping);

      } else {

         if (block_hierarchy->levelExists(level_number)) {
            
            old_level[i] =
               block_hierarchy->getPatchLevel(level_number);
            
            block_hierarchy->removePatchLevel(level_number);
            
            block_hierarchy->makeNewPatchLevel(level_number,
                                               ratio,
                                               level_boxes,
                                               mapping);
         }

      }


   }

   bool initial_time;
   tbox::Pointer< hier::MultiblockPatchLevel<DIM> > old_mb_level;
   if (!level_zero_exists) {
      initial_time = true;
      old_mb_level = (tbox::Pointer< hier::MultiblockPatchLevel<DIM> >) NULL;
   } else {
      initial_time = false;
      old_mb_level = new hier::MultiblockPatchLevel<DIM>(old_level);
   }

   d_tag_init_strategy->initializeLevelData(multiblock, level_number,
                                            level_time,
                                            levelCanBeRefined(level_number),
                                            initial_time,
                                            old_mb_level);

   if (level_zero_exists) {
      old_mb_level.setNull();
   }

   d_tag_init_strategy->resetHierarchyConfiguration(multiblock,
                                                 level_number, level_number);

   multiblock->adjustMultiblockPatchLevelBoundaries(
      multiblock->getPatchLevel(0)); 


   t_make_coarsest->stop();

}

/*
*************************************************************************
*                                                                       *
* Write out class version number and data members to database.          *
*                                                                       *
*************************************************************************
*/

template<int DIM>
void MultiblockGriddingAlgorithm<DIM>::putToDatabase(
   tbox::Pointer<tbox::Database> db)
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

template<int DIM>
int MultiblockGriddingAlgorithm<DIM>::getMaxLevels() const
{
   return(d_max_levels);
}

template<int DIM>
bool MultiblockGriddingAlgorithm<DIM>::levelCanBeRefined(
   const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(level_number >= 0);
#endif
   return(level_number < d_max_levels-1);
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

template<int DIM>
void MultiblockGriddingAlgorithm<DIM>::getGriddingParameters(
   hier::IntVector<DIM>& smallest_patch,
   hier::IntVector<DIM>& smallest_box_to_refine,
   hier::IntVector<DIM>& largest_patch,
   hier::IntVector<DIM>& extend_ghosts,
   const tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
   const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!(hierarchy.isNull()));
   TBOX_ASSERT( (level_number >= 0) && (level_number < d_max_levels) );
#endif

   int fine_level_number = level_number + 1;

   /*
    * Determine maximum ghost cell width needed over all variables
    * currently known to the patch descriptor.
    */

   hier::IntVector<DIM> max_ghosts =
      hierarchy->getPatchDescriptor()->getMaxGhostWidth();

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

   largest_patch = getLargestPatchSize(level_number);
   largest_patch.max(smallest_patch);

   /*
    * Set number of cells that coarsened patches must accomodate
    * to meet ghost cell needs of variables.
    */

   if (level_number < d_max_levels-1) {

      for (int i = 0; i < DIM; i++) {
         int den = d_ratio_to_coarser[fine_level_number](i)
                   / d_tag_init_strategy->getErrorCoarsenRatio();
         int sz = max_ghosts(i)/den;
         if ( max_ghosts(i) - sz * den ) sz++;

         smallest_box_to_refine(i) = 
            tbox::MathUtilities<int>::Max( sz, 
               smallest_patch(i)/d_ratio_to_coarser[fine_level_number](i) );
      }

   } else {
      smallest_box_to_refine = smallest_patch;
   }

   /*
    * Determine number of cells box may be extended to physical
    * domain boundary to accomodate ghost cells.
    */
   extend_ghosts = max_ghosts
                   * d_tag_init_strategy->getErrorCoarsenRatio();

}

template<int DIM>
const hier::IntVector<DIM>&
MultiblockGriddingAlgorithm<DIM>::getSmallestPatchSize(
   const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( (level_number >= 0) && (level_number < d_max_levels) );
#endif
   int size = d_smallest_patch_size.getSize();
   return( (level_number < size)
           ? d_smallest_patch_size[level_number]
           : d_smallest_patch_size[size-1] );
}

template<int DIM>
const hier::IntVector<DIM>&
MultiblockGriddingAlgorithm<DIM>::getLargestPatchSize(
   const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( (level_number >= 0) && (level_number < d_max_levels) );
#endif
   int size = d_largest_patch_size.getSize();
   return( (level_number < size)
           ? d_largest_patch_size[level_number]
           : d_largest_patch_size[size-1] );
}

template<int DIM>
void MultiblockGriddingAlgorithm<DIM>::getFromInput(
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

   for (ln = ln_lo; ln <= ln_hi; ln++) {
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
                     << "' in ratio_to_coarser input has entries less than 1.");      }

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
            smallest_patch_db->getIntegerArray(smallest_patch_keys[ln],
                                               temp_array, DIM);
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
                    << "\nGriddingAlgorithm<DIM>, you must specify the"
                    << "\n'regrid_boxes_filename' in the input file.");
      }

   }

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

template<int DIM>
void MultiblockGriddingAlgorithm<DIM>::regridAllFinerLevels(
   tbox::Pointer< hier::BasePatchHierarchy<DIM> > hierarchy,
   const int level_number,
   const double regrid_time,
   const tbox::Array<int>& tag_buffer,
   tbox::Array<double> regrid_start_time,
   const bool level_is_coarsest_sync_level)
{
   tbox::Pointer< hier::MultiblockPatchHierarchy<DIM> > multiblock = hierarchy;

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!(multiblock.isNull()));
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

      int nblocks = multiblock->getNumberOfBlocks();
      for (int n = 0; n < nblocks; n++) {
         tbox::Pointer< hier::PatchHierarchy<DIM> > this_hierarchy =
            multiblock->getHierarchy(n);

         if (this_hierarchy->levelExists(level_number)) {
            const tbox::Pointer< hier::PatchLevel<DIM> > patch_level =
               this_hierarchy->getPatchLevel(level_number);

            if (!patch_level.isNull()) {
               hier::BoxList<DIM> complement(patch_level->getPhysicalDomain());
               hier::BoxList<DIM> level_boxes(patch_level->getBoxes());

               t_remove_intersections_regrid_all->start();
               complement.removeIntersections(level_boxes);
               t_remove_intersections_regrid_all->stop();

               t_find_proper_nesting->start();
               findProperNestingBoxes(this_hierarchy, level_number, n, complement);
               t_find_proper_nesting->stop();

               /*
                * Perform pre-processing of error estimation data, if
                * appropriate.
                */
               if (errorEstimationUsesTimeIntegration() ) {
                  for (int ln = level_number;
                       ln <= this_hierarchy->getFinestLevelNumber(); ln++) {
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
                           preprocessErrorEstimation(this_hierarchy,
                                                     ln,
                                                     regrid_time,
                                                     level_regrid_start_time,
                                                     initial_time);
                     }
                  }
               }
            } // !level.isNull()
         } // level exists
      } // loop over blocks
      

      /*
       * Recursively regrid each finer level.
       */
      const int finest_level_not_regridded = level_number;
      regridFinerLevel(multiblock,
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

      if ( multiblock->getFinestLevelNumber() >= (level_number+1) ) {
         d_tag_init_strategy->
            resetHierarchyConfiguration(multiblock,
                                        level_number+1,
                                        multiblock->
                                           getFinestLevelNumber());

      }

   } //  if level cannot be refined, the routine drops through...

   t_regrid_all_finer->stop();
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

template<int DIM>
void MultiblockGriddingAlgorithm<DIM>::makeFinerLevel(
   tbox::Pointer< hier::BasePatchHierarchy<DIM> > hierarchy,
   const double level_time,
   const bool initial_time,
   const int tag_buffer,
   const double regrid_start_time)
{
   tbox::Pointer< hier::MultiblockPatchHierarchy<DIM> > multiblock = hierarchy;

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!(multiblock.isNull()));
   TBOX_ASSERT(!(multiblock->getPatchLevel(
                       multiblock->getFinestLevelNumber()).isNull()));
   TBOX_ASSERT(tag_buffer >= 0);
#endif

   t_make_finer->start();

   const int level_number = multiblock->getFinestLevelNumber();
   const int fine_level_number = level_number + 1;
   const int nblocks = multiblock->getNumberOfBlocks();

   if ( levelCanBeRefined(level_number) ) {

      tbox::Pointer< hier::MultiblockPatchLevel<DIM> > mb_level =
         multiblock->getPatchLevel(level_number);

      tbox::Array<hier::BoxArray<DIM> > fine_boxes(nblocks);
      hier::ProcessorMapping* mapping;
      mapping = new hier::ProcessorMapping[nblocks];

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
          d_tag_init_strategy->refineUserBoxInputOnly()) {
         do_tagging = false;
      }

      /*
       * Tag cells, determine refine boxes from tagged regions, and
       * load balance in preparation for constructing new refined level.
       */
      if (do_tagging) {

         /*
          * Determine proper nesting box array for specified level.  If the
          * level is the coarsest, the proper nesting box array covers same
          * region as the physical domain.  Otherwise, the nesting boxes are
          * computed from a proper interior of the level.
          */
         
         for (int nb = 0; nb < nblocks; nb++) {
            tbox::Pointer< hier::PatchHierarchy<DIM> > block_hierarchy =
               multiblock->getHierarchy(nb);

            tbox::Pointer< hier::PatchLevel<DIM> > patch_level =
               mb_level->getPatchLevelForBlock(nb); 

            if (!patch_level.isNull()) {
 
               if ( level_number == 0 ) {
                  d_proper_nesting_boxes[nb][0] =
                     hier::BoxList<DIM>(patch_level->getBoxes());

               } else {

                  hier::BoxList<DIM> complement(
                     patch_level->getPhysicalDomain());
                  hier::BoxList<DIM> level_boxes(patch_level->getBoxes());
                  t_remove_intersections_make_finer->start();
                  complement.removeIntersections(level_boxes);
                  t_remove_intersections_make_finer->stop();

                  t_find_proper_nesting->start();
                  findProperNestingBoxes(block_hierarchy,
                                         level_number,
                                         nb,
                                         complement);
                  t_find_proper_nesting->stop();
               }
            }
         }

         t_bdry_fill_tags_create->start();
         d_bdry_sched_tags[level_number] =
            d_bdry_fill_tags->createSchedule(mb_level, d_mb_tagger_strategy);
         t_bdry_fill_tags_create->stop();

         /*
          * Initialize integer tag arrays on level to false.
          */

         mb_level->allocatePatchData(d_tag_indx);
         for (int st = 0; st < nblocks; st++) {
            tbox::Pointer< hier::PatchLevel<DIM> > patch_level = 
               mb_level->getPatchLevelForBlock(st); 
            if (!patch_level.isNull()) {
               setTagsOnLevel(d_false_tag,
                              patch_level,
                              d_tag_indx,
                              patch_level->getBoxes());
            }
         }

         /*
          * Perform pre-processing of error estimation data, if appropriate.
          */
         if (errorEstimationUsesTimeIntegration() ) {
            d_tag_init_strategy->
                preprocessErrorEstimation(multiblock,
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
         d_tag_init_strategy->
            tagCellsForRefinement(multiblock,
                                  level_number,
                                  level_time,
                                  d_tag_indx,
                                  initial_time,
                                  coarsest_sync_level,
                                  levelCanBeRefined(level_number),
                                  regrid_start_time);

         /*
          * Buffer true tagged cells by specified amount which should be
          * sufficient to keep disturbance on refined region until next
          * regrid of the level occurs.
          */
         bufferTagsOnLevel(d_true_tag, mb_level, tag_buffer);

         /*
          * Determine box array and processor mapping for new fine level.
          */
         for (int fr = 0; fr < nblocks; fr++) {
            if (multiblock->getHierarchy(fr)->getFinestLevelNumber() >=
                level_number) {
               findRefinementBoxes(fine_boxes[fr],
                                   mapping[fr],
                                   multiblock->getHierarchy(fr),
                                   level_number,
                                   d_proper_nesting_boxes[fr][level_number]);
            }
         }

         /*
          * Deallocate tag arrays and schedule -- no longer needed.
          */
         mb_level->deallocatePatchData(d_tag_indx);
         d_bdry_sched_tags[level_number].setNull();

         /*
          * Write refine boxes to file if requested.
          */
         if (d_write_dumped_level_boxes) {
            for (int pl = 0; pl < nblocks; pl++) {
               d_regrid_box_utility->putLevelBoxes(fine_boxes[pl],
                                      level_number+1,
                                      d_regrid_box_counter[pl][level_number+1]);
               d_regrid_box_counter[pl][level_number+1]++;
            } 
         }
     
      } else {

         for (int nt = 0; nt < nblocks; nt++) {

            tbox::Pointer< hier::PatchHierarchy<DIM> > block_hierarchy =
               multiblock->getHierarchy(nt);

            /*
             * If tagging is not necessary (do_tagging = false) we simply
             * need to access the level boxes, either from a dumped file or
             * from user-supplied refine boxes, and load balance them before
             * constructing the finer level.
             */
            hier::BoxArray<DIM> boxes_to_refine;

            if (d_read_dumped_level_boxes) {

               d_regrid_box_utility->getLevelBoxes(boxes_to_refine,
                  fine_level_number,
                  d_regrid_box_counter[nt][fine_level_number]);
               d_regrid_box_counter[nt][fine_level_number]++;

            }

            /*
             * Access the user supplied refine boxes.
             */
            if (d_tag_init_strategy->refineUserBoxInputOnly()) {

               (void) d_tag_init_strategy->
                  getUserSuppliedRefineBoxes(boxes_to_refine,
                                             level_number,
                                             level_time);

               boxes_to_refine.refine(d_ratio_to_coarser[fine_level_number]);

            }

            /*
             * Avoid further work if pre-defined refine boxes are empty.  In
             * this case, the "fine_boxes" boxarray will remain empty and
             * no further levels will be constructed.
             */
            if (!(boxes_to_refine.getNumberOfBoxes() == 0)) {

               hier::BoxList<DIM> fine_box_list(boxes_to_refine);

               hier::IntVector<DIM> ratio_to_level_zero =
                  mb_level->getRatio() *
                     d_ratio_to_coarser[fine_level_number];

               hier::BoxArray<DIM> physical_domain;
               block_hierarchy->getGridGeometry()->
                  computePhysicalDomain(physical_domain, ratio_to_level_zero);

               hier::IntVector<DIM>
                  patch_cut_factor(d_ratio_to_coarser[fine_level_number]);

               hier::IntVector<DIM> smallest_patch;
               hier::IntVector<DIM> smallest_box_to_refine;
               hier::IntVector<DIM> largest_patch;
               hier::IntVector<DIM> extend_ghosts;
               getGriddingParameters(smallest_patch,
                                     smallest_box_to_refine,
                                     largest_patch,
                                     extend_ghosts,
                                     block_hierarchy, fine_level_number);

               d_load_balancer->loadBalanceBoxes(
                  fine_boxes[nt], mapping[nt],
                  fine_box_list, block_hierarchy, fine_level_number,
                  physical_domain,
                  ratio_to_level_zero,
                  smallest_patch,
                  getLargestPatchSize(fine_level_number),
                  patch_cut_factor, extend_ghosts);

            }
         }
      }

      /*
       * Make new finer level (level_number+1), if appropriate.
       */
      bool new_level_created = false;
      for (int fl = 0; fl < nblocks; fl++) {
         tbox::Pointer< hier::PatchHierarchy<DIM> > block_hierarchy =
            multiblock->getHierarchy(fl);
         tbox::Pointer< hier::PatchLevel<DIM> > patch_level =
            mb_level->getPatchLevelForBlock(fl);
         if (!patch_level.isNull()) {
            if ( !(fine_boxes[fl].getNumberOfBoxes() == 0) ) {
               hier::IntVector<DIM> ratio =
                  patch_level->getRatio() * getRatioToCoarserLevel(level_number+1);
               block_hierarchy->makeNewPatchLevel(level_number+1,
                                                  ratio, fine_boxes[fl],
                                                  mapping[fl]);
               new_level_created = true;
            }
         }
      }
      if (new_level_created) {
         d_tag_init_strategy->initializeLevelData(multiblock,
                                                  level_number+1,
                                                  level_time,
                                                  levelCanBeRefined(
                                                     level_number+1),
                                                  initial_time);

         d_tag_init_strategy->resetHierarchyConfiguration(multiblock,
                                                          level_number+1,
                                                          level_number+1);

      }

      multiblock->adjustMultiblockPatchLevelBoundaries(
         multiblock->getMultiblockPatchLevel(level_number+1));

      delete[] mapping;
   }

   t_make_finer->stop();

}

/*
*************************************************************************
*                                                                       *
* Set each integer value in specified tag array to tag_value where      *
* patch level intersects given box array.                               *
*                                                                       *
*************************************************************************
*/

template<int DIM>
void MultiblockGriddingAlgorithm<DIM>::setTagsOnLevel(
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

template<int DIM>
bool MultiblockGriddingAlgorithm<DIM>::errorEstimationUsesTimeIntegration()
const
{
   return(d_tag_init_strategy->usesTimeIntegration());
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

template<int DIM>
void MultiblockGriddingAlgorithm<DIM>::bufferTagsOnLevel(
   const int tag_value,
   const tbox::Pointer< hier::MultiblockPatchLevel<DIM> > level,
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

   level->allocatePatchData(d_buf_tag_indx);
   level->allocatePatchData(d_buf_tag_src_indx);

   /*
    * Set temporary buffered tags based on buffer width and
    * distance from actual tags.
    */
   const int not_tag = ((tag_value == d_true_tag) ? d_false_tag : d_true_tag);
   for (int nb = 0; nb < level->getNumberOfBlocks(); nb++) {
      tbox::Pointer< hier::PatchLevel<DIM> > block_level =
         level->getPatchLevelForBlock(nb);

      if (!block_level.isNull()) {
         for (typename hier::PatchLevel<DIM>::Iterator ip1(block_level); ip1; ip1++) {
            tbox::Pointer< hier::Patch<DIM> > patch = block_level->getPatch(ip1());

            tbox::Pointer< pdat::CellData<DIM,int> >
               buf_tag_src_data = patch->getPatchData(d_buf_tag_src_indx);
            tbox::Pointer< pdat::CellData<DIM,int> >
               buf_tag_data = patch->getPatchData(d_buf_tag_indx);
            tbox::Pointer< pdat::CellData<DIM,int> >
               tag_data = patch->getPatchData(d_tag_indx);

            buf_tag_data->fillAll(not_tag);
            buf_tag_src_data->fillAll(not_tag);

            hier::Box<DIM> interior = patch->getBox();

            for (int bc = buffer_size; bc >= 0; bc--) {

               int fill_val = buffer_size - bc + 1;

               for (typename pdat::CellData<DIM,int>::Iterator ic(interior); ic; ic++) {
                  if ( (*tag_data)(ic()) == tag_value ) {
                     hier::Box<DIM> buf_box(ic()-bc, ic()+bc);
                     buf_tag_src_data->fill(fill_val, buf_box * interior);
                  }
               }

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
   for (int bc = 0; bc < level->getNumberOfBlocks(); bc++) {
      tbox::Pointer< hier::PatchLevel<DIM> > block_level =
         level->getPatchLevelForBlock(bc);
      if (!block_level.isNull()) {
         for (typename hier::PatchLevel<DIM>::Iterator ip2(block_level); ip2; ip2++) {
            tbox::Pointer< hier::Patch<DIM> > patch = block_level->getPatch(ip2());

            tbox::Pointer< pdat::CellData<DIM,int> > buf_tag_data =
                                       patch->getPatchData(d_buf_tag_indx);
            tbox::Pointer< pdat::CellData<DIM,int> > tag_data =
                                       patch->getPatchData(d_tag_indx);

            hier::Box<DIM> buf_tag_box = buf_tag_data->getGhostBox();
            hier::Box<DIM> tag_box = tag_data->getBox();

            /*
             * Set all interior tags to tag value where buffer tags non-zero.
             */
            for (typename pdat::CellData<DIM,int>::Iterator ic(tag_box); ic; ic++) {
               (*tag_data)(ic()) = ( (*buf_tag_data)(ic()) ? tag_value : not_tag );
            }

            /*
             * Set all interior tags in buffers around tags in ghosts.
             */
            for (typename pdat::CellData<DIM,int>::Iterator ic2(buf_tag_box); ic2; ic2++) {
               int tval = (*buf_tag_data)(ic2());
               if ( tval > 1 ) {
                  int buf_size = tval - 1;
                  hier::Box<DIM> buf_box(ic2()-buf_size, ic2()+buf_size);
                  tag_data->fill(tag_value, buf_box);
               }
            }
         }
      }
   }

   level->deallocatePatchData(d_buf_tag_indx);
   level->deallocatePatchData(d_buf_tag_src_indx);

   t_buffer_tags->stop();

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

template<int DIM>
void MultiblockGriddingAlgorithm<DIM>::findRefinementBoxes(
   hier::BoxArray<DIM>& fine_boxes,
   hier::ProcessorMapping& mapping,
   const tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
   const int coarse_level_number,
   const hier::BoxList<DIM>& nesting_boxes) const
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
    * Determine single smallest bounding box for all nesting boxes.
    */
   typename hier::BoxList<DIM>::Iterator lb(nesting_boxes);
   hier::Box<DIM> bounding_box(lb());
   for (; lb; lb++) {
      bounding_box += lb();
   }

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
   getGriddingParameters(smallest_patch,
                         smallest_box_to_refine,
                         largest_patch,
                         extend_ghosts,
                         hierarchy, coarse_level_number);

   hier::BoxList<DIM> box_list;
   tbox::Pointer< hier::PatchLevel<DIM> > level =
      hierarchy->getPatchLevel(coarse_level_number);

   t_find_boxes_containing_tags->start();
   d_box_generator->findBoxesContainingTags(
         box_list, level, d_tag_indx, d_true_tag, bounding_box,
         smallest_box_to_refine,
         getEfficiencyTolerance(coarse_level_number),
         getCombineEfficiency(coarse_level_number));
   t_find_boxes_containing_tags->stop();

   if (!box_list.isEmpty()) {

      t_intersect_boxes_find_refinement->start();
      box_list.intersectBoxes(nesting_boxes);
      t_intersect_boxes_find_refinement->stop();

      /*
       * Before constructing a new patch level for the hierarchy, the boxes
       * that cover the tags are massaged further.  An attempt is made to
       * produce a collection of boxes, from which the new patch level will
       * be constructed, that obeys several constraints.  In particular, the
       * boxes should not overlap, and each box when grown by the appropriate
       * ghost cell width should not intersect the physical domain boundary
       * improperly.  Also, each box should be within a certain size range.
       * The box generator will not produce boxes smaller than the minimum
       * size.  However, intersecting against the nesting boxes may.  Thus,
       * we may need to grow some on these boxes within the nesting domain,
       * which may introduce overlaps.
       *
       * If overlaps are to be avoided, at the possible expense of constructing
       * boxes smaller than the minimum patch size, set the
       * "grow_after_nesting" argument to false in input.  By default, this
       * argument is true.  If set false, the boxes will not be grown within
       * the nesting domain and no overlaps will be introduced.
       */

      hier::BoxList<DIM> level_domain(level->getPhysicalDomain());
      for (int i = 0; i < DIM; i++) {
         int extend = extend_ghosts(i)
                      / d_ratio_to_coarser[fine_level_number](i);
         int modulo = extend_ghosts(i)
                      % d_ratio_to_coarser[fine_level_number](i);
         if (modulo) {
            extend++;
         }
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

      box_list.coalesceBoxes();

      if (!d_allow_patches_smaller_than_minimum_size_to_prevent_overlaps) {
         hier::BoxUtilities<DIM>::growBoxesWithinDomain(box_list,
                                                   nesting_boxes,
                                                   smallest_box_to_refine);
      }  else {
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
            hier::BoxUtilities<DIM>::growBoxesWithinDomain(box_list,
                                                           nesting_boxes,
                                                           min_size);
         }
      }

      /*
       * Refine boxes to index space of next finer level and load balance.
       */

      box_list.refine(d_ratio_to_coarser[fine_level_number]);

      hier::IntVector<DIM> ratio_to_level_zero =
         level->getRatio() * d_ratio_to_coarser[fine_level_number];

      hier::BoxArray<DIM> physical_domain;
      hierarchy->getGridGeometry()->computePhysicalDomain(physical_domain,
                                                          ratio_to_level_zero);

      hier::IntVector<DIM> patch_cut_factor(d_ratio_to_coarser[fine_level_number]);
      getGriddingParameters(smallest_patch,
                            smallest_box_to_refine,
                            largest_patch,
                            extend_ghosts,
                            hierarchy, fine_level_number);

      t_load_balance_boxes->start();
      d_load_balancer->loadBalanceBoxes(
         fine_boxes, mapping,
         box_list, hierarchy, fine_level_number, physical_domain,
         ratio_to_level_zero,
         smallest_patch,
         getLargestPatchSize(fine_level_number),
         patch_cut_factor, extend_ghosts);
      t_load_balance_boxes->stop();

   } else {

      fine_boxes = hier::BoxArray<DIM>(box_list);

   }

   t_find_refinement->stop();
}

template<int DIM>
double MultiblockGriddingAlgorithm<DIM>::getCombineEfficiency(
   const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( (level_number >= 0) && (level_number < d_max_levels - 1) );
#endif
   int size = d_combine_efficiency.getSize();
   return( (level_number < size)
           ? d_combine_efficiency[level_number]
           : d_combine_efficiency[size-1] );
}

template<int DIM>
double MultiblockGriddingAlgorithm<DIM>::getEfficiencyTolerance(
   const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( (level_number >= 0) && (level_number < d_max_levels - 1) );
#endif
   int size = d_efficiency_tolerance.getSize();
   return( (level_number < size)
           ? d_efficiency_tolerance[level_number]
           : d_efficiency_tolerance[size-1] );
}

template<int DIM>
const hier::IntVector<DIM>&
MultiblockGriddingAlgorithm<DIM>::getRatioToCoarserLevel(
   const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( (level_number >= 0) && (level_number < d_max_levels) );
#endif
   return(d_ratio_to_coarser[level_number]);
}

/*
*************************************************************************
*                                                                       *
* Compute box array comprising the proper nesting region for given      *
* level.  The proper nesting region represents the allowable extent     *
* of the next finer level within the current level.  This box array     *
* must be properly contained within the complement of the complement    *
* box list.  Then, repeat this process recursively on successively      *
* finer levels so that each finer proper nesting array is properly      *
* contained within the proper nesting region of next coarser level.     *
*                                                                       *
*************************************************************************
*/

template<int DIM>
void MultiblockGriddingAlgorithm<DIM>::findProperNestingBoxes(
   const tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
   const int level_number,
   const int block_number,
   hier::BoxList<DIM>& complement)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!(hierarchy.isNull()));
   TBOX_ASSERT( (level_number >= 0)
            && (level_number <= hierarchy->getFinestLevelNumber()) );
   TBOX_ASSERT(!(hierarchy->getPatchLevel(level_number).isNull()));
#endif

   int fine_level_number = level_number + 1;


   if (hierarchy->levelExists(level_number)) {
      tbox::Pointer< hier::PatchLevel<DIM> > patch_level =
         hierarchy->getPatchLevel(level_number);
      d_proper_nesting_boxes[block_number][level_number] =
         hier::BoxList<DIM>(patch_level->getPhysicalDomain());

      /*
       * Grow boxes in complement and throw away pieces outside of physical
       * domain.  Then, proper nesting boxes are complement of the complement.
       */
      complement.grow(getProperNestingBuffer(level_number));
      t_intersect_boxes_find_proper->start();
      complement.intersectBoxes(
         d_proper_nesting_boxes[block_number][level_number]);
      t_intersect_boxes_find_proper->stop();

      t_remove_intersections_find_proper->start();
      d_proper_nesting_boxes[block_number][level_number].removeIntersections(
         complement);
      t_remove_intersections_find_proper->stop();

      /*
       * Recurse to finer level, if error estimation may be performed
       * on that level.  Note that the complement on the next finer
       * level will contain the complement on the current level.
       */
      if ( hierarchy->finerLevelExists(level_number)
           && levelCanBeRefined(fine_level_number) ) {
         hier::BoxList<DIM> coarse_complement(complement);
         coarse_complement.refine(d_ratio_to_coarser[fine_level_number]);
         findProperNestingBoxes(hierarchy, fine_level_number, block_number,
                                coarse_complement);
   
      }
   } // level exists

}

template<int DIM>
int MultiblockGriddingAlgorithm<DIM>::getProperNestingBuffer(
   const int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( (level_number >= 0) && (level_number < d_max_levels - 1) );
#endif
   int size = d_proper_nesting_buffer.getSize();
   return( (level_number < size)
           ? d_proper_nesting_buffer[level_number]
           : d_proper_nesting_buffer[size-1] );
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

template<int DIM>
void MultiblockGriddingAlgorithm<DIM>::regridFinerLevel(
   tbox::Pointer< hier::MultiblockPatchHierarchy<DIM> > multiblock,
   const int level_number,
   const double regrid_time,
   const int finest_level_not_regridded,
   const bool level_is_coarsest_sync_level,
   const tbox::Array<int>& tag_buffer,
   const tbox::Array<double>& regrid_start_time)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!(multiblock.isNull()));
   TBOX_ASSERT( (level_number >= 0)
          && (level_number <= multiblock->getFinestLevelNumber()) );
   TBOX_ASSERT(!(multiblock->getPatchLevel(level_number).isNull()));
   TBOX_ASSERT(finest_level_not_regridded >= 0
          && finest_level_not_regridded <= level_number);
   TBOX_ASSERT(tag_buffer.getSize() >= level_number+1);
   for (int i = 0; i < tag_buffer.getSize(); i++) {
      TBOX_ASSERT(tag_buffer[i] >= 0);
   }
#endif

   if ( levelCanBeRefined(level_number) ) {

      int fine_level_number = level_number+1;

      tbox::Pointer< hier::MultiblockPatchLevel<DIM> >
         mb_level = multiblock->getPatchLevel(level_number);

      tbox::Array<hier::BoxArray<DIM> > fine_boxes(multiblock->getNumberOfBlocks());
      hier::ProcessorMapping* mapping;
      mapping = new hier::ProcessorMapping[multiblock->getNumberOfBlocks()];

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

         d_bdry_sched_tags[level_number] =
           d_bdry_fill_tags->createSchedule(mb_level, d_mb_tagger_strategy);

         mb_level->allocatePatchData(d_tag_indx);

         for (int st = 0; st < mb_level->getNumberOfBlocks(); st++) {
            tbox::Pointer< hier::PatchLevel<DIM> > patch_level =
               mb_level->getPatchLevelForBlock(st);
            if (!patch_level.isNull()) {
               setTagsOnLevel(d_false_tag,
                              patch_level, d_tag_indx, patch_level->getBoxes());
            }
         }

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

         if ( multiblock->finerLevelExists(level_number) ) {
            tbox::Pointer< hier::MultiblockPatchLevel<DIM> > fine_level =
               multiblock->getPatchLevel(fine_level_number);
            for (int fb = 0; fb < mb_level->getNumberOfBlocks(); fb++) {
               tbox::Pointer< hier::PatchLevel<DIM> > fine_block_level =
                  fine_level->getPatchLevelForBlock(fb);
               if (!fine_block_level.isNull()) {
                  fine_boxes[fb] = fine_block_level->getBoxes();
                  fine_boxes[fb].coarsen(d_ratio_to_coarser[fine_level_number]);
                  setTagsOnLevel(d_true_tag,
                                 mb_level->getPatchLevelForBlock(fb),
                                 d_tag_indx,
                                 fine_boxes[fb]);
               }
            }
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
            tagCellsForRefinement(multiblock,
                                  level_number,
                                  regrid_time,
                                  d_tag_indx,
                                  initial_time,
                                  coarsest_sync_level,
                                  levelCanBeRefined(level_number),
                                  level_regrid_start_time);

         /*
          * Perform regridding recursively on finer levels, if appropriate.
          */
         if ( multiblock->finerLevelExists(level_number)
              && levelCanBeRefined(fine_level_number) ) {
            regridFinerLevel(multiblock,
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
 
         bufferTagsOnLevel(d_true_tag, mb_level, tag_buffer[level_number]);

         /*
          * Determine box array containing cells on level with a true tag
          * value.  The box array must be contained in array of proper
          * nesting boxes.
          */
         for (int fr = 0; fr < multiblock->getNumberOfBlocks(); fr++) {
            if (multiblock->getHierarchy(fr)->getFinestLevelNumber() >=
                level_number) {
               findRefinementBoxes(fine_boxes[fr],
                                   mapping[fr],
                                   multiblock->getHierarchy(fr),
                                   level_number,
                                   d_proper_nesting_boxes[fr][level_number]);
            }
         }

         /*
          * Deallocate tag arrays and schedule; no longer needed on current
          * level.
          */

         mb_level->deallocatePatchData(d_tag_indx);
         d_bdry_sched_tags[level_number].setNull();

         /*
          * Now that we have determined refine boxes, write those to
          * file.
          */

         if (d_write_dumped_level_boxes) {
            for (int pl = 0; pl < multiblock->getNumberOfBlocks(); pl++) {
               d_regrid_box_utility->putLevelBoxes(
                  fine_boxes[pl],
                  fine_level_number, 
                     d_regrid_box_counter[pl][fine_level_number]);
               d_regrid_box_counter[pl][fine_level_number]++;
            }
         }

      } else {

         if ( multiblock->finerLevelExists(level_number)
             && levelCanBeRefined(fine_level_number) ) {
            regridFinerLevel(multiblock,
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

         for (int br = 0; br < multiblock->getNumberOfBlocks(); br++) {

            tbox::Pointer< hier::PatchHierarchy<DIM> > block_hierarchy =
               multiblock->getHierarchy(br);

            hier::BoxArray<DIM> boxes_to_refine;

            if (d_read_dumped_level_boxes) {

               d_regrid_box_utility->getLevelBoxes(
                  boxes_to_refine,
                  level_number+1,
                  d_regrid_box_counter[br][level_number+1]);

               d_regrid_box_counter[br][level_number+1]++;

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

               boxes_to_refine.refine(d_ratio_to_coarser[fine_level_number]);

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

               hier::IntVector<DIM> ratio_to_level_zero =
                  mb_level->getRatio() *
                     d_ratio_to_coarser[fine_level_number];

               hier::BoxArray<DIM> physical_domain;
               block_hierarchy->getGridGeometry()->
                  computePhysicalDomain(physical_domain, ratio_to_level_zero);

               hier::IntVector<DIM>
                  patch_cut_factor(d_ratio_to_coarser[fine_level_number]);

               hier::IntVector<DIM> smallest_patch;
               hier::IntVector<DIM> smallest_box_to_refine;
               hier::IntVector<DIM> largest_patch;
               hier::IntVector<DIM> extend_ghosts;
               getGriddingParameters(smallest_patch,
                                     smallest_box_to_refine,
                                     largest_patch,
                                     extend_ghosts,
                                     block_hierarchy, fine_level_number);

               t_load_balance_boxes->start();
               d_load_balancer->loadBalanceBoxes(
                  fine_boxes[br], mapping[br],
                  fine_box_list, block_hierarchy,
                  fine_level_number, physical_domain,
                  ratio_to_level_zero,
                  smallest_patch,
                  getLargestPatchSize(fine_level_number),
                  patch_cut_factor, extend_ghosts);
               t_load_balance_boxes->stop();
            }
         }
      }

      /*
       * Make new finer level (fine_level_number) if necessary, or remove
       * next finer level if it is no longer needed.
       */

      bool fine_boxes_exist = false;
      for (int fb = 0; fb < multiblock->getNumberOfBlocks(); fb++) {
         if ( (fine_boxes[fb].getNumberOfBoxes() != 0) ) {
            fine_boxes_exist = true;
            break;
         }
      }
      if (fine_boxes_exist) {

         /*
          * Either remove pre-existing fine level from hierarchy and make
          * a new level, or just make a new fine level for hierarchy.
          */

         tbox::Pointer< hier::MultiblockPatchLevel<DIM> > old_fine_level;

         hier::IntVector<DIM> ratio(mb_level->getRatio()
                            * getRatioToCoarserLevel(fine_level_number));

         if ( multiblock->finerLevelExists(level_number) ) {
            old_fine_level = multiblock->getPatchLevel(fine_level_number);
            for (int rp = 0; rp < multiblock->getNumberOfBlocks(); rp++) {
               tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy =
                  multiblock->getHierarchy(rp);
               if (hierarchy->finerLevelExists(level_number)) {
                  hierarchy->removePatchLevel(fine_level_number);
               }
            }
            ratio = old_fine_level->getRatio();
         }

         for (int nb = 0; nb < multiblock->getNumberOfBlocks(); nb++) {
            if (fine_boxes[nb].getNumberOfBoxes() != 0) {
               multiblock->getHierarchy(nb)->
                  makeNewPatchLevel(fine_level_number, ratio,
                                    fine_boxes[nb], mapping[nb]);
            }
         }

         multiblock->adjustMultiblockPatchLevelBoundaries(
            multiblock->getMultiblockPatchLevel(fine_level_number));

         // "false" argument": const bool initial_time = false;
         d_tag_init_strategy->initializeLevelData(multiblock,
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

         if ( do_tagging ) {

            /*
             * If current level (level_number) is coarser than the finest
             * level subject to regridding, then adjust tags for next coarser
             * level (i.e., level_number-1) so that the next finer level
             * (fine_level_number) will be nested within the new level
             * (level_number) on return from recursion.
             */

            if (level_number > finest_level_not_regridded) {
               for (int bc = 0; bc < multiblock->getNumberOfBlocks(); bc++) {
                  fine_boxes[bc].coarsen(d_ratio_to_coarser[fine_level_number]
                                         * d_ratio_to_coarser[level_number]);

                  tbox::Pointer< hier::PatchHierarchy<DIM> > block_hierarchy =
                     multiblock->getHierarchy(bc);

                  if (block_hierarchy->getFinestLevelNumber() >=
                      level_number-1) {

                     tbox::Pointer< hier::PatchLevel<DIM> > block_level =
                        block_hierarchy->getPatchLevel(level_number-1);

                     if (!block_level.isNull()) {
                        
                        setTagsOnLevel(d_true_tag,
                                       block_level,
                                       d_tag_indx,
                                       fine_boxes[bc]);

                     }
                     
                  }
               }
            }
         }

      } else {

         /*
          * If there are no boxes for the new fine level, remove the
          * pre-existing fine level if it existed.
          */

         if ( multiblock->finerLevelExists(level_number)
              && remove_old_fine_level) {
            for (int rl = 0; rl < multiblock->getNumberOfBlocks(); rl++) {
               if (multiblock->getHierarchy(rl)->
                   finerLevelExists(level_number)) {
                  multiblock->getHierarchy(rl)->
                     removePatchLevel(fine_level_number);
               }
            }
         }

      } // if we are not re-gridding the level

      delete[] mapping;

   } //  if level cannot be refined, the routine drops through...

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

template<int DIM>
void MultiblockGriddingAlgorithm<DIM>::getFromRestart()
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

template<int DIM>
int MultiblockGriddingAlgorithm<DIM>::getErrorCoarsenRatio()
const
{
   return(d_tag_init_strategy->getErrorCoarsenRatio());
}

template<int DIM>
tbox::Pointer<mesh::TagAndInitializeStrategy<DIM> >
MultiblockGriddingAlgorithm<DIM>::getTagAndInitializeStrategy() const
{
   return(d_tag_init_strategy);
}


}
}
#endif
