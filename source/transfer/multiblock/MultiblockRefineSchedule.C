//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-1/source/transfer/multiblock/MultiblockRefineSchedule.C $
// Package:     SAMRAI multiblock package
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1846 $
// Modified:    $LastChangedDate: 2008-01-11 09:51:05 -0800 (Fri, 11 Jan 2008) $
// Description: Base class for geometry management on patches
//

#ifndef included_xfer_MultiblockRefineSchedule_C
#define included_xfer_MultiblockRefineSchedule_C

#include "MultiblockRefineSchedule.h"

#include "BoxUtilities.h"
#include "MBUtilities.h"
#include "MultiblockRefinePatchStrategy.h"
#include "MultiblockRefineAlgorithm.h"
#include "StandardRefineTransactionFactory.h"


namespace SAMRAI {
    namespace xfer {

/*
 * ************************************************************************
 *                                                                        *
 * Create a refine schedule that copies data from the source level into   *
 * the destination level on the components represented by the refine      *
 * It is assumed that the index spaces of the source and destination      *
 * levels represent the same grid resolution.                             *
 *                                                                        *
 * ************************************************************************
 */
//
template<int DIM>
MultiblockRefineSchedule<DIM>::MultiblockRefineSchedule(
   const std::string& fill_pattern,
   tbox::Pointer< hier::MultiblockPatchLevel<DIM> > dst_level, 
   tbox::Pointer< hier::MultiblockPatchLevel<DIM> > src_level, 
   tbox::Pointer< hier::MultiblockPatchHierarchy<DIM> > multiblock,
   tbox::Pointer< xfer::RefineAlgorithm<DIM> > refine_alg,
   tbox::Pointer< xfer::RefineTransactionFactory<DIM> > transaction_factory,
   MultiblockRefinePatchStrategy<DIM>* strategy,
   bool use_time_refinement)
{
   d_multiblock_hierarchy = multiblock;
   d_single_block_refine_alg = refine_alg;

   d_fill_pattern = fill_pattern;
   d_transaction_factory = transaction_factory;

   d_using_standard_transaction = false;
   tbox::Pointer< xfer::StandardRefineTransactionFactory<DIM> > t_factory =
      d_transaction_factory;

   if (!(t_factory.isNull())) {
      d_using_standard_transaction = true;
   }

   d_single_block_scratch_refine_alg.setNull();

   d_multiblock_dst_level = dst_level;

   d_single_block_fill_local.resizeArray(multiblock->getNumberBlocks());
   d_local_fill_only.resizeArray(multiblock->getNumberBlocks());
   d_unfilled_boxes.resizeArray(multiblock->getNumberBlocks());
   d_multiblock_coarse_scratch_level.resizeArray(multiblock->getNumberBlocks());
   d_multiblock_coarse_schedule.resizeArray(multiblock->getNumberBlocks());
   d_coarse_selector.resizeArray(multiblock->getNumberBlocks());

   d_multiblock_strategy = strategy;

   d_neighbor_single_block_refine_schedule.resizeArray(multiblock->getNumberBlocks());
   d_neighbor_ghost_level.resizeArray(multiblock->getNumberBlocks());
   d_finalize_ghost_level.resizeArray(multiblock->getNumberBlocks());
   d_finalize_ghost_patch_numbers.resizeArray(multiblock->getNumberBlocks());
   d_finalize_ghost_num_src_patches.resizeArray(multiblock->getNumberBlocks());

   d_neighbor_copy_only.resizeArray(multiblock->getNumberBlocks());
   d_neighbor_unfilled_boxes.resizeArray(multiblock->getNumberBlocks());
   d_neighbor_multiblock_coarse_level.resizeArray(multiblock->getNumberBlocks());
   d_neighbor_multiblock_coarse_schedule.resizeArray(multiblock->getNumberBlocks());

   for (int nb = 0; nb < d_multiblock_hierarchy->getNumberBlocks(); nb++) {
      d_coarse_selector[nb] = new hier::ComponentSelector();
      d_coarse_selector[nb]->clrAllFlags();

      tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy =
         d_multiblock_hierarchy->getHierarchy(nb);

      tbox::Pointer< hier::PatchLevel<DIM> > dst_patch_level =
         dst_level->getPatchLevelForBlock(nb);

      if (!dst_patch_level.isNull()) {
         tbox::Pointer< hier::PatchLevel<DIM> > src_patch_level;
         if (!(src_level.isNull())) {
            src_patch_level = src_level->getPatchLevelForBlock(nb);
         } else {
            src_patch_level.setNull();
         }
         if (!src_patch_level.isNull()) {
            d_single_block_fill_local[nb] =
               d_single_block_refine_alg->createSchedule(
                  d_fill_pattern,
                  dst_patch_level,
                  src_patch_level,
                  (xfer::RefinePatchStrategy<DIM>*) strategy,
                  use_time_refinement,
                  d_transaction_factory);
         } else {
            d_single_block_fill_local[nb].setNull();
         }
      } else {

         d_single_block_fill_local[nb].setNull();

      }

      d_local_fill_only[nb] = true;
   }

   if (strategy != NULL) {
      createInterblockSchedules(dst_level,
                                src_level,
                                (xfer::RefinePatchStrategy<DIM>*)strategy,
                                MULTIBLOCK_FAKE_LEVEL_NUMBER,
                                use_time_refinement);
   }
}

/*
 * ************************************************************************
 *                                                                        *
 * Create a refine schedule that copies data from the source level into   *
 * the destination level on the components represented by the refine      *
 * list.  If portions of the destination level remain unfilled, then      *
 * the algorithm recursively fills those unfilled portions from coarser   *
 * levels in the hierarchies of the multiblock object.  It is assumed     *
 * that the index spaces of the source and destination levels represent   *
 * the same grid resolution.  Also, the next coarser level integer        *
 * argument must be the number of level in the multiblock hierarchies     *
 * representing the next coarser level of mesh resolution to the          *
 * destination level.                                                     *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM> 
MultiblockRefineSchedule<DIM>::MultiblockRefineSchedule(
   const std::string& fill_pattern,
   tbox::Pointer< hier::MultiblockPatchLevel<DIM> > dst_level, 
   tbox::Pointer< hier::MultiblockPatchLevel<DIM> > src_level,
   const int next_coarser_level,
   tbox::Pointer< hier::MultiblockPatchHierarchy<DIM> > multiblock,
   tbox::Pointer< xfer::RefineAlgorithm<DIM> > refine_alg,
   tbox::Pointer< xfer::RefineTransactionFactory<DIM> > transaction_factory,
   MultiblockRefinePatchStrategy<DIM>* strategy,
   bool use_time_refinement)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(strategy != NULL);
#endif

   d_multiblock_hierarchy = multiblock;
   d_single_block_refine_alg = refine_alg;

   d_fill_pattern = fill_pattern;
   d_transaction_factory = transaction_factory;

   d_using_standard_transaction = false;
   tbox::Pointer< xfer::StandardRefineTransactionFactory<DIM> > t_factory =
      d_transaction_factory;

   if (!(t_factory.isNull())) {
      d_using_standard_transaction = true;
   }

   d_single_block_scratch_refine_alg.setNull();

   constructScratchRefineAlgorithm();

   d_multiblock_dst_level = dst_level;

   d_multiblock_strategy = strategy;

   d_single_block_fill_local.resizeArray(multiblock->getNumberBlocks());
   d_local_fill_only.resizeArray(multiblock->getNumberBlocks());
   d_multiblock_coarse_schedule.resizeArray(multiblock->getNumberBlocks());
   d_coarse_selector.resizeArray(multiblock->getNumberBlocks());
   d_unfilled_boxes.resizeArray(multiblock->getNumberBlocks());
   d_multiblock_coarse_scratch_level.resizeArray(multiblock->getNumberBlocks());

   d_neighbor_single_block_refine_schedule.resizeArray(multiblock->getNumberBlocks());
   d_neighbor_ghost_level.resizeArray(multiblock->getNumberBlocks());
   d_finalize_ghost_level.resizeArray(multiblock->getNumberBlocks());
   d_finalize_ghost_patch_numbers.resizeArray(multiblock->getNumberBlocks());
   d_finalize_ghost_num_src_patches.resizeArray(multiblock->getNumberBlocks());

   d_neighbor_copy_only.resizeArray(multiblock->getNumberBlocks());
   d_neighbor_unfilled_boxes.resizeArray(multiblock->getNumberBlocks());
   d_neighbor_multiblock_coarse_level.resizeArray(multiblock->getNumberBlocks());
   d_neighbor_multiblock_coarse_schedule.resizeArray(multiblock->getNumberBlocks());

   for (int nb = 0; nb < d_multiblock_hierarchy->getNumberBlocks(); nb++) {

      d_coarse_selector[nb] = new hier::ComponentSelector();
      d_coarse_selector[nb]->clrAllFlags();

      d_local_fill_only[nb] = true;

      tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy =
         d_multiblock_hierarchy->getHierarchy(nb);

      tbox::Pointer< hier::PatchLevel<DIM> > dst_patch_level =
         dst_level->getPatchLevelForBlock(nb);

      if (!dst_patch_level.isNull()) {

         tbox::Pointer< hier::PatchLevel<DIM> > src_patch_level;
         if (!(src_level.isNull())) {
            src_patch_level = src_level->getPatchLevelForBlock(nb);
         } else {
            src_patch_level.setNull();
         }

         bool need_other_source_blocks = true;

         tbox::Pointer< hier::PatchDescriptor<DIM> > descriptor =
            dst_patch_level->getPatchDescriptor();

         hier::IntVector<DIM> gcw(descriptor->getMaxGhostWidth());

         if (!src_patch_level.isNull()) {

            hier::BoxList<DIM> dst_boxes(dst_patch_level->getBoxes());
            hier::BoxList<DIM> src_boxes(src_patch_level->getBoxes());

            dst_boxes.grow(gcw);
            dst_boxes.intersectBoxes(dst_patch_level->getPhysicalDomain());

            hier::BoxList<DIM> domain_outside_block;
            d_multiblock_hierarchy->getDomainOutsideBlock(domain_outside_block, nb);
            domain_outside_block.refine(dst_patch_level->getRatio());

            need_other_source_blocks =
               needOtherSourceBlocks(dst_boxes,
                                     src_boxes,
                                     domain_outside_block);
         }

         if (need_other_source_blocks) {

            d_local_fill_only[nb] = false;

            d_unfilled_boxes[nb] = dst_patch_level->getBoxes();
            d_unfilled_boxes[nb].grow(gcw);
            d_unfilled_boxes[nb].intersectBoxes(dst_patch_level->
                                                    getPhysicalDomain());
            if (!src_patch_level.isNull()) {
               d_unfilled_boxes[nb].removeIntersections(
                  src_patch_level->getBoxes());
            }

            d_unfilled_boxes[nb].coalesceBoxes();

            hier::IntVector<DIM> ratio(dst_patch_level->getRatio());

            hier::IntVector<DIM> extend_ghosts;

            for (int d = 0; d < DIM; d++) {
               extend_ghosts(d) = gcw(d) * ratio(d);
            }

            hier::BoxArray<DIM> tc_boxes(dst_patch_level->getNumberOfPatches());
            tbox::Array<int> map_array(dst_patch_level->getNumberOfPatches());
            int num_tc_boxes = 0;
            for (int j = 0; j < dst_patch_level->getNumberOfPatches(); j++) {
               
               hier::Box<DIM> dst_box(dst_patch_level->getBoxes()(j));
               dst_box.grow(gcw);
               hier::BoxList<DIM> dst_box_list;
               dst_box_list.unionBoxes(dst_box);
               dst_box_list.intersectBoxes(
                  dst_patch_level->getPhysicalDomain());

               for (typename hier::BoxList<DIM>::Iterator
                    db(dst_box_list); db; db++) {
                  if (num_tc_boxes == map_array.size()) {
                     tc_boxes.resizeBoxArray(num_tc_boxes*2);
                     map_array.resizeArray(num_tc_boxes*2);
                  }
                  tc_boxes(num_tc_boxes) = db();
                  map_array[num_tc_boxes] =
                     dst_patch_level->getMappingForPatch(j); 
                     num_tc_boxes++;
               } 
            }
            tc_boxes.resizeBoxArray(num_tc_boxes);
            map_array.resizeArray(num_tc_boxes);

            hier::ProcessorMapping coarse_map(map_array);

            hier::BoxList<DIM> tmp_crs_boxes(tc_boxes);

            (void) hier::BoxUtilities<DIM>::extendBoxesToDomainBoundary(
               tmp_crs_boxes,
               dst_patch_level->getPhysicalDomain(),
               extend_ghosts);

            hier::IntVector<DIM> coarse_ratio =
               dst_patch_level->getRatioToCoarserLevel();
            if (coarse_ratio == hier::IntVector<DIM>(0)) {
               tbox::Pointer< hier::MultiblockPatchLevel<DIM> > example_level =
                   d_multiblock_hierarchy->getPatchLevel(next_coarser_level);
               coarse_ratio = dst_patch_level->getRatio() /
                              example_level->getRatio();
            }

            tmp_crs_boxes.coarsen(coarse_ratio);

            tbox::Pointer< hier::PatchLevel<DIM> > coarse_level =
               new hier::PatchLevel<DIM>(
                  tmp_crs_boxes,
                  coarse_map,
                  d_multiblock_hierarchy->getPatchLevel(next_coarser_level)
                              ->getRatio(),
                  dst_patch_level->getGridGeometry(),
                  dst_patch_level->getPatchDescriptor(),
                  dst_patch_level->getPatchFactory());

            createCoarseSchedule(coarse_level, next_coarser_level,
                                 dst_patch_level->getRatioToCoarserLevel(),
                                 hierarchy, nb);

         } else {

            d_local_fill_only[nb] = true;

         }
         if (d_local_fill_only[nb]) {
            d_single_block_fill_local[nb] =
               d_single_block_refine_alg->createSchedule(
                  d_fill_pattern,
                  dst_patch_level,
                  src_patch_level,
                  next_coarser_level,
                  hierarchy,
                  (xfer::RefinePatchStrategy<DIM>*) strategy,
                  use_time_refinement,
                  d_transaction_factory);
         } else {
            if (!src_patch_level.isNull()) {
               d_single_block_fill_local[nb] =
                  d_single_block_refine_alg->createSchedule(
                     d_fill_pattern,
                     dst_patch_level,
                     src_patch_level,
                     (xfer::RefinePatchStrategy<DIM>*) NULL,
                     use_time_refinement,
                     d_transaction_factory);
            } else {
               d_single_block_fill_local[nb].setNull();
            }
         }
      } else {
         d_single_block_fill_local[nb].setNull();
      }
   }

   createInterblockSchedules(dst_level,
                             src_level,
                             (xfer::RefinePatchStrategy<DIM>*) strategy,
                             next_coarser_level+1,
                             use_time_refinement);
}

/*
 * ************************************************************************
 *                                                                        *
 * The destructor implicitly deallocates all of the data associated with  *
 * the communication schedule.                                            *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM> 
MultiblockRefineSchedule<DIM>::~MultiblockRefineSchedule()
{
}

/*
 * ************************************************************************
 *                                                                        *
 * Create the xfer::RefineSchedules that will transfer data across block   *
 * boundaries.                                                            *
 *                                                                        *
 * ************************************************************************
 */
template<int DIM>
void MultiblockRefineSchedule<DIM>::createInterblockSchedules(
   tbox::Pointer< hier::MultiblockPatchLevel<DIM> > dst_level,
   tbox::Pointer< hier::MultiblockPatchLevel<DIM> > src_level,
   xfer::RefinePatchStrategy<DIM>* refine_strategy,
   int level_number,
   bool use_time_refinement)
{

   for (int nb = 0; nb < d_multiblock_hierarchy->getNumberBlocks(); nb++) {

      const int ln = level_number;

      tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy =
         d_multiblock_hierarchy->getHierarchy(nb);

      tbox::Pointer< hier::PatchLevel<DIM> > level = 
         dst_level->getPatchLevelForBlock(nb);

      if (!level.isNull()) {
         int num_neighbors = d_multiblock_hierarchy->getNumberNeighbors(nb);

         /*
          * resize arrays to proper size if needed
          */
         if (d_neighbor_single_block_refine_schedule[nb].getSize() < num_neighbors) {
            d_neighbor_single_block_refine_schedule[nb].resizeArray(num_neighbors);
         }

         if (d_neighbor_ghost_level[nb].getSize() < num_neighbors) {
            d_neighbor_ghost_level[nb].resizeArray(num_neighbors);
         }
         if (d_finalize_ghost_level[nb].getSize() < num_neighbors) {
            d_finalize_ghost_level[nb].resizeArray(num_neighbors);
         }
         if (d_finalize_ghost_patch_numbers[nb].getSize() < num_neighbors) {
            d_finalize_ghost_patch_numbers[nb].resizeArray(num_neighbors);
         }
         if (d_finalize_ghost_num_src_patches[nb].getSize() < num_neighbors) {
            d_finalize_ghost_num_src_patches[nb].resizeArray(num_neighbors);
         }

         if (d_neighbor_copy_only[nb].getSize() < num_neighbors) {
            d_neighbor_copy_only[nb].resizeArray(num_neighbors);
         }
         if (d_neighbor_multiblock_coarse_level[nb].getSize() < num_neighbors) {
            d_neighbor_multiblock_coarse_level[nb].resizeArray(num_neighbors);
         }
         if (d_neighbor_multiblock_coarse_schedule[nb].getSize() < num_neighbors) {
            d_neighbor_multiblock_coarse_schedule[nb].resizeArray(num_neighbors);
         }
         if (d_neighbor_unfilled_boxes[nb].getSize() < num_neighbors) {
            d_neighbor_unfilled_boxes[nb].resizeArray(num_neighbors);
         }

         for (int nn = 0; nn < num_neighbors; nn++) {
            d_neighbor_copy_only[nb][nn] = true;
         }

         int nc = 0;
         for (typename
              tbox::List<typename hier::MultiblockPatchHierarchy<DIM>::Neighbor>::
              Iterator ni(d_multiblock_hierarchy->getNeighbors(nb));
              ni; ni++) {
            int id = ni().d_id;
            hier::IntVector<DIM> shift (ni().d_shift);

            hier::IntVector<DIM> ghost_width =
               hierarchy->getGridGeometry()->computeMaxGhostWidth(
                  hierarchy->getPatchDescriptor());

            hier::BoxList<DIM> dst_ghost_boxes(level->getBoxes());
            dst_ghost_boxes.grow(ghost_width);
            dst_ghost_boxes.removeIntersections(level->getBoxes());

            /*
             * get intersection with translated neighbor_boxes and
             * dst_ghost_boxes
             */
            hier::BoxArray<DIM> trans_neighbor_boxes;
            d_multiblock_hierarchy->getTranslatedBlock(trans_neighbor_boxes, nb, id);
            trans_neighbor_boxes.refine(level->getRatio());

            hier::BoxList<DIM> trans_neighbor_list(trans_neighbor_boxes);
            trans_neighbor_list.intersectBoxes(dst_ghost_boxes);

            /*
             * trans_neighbor_list is a list that contains the boxes for the
             * temporary dst level.  We need to create an equivalent array such
             * that each box in the array corresponds to a box on the
             * destination level, and is mapped to the same processor as the
             * corresponding box.
             */
            hier::BoxArray<DIM> finalize_box_array(level->getNumberOfPatches());
            tbox::Array<int> finalize_map_array(level->getNumberOfPatches());
            d_finalize_ghost_patch_numbers[nb][nc].resizeArray(
               level->getNumberOfPatches());
            d_finalize_ghost_num_src_patches[nb][nc].resizeArray(
               level->getNumberOfPatches());
            int num_finalize_boxes = 0;

            for (int p = 0; p < level->getNumberOfPatches(); p++) {
               if (level->patchTouchesRegularBoundary(p)) {
                  hier::Box<DIM> dst_grow_box(level->getBoxes()(p));
                  dst_grow_box.grow(ghost_width);

                  hier::BoxList<DIM> tmp_list(trans_neighbor_list);
                  tmp_list.intersectBoxes(dst_grow_box);

                  if (ln == MULTIBLOCK_FAKE_LEVEL_NUMBER && tmp_list.size() > 0) {

                     tbox::Pointer< hier::PatchLevel<DIM> > neighbor_level =
                        src_level->getPatchLevelForBlock(id);

                     if (!neighbor_level.isNull()) {

                        hier::BoxArray<DIM> neighbor_boxes(neighbor_level->
                                                         getBoxes());

                        d_multiblock_hierarchy->translateBoxArray(
                           neighbor_boxes,
                           neighbor_level->getRatio(),
                           nb, id);

                        tmp_list.intersectBoxes(neighbor_boxes);

                     }
                  }

                  if (tmp_list.size() > 0) {

                     tmp_list.coalesceBoxes();

                     if (tmp_list.size() > 1) {
                        hier::Box<DIM> bound_box(tmp_list.getBoundingBox());
                        hier::BoxList<DIM> bound_list;
                        bound_list.addItem(bound_box);
                        bound_list.removeIntersections(tmp_list);
                        if (bound_list.size() == 0) {
                           tmp_list.clearItems();
                           tmp_list.addItem(bound_box);
                        }                  
                     }
                     int map_num = level->getMappingForPatch(p);

                     d_finalize_ghost_patch_numbers[nb][nc][p] =
                        num_finalize_boxes;
                     d_finalize_ghost_num_src_patches[nb][nc][p] =
                        tmp_list.size();
                     for (typename hier::BoxList<DIM>::Iterator bli(tmp_list);
                          bli; bli++) {
                        if (num_finalize_boxes ==
                            finalize_box_array.getNumberOfBoxes()) {
   
                           finalize_box_array.resizeBoxArray(
                              2*finalize_box_array.getNumberOfBoxes());
                           finalize_map_array.resizeArray(
                              2*finalize_box_array.getNumberOfBoxes());
                        }
                        finalize_box_array(num_finalize_boxes) = bli();
                        finalize_map_array[num_finalize_boxes] = map_num;
                        num_finalize_boxes++;
                     }
                  } else {
                     d_finalize_ghost_patch_numbers[nb][nc][p] = -1;
                     d_finalize_ghost_num_src_patches[nb][nc][p] = 0;
                  }
               } else {
                  d_finalize_ghost_patch_numbers[nb][nc][p] = -1;
                     d_finalize_ghost_num_src_patches[nb][nc][p] = 0;
               }
            }
            finalize_box_array.resizeBoxArray(num_finalize_boxes);
            finalize_map_array.resizeArray(num_finalize_boxes);

            hier::ProcessorMapping fin_mapping(finalize_map_array);

            if (num_finalize_boxes) {
               d_finalize_ghost_level[nb][nc] =
                  new hier::PatchLevel<DIM>(finalize_box_array,
                                       fin_mapping,
                                       level->getRatio(),
                                       hierarchy->getGridGeometry(),
                                       hierarchy->getPatchDescriptor());

               /*
                * From here take the finalize_box_array and translate each box
                * to the neighbor's index space.  Create a level with the same
                * mapping.  Then create only one refine schedule to copy from
                * the neighbor source level to the neighbor temp level.
                */
               hier::BoxArray<DIM> neighbor_ghost_array(finalize_box_array);

               d_multiblock_hierarchy->translateBoxArray(
                  neighbor_ghost_array, level->getRatio(),
                  ni().d_id,
                  nb);

               tbox::Pointer< hier::PatchHierarchy<DIM> > neighbor_hierarchy =
                  d_multiblock_hierarchy->getHierarchy(ni().d_id);

               d_neighbor_ghost_level[nb][nc] =
                  new hier::PatchLevel<DIM>(neighbor_ghost_array,
                                       fin_mapping,
                                       level->getRatio(),
                                       neighbor_hierarchy->
                                          getGridGeometry(),
                                       neighbor_hierarchy->
                                          getPatchDescriptor());

               if (ln != MULTIBLOCK_FAKE_LEVEL_NUMBER) {

                  tbox::Pointer< hier::PatchLevel<DIM> > neighbor_level;

                  if (neighbor_hierarchy->getNumberLevels() > ln) {
                     neighbor_level = neighbor_hierarchy->getPatchLevel(ln);
                  } else {
                     neighbor_level.setNull();
                  }

                  bool need_other_source_blocks = true;

                  if (!neighbor_level.isNull()) {
                     hier::BoxList<DIM> dst_boxes(
                        d_neighbor_ghost_level[nb][nc]->getBoxes());
                     hier::BoxList<DIM> src_boxes(neighbor_level->getBoxes());

                     hier::BoxList<DIM> domain_outside_block;
                     d_multiblock_hierarchy->getDomainOutsideBlock(domain_outside_block,
                                                         ni().d_id);
                     domain_outside_block.refine(
                        d_neighbor_ghost_level[nb][nc]->getRatio());

                     need_other_source_blocks =
                        needOtherSourceBlocks(dst_boxes,
                                           src_boxes,
                                           domain_outside_block);
                  }

                  if (need_other_source_blocks) {

                     tbox::Pointer< hier::PatchLevel<DIM> > neighbor_dst_level =
                        d_neighbor_ghost_level[nb][nc];

                     d_neighbor_copy_only[nb][nc] = false;

                     d_neighbor_unfilled_boxes[nb][nc] =
                        neighbor_dst_level->getBoxes();
                     if (!neighbor_level.isNull()) {
                        d_neighbor_unfilled_boxes[nb][nc].removeIntersections(
                           neighbor_level->getBoxes());
                     }

                     d_neighbor_unfilled_boxes[nb][nc].coalesceBoxes();

                     hier::IntVector<DIM> gcw(neighbor_dst_level->
                                            getPatchDescriptor()->
                                               getMaxGhostWidth());

                     hier::IntVector<DIM> ratio(neighbor_dst_level->getRatio());

                     hier::IntVector<DIM> extend_ghosts;

                     for (int d = 0; d < DIM; d++) {
                        extend_ghosts(d) = gcw(d) * ratio(d);
                     }

                     hier::BoxList<DIM> tmp_src_boxes(
                        neighbor_dst_level->getBoxes());

                     (void)hier::BoxUtilities<DIM>::extendBoxesToDomainBoundary(
                        tmp_src_boxes,
                        neighbor_dst_level->getPhysicalDomain(),
                        extend_ghosts);

                     hier::IntVector<DIM> coarse_ratio =
                        level->getRatioToCoarserLevel();
                     if (coarse_ratio == hier::IntVector<DIM>(0)) {
                        tbox::Pointer< hier::MultiblockPatchLevel<DIM> >
                           example_level = d_multiblock_hierarchy->getPatchLevel(ln-1);
                        coarse_ratio = neighbor_dst_level->getRatio() /
                                       example_level->getRatio();
                     }
                     tmp_src_boxes.coarsen(coarse_ratio);

                     tbox::Pointer< hier::PatchLevel<DIM> >
                        neighbor_coarse_level =
                           new hier::PatchLevel<DIM>(
                              tmp_src_boxes,
                              neighbor_dst_level->getProcessorMapping(),
                              d_multiblock_hierarchy->getPatchLevel(ln-1)->getRatio(),
                              neighbor_dst_level->getGridGeometry(),
                              neighbor_dst_level->getPatchDescriptor(),
                              neighbor_dst_level->getPatchFactory());

                     createNeighborCoarseSchedule(
                        neighbor_coarse_level, ln-1,
                        level->getRatioToCoarserLevel(),
                        hierarchy, ni().d_id, nb, nc);

                  } else {

                     d_neighbor_copy_only[nb][nc] = true;

                  }

                  if (d_neighbor_copy_only[nb][nc]) {
                     d_neighbor_single_block_refine_schedule[nb][nc] =
                        d_single_block_refine_alg->createSchedule(
                           d_fill_pattern,
                           d_neighbor_ghost_level[nb][nc],
                           neighbor_level,
                           ln-1,
                           neighbor_hierarchy,
                           refine_strategy,
                           use_time_refinement,
                           d_transaction_factory);
                  } else {
                     if (!neighbor_level.isNull()) {
                        d_neighbor_single_block_refine_schedule[nb][nc] =
                           d_single_block_refine_alg->createSchedule(
                              d_fill_pattern,
                              d_neighbor_ghost_level[nb][nc],
                              neighbor_level,
                              NULL,
                              use_time_refinement,
                              d_transaction_factory);
                     } else {
                        d_neighbor_single_block_refine_schedule[nb][nc].
                           setNull();
                     }
                  }
               } else {

                  tbox::Pointer< hier::PatchLevel<DIM> > neighbor_level =
                     src_level->getPatchLevelForBlock(ni().d_id);

                  if (!neighbor_level.isNull()) { 
                     d_neighbor_single_block_refine_schedule[nb][nc] =
                        d_single_block_refine_alg->createSchedule(
                           d_fill_pattern,
                           d_neighbor_ghost_level[nb][nc],
                           neighbor_level,
                           refine_strategy,
                           use_time_refinement,
                           d_transaction_factory);
                  } else {
                     d_neighbor_single_block_refine_schedule[nb][nc].setNull();
                  }
               }
            } else {
               d_neighbor_single_block_refine_schedule[nb][nc].setNull();
            }

            nc++; 
         }
      }
   }
}

/*
 * ************************************************************************
 *                                                                        *
 * Execute the communication schedules that copy data into the            *
 * destination component of the destination level.                        *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM>
void MultiblockRefineSchedule<DIM>::fillData(double fill_time,
                                         bool do_physical_boundary_fill) const
{
   fillData(fill_time,
            do_physical_boundary_fill,
            false,
            false);
}

template<int DIM>
void MultiblockRefineSchedule<DIM>::fillData(double fill_time,
                                         bool do_physical_boundary_fill,
                                         bool filling_coarse_scratch,
                                         bool filling_crse_scr_recursive) const
{
   bool is_scr_recursive = false;
   tbox::Array<hier::ComponentSelector> fill_crse_scr_selector(0);

   if (filling_coarse_scratch) {
      if (d_multiblock_strategy != NULL) {
         d_multiblock_strategy->setFillingCoarseScratch(true);
         is_scr_recursive = true;
      }
      fill_crse_scr_selector.resizeArray(d_multiblock_hierarchy->getNumberBlocks());
   }

   tbox::Array<hier::ComponentSelector>
      allocate_scr_vector(d_multiblock_hierarchy->getNumberBlocks());

   int nb;

   if (!filling_coarse_scratch) {
      for (nb = 0; nb < d_multiblock_hierarchy->getNumberBlocks(); nb++) {
         tbox::Pointer< hier::PatchLevel<DIM> > dst_patch_level=
            d_multiblock_dst_level->getPatchLevelForBlock(nb);

         if (!dst_patch_level.isNull()) {
            allocateScratchSpace(
               dst_patch_level,
               fill_time,
               allocate_scr_vector[nb]);
         }
      }
   }

   for (nb = 0; nb < d_multiblock_hierarchy->getNumberBlocks(); nb++) {

      if (d_local_fill_only[nb]) {
         if (!d_single_block_fill_local[nb].isNull()) {
            if (filling_coarse_scratch) {
               d_single_block_fill_local[nb]->allocateDestinationSpace(
                  fill_time, *d_coarse_selector[nb]);
            }
            if (d_multiblock_strategy != NULL) {
               d_multiblock_strategy->setBlockNumber(nb);
            }
            d_single_block_fill_local[nb]->fillData(fill_time, false);
            if (d_multiblock_strategy != NULL) {
               d_multiblock_strategy->clearBlockNumber();
            }
         }
      } else {
         if (d_multiblock_strategy != NULL) {
            d_multiblock_strategy->setBlockNumber(nb);
         }
         d_multiblock_coarse_schedule[nb]->fillData(fill_time,
                                        true,
                                        true,
                                        is_scr_recursive);
         if (d_multiblock_strategy != NULL) {
            d_multiblock_strategy->clearBlockNumber();
         }
         tbox::Pointer< hier::PatchLevel<DIM> > dst_patch_level=
            d_multiblock_dst_level->getPatchLevelForBlock(nb);

         if (filling_coarse_scratch) {
            if (!dst_patch_level.isNull()) {
               d_multiblock_coarse_schedule[nb]->allocateScratchSpace(
                  dst_patch_level,
                  fill_time, *d_coarse_selector[nb]);
            }
         }
 
         if (!d_single_block_fill_local[nb].isNull()) {
            if (d_multiblock_strategy != NULL) {
               d_multiblock_strategy->setBlockNumber(nb);
            }
            d_single_block_fill_local[nb]->fillData(fill_time, false);
            if (d_multiblock_strategy != NULL) {
               d_multiblock_strategy->clearBlockNumber();
            }
         }

         if (!dst_patch_level.isNull()) {
            if (d_multiblock_strategy != NULL) {
               d_multiblock_strategy->setBlockNumber(nb);
            }
            refineScratchData(d_multiblock_coarse_scratch_level[nb],
                              dst_patch_level,
                              d_unfilled_boxes[nb],
                              nb,
                              do_physical_boundary_fill);

            copyScratchToDestination(dst_patch_level,
                                     d_unfilled_boxes[nb],
                                     d_single_block_refine_alg->getEquivalenceClasses()); 
            if (d_multiblock_strategy != NULL) {
               d_multiblock_strategy->clearBlockNumber();
            }
         }

         tbox::Pointer<hier::ComponentSelector> crs_scr_vector =
            d_multiblock_coarse_schedule[nb]->getCoarseScratchVector(nb);
         d_multiblock_coarse_scratch_level[nb]->
            deallocatePatchData(*crs_scr_vector);
      }
   }

   if (d_multiblock_strategy != NULL) {
      for (nb = 0; nb < d_multiblock_hierarchy->getNumberBlocks(); nb++) {

         tbox::Pointer< hier::PatchLevel<DIM> > level =
            d_multiblock_dst_level->getPatchLevelForBlock(nb);

         d_multiblock_strategy->setBlockNumber(nb);
         
         if (!level.isNull()) {
            tbox::Array< tbox::List<SingularityPatch> >
               singularity_patches(level->getNumberOfPatches());

            bool singularity_to_fill = false;
            if (d_multiblock_hierarchy->reducedConnectivityExists(nb)) {
               singularity_to_fill = true;
            }

            tbox::Array<hier::ComponentSelector>
               local_selector(level->getNumberOfPatches());

            int nc = 0;
            for (typename
                 tbox::List<typename hier::MultiblockPatchHierarchy<DIM>::Neighbor>::
                 Iterator ni(d_multiblock_hierarchy->getNeighbors(nb));
                 ni; ni++) {
               int id = ni().d_id;

               hier::ComponentSelector dst_vector;
               if (!(d_neighbor_single_block_refine_schedule[nb][nc].isNull()) ||
                   !(d_neighbor_multiblock_coarse_schedule[nb][nc].isNull())) {

                  if (!d_neighbor_single_block_refine_schedule[nb][nc].isNull()) {
                     d_neighbor_single_block_refine_schedule[nb][nc]->
                        allocateDestinationSpace(fill_time, dst_vector);
                  } else {
                     if (!d_single_block_fill_local[nb].isNull()) {
                        d_single_block_fill_local[nb]->initializeDestinationVector(
                           dst_vector);
                     } else {
                        initializeDestinationVector(dst_vector);
                     }

                     const int ncomponents = d_neighbor_ghost_level[nb][nc]->
                                                getPatchDescriptor()->
                                                   getMaxNumberRegisteredComponents();

                     for (int di = 0; di < ncomponents; di++) {
                        if (dst_vector.isSet(di)) {
                           if (d_neighbor_ghost_level[nb][nc]->
                                  checkAllocated(di)) {
                              dst_vector.clrFlag(di);
                           }
                        }
                     }

                     d_neighbor_ghost_level[nb][nc]->allocatePatchData(
                        dst_vector, fill_time);
 
                  }

                  if (d_neighbor_copy_only[nb][nc]) {

                     d_multiblock_strategy->setBlockNumber(nb);
                     d_neighbor_single_block_refine_schedule[nb][nc]->
                        fillData(fill_time, false);
                     d_multiblock_strategy->clearBlockNumber();
                  } else {
                     d_multiblock_strategy->setBlockNumber(nb);
                     d_neighbor_multiblock_coarse_schedule[nb][nc]->
                        fillData(fill_time, true, true, is_scr_recursive);
                     d_multiblock_strategy->clearBlockNumber();

                     if (!d_neighbor_single_block_refine_schedule[nb][nc].isNull()) {
                        d_multiblock_strategy->setBlockNumber(nb);
                        d_neighbor_single_block_refine_schedule[nb][nc]->
                           fillData(fill_time, false);
                        d_multiblock_strategy->clearBlockNumber();
                     }

                     hier::ComponentSelector scr_vector;

                     d_neighbor_multiblock_coarse_schedule[nb][nc]->allocateScratchSpace(
                        d_neighbor_ghost_level[nb][nc],
                        fill_time,
                        scr_vector);
                        
                     d_multiblock_strategy->setBlockNumber(nb);
                     refineScratchData(d_neighbor_multiblock_coarse_level[nb][nc],
                                       d_neighbor_ghost_level[nb][nc],
                                       d_neighbor_unfilled_boxes[nb][nc],
                                       id,
                                       false);
                     d_multiblock_strategy->clearBlockNumber();

                     copyScratchToDestination(
                        d_neighbor_ghost_level[nb][nc],
                        d_neighbor_unfilled_boxes[nb][nc],
                        d_single_block_refine_alg->getEquivalenceClasses());

                     d_neighbor_ghost_level[nb][nc]->
                        deallocatePatchData(scr_vector);

                     tbox::Pointer<hier::ComponentSelector> crs_scr_vector =
                        d_neighbor_multiblock_coarse_schedule[nb][nc]->
                           getCoarseScratchVector(id);
                     d_neighbor_multiblock_coarse_level[nb][nc]->
                        deallocatePatchData(*crs_scr_vector);
                  }

                  /*
                   * get a src_vector to allocate on the finalize_ghost_level,
                   * then copy between blocks.
                   */

                  hier::ComponentSelector src_vector;

                  if (d_neighbor_single_block_refine_schedule[nb][nc].isNull()) {
                     d_neighbor_multiblock_coarse_schedule[nb][nc]->initializeSourceVector(
                                                           src_vector);
                  } else {
                     d_neighbor_single_block_refine_schedule[nb][nc]->initializeSourceVector(
                                                           src_vector);
                  }

                  d_finalize_ghost_level[nb][nc]->allocatePatchData(src_vector);

                  const tbox::Pointer< xfer::RefineClasses<DIM> >
                     equiv_classes = d_single_block_refine_alg->getEquivalenceClasses();
                  int num_classes =
                     equiv_classes->getNumberEquivalenceClasses();

                  if (d_using_standard_transaction) {
                     copyBetweenBlocks(d_finalize_ghost_level[nb][nc],
                                       d_neighbor_ghost_level[nb][nc],
                                       (ni().d_shift)*(level->getRatio()),
                                       ni().d_rotation, equiv_classes);
                  } else {
                     fillBetweenBlocks(d_finalize_ghost_level[nb][nc],
                                       d_neighbor_ghost_level[nb][nc],
                                       (ni().d_shift)*(level->getRatio()),
                                       ni().d_rotation, equiv_classes);
                  }

                  for (typename hier::PatchLevel<DIM>::Iterator p(level);
                       p; p++) {
                     if (d_finalize_ghost_patch_numbers[nb][nc][p()] >= 0) {
                        tbox::Pointer< hier::Patch<DIM> > dst_patch =
                           level->getPatch(p());

                        const int num_src_patches =
                           d_finalize_ghost_num_src_patches[nb][nc][p()];
                        for (int ns = 0; ns < num_src_patches; ns++) {
                           tbox::Pointer< hier::Patch<DIM> > src_patch =
                              d_finalize_ghost_level[nb][nc]->getPatch(
                                 d_finalize_ghost_patch_numbers[nb][nc][p()]+ns);

                           if (!ni().d_is_singularity) {
                              for (int ne = 0; ne < num_classes; ne++) {

                                 tbox::Pointer< hier::BoxOverlap<DIM> >
                                    overlap = calculateOverlap(*dst_patch,
                                                               *src_patch,
                                                               *equiv_classes,
                                                               ne); 

                                 for (typename
                                      tbox::List<typename
                                      xfer::RefineClasses<DIM>::Data>::Iterator
                                      l(equiv_classes->getIterator(ne));
                                      l; l++) {
                                    const int dst = l().d_dst;
                                    const int src = l().d_src;

                                    dst_patch->getPatchData(dst)->
                                       copy(*(src_patch->getPatchData(src)),
                                            *(overlap));
                                 }
                              }
                           } else {
                              singularity_to_fill = true;
                              tbox::Pointer< hier::Patch<DIM> > patch
                                = new hier::Patch<DIM>(
                                   src_patch->getBox(),
                                   src_patch->getPatchDescriptor());
                              for (int ne = 0; ne < num_classes; ne++) {
                                 for (typename
                                      tbox::List<typename
                                      xfer::RefineClasses<DIM>::Data>::Iterator
                                      l(equiv_classes->getIterator(ne));
                                      l; l++) {
                                    const int dst = l().d_dst;
                                    const int src = l().d_src;

                                    patch->allocatePatchData(dst);
                                    local_selector[p()].setFlag(dst);

                                    patch->getPatchData(dst)->
                                       copy(*(src_patch->getPatchData(src)));
                                 }
                              }

                              SingularityPatch sing_patch;
                              sing_patch.d_patch = patch;
                              sing_patch.d_id = id;
                              singularity_patches[p()].addItem(sing_patch); 
                           }
                        }
                     }
                  }

                  d_finalize_ghost_level[nb][nc]->
                     deallocatePatchData(src_vector);
                  d_neighbor_ghost_level[nb][nc]->
                     deallocatePatchData(dst_vector);
               }
               nc++;
            }

            if (singularity_to_fill) {
               fillSingularityBoundary(level, singularity_patches,
                                       nb, fill_time);
            }

            if (do_physical_boundary_fill) {
               const hier::IntVector<DIM> gcw = getBoundaryFillGhostWidth();
               for (typename hier::PatchLevel<DIM>::Iterator p(level);
                    p; p++) {
                  tbox::Pointer< hier::Patch<DIM> > patch =
                     level->getPatch(p());
                  if (patch->getPatchGeometry()->
                         intersectsPhysicalBoundary()) {
                     d_multiblock_strategy->setBlockNumber(nb);
                     d_multiblock_strategy->
                        setPhysicalBoundaryConditions(*patch,
                                                      fill_time,
                                                      gcw);
                     d_multiblock_strategy->clearBlockNumber();
                  }
               }
            }

            for (int j = 0; j < level->getNumberOfPatches(); j++) {
               for (typename tbox::List<SingularityPatch>::Iterator
                    sp(singularity_patches[j]);
                    sp; sp++) {
                  sp().d_patch->deallocatePatchData(local_selector[j]);
               }
            }
         }

         d_multiblock_strategy->clearBlockNumber();
      }
   }

   if (!filling_coarse_scratch) {
      for (nb = 0; nb < d_multiblock_hierarchy->getNumberBlocks(); nb++) {

         tbox::Pointer< hier::PatchLevel<DIM> > dst_patch_level=
            d_multiblock_dst_level->getPatchLevelForBlock(nb);

         if (!dst_patch_level.isNull()) {
            dst_patch_level->deallocatePatchData(allocate_scr_vector[nb]);
         }
      }
   }


   if (filling_coarse_scratch) {
      if (!filling_crse_scr_recursive) {
         if (d_multiblock_strategy != NULL) {
            d_multiblock_strategy->setFillingCoarseScratch(false);
         }
      }

   }
}

/*
 * ************************************************************************
 *                                                                        *
 * Fill the boundaries of the specified level at areas of singularity.    *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM>
void MultiblockRefineSchedule<DIM>::fillSingularityBoundary(
      tbox::Pointer< hier::PatchLevel<DIM> >& level,
      tbox::Array< tbox::List<SingularityPatch> >& singularity_patches,
      const int block_number,
      const double fill_time) const
{

   hier::IntVector<DIM> ratio = level->getRatio();

   for (typename hier::BoxList<DIM>::Iterator sb(
        d_multiblock_hierarchy->getSingularityBoxList(block_number)); sb; sb++) {
      hier::Box<DIM> singularity(sb());
      if (level->getLevelNumber() != 0) {
         singularity.refine(ratio);
      }

      hier::IntVector<DIM> gcw(
         level->getPatchDescriptor()->getMaxGhostWidth());

      if (d_multiblock_strategy != NULL) {
         for (typename hier::PatchLevel<DIM>::Iterator p(level); p; p++) {
            tbox::Pointer< hier::Patch<DIM> > patch = level->getPatch(p());
            tbox::Pointer< hier::PatchGeometry<DIM> > pgeom =
               patch->getPatchGeometry();

            tbox::Array<hier::BoundaryBox<DIM> > nboxes =
               pgeom->getNodeBoundary();

            if (nboxes.getSize()) {
               for (int bb = 0; bb < nboxes.getSize(); bb++) {
                  hier::Box<DIM> intersection = (nboxes[bb].getBox()) *
                                                singularity;
                  if (!(intersection.empty())) {
                     hier::Box<DIM> fill_box =
                        pgeom->getBoundaryFillBox(nboxes[bb],
                                                  patch->getBox(),
                                                  gcw);
                     d_multiblock_strategy->fillSingularityBoundaryConditions(
                        *patch, singularity_patches[p()],
                        fill_time, fill_box, nboxes[bb]);
                  }
               }
            }

            if (DIM == 3) {
               tbox::Array<hier::BoundaryBox<DIM> > eboxes =
                  pgeom->getEdgeBoundary();

               if (eboxes.getSize()) {
                  for (int bb = 0; bb < eboxes.getSize(); bb++) {
                     hier::Box<DIM> intersection =
                        (eboxes[bb].getBox()) * singularity;
                     if (!(intersection.empty())) {
                        hier::Box<DIM> fill_box =
                           pgeom->getBoundaryFillBox(eboxes[bb],
                                                     patch->getBox(),
                                                     gcw);
                        d_multiblock_strategy->fillSingularityBoundaryConditions(
                           *patch, singularity_patches[p()],
                           fill_time, fill_box, eboxes[bb]);
                     }
                  }
               }
            }
         }
      }
   }

}

/*
 * ************************************************************************
 *                                                                        *
 * Call the routines that will execute a copy from one block to another.  *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM>
void MultiblockRefineSchedule<DIM>::copyBetweenBlocks(
   tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
   const tbox::Pointer< hier::PatchLevel<DIM> > src_level,
   const hier::IntVector<DIM>& shift,
   const typename hier::MultiblockPatchHierarchy<DIM>::RotationIdentifier rotate,
   const tbox::Pointer< xfer::RefineClasses<DIM> > refine_classes) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst_level.isNull());
   TBOX_ASSERT(!src_level.isNull());
#endif

   const int num_equiv_classes =
      refine_classes->getNumberEquivalenceClasses();

   for (typename hier::PatchLevel<DIM>::Iterator p(src_level); p; p++) {
      tbox::Pointer< hier::Patch<DIM> > src_patch = src_level->getPatch(p());
      tbox::Pointer< hier::Patch<DIM> > dst_patch = dst_level->getPatch(p());
      for (int nc = 0; nc < num_equiv_classes; nc++) {
         for (typename
              tbox::List<typename xfer::RefineClasses<DIM>::Data>::
              Iterator l(refine_classes->getIterator(nc)); l; l++) {
            const int src = l().d_dst; //dst for src_level
            const int dst = l().d_src; //src for dst_level

            hier::MBUtilities<DIM>::translateAndCopyData(
               *dst_patch,
               dst,
               *src_patch,
               src,
               shift,
               rotate);

         }
      }
   }
}

/*
 * ************************************************************************
 *                                                                        *
 * Call the routines that will fill data from one block to another.       *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM>
void MultiblockRefineSchedule<DIM>::fillBetweenBlocks(
   tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
   const tbox::Pointer< hier::PatchLevel<DIM> > src_level,
   const hier::IntVector<DIM>& shift,
   const typename hier::MultiblockPatchHierarchy<DIM>::RotationIdentifier rotate,
   const tbox::Pointer< xfer::RefineClasses<DIM> > refine_classes) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst_level.isNull());
   TBOX_ASSERT(!src_level.isNull());
#endif

   const int num_equiv_classes =
      refine_classes->getNumberEquivalenceClasses();

   for (typename hier::PatchLevel<DIM>::Iterator p(src_level); p; p++) {
      tbox::Pointer< hier::Patch<DIM> > src_patch = src_level->getPatch(p());
      tbox::Pointer< hier::Patch<DIM> > dst_patch = dst_level->getPatch(p());
      for (int nc = 0; nc < num_equiv_classes; nc++) {
         for (typename
              tbox::List<typename xfer::RefineClasses<DIM>::Data>::
              Iterator l(refine_classes->getIterator(nc)); l; l++) {
            const int src = l().d_dst; //dst for src_level
            const int dst = l().d_src; //src for dst_level

            hier::MBUtilities<DIM>::translateAndFillData(
               *dst_patch,
               dst,
               *src_patch,
               src,
               shift,
               rotate);
         }
      }
   }
}


/*
 * ************************************************************************
 *                                                                        *
 * Initialize a component selector to contain the source data components  *
 * for this schedule.                                                     *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM>
void MultiblockRefineSchedule<DIM>::initializeSourceVector(
   hier::ComponentSelector& allocate_vector) const
{
   int nb;

   for (nb = 0; nb < d_multiblock_hierarchy->getNumberBlocks(); nb++) {

      if (!d_single_block_fill_local[nb].isNull()) {
         d_single_block_fill_local[nb]->initializeSourceVector(allocate_vector);
         return;
      }
   }

   for (nb = 0; nb < d_multiblock_hierarchy->getNumberBlocks(); nb++) {

      if (!d_multiblock_coarse_schedule[nb].isNull()) {
         d_multiblock_coarse_schedule[nb]->initializeSourceVector(allocate_vector);
         return;
      }
   }

   /*
    * Error results if this line is reached.
    */
   TBOX_ERROR("Schedules not properly constructed");
}

/*
 * ************************************************************************
 *                                                                        *
 * Get the equivalence classes from the refine algorithm.                 *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM>
const tbox::Pointer< xfer::RefineClasses<DIM> >&
MultiblockRefineSchedule<DIM>::getEquivalenceClasses() const
{
   return (d_single_block_refine_alg->getEquivalenceClasses());
}

/*
 * ************************************************************************
 *                                                                        *
 * Determine if the fillData operation requires source data from blocks   *
 * other than the one containing the data being filled                    *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM> bool
MultiblockRefineSchedule<DIM>::needOtherSourceBlocks(
   hier::BoxList<DIM>& dst_boxes,
   hier::BoxList<DIM>& src_boxes,
   hier::BoxList<DIM>& domain_outside_block) const
{
   bool needed;

   hier::BoxList<DIM> dst_grow_list(dst_boxes);
   dst_grow_list.grow(hier::IntVector<DIM>(1));

   dst_grow_list.intersectBoxes(domain_outside_block);

   if (dst_grow_list.size() == 0) {
      needed = false;
   } else {
 
      if (dst_boxes.getTotalSizeOfBoxes() > src_boxes.getTotalSizeOfBoxes()) {

         hier::BoxList<DIM> dst_remove_list(dst_boxes);
         dst_remove_list.removeIntersections(src_boxes);

         dst_remove_list.grow(hier::IntVector<DIM>(1));

         dst_remove_list.intersectBoxes(domain_outside_block);

         if (dst_remove_list.size() == 0) {

            needed = false; 

         } else {

            needed = true;

         } 
      } else {

         hier::BoxList<DIM> intersection(src_boxes);
         intersection.intersectBoxes(dst_boxes);

         if (dst_boxes.getTotalSizeOfBoxes() >
             intersection.getTotalSizeOfBoxes()) {

            hier::BoxList<DIM> dst_remove_list(dst_boxes);
            dst_remove_list.removeIntersections(src_boxes);

            dst_remove_list.grow(hier::IntVector<DIM>(1));

            dst_remove_list.intersectBoxes(domain_outside_block);

            if (dst_remove_list.size() == 0) {

               needed = false;

            } else {

               needed = true;

            }

         } else {

            dst_boxes.removeIntersections(intersection);

            if (dst_boxes.size() > 0) {

               needed = true;

            } else {

               needed = false;

            }
         }
      }
   }

   if (needed) {
      return(true);
   } else {
      return(false);
   }
}



/*
 * ************************************************************************
 *                                                                        *
 * Create a schedule that will fill data on a temporary coarse level on a *
 * neighboring block.                                                     *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM>
void MultiblockRefineSchedule<DIM>::createNeighborCoarseSchedule(
   tbox::Pointer< hier::PatchLevel<DIM> >& neighbor_coarse_level,
   int next_coarser_level,
   const hier::IntVector<DIM>& ratio_to_coarser,
   tbox::Pointer< hier::PatchHierarchy<DIM> >& hierarchy,
   int neighbor_block_number,
   int dst_block_number,
   int neighbor_counter)
{
   NULL_USE(ratio_to_coarser);
   NULL_USE(hierarchy);

   tbox::Pointer< hier::PatchLevel<DIM> > coarse_level = neighbor_coarse_level;

   hier::BoxList<DIM> domain_outside_block;
   tbox::List<typename hier::MultiblockPatchHierarchy<DIM>::Neighbor> neighbors =
      d_multiblock_hierarchy->getNeighbors(neighbor_block_number);

   for (typename tbox::List<typename hier::MultiblockPatchHierarchy<DIM>::Neighbor>::
        Iterator nei(neighbors); nei; nei++) {
      domain_outside_block.unionBoxes(nei().d_translated_domain);
   }

   domain_outside_block.refine(coarse_level->getRatio());

   hier::IntVector<DIM> gcw =
      coarse_level->getPatchDescriptor()->getMaxGhostWidth();

   tbox::Array< tbox::Pointer< hier::PatchLevel<DIM> > >
      coarse_level_array(d_multiblock_hierarchy->getNumberBlocks());

   for (int nb = 0; nb < d_multiblock_hierarchy->getNumberBlocks(); nb++) {
      coarse_level_array[nb].setNull();
   }

   coarse_level_array[neighbor_block_number] = coarse_level;

   d_neighbor_multiblock_coarse_level[dst_block_number][neighbor_counter] =
      new hier::MultiblockPatchLevel<DIM>(coarse_level_array);

   d_multiblock_hierarchy->adjustMultiblockPatchLevelBoundaries(
      d_neighbor_multiblock_coarse_level[dst_block_number][neighbor_counter]);

   tbox::Pointer< hier::MultiblockPatchLevel<DIM> > coarse_hierarchy_level =
      d_multiblock_hierarchy->getPatchLevel(next_coarser_level);

   hier::BoxList<DIM> pseudo_domain(domain_outside_block);
   pseudo_domain.unionBoxes(coarse_level->getPhysicalDomain());

   hier::BoxList<DIM> unfilled_boxes(coarse_level->getBoxes());

   findUnfilledBoxes(unfilled_boxes,
                     neighbor_block_number,
                     coarse_hierarchy_level,
                     pseudo_domain,
                     gcw);

   if (unfilled_boxes.size() == 0) {
      d_neighbor_multiblock_coarse_schedule[dst_block_number][neighbor_counter] =
         new MultiblockRefineSchedule<DIM>(
            d_fill_pattern,
            d_neighbor_multiblock_coarse_level[dst_block_number][neighbor_counter],
            coarse_hierarchy_level,
            d_multiblock_hierarchy,
            d_single_block_scratch_refine_alg,
            d_transaction_factory,
            d_multiblock_strategy,
            true);
   } else {
      d_neighbor_multiblock_coarse_schedule[dst_block_number][neighbor_counter] =
         new MultiblockRefineSchedule<DIM>(
            d_fill_pattern,
            d_neighbor_multiblock_coarse_level[dst_block_number][neighbor_counter],
            coarse_hierarchy_level,
            next_coarser_level-1,
            d_multiblock_hierarchy,
            d_single_block_scratch_refine_alg,
            d_transaction_factory,
            d_multiblock_strategy,
            true);

   }
}



/*
 * ************************************************************************
 *                                                                        *
 * Create a schedule to fill data on a temporary coarse level             *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM>
void MultiblockRefineSchedule<DIM>::createCoarseSchedule(
   tbox::Pointer< hier::PatchLevel<DIM> >& coarse_level,
   int next_coarser_level,
   const hier::IntVector<DIM>& ratio_to_coarser,
   tbox::Pointer< hier::PatchHierarchy<DIM> >& hierarchy,
   int block_number)
{
   NULL_USE(ratio_to_coarser);
   NULL_USE(hierarchy);

   hier::BoxList<DIM> domain_outside_block;
   tbox::List<typename hier::MultiblockPatchHierarchy<DIM>::Neighbor> neighbors =
      d_multiblock_hierarchy->getNeighbors(block_number);

   for (typename tbox::List<typename hier::MultiblockPatchHierarchy<DIM>::Neighbor>::
        Iterator nei(neighbors); nei; nei++) {
      domain_outside_block.unionBoxes(nei().d_translated_domain);
   }

   domain_outside_block.refine(coarse_level->getRatio());

   hier::IntVector<DIM> gcw =
      coarse_level->getPatchDescriptor()->getMaxGhostWidth();

   tbox::Array< tbox::Pointer< hier::PatchLevel<DIM> > >
      coarse_level_array(d_multiblock_hierarchy->getNumberBlocks());

   for (int nb = 0; nb < d_multiblock_hierarchy->getNumberBlocks(); nb++) {
      coarse_level_array[nb].setNull();
   }

   coarse_level_array[block_number] = coarse_level;

   d_multiblock_coarse_scratch_level[block_number] =
      new hier::MultiblockPatchLevel<DIM>(coarse_level_array); 

   d_multiblock_hierarchy->adjustMultiblockPatchLevelBoundaries(
      d_multiblock_coarse_scratch_level[block_number]);

   tbox::Pointer< hier::MultiblockPatchLevel<DIM> > coarse_hierarchy_level =
      d_multiblock_hierarchy->getPatchLevel(next_coarser_level);

   hier::BoxList<DIM> pseudo_domain(domain_outside_block);
   pseudo_domain.unionBoxes(coarse_level->getPhysicalDomain()); 

   hier::BoxList<DIM> unfilled_boxes(coarse_level->getBoxes());

   findUnfilledBoxes(unfilled_boxes,
                     block_number,
                     coarse_hierarchy_level,
                     pseudo_domain,
                     gcw);

   if (unfilled_boxes.size() == 0) {
      d_multiblock_coarse_schedule[block_number] =
         new MultiblockRefineSchedule<DIM>(
            d_fill_pattern,
            d_multiblock_coarse_scratch_level[block_number],
            coarse_hierarchy_level,
            d_multiblock_hierarchy,
            d_single_block_scratch_refine_alg,
            d_transaction_factory,
            d_multiblock_strategy,
            true);
   } else {
      d_multiblock_coarse_schedule[block_number] =
         new MultiblockRefineSchedule<DIM>(
            d_fill_pattern,
            d_multiblock_coarse_scratch_level[block_number],
            coarse_hierarchy_level,
            next_coarser_level-1,
            d_multiblock_hierarchy,
            d_single_block_scratch_refine_alg,
            d_transaction_factory,
            d_multiblock_strategy,
            true);
   }
}


/*
 * ************************************************************************
 *                                                                        *
 * Find the boxes that are unfilled by copying and communicating data     *
 * from other patches of the same resolution as the destination level     *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM>
void MultiblockRefineSchedule<DIM>::findUnfilledBoxes(
   hier::BoxList<DIM>& unfilled_boxes,
   const int block_number,
   tbox::Pointer< hier::MultiblockPatchLevel<DIM> > coarse_hierarchy_level,
   const hier::BoxList<DIM>& pseudo_domain,
   const hier::IntVector<DIM>& gcw)
{
   unfilled_boxes.grow(gcw);
   unfilled_boxes.intersectBoxes(pseudo_domain);

   for (int nb = 0; nb < d_multiblock_hierarchy->getNumberBlocks(); nb++) {

      tbox::Pointer< hier::PatchLevel<DIM> > patch_level=
         coarse_hierarchy_level->getPatchLevelForBlock(nb);

      if (nb == block_number) {

         if (!patch_level.isNull()) {

            unfilled_boxes.removeIntersections(patch_level->getBoxes());

         }

      } else {

         if (!patch_level.isNull()) {

            hier::BoxArray<DIM> level_boxes(patch_level->getBoxes());

            d_multiblock_hierarchy->translateBoxArray(level_boxes,
                                            patch_level->getRatio(),
                                            block_number,
                                            nb);

            unfilled_boxes.removeIntersections(level_boxes);
         }

      }

      if (unfilled_boxes.size() == 0) {
         break;
      }

   }

}


/*
 * ************************************************************************
 *                                                                        *
 * Construct a refine algorithm whose destination components are the      *
 * scratch components of the algorithm that created this                  *
 * MultiblockRefineSchedule.                                              *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM>
void MultiblockRefineSchedule<DIM>::constructScratchRefineAlgorithm()
{
   d_single_block_scratch_refine_alg = new xfer::RefineAlgorithm<DIM>();

   const tbox::Pointer< xfer::RefineClasses<DIM> > refine_classes =
      d_single_block_refine_alg->getEquivalenceClasses();

   tbox::Pointer< xfer::RefineClasses<DIM> > scratch_refine_classes =
      new xfer::RefineClasses<DIM>();

   int num_classes = refine_classes->getNumberEquivalenceClasses();

   for (int ne = 0; ne < num_classes; ne++) {
      for (typename tbox::List<typename
           xfer::RefineClasses<DIM>::Data>::
           Iterator rc(refine_classes->getIterator(ne)); rc; rc++) {

         typename xfer::RefineClasses<DIM>::Data data = rc();

         data.d_dst = data.d_scratch;

         scratch_refine_classes->insertEquivalenceClassItem(data);

      }
   }

   d_single_block_scratch_refine_alg->setEquivalenceClasses(scratch_refine_classes);

}


/*
 * ************************************************************************
 *                                                                        *
 * Copy data from scratch to destination where there is overlap between   *
 * the data.  Nothing is done if scratch and destination are the same.    *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM>
void
MultiblockRefineSchedule<DIM>::copyScratchToDestination(
   tbox::Pointer< hier::PatchLevel<DIM> > level,
   const hier::BoxList<DIM>& unfilled_boxes,
   tbox::Pointer< xfer::RefineClasses<DIM> > refine_classes) const
{
   const int num_classes = refine_classes->getNumberEquivalenceClasses();

   tbox::Pointer< hier::PatchDescriptor<DIM> > descriptor =
      level->getPatchDescriptor();

   bool copy_needed = false;
   for (int nc = 0; nc < num_classes; nc++) {
      for (typename tbox::List<typename
           xfer::RefineClasses<DIM>::Data>::
           Iterator rc(refine_classes->getIterator(nc)); rc; rc++) {
         const int dst = rc().d_dst;
         const int src = rc().d_scratch;

         if (dst != src) {
            copy_needed = true;
            break;
         }
      }
      if (copy_needed) {
         break;
      }
   }

   if (copy_needed) {
      for (typename hier::PatchLevel<DIM>::Iterator p(level); p; p++) {

         tbox::Pointer< hier::Patch<DIM> > patch = level->getPatch(p());
         hier::Box<DIM> patch_box(patch->getBox());

         for (int ne = 0; ne < num_classes; ne++) {

            const typename xfer::RefineClasses<DIM>::Data& rep_item =
               refine_classes->getClassRepresentative(ne); 

            const int rep_item_src_id = rep_item.d_scratch;
            const int rep_item_dst_id = rep_item.d_scratch;

            tbox::Pointer< hier::PatchDataFactory<DIM> > src_pdf =
               descriptor->getPatchDataFactory(rep_item_src_id);
            tbox::Pointer< hier::PatchDataFactory<DIM> > dst_pdf =
               descriptor->getPatchDataFactory(rep_item_dst_id); 

            for (typename hier::BoxList<DIM>::Iterator b(unfilled_boxes);
                 b; b++) {

               const hier::Box<DIM> fill_box(b());

               const hier::Box<DIM> src_mask(fill_box);

               if (!src_mask.empty()) {
                  tbox::Pointer< hier::BoxOverlap<DIM> > overlap =
                     dst_pdf->getBoxGeometry(fill_box)
                            ->calculateOverlap(
                               *src_pdf->getBoxGeometry(patch_box),
                               src_mask,
                               true, hier::IntVector<DIM>(0));

                  for (typename
                       tbox::List<typename xfer::RefineClasses<DIM>::Data>::
                       Iterator l(refine_classes->getIterator(ne)); l; l++) {
                     const int dst = l().d_dst;
                     const int src = l().d_scratch; 

                     if (dst != src) {

                        patch->getPatchData(dst)->
                           copy(*(patch->getPatchData(src)), *overlap);
                     }
                  }
               }
            }
         }
      }
   }
}


/*
 * ************************************************************************
 *                                                                        *
 * Execute refinement operations on scracth data using the preprocess and *
 * postprocess refine routines.                                           *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM>
void MultiblockRefineSchedule<DIM>::refineScratchData(
   tbox::Pointer< hier::MultiblockPatchLevel<DIM> > coarse_mb_level,
   tbox::Pointer< hier::PatchLevel<DIM> > fine_level,
   const hier::BoxList<DIM>& unfilled_boxes,
   const int block_number,
   const bool add_fine_gcw) const
{
   tbox::Pointer< hier::PatchLevel<DIM> > coarse_level =
      coarse_mb_level->getPatchLevelForBlock(block_number);

   const hier::IntVector<DIM> ratio =
      fine_level->getRatio() / coarse_level->getRatio();

   hier::IntVector<DIM> gcw(0);

   if (add_fine_gcw) {
      gcw = fine_level->getPatchDescriptor()->getMaxGhostWidth();
   }

   for (typename hier::PatchLevel<DIM>::Iterator p(coarse_level); p; p++) {
      tbox::Pointer< hier::Patch<DIM> > fine_patch =
         fine_level->getPatch(p());
      tbox::Pointer< hier::Patch<DIM> > crse_patch =
         coarse_level->getPatch(p());

      hier::Box<DIM> fine_box(fine_patch->getBox());
      fine_box.grow(gcw);

      hier::BoxList<DIM> fill_boxes(unfilled_boxes);
      fill_boxes.intersectBoxes(fine_box);

      fill_boxes.coalesceBoxes();

      if (d_multiblock_strategy != NULL) {
         d_multiblock_strategy->preprocessRefineBoxes(*fine_patch,
                                           *crse_patch,
                                           fill_boxes,
                                           ratio);
      }

      if (d_multiblock_strategy != NULL) {
         d_multiblock_strategy->postprocessRefineBoxes(*fine_patch,
                                            *crse_patch,
                                            fill_boxes,
                                            ratio);

      }
   }
}


/*
 * ************************************************************************
 *                                                                        *
 * Calculate the overlap between data on two patches for a single         *
 * equivalence class.                                                     *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM>
tbox::Pointer< hier::BoxOverlap<DIM> >
MultiblockRefineSchedule<DIM>::calculateOverlap(
   const hier::Patch<DIM>& dst_patch,
   const hier::Patch<DIM>& src_patch,
   const xfer::RefineClasses<DIM>& refine_classes,
   const int refine_class_id) const
{
   tbox::Pointer< hier::PatchDescriptor<DIM> > dst_patch_descriptor =
      dst_patch.getPatchDescriptor();
   tbox::Pointer< hier::PatchDescriptor<DIM> > src_patch_descriptor =
      src_patch.getPatchDescriptor();

   const typename xfer::RefineClasses<DIM>::Data& rep_item =
      refine_classes.getClassRepresentative(refine_class_id);

   const int rep_item_dst_id = rep_item.d_dst;
   const int rep_item_src_id = rep_item.d_src;

   tbox::Pointer< hier::PatchDataFactory<DIM> > src_pdf =
      src_patch_descriptor->getPatchDataFactory(rep_item_src_id);
   tbox::Pointer< hier::PatchDataFactory<DIM> > dst_pdf =
      dst_patch_descriptor->getPatchDataFactory(rep_item_dst_id);

   hier::Box<DIM> dst_fill_box(dst_patch.getPatchData(rep_item_dst_id)->
                          getGhostBox() * src_patch.getBox());

   return (dst_pdf->getBoxGeometry(dst_fill_box)
                  ->calculateOverlap(
                     *src_pdf->getBoxGeometry(src_patch.getBox()),
                     dst_fill_box,
                     true, hier::IntVector<DIM>(0)));

}




/*
 * ************************************************************************
 *                                                                        *
 * Allocate scratch space on a patch level.                               *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM>
void MultiblockRefineSchedule<DIM>::allocateScratchSpace(
   tbox::Pointer< hier::PatchLevel<DIM> > level,
   double fill_time,
   hier::ComponentSelector& allocate_vector) const
{

   int nb;

   for (nb = 0; nb < d_multiblock_hierarchy->getNumberBlocks(); nb++) {

      if (!d_single_block_fill_local[nb].isNull()) {
         d_single_block_fill_local[nb]->allocateScratchSpace(
            level, fill_time, allocate_vector);
         return;
      }
   }

   for (nb = 0; nb < d_multiblock_hierarchy->getNumberBlocks(); nb++) {

      if (!d_multiblock_coarse_schedule[nb].isNull()) {
         d_multiblock_coarse_schedule[nb]->allocateScratchSpace(
            level, fill_time, allocate_vector);
         return;
      }
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(0);
#endif

   return;
}


/*
 * ************************************************************************
 *                                                                        *
 * Get the ghost width needed for boundary filling.                       *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM>
hier::IntVector<DIM>
MultiblockRefineSchedule<DIM>::getBoundaryFillGhostWidth() const
{

   int nb; 

   for (nb = 0; nb < d_multiblock_hierarchy->getNumberBlocks(); nb++) {
      if (!d_single_block_fill_local[nb].isNull()) {
         return (d_single_block_fill_local[nb]->getBoundaryFillGhostWidth());
      }
   }

   for (nb = 0; nb < d_multiblock_hierarchy->getNumberBlocks(); nb++) {
      if (!d_multiblock_coarse_schedule[nb].isNull()) {
         return (d_multiblock_coarse_schedule[nb]->getBoundaryFillGhostWidth());
      }
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(0);
#endif

   return (hier::IntVector<DIM>(-1));
}

/*
 * ************************************************************************
 *                                                                        *
 * Get a pointer a ComponentSelector for scratch components on the coarse *
 * level.                                                                 *
 *                                                                        *
 * ************************************************************************
 */
template<int DIM>
tbox::Pointer<hier::ComponentSelector>&
MultiblockRefineSchedule<DIM>::getCoarseScratchVector(
   const int block_num)
{
   return(d_coarse_selector[block_num]);
}

/*
 * ************************************************************************
 *                                                                        *
 * Initialize a ComponentSelector to store all destination components for *
 * this schedule.                                                         * 
 *                                                                        *
 * ************************************************************************
 */

template<int DIM>
void MultiblockRefineSchedule<DIM>::initializeDestinationVector(
   hier::ComponentSelector& dst_vector) const
{
   dst_vector.clrAllFlags();

   const tbox::Pointer< xfer::RefineClasses<DIM> > refine_classes =
      d_single_block_refine_alg->getEquivalenceClasses();

   int num_classes = refine_classes->getNumberEquivalenceClasses();

   for (int ne = 0; ne < num_classes; ne++) {
      for (typename
           tbox::List<typename xfer::RefineClasses<DIM>::Data>::
           Iterator rc(refine_classes->getIterator(ne)); rc; rc++) {

         typename xfer::RefineClasses<DIM>::Data data = rc();

	 // Intel compiler was issuing a warning about data not being used.
	 // Which is odd.
	 NULL_USE(data);

         dst_vector.setFlag(data.d_dst);

      }
   }
}


}
}
#endif
