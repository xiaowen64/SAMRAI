/************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Base class for geometry management on patches 
 *
 ************************************************************************/

#ifndef included_xfer_MultiblockRefineSchedule_C
#define included_xfer_MultiblockRefineSchedule_C

#include "SAMRAI/xfer/MultiblockRefineSchedule.h"

#include "SAMRAI/hier/MappedBoxSetSingleBlockIterator.h"
#include "SAMRAI/hier/OverlapConnectorAlgorithm.h"
#include "SAMRAI/hier/PatchGeometry.h"
#include "SAMRAI/hier/RealMappedBoxConstIterator.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/xfer/MultiblockRefinePatchStrategy.h"
#include "SAMRAI/xfer/PatchLevelFullFillPattern.h"
#include "SAMRAI/xfer/StandardRefineTransactionFactory.h"

namespace SAMRAI {
namespace xfer {

const int MultiblockRefineSchedule::MULTIBLOCK_FAKE_LEVEL_NUMBER = -24;

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

MultiblockRefineSchedule::MultiblockRefineSchedule(
   tbox::Pointer<PatchLevelFillPattern> fill_pattern,
   tbox::Pointer<hier::PatchLevel> dst_level,
   tbox::Pointer<hier::PatchLevel> src_level,
   tbox::Pointer<hier::PatchHierarchy> hierarchy,
   tbox::Pointer<xfer::RefineAlgorithm> refine_alg,
   tbox::Pointer<xfer::RefineTransactionFactory> transaction_factory,
   MultiblockRefinePatchStrategy* strategy,
   bool use_time_refinement,
   bool enable_singularity_patches):
   d_multiblock_hierarchy(hierarchy),
   d_fill_pattern(fill_pattern),
   d_data_fill_gcw(dst_level->getDim()),
   d_enable_singularity_patches(enable_singularity_patches)
{
#ifdef DEBUG_CHECK_DIM_ASSERTIONS
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst_level, *src_level, *hierarchy);
   if (refine_alg) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*dst_level, *refine_alg);
   }
   if (strategy) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*dst_level, *strategy);
   }
#endif

   const int num_blocks = hierarchy->getGridGeometry()->getNumberBlocks();

   d_single_block_refine_alg = refine_alg;

   d_transaction_factory = transaction_factory;

   d_using_standard_transaction = false;
   tbox::Pointer<xfer::StandardRefineTransactionFactory> t_factory =
      d_transaction_factory;

   if (!(t_factory.isNull())) {
      d_using_standard_transaction = true;
   }

   d_single_block_scratch_refine_alg.setNull();

   d_multiblock_dst_level = dst_level;

   d_unfilled_boxes.resizeArray(num_blocks);
   d_local_refine_overlaps.resize(num_blocks);

   d_multiblock_strategy = strategy;

   d_finalize_ghost_patch_numbers.resizeArray(num_blocks);
   d_finalize_ghost_num_src_patches.resizeArray(num_blocks);

   d_neighbor_unfilled_boxes.resizeArray(num_blocks);
   d_neighbor_copy_overlaps.resizeArray(num_blocks);
   d_neighbor_refine_overlaps.resize(num_blocks);

   d_only_copy_source_intersections = true;

   const tbox::Pointer<RefineClasses>& refine_classes =
      d_single_block_refine_alg->getEquivalenceClasses();
   const int num_classes = refine_classes->getNumberOfEquivalenceClasses();

   tbox::Pointer<hier::PatchDescriptor> descriptor =
      hier::VariableDatabase::getDatabase()->getPatchDescriptor();

   const tbox::Dimension& dim = dst_level->getDim();
   hier::IntVector desc_gcw(descriptor->getMaxGhostWidth(dim));
   hier::IntVector gcw(hier::IntVector::getOne(dim));

   for (int ne = 0; ne < num_classes; ne++) {

      const RefineClasses::Data& rep_item =
         refine_classes->getClassRepresentative(ne);

      const tbox::Pointer<VariableFillPattern>& var_fill_pattern =
         rep_item.d_var_fill_pattern;

      if (var_fill_pattern->getPatternName() !=
          "BOX_GEOMETRY_FILL_PATTERN") {
         gcw = hier::IntVector::max(gcw,
               var_fill_pattern->getStencilWidth());
      } else {
         gcw = hier::IntVector::max(gcw, desc_gcw);
      }
   }

   d_data_fill_gcw = gcw;

   bool encon_only = d_fill_pattern->fillingEnhancedConnectivityOnly();

   dst_level->getBoxes();
   const hier::MappedBoxSet& dst_global_mapped_boxes =
      dst_level->getMappedBoxLevel()->getGlobalizedVersion().getGlobalMappedBoxes();
   hier::MappedBoxSet src_mapped_boxes;
   hier::MappedBoxSet src_global_mapped_boxes;

   if (!(src_level.isNull())) {
      src_level->getBoxes();
      src_mapped_boxes = src_level->getMappedBoxLevel()->getMappedBoxes();
      src_global_mapped_boxes = src_level->getMappedBoxLevel()->getGlobalizedVersion().getGlobalMappedBoxes();
   }

   d_coarse_selector = new hier::ComponentSelector();
   d_coarse_selector->clrAllFlags();

   if (!encon_only && dst_global_mapped_boxes.size()) {

      if (src_global_mapped_boxes.size()) {

         /*
          * This private constructor is written for same src and
          * dst refinement ratios.  It is not meant to be used when
          * they are not the same.
          */
         TBOX_ASSERT(
            src_level->getRatioToLevelZero() ==
            dst_level->getRatioToLevelZero());

         /*
          * dst_patch_level does not have a real level number, but
          * note that it is at the same resolution as the
          * src_patch_level.
          */

         d_single_block_fill_local =
            d_single_block_refine_alg->createSchedule(
               d_fill_pattern,
               dst_level,
               src_level,
               (xfer::RefinePatchStrategy *)strategy,
               hier::BlockId::invalidId(),
               use_time_refinement,
               d_transaction_factory);
      } else {
         d_single_block_fill_local.setNull();
      }
   } else {
      d_single_block_fill_local.setNull();
   }
   d_local_fill_only = true;

   if (strategy != NULL) {
      createInterblockSchedules(dst_level,
         src_level,
         (xfer::RefinePatchStrategy *)strategy,
         MULTIBLOCK_FAKE_LEVEL_NUMBER,
         use_time_refinement);

      createOverlapsForRefine(dst_level, src_level);
      createOverlapsForCopy(dst_level, src_level);
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

MultiblockRefineSchedule::MultiblockRefineSchedule(
   tbox::Pointer<PatchLevelFillPattern> fill_pattern,
   tbox::Pointer<hier::PatchLevel> dst_level,
   tbox::Pointer<hier::PatchLevel> src_level,
   const int next_coarser_level,
   tbox::Pointer<hier::PatchHierarchy> hierarchy,
   tbox::Pointer<xfer::RefineAlgorithm> refine_alg,
   tbox::Pointer<xfer::RefineTransactionFactory> transaction_factory,
   MultiblockRefinePatchStrategy* strategy,
   bool use_time_refinement,
   bool enable_singularity_patches):
   d_multiblock_hierarchy(hierarchy),
   d_fill_pattern(fill_pattern),
   d_data_fill_gcw(dst_level->getDim()),
   d_enable_singularity_patches(enable_singularity_patches)
{
#ifdef DEBUG_CHECK_DIM_ASSERTIONS
   TBOX_DIM_ASSERT_CHECK_ARGS2(*dst_level, *strategy);
   if (!src_level.isNull()) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*dst_level, *src_level);
   }
   if (!hierarchy.isNull()) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*dst_level, *hierarchy);
   }
   if (refine_alg) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*dst_level, *refine_alg);
   }
#endif

   const tbox::SAMRAI_MPI &mpi(d_multiblock_hierarchy->getMPI());

   const int num_blocks = hierarchy->getGridGeometry()->getNumberBlocks();

   if (strategy == NULL && num_blocks > 1) {
      TBOX_ERROR("A non-null MultiblockRefinePatchStrategy is required" << std::endl
                 << "when the number of blocks is greater than 1.");
   }

   const tbox::Dimension& dim(d_multiblock_hierarchy->getDim());

   d_single_block_refine_alg = refine_alg;

   d_transaction_factory = transaction_factory;

   d_using_standard_transaction = false;
   tbox::Pointer<xfer::StandardRefineTransactionFactory> t_factory =
      d_transaction_factory;

   if (!(t_factory.isNull())) {
      d_using_standard_transaction = true;
   }

   d_single_block_scratch_refine_alg.setNull();

   constructScratchRefineAlgorithm();

   d_multiblock_dst_level = dst_level;

   d_multiblock_strategy = strategy;

   d_unfilled_boxes.resizeArray(num_blocks);
   d_local_refine_overlaps.resize(num_blocks);

   d_finalize_ghost_patch_numbers.resizeArray(num_blocks);
   d_finalize_ghost_num_src_patches.resizeArray(num_blocks);

   d_neighbor_unfilled_boxes.resizeArray(num_blocks);
   d_neighbor_copy_overlaps.resizeArray(num_blocks);
   d_neighbor_refine_overlaps.resize(num_blocks);

   d_only_copy_source_intersections = false;

   tbox::Pointer<hier::GridGeometry> grid_geometry(
      d_multiblock_hierarchy->getGridGeometry());

   const tbox::Pointer<RefineClasses>& refine_classes =
      d_single_block_refine_alg->getEquivalenceClasses();
   const int num_classes = refine_classes->getNumberOfEquivalenceClasses();

   tbox::Pointer<hier::PatchDescriptor> descriptor =
      hier::VariableDatabase::getDatabase()->getPatchDescriptor();

   hier::IntVector desc_gcw(descriptor->getMaxGhostWidth(dim));
   hier::IntVector gcw(hier::IntVector::getOne(dim));

   for (int ne = 0; ne < num_classes; ne++) {

      const RefineClasses::Data& rep_item =
         refine_classes->getClassRepresentative(ne);

      const tbox::Pointer<VariableFillPattern>& var_fill_pattern =
         rep_item.d_var_fill_pattern;

      if (var_fill_pattern->getPatternName() !=
          "BOX_GEOMETRY_FILL_PATTERN") {
         gcw = hier::IntVector::max(gcw,
               var_fill_pattern->getStencilWidth());
      } else {
         gcw = hier::IntVector::max(gcw, desc_gcw);
      }
   }

   d_data_fill_gcw = gcw;

   dst_level->getBoxes();
   const hier::MappedBoxSet& dst_global_mapped_boxes =
      dst_level->getMappedBoxLevel()->getGlobalizedVersion().getGlobalMappedBoxes();
   hier::MappedBoxSet src_mapped_boxes;
   hier::MappedBoxSet src_global_mapped_boxes;

   if (!(src_level.isNull())) {
      src_level->getBoxes();
      src_mapped_boxes = src_level->getMappedBoxLevel()->getMappedBoxes();
      src_global_mapped_boxes = src_level->getMappedBoxLevel()->getGlobalizedVersion().getGlobalMappedBoxes();
   }

   bool encon_only = d_fill_pattern->fillingEnhancedConnectivityOnly();

   d_coarse_selector = new hier::ComponentSelector();
   d_coarse_selector->clrAllFlags();

   d_local_fill_only = true;

   if (!encon_only && dst_global_mapped_boxes.size()) {

      if ((num_blocks > 1 && next_coarser_level >= 0)) {

         tbox::Array<hier::BoxList> dst_level_physdomain_boxlist(num_blocks);
         for (int nb = 0; nb < num_blocks; nb++) {
            hier::BlockId block_id(nb);

            dst_level_physdomain_boxlist[nb] =
               hier::BoxList(dst_level->getPhysicalDomain(block_id));
         }

         hier::OverlapConnectorAlgorithm connect_util;

         d_local_fill_only = false;

         for (int nb = 0; nb < num_blocks; nb++) {
            d_unfilled_boxes[nb].clearItems();
         }

         if (src_global_mapped_boxes.size()) {

            const hier::IntVector connector_gcw =
               hierarchy->getRequiredConnectorWidth(src_level->getLevelNumber(),
                  src_level->getLevelNumber());
            const hier::Connector& dst_to_dst =
               dst_level->getMappedBoxLevel()->
               getPersistentOverlapConnectors().findOrCreateConnector(
                  *dst_level->getMappedBoxLevel(),
                  connector_gcw);

            const hier::Connector& dst_to_src =
               dst_level->getMappedBoxLevel()->
               getPersistentOverlapConnectors().findOrCreateConnector(
                  *src_level->getMappedBoxLevel(),
                  connector_gcw);

            const hier::Connector& src_to_dst =
               src_level->getMappedBoxLevel()->
               getPersistentOverlapConnectors().findOrCreateConnector(
                  *dst_level->getMappedBoxLevel(),
                  connector_gcw);

            hier::MappedBoxSet fill_mapped_boxes;
            hier::NeighborhoodSet dst_to_fill_edges;
            d_fill_pattern->computeFillMappedBoxesAndNeighborhoodSets(
               fill_mapped_boxes,
               dst_to_fill_edges,
               *(dst_level->getMappedBoxLevel()),
               dst_to_dst,
               dst_to_src,
               src_to_dst,
               gcw);

            tbox::Array<hier::BoxList> local_fill_boxes(num_blocks);
            tbox::Array<hier::BoxList> nabr_src_boxes(num_blocks);
            for (hier::MappedBoxSet::const_iterator ni =
                    dst_level->getMappedBoxLevel()->getMappedBoxes().begin();
                 ni != dst_level->getMappedBoxLevel()->getMappedBoxes().end();
                 ++ni) {
               const hier::MappedBoxId& mbid = ni->getId();
               const hier::BlockId& block_id = mbid.getBlockId();
               const int block_num = block_id.getBlockValue();

               const hier::MappedBoxSet& dst_to_fill_nabrs =
                  dst_to_fill_edges[mbid];

               for (hier::RealMappedBoxConstIterator di(dst_to_fill_nabrs);
                    di.isValid(); ++di) {
                  local_fill_boxes[block_num].appendItem((*di).getBox());
               }

               if (dst_to_src.hasNeighborSet(mbid)) {
                  const hier::MappedBoxSet& src_nabrs =
                     dst_to_src.getNeighborSet(mbid);
                  for (hier::MappedBoxSet::const_iterator si = src_nabrs.begin();
                       si != src_nabrs.end(); ++si) {
                     if ((*si).getBlockId() == block_id) {  
                        nabr_src_boxes[block_num].appendItem((*si).getBox());
                     }
                  }
               }
            }

            for (int nb = 0; nb < num_blocks; nb++) {
               local_fill_boxes[nb].removeIntersections(nabr_src_boxes[nb]);

               d_unfilled_boxes[nb].unionBoxes(local_fill_boxes[nb]);
            }

         } else {

            hier::IntVector connector_gcw(gcw);
            if (dst_level->getLevelNumber() >= 0) {
               connector_gcw =
                  hierarchy->getRequiredConnectorWidth(
                     dst_level->getLevelNumber(),
                     dst_level->getLevelNumber());
            }

            const hier::Connector& dst_to_dst =
               dst_level->getMappedBoxLevel()->
               getPersistentOverlapConnectors().findOrCreateConnector(
                  *dst_level->getMappedBoxLevel(),
                  connector_gcw);

            // Empty connectors since there is no source level
            hier::Connector src_to_dst;
            hier::Connector dst_to_src;

            hier::MappedBoxSet fill_mapped_boxes;
            hier::NeighborhoodSet dst_to_fill_edges;
            d_fill_pattern->computeFillMappedBoxesAndNeighborhoodSets(
               fill_mapped_boxes,
               dst_to_fill_edges,
               *(dst_level->getMappedBoxLevel()),
               dst_to_dst,
               dst_to_src,
               src_to_dst,
               gcw);

            for (hier::MappedBoxSet::const_iterator ni =
                    dst_level->getMappedBoxLevel()->getMappedBoxes().begin();
                 ni != dst_level->getMappedBoxLevel()->getMappedBoxes().end();
                 ++ni) {
               const hier::MappedBoxId& mbid = ni->getId();
               const hier::BlockId& block_id = mbid.getBlockId();
               const int block_num = block_id.getBlockValue();

               const hier::MappedBoxSet& dst_to_fill_nabrs =
                  dst_to_fill_edges[mbid];

               for (hier::RealMappedBoxConstIterator di(dst_to_fill_nabrs);
                    di.isValid(); ++di) {
                  d_unfilled_boxes[block_num].appendItem((*di).getBox());
               }
            }
         }

         for (int nb = 0; nb < num_blocks; nb++) {
            d_unfilled_boxes[nb].intersectBoxes(dst_level_physdomain_boxlist[nb]);
            d_unfilled_boxes[nb].coalesceBoxes();
         }

         std::vector<hier::MappedBox> dst_mapped_boxes;
         std::vector<int> dst_box_for_crs_box;

         const int rank = dst_level->getMappedBoxLevel()->getMPI().getRank();

         tbox::Array<hier::BoxList> tc_boxes(num_blocks);
         hier::MappedBoxSet tmp_crs_mapped_box_set;
         int num_local_dst = 0;

         for (int nb = 0; nb < num_blocks; nb++) {
            hier::BlockId block_id(nb);

            if (d_unfilled_boxes[nb].size()) {
               hier::BlockId block_id(nb);

               for (hier::MappedBoxSetSingleBlockIterator ni(dst_level->getMappedBoxLevel()->
                       getMappedBoxes(), block_id);
                    ni.isValid(); ++ni) {
                  const hier::MappedBoxId& mapped_box_id = ni->getId();
                  tbox::Pointer<hier::Patch> patch(
                     dst_level->getPatch(mapped_box_id));

                  hier::Box dst_box(patch->getBox());
                  dst_box.grow(gcw);
                  hier::BoxList dst_box_list(dst_box);
                  dst_box_list.intersectBoxes(dst_level_physdomain_boxlist[nb]);
                  dst_box_list.intersectBoxes(d_unfilled_boxes[nb]);
                  dst_box_list.coalesceBoxes();

                  dst_mapped_boxes.push_back(patch->getMappedBox());

                  for (hier::BoxList::Iterator
                       db(dst_box_list); db; db++) {
                     tc_boxes[nb].appendItem(db());
                     dst_box_for_crs_box.push_back(num_local_dst);
                  }
                  num_local_dst++;
               }
            }

            if (tc_boxes[nb].size()) {
               hier::BoxList tmp_crs_boxes(tc_boxes[nb]);

               hier::IntVector coarse_ratio(
                  dst_level->getRatioToCoarserLevel());
               tbox::Pointer<hier::PatchLevel> coarser_mblk_level =
                  d_multiblock_hierarchy->getPatchLevel(next_coarser_level);
               if (coarse_ratio == hier::IntVector::getZero(dim)) {
                  coarse_ratio = dst_level->getRatioToLevelZero()
                     / coarser_mblk_level->getRatioToLevelZero();
               }

               tmp_crs_boxes.coarsen(coarse_ratio);
               hier::BoxList tmp_crs_list(tmp_crs_boxes);
   
               hier::LocalId coarse_index(0);
               for (hier::BoxList::Iterator ts(tmp_crs_list); ts; ts++) {
                  hier::MappedBox local_mapped_box(*ts,
                                                   coarse_index++,
                                                   rank,
                                                   block_id);

                  tmp_crs_mapped_box_set.insert(tmp_crs_mapped_box_set.end(),
                     local_mapped_box);
               }
            }
         }

         int num_local_coarse_boxes = static_cast<int>(tmp_crs_mapped_box_set.size());
         int num_global_coarse_boxes = num_local_coarse_boxes;
         if (mpi.getSize() > 1) {
            mpi.Allreduce(&num_local_coarse_boxes, &num_global_coarse_boxes,
                          1, MPI_INT, MPI_SUM);
         }

         if (num_global_coarse_boxes) { 
            tbox::ConstPointer<hier::MappedBoxLevel> coarse_mapped_box_level =
               new hier::MappedBoxLevel(
                  tmp_crs_mapped_box_set,
                  d_multiblock_hierarchy->
                  getPatchLevel(next_coarser_level)->getRatioToLevelZero(),
                  grid_geometry);

            tbox::Pointer<hier::PatchLevel> coarse_level(
               new hier::PatchLevel(
                  *coarse_mapped_box_level,
                  grid_geometry,
                  dst_level->getPatchDescriptor(),
                  dst_level->getPatchFactory()));
            coarse_level->setLevelNumber( next_coarser_level );

            hier::NeighborhoodSet crs_to_dst_edges;

            int crs_box_counter = 0;
            for (hier::MappedBoxSet::const_iterator ni =
                    tmp_crs_mapped_box_set.begin();
                 ni != tmp_crs_mapped_box_set.end(); ++ni) {

               crs_to_dst_edges[(*ni).getId()].insert(
                  dst_mapped_boxes[dst_box_for_crs_box[crs_box_counter]]);
               ++crs_box_counter;
            }

            d_crs_to_dst_connector.initialize(
               *coarse_mapped_box_level,
               *(dst_level->getMappedBoxLevel()),
               hier::IntVector::getZero(dim),
               crs_to_dst_edges);

            createCoarseSchedule(coarse_level, next_coarser_level);
         } else {
            d_local_fill_only = true;
         } 
      } else {

         d_local_fill_only = true;

      }
      if (d_local_fill_only) {

         hier::OverlapConnectorAlgorithm oca;
         const int dst_ln = dst_level->getLevelNumber();
         const hier::IntVector dst_to_dst_gcw =
            hierarchy->getRequiredConnectorWidth(
               dst_ln,
               dst_ln);
         hierarchy->getRequiredConnectorWidth(dst_ln, dst_ln);
         dst_level->getMappedBoxLevel()->
         getPersistentOverlapConnectors().findOrCreateConnector(
            *dst_level->getMappedBoxLevel(),
            dst_to_dst_gcw);

         if (!src_level.isNull()) {
            const int src_ln = src_level->getLevelNumber();
            TBOX_ASSERT(src_ln == dst_ln);

            if (src_global_mapped_boxes.size()) {
               const hier::IntVector dst_to_src_gcw =
                  hierarchy->getRequiredConnectorWidth(dst_ln, src_ln);
               const hier::IntVector src_to_dst_gcw =
                  hierarchy->getRequiredConnectorWidth(src_ln, dst_ln);
               dst_level->getMappedBoxLevel()->
               getPersistentOverlapConnectors().findOrCreateConnector(
                  *src_level->getMappedBoxLevel(),
                  dst_to_src_gcw);
               src_level->getMappedBoxLevel()->
               getPersistentOverlapConnectors().findOrCreateConnector(
                  *dst_level->getMappedBoxLevel(),
                  src_to_dst_gcw);
            }
         }

         if (next_coarser_level >= 0) {
            const hier::IntVector dst_to_hiercoarse_width(
               hierarchy->getRequiredConnectorWidth(dst_ln, next_coarser_level));
            const hier::IntVector hiercoarse_to_dst_width(
               hierarchy->getRequiredConnectorWidth(next_coarser_level, dst_ln));

            tbox::Pointer<hier::PatchLevel> hiercoarse_level =
               hierarchy->getPatchLevel(next_coarser_level);
            const hier::MappedBoxLevel& hiercoarse_mapped_box_level =
               *(hiercoarse_level->getMappedBoxLevel());

            dst_level->getMappedBoxLevel()->
            getPersistentOverlapConnectors().findOrCreateConnector(
               hiercoarse_mapped_box_level,
               dst_to_hiercoarse_width);

            hiercoarse_mapped_box_level.getPersistentOverlapConnectors().findOrCreateConnector(
               *(dst_level->getMappedBoxLevel()),
               hiercoarse_to_dst_width);
         }

         hier::BlockId sched_block_id(hier::BlockId::invalidId());
         if (num_blocks == 1) {
            sched_block_id = hier::BlockId(0);
         }

         d_single_block_fill_local =
            d_single_block_refine_alg->createSchedule(
               d_fill_pattern,
               dst_level,
               src_level,
               next_coarser_level,
               hierarchy,
               (xfer::RefinePatchStrategy *)strategy,
               sched_block_id,
               use_time_refinement,
               d_transaction_factory);
      } else {
         if (src_global_mapped_boxes.size()) {

            hier::BlockId sched_block_id(hier::BlockId::invalidId());
            if (num_blocks == 1) {
               sched_block_id = hier::BlockId(0);
            }

            d_single_block_fill_local =
               d_single_block_refine_alg->createSchedule(
                  d_fill_pattern,
                  dst_level,
                  src_level,
                  (xfer::RefinePatchStrategy *)NULL,
                  sched_block_id,
                  use_time_refinement,
                  d_transaction_factory);
         } else {
            d_single_block_fill_local.setNull();
         }
      }
   } else {
      d_single_block_fill_local.setNull();
   }

   if (num_blocks > 1) {
      createInterblockSchedules(dst_level,
         src_level,
         (xfer::RefinePatchStrategy *)strategy,
         next_coarser_level + 1,
         use_time_refinement);
   }

   createOverlapsForRefine(dst_level, src_level);
   createOverlapsForCopy(dst_level, src_level);

}

/*
 * ************************************************************************
 *                                                                        *
 * The destructor implicitly deallocates all of the data associated with  *
 * the communication schedule.                                            *
 *                                                                        *
 * ************************************************************************
 */

MultiblockRefineSchedule::~MultiblockRefineSchedule()
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
void MultiblockRefineSchedule::createInterblockSchedules(
   tbox::Pointer<hier::PatchLevel> dst_level,
   tbox::Pointer<hier::PatchLevel> src_level,
   xfer::RefinePatchStrategy* refine_strategy,
   int level_number,
   bool use_time_refinement)
{
   const tbox::Dimension& dim(d_multiblock_hierarchy->getDim());

   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(dim, *dst_level);
   if (!src_level.isNull()) {
      TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(dim, *src_level);
   }

   tbox::Pointer<hier::GridGeometry> grid_geometry(
      d_multiblock_hierarchy->getGridGeometry());

   bool copy_all = d_fill_pattern->doesSourceLevelCommunicateToDestination();
   bool do_fill_fine_gcw = d_fill_pattern->fillingCoarseFineGhosts();

   bool encon_only = d_fill_pattern->fillingEnhancedConnectivityOnly();

   const hier::IntVector& constant_zero_intvector(hier::IntVector::getZero(dim));

   const int num_blocks =
      d_multiblock_hierarchy->getGridGeometry()->getNumberBlocks();

   const hier::MappedBoxSet& dst_mapped_boxes =
      dst_level->getMappedBoxLevel()->getMappedBoxes();
   const hier::MappedBoxSet& dst_global_mapped_boxes =
      dst_level->getMappedBoxLevel()->getGlobalizedVersion().getGlobalMappedBoxes();
   hier::MappedBoxSet src_mapped_boxes;
   hier::MappedBoxSet src_global_mapped_boxes;

   if (!(src_level.isNull())) {
      src_mapped_boxes = src_level->getMappedBoxLevel()->getMappedBoxes();
      src_global_mapped_boxes = src_level->getMappedBoxLevel()->getGlobalizedVersion().getGlobalMappedBoxes();
   }

   const tbox::SAMRAI_MPI &mpi(d_multiblock_hierarchy->getMPI());

   if (dst_global_mapped_boxes.size()) {


      hier::MappedBoxSet finalize_mapped_box_set;
      tbox::Array<tbox::Array<hier::BoxList> > finalize_boxes(num_blocks);
      int num_finalize_boxes = 0;

      for (int nb = 0; nb < num_blocks; nb++) {
         int num_neighbors = grid_geometry->getNumberOfNeighbors(nb);
         d_neighbor_copy_overlaps[nb].resizeArray(num_neighbors);
         d_neighbor_refine_overlaps[nb].resize(num_neighbors);
         d_finalize_ghost_patch_numbers[nb].resizeArray(num_neighbors);
         d_finalize_ghost_num_src_patches[nb].resizeArray(num_neighbors);
         d_neighbor_unfilled_boxes[nb].resizeArray(num_neighbors);

         finalize_boxes[nb].resizeArray(num_neighbors, hier::BoxList(dim));
         hier::BlockId block_id(nb);

         hier::BoxList dst_level_boxes;
         dst_level->getBoxes(dst_level_boxes, block_id);

         hier::BoxList dst_ghost_boxes(dst_level_boxes);
         dst_ghost_boxes.grow(d_data_fill_gcw);
         dst_ghost_boxes.removeIntersections(dst_level_boxes);

         int nc = 0;
         for (tbox::List<hier::GridGeometry::Neighbor>::
              Iterator ni(grid_geometry->getNeighbors(nb));
              ni; ni++) {

            if (encon_only) { 
               if (!ni().isSingularity()) {
                  nc++;
                  continue;
               }
            }

            hier::BlockId nbr_id(ni().getBlockNumber());
            hier::IntVector shift(ni().getShift());

            /*
             * get intersection with translated neighbor_boxes and
             * dst_ghost_boxes
             */
            hier::BoxList trans_neighbor_list(dim);
            grid_geometry->getTranslatedBlock(trans_neighbor_list,
               nb,
               nbr_id.getBlockValue());
            trans_neighbor_list.refine(dst_level->getRatioToLevelZero());
            trans_neighbor_list.intersectBoxes(dst_ghost_boxes);

            if (trans_neighbor_list.size() > 0) {

               /*
                * trans_neighbor_list is a list that contains the boxes for the
                * temporary dst level.  We need to create an equivalent array
                * such that each box in the array corresponds to a box on the
                * destination level, and is mapped to the same processor as the
                * corresponding box.
                */

               hier::BoxList neighbor_boxes(dim);

               if (dst_level->getLocalNumberOfPatches()) {
                  if (!src_level.isNull()) {
                     hier::BoxList neighbor_boxlist;
                     src_level->getBoxes(neighbor_boxlist, nbr_id);

                     neighbor_boxes = neighbor_boxlist;
                     grid_geometry->translateBoxList(
                        neighbor_boxes,
                        src_level->getRatioToLevelZero(),
                        block_id, nbr_id);
                  }

                  for (hier::MappedBoxSetSingleBlockIterator
                       dmb(dst_mapped_boxes, block_id); dmb.isValid(); dmb++) {

                     const hier::MappedBoxId& mapped_box_id = dmb->getId();
                     tbox::Pointer<hier::Patch> patch(
                        dst_level->getPatch(mapped_box_id));
 
                     if (patch->getPatchGeometry()->getTouchesRegularBoundary()) {
                        hier::Box dst_grow_box(patch->getBox());
                        dst_grow_box.grow(d_data_fill_gcw);

                        hier::BoxList tmp_list(trans_neighbor_list);
                        tmp_list.intersectBoxes(dst_grow_box);

                        if (tmp_list.size() > 0) {
                           if (level_number == MULTIBLOCK_FAKE_LEVEL_NUMBER) {
                              if (neighbor_boxes.size()) {
                                 tmp_list.intersectBoxes(neighbor_boxes);
                              }
                           } else {
                              if (!copy_all && neighbor_boxes.size()) {
                                 tmp_list.removeIntersections(neighbor_boxes);
                              }
                           }
                        }

                        if (tmp_list.size() > 0) {

                           tmp_list.coalesceBoxes();

                           if (tmp_list.size() > 1) {
                              hier::Box bound_box(tmp_list.getBoundingBox());
                              hier::BoxList bound_list;
                              bound_list.addItem(bound_box);
                              bound_list.removeIntersections(tmp_list);
                              if (bound_list.size() == 0) {
                                 tmp_list.clearItems();
                                 tmp_list.addItem(bound_box);
                              }
                           }
                           int rank = mpi.getRank();

                           d_finalize_ghost_patch_numbers[nb][nc].push_back(
                              num_finalize_boxes);
                           d_finalize_ghost_num_src_patches[nb][nc].push_back(
                              tmp_list.size());
                           for (hier::BoxList::Iterator bli(tmp_list);
                                bli; bli++) {
                              finalize_boxes[nb][nc].appendItem(bli());

                              hier::MappedBox finalize_mapped_box(bli(),
                                                                  hier::LocalId(num_finalize_boxes),
                                                                  rank,
                                                                  block_id);

                              finalize_mapped_box_set.insert(
                                 finalize_mapped_box_set.end(),
                                 finalize_mapped_box);

                              num_finalize_boxes++;
                           }
                        } else {
                           d_finalize_ghost_patch_numbers[nb][nc].push_back(-1);
                           d_finalize_ghost_num_src_patches[nb][nc].push_back(0);
                        }
                     } else {
                        d_finalize_ghost_patch_numbers[nb][nc].push_back(-1);
                        d_finalize_ghost_num_src_patches[nb][nc].push_back(0);
                     }
                  }
               }
            }
            nc++;
         }
      }
 
      int global_finalize_boxes = num_finalize_boxes;
      if (mpi.getSize() > 1) {
         mpi.Allreduce(&num_finalize_boxes, &global_finalize_boxes, 1, MPI_INT, MPI_SUM);
      }

      if (global_finalize_boxes) {

         const int rank = mpi.getRank();

         tbox::ConstPointer<hier::MappedBoxLevel>
         finalize_mapped_box_level =
            new hier::MappedBoxLevel(finalize_mapped_box_set,
               dst_level->getRatioToLevelZero(),
               dst_level->getGridGeometry());

         d_finalize_ghost_level =
            new hier::PatchLevel(*finalize_mapped_box_level,
               d_multiblock_hierarchy->getGridGeometry(),
               d_multiblock_hierarchy->getPatchDescriptor());
         d_finalize_ghost_level->setLevelNumber(
            dst_level->getLevelNumber() );

         hier::MappedBoxSet neighbor_mapped_box_set;
         hier::LocalId neighbor_index(0);
         for (int nb = 0; nb < num_blocks; nb++) {
            hier::BlockId block_id(nb);

            int nc = 0;
            for (tbox::List<hier::GridGeometry::Neighbor>::
                 Iterator ni(grid_geometry->getNeighbors(nb));
                 ni; ni++) {

               if (encon_only) {
                  if (!ni().isSingularity()) {
                     nc++;
                     continue;
                  }
               }

               hier::BlockId nbr_id(ni().getBlockNumber());

               /*
                * From here take the finalize_box_array and translate each box
                * to the neighbor's index space.  Create a level with the same
                * mapping.  Then create only one refine schedule to copy from
                * the neighbor source level to the neighbor temp level.
                */
               hier::BoxList& neighbor_ghost_list = finalize_boxes[nb][nc];

               grid_geometry->translateBoxList(
                  neighbor_ghost_list, dst_level->getRatioToLevelZero(),
                  hier::BlockId(ni().getBlockNumber()),
                  block_id);

               for (hier::BoxList::Iterator nf(neighbor_ghost_list); nf; nf++) {

                  hier::MappedBox local_mapped_box(*nf,
                                                   neighbor_index++,
                                                   rank,
                                                   nbr_id);

                  neighbor_mapped_box_set.insert(neighbor_mapped_box_set.end(),
                     local_mapped_box);

               }
               nc++;
            }
         }

         tbox::ConstPointer<hier::MappedBoxLevel>
         neighbor_mapped_box_level =
            new hier::MappedBoxLevel(neighbor_mapped_box_set,
               dst_level->getRatioToLevelZero(),
               dst_level->getGridGeometry());

         d_neighbor_ghost_level =
            new hier::PatchLevel(*neighbor_mapped_box_level,
               d_multiblock_hierarchy->getGridGeometry(),
               d_multiblock_hierarchy->getPatchDescriptor());
         d_neighbor_ghost_level->setLevelNumber(
            dst_level->getLevelNumber() );

         tbox::Pointer<PatchLevelFullFillPattern> fill_pattern(
            new PatchLevelFullFillPattern());

         bool unfilled_boxes_exist = false;
         if (level_number != MULTIBLOCK_FAKE_LEVEL_NUMBER) {

            const hier::MappedBoxSet& neighbor_global_mapped_boxes =
               d_multiblock_hierarchy->getPatchLevel(level_number)->
               getMappedBoxLevel()->getGlobalizedVersion().
               getGlobalMappedBoxes();
                        
            tbox::Pointer<hier::PatchLevel> neighbor_dst_level =
               d_neighbor_ghost_level;

            hier::MappedBoxSet neighbor_coarse_mapped_box_set;
            hier::LocalId neighbor_coarse_index(0);
            hier::NeighborhoodSet crs_to_dst_edges;

            for (int nb = 0; nb < num_blocks; nb++) {

               int nc = 0;
               for (tbox::List<hier::GridGeometry::Neighbor>::
                    Iterator ni(grid_geometry->getNeighbors(nb));
                    ni; ni++) {

                  if (encon_only) {
                     if (!ni().isSingularity()) {
                        nc++;
                        continue;
                     }
                  }

                  hier::BlockId nbr_id(ni().getBlockNumber());

                  neighbor_dst_level->getBoxes(
                     d_neighbor_unfilled_boxes[nb][nc],
                     nbr_id);

                  if (neighbor_global_mapped_boxes.size()) {
                     hier::BoxList neighbor_hier_boxes;
                     d_multiblock_hierarchy->getPatchLevel(level_number)->getBoxes(
                        neighbor_hier_boxes, nbr_id);
                     d_neighbor_unfilled_boxes[nb][nc].
                        removeIntersections(neighbor_hier_boxes);
                  }

                  if (!unfilled_boxes_exist &&
                      d_neighbor_unfilled_boxes[nb][nc].size()) {
                     unfilled_boxes_exist = true;
                  }

                  if ((copy_all ||
                      d_neighbor_unfilled_boxes[nb][nc].size()) &&
                      do_fill_fine_gcw) {

                     d_neighbor_unfilled_boxes[nb][nc].coalesceBoxes();

                     const hier::MappedBoxSet& nbor_dst_mapped_boxes =
                        neighbor_dst_level->getMappedBoxLevel()->
                           getMappedBoxes();

                     tbox::List<hier::MappedBox> dst_mapped_box_list;
                     hier::BoxList needed_src_boxes;
                     hier::BoxList tmp_src_boxes;
                     for (hier::MappedBoxSetSingleBlockIterator
                          dn(nbor_dst_mapped_boxes, nbr_id);
                          dn.isValid(); ++dn) {

                        needed_src_boxes.appendItem((*dn).getBox());
                        needed_src_boxes.intersectBoxes(
                           d_neighbor_unfilled_boxes[nb][nc]);

                        if (needed_src_boxes.size()) {
                           for (hier::BoxList::Iterator ns(needed_src_boxes);
                                ns; ns++) {

                              tmp_src_boxes.appendItem(ns());
                              dst_mapped_box_list.appendItem(*dn); 
                           }
                        }

                        needed_src_boxes.clearItems();
                     }

                     if (tmp_src_boxes.size()) {
                        hier::IntVector coarse_ratio =
                           dst_level->getRatioToCoarserLevel();
                        if (coarse_ratio == constant_zero_intvector) {
                           tbox::Pointer<hier::PatchLevel>
                           example_level = d_multiblock_hierarchy->getPatchLevel(
                                 level_number - 1);
                           coarse_ratio =
                              neighbor_dst_level->getRatioToLevelZero()
                              / example_level->getRatioToLevelZero();
                        }
                        tmp_src_boxes.coarsen(coarse_ratio);

                        tbox::List<hier::MappedBox>::Iterator
                           db(dst_mapped_box_list);
                        for (hier::BoxList::Iterator ts(tmp_src_boxes);
                             ts; ts++) {
   
                           hier::MappedBox local_mapped_box(
                              ts(),
                              neighbor_coarse_index++,
                              rank,
                              nbr_id);

                           neighbor_coarse_mapped_box_set.insert(
                              neighbor_coarse_mapped_box_set.end(),
                              local_mapped_box);

                           crs_to_dst_edges[local_mapped_box.getId()].
                              insert(db());
                           db++;
                        }
                     }
                  }
                  nc++;
               }
            }

            int num_local_coarse_boxes = static_cast<int>(neighbor_coarse_mapped_box_set.size());
            int num_global_coarse_boxes = num_local_coarse_boxes;
            if (mpi.getSize() > 1) {
               mpi.Allreduce(&num_local_coarse_boxes, &num_global_coarse_boxes,
                             1, MPI_INT, MPI_SUM);
            }

            if (num_global_coarse_boxes) {
               tbox::ConstPointer<hier::MappedBoxLevel>
               neighbor_coarse_mapped_box_level =
                  new hier::MappedBoxLevel(
                     neighbor_coarse_mapped_box_set,
                     d_multiblock_hierarchy->
                     getPatchLevel(level_number - 1)->getRatioToLevelZero(),
                     grid_geometry);

               tbox::Pointer<hier::PatchLevel>
               neighbor_coarse_level(
                  new hier::PatchLevel(
                     *neighbor_coarse_mapped_box_level,
                     neighbor_dst_level->getGridGeometry(),
                     neighbor_dst_level->getPatchDescriptor(),
                     neighbor_dst_level->getPatchFactory()));
               neighbor_coarse_level->setLevelNumber(level_number - 1);

               d_neighbor_crs_to_dst_connector.initialize(
                  *neighbor_coarse_mapped_box_level,
                  *(neighbor_dst_level->getMappedBoxLevel()),
                  hier::IntVector::getZero(dim),
                  crs_to_dst_edges);

               createNeighborCoarseSchedule(
                  neighbor_coarse_level, level_number - 1,
                  dst_level->getRatioToCoarserLevel(),
                  d_multiblock_hierarchy);

            }

            if (neighbor_global_mapped_boxes.size()) {
               if ((copy_all || unfilled_boxes_exist)
                   && do_fill_fine_gcw) {
                  const hier::IntVector neighbor_connector_width(
                     d_multiblock_hierarchy->getRequiredConnectorWidth(
                        level_number,
                        level_number));

                  d_neighbor_ghost_level->getMappedBoxLevel()->
                  getPersistentOverlapConnectors().createConnector(
                     *(d_multiblock_hierarchy->getPatchLevel(level_number)->getMappedBoxLevel()),
                     neighbor_connector_width);

                  d_multiblock_hierarchy->getPatchLevel(level_number)->getMappedBoxLevel()->
                  getPersistentOverlapConnectors().createConnector(
                     *(d_neighbor_ghost_level->getMappedBoxLevel()),
                     neighbor_connector_width);

                  d_neighbor_single_block_refine_schedule =
                     d_single_block_refine_alg->createSchedule(
                        fill_pattern,
                        d_neighbor_ghost_level,
                        d_multiblock_hierarchy->getPatchLevel(level_number),
                        NULL,
                        hier::BlockId::invalidId(),
                        use_time_refinement,
                        d_transaction_factory);
               }
            } else {
               d_neighbor_single_block_refine_schedule.setNull();
            }
         } else {

            if (!src_level.isNull()) {

               if (src_global_mapped_boxes.size()) {
                  if ((copy_all || unfilled_boxes_exist)
                      && do_fill_fine_gcw) {

                     const int neighbor_ln =
                        src_level->getLevelNumber();
                     const hier::IntVector neighbor_connector_width(
                        d_multiblock_hierarchy->getRequiredConnectorWidth(
                           neighbor_ln,
                           neighbor_ln));

                     d_neighbor_ghost_level->getMappedBoxLevel()->
                     getPersistentOverlapConnectors().createConnector(
                        *(d_neighbor_ghost_level->getMappedBoxLevel()),
                                 neighbor_connector_width);

                     d_neighbor_ghost_level->getMappedBoxLevel()->
                     getPersistentOverlapConnectors().createConnector(
                        *(src_level->getMappedBoxLevel()),
                        neighbor_connector_width);

                     src_level->getMappedBoxLevel()->getPersistentOverlapConnectors().
                     createConnector(
                        *(d_neighbor_ghost_level->getMappedBoxLevel()),
                        neighbor_connector_width);

                     d_neighbor_single_block_refine_schedule =
                        d_single_block_refine_alg->createSchedule(
                           fill_pattern,
                           d_neighbor_ghost_level,
                           src_level,
                           refine_strategy,
                           hier::BlockId::invalidId(),
                           use_time_refinement,
                           d_transaction_factory);
                  }
               }
            } else {
               d_neighbor_single_block_refine_schedule.setNull();
            }
         }
      } else {
         d_neighbor_single_block_refine_schedule.setNull();
      }
   }
}

/*
 * ************************************************************************
 *
 * Create the BoxOverlaps that will be used in refineScratchData.         *
 *                                                                        *
 * ************************************************************************
 */
void MultiblockRefineSchedule::createOverlapsForRefine(
   tbox::Pointer<hier::PatchLevel>& dst_level,
   tbox::Pointer<hier::PatchLevel>& src_level)
{
   NULL_USE(src_level); 
   TBOX_DIM_ASSERT_CHECK_ARGS2(*d_multiblock_hierarchy, *dst_level);

   tbox::Pointer<hier::GridGeometry> grid_geometry(
      d_multiblock_hierarchy->getGridGeometry());

   bool fill_fine_gcw = d_fill_pattern->fillingCoarseFineGhosts();

   const hier::MappedBoxSet& dst_global_mapped_boxes =
      dst_level->getMappedBoxLevel()->getGlobalizedVersion().getGlobalMappedBoxes(); 

   const int nblocks = d_multiblock_hierarchy->getGridGeometry()->getNumberBlocks();

   if (!d_multiblock_coarse_scratch_level.isNull()) {

      for (int nb = 0; nb < nblocks; nb++) {
         const bool filling_neighbor = false;
         createRefineOverlaps(d_local_refine_overlaps[nb],
            d_multiblock_coarse_scratch_level,
            dst_level,
            d_crs_to_dst_connector,
            d_unfilled_boxes[nb],
            nb,
            fill_fine_gcw,
            filling_neighbor);
      }
   }

   for (int nb = 0; nb < nblocks; nb++) {

      hier::MappedBoxSetSingleBlockIterator dst_iter(
         dst_global_mapped_boxes, hier::BlockId(nb));

      if (dst_iter.isValid()) {

         int nc = 0;
         for (tbox::List<hier::GridGeometry::Neighbor>::
              Iterator ni(grid_geometry->getNeighbors(nb));
              ni; ni++) {
            int id = ni().getBlockNumber();

            if (!d_neighbor_multiblock_coarse_level.isNull()) {

               const bool fill_gcw = false;
               const bool filling_neighbor = true;
               createRefineOverlaps(d_neighbor_refine_overlaps[nb][nc],
                  d_neighbor_multiblock_coarse_level,
                  d_neighbor_ghost_level,
                  d_neighbor_crs_to_dst_connector,
                  d_neighbor_unfilled_boxes[nb][nc],
                  id,
                  fill_gcw,
                  filling_neighbor);
            }
            nc++;
         }
      }
   }
}

/*
 * ************************************************************************
 *
 * Create the BoxOverlaps that will be used in fillData.                  *
 *                                                                        *
 * ************************************************************************
 */
void MultiblockRefineSchedule::createOverlapsForCopy(
   tbox::Pointer<hier::PatchLevel>& dst_level,
   tbox::Pointer<hier::PatchLevel>& src_level)
{

   tbox::Pointer<hier::GridGeometry> grid_geometry(
      d_multiblock_hierarchy->getGridGeometry());

   const hier::MappedBoxSet& dst_global_mapped_boxes =
      dst_level->getMappedBoxLevel()->getGlobalizedVersion().getGlobalMappedBoxes();

   const int nblocks = d_multiblock_hierarchy->getGridGeometry()->getNumberBlocks();
   for (int nb = 0; nb < nblocks; nb++) {

      hier::MappedBoxSetSingleBlockIterator dst_iter(
         dst_global_mapped_boxes, hier::BlockId(nb));

      if (dst_iter.isValid()) {

         int nc = 0;
         for (tbox::List<hier::GridGeometry::Neighbor>::
              Iterator ni(grid_geometry->getNeighbors(nb));
              ni; ni++) {

            if (!(d_neighbor_single_block_refine_schedule.isNull())
                ||
                !(d_neighbor_multiblock_coarse_schedule.isNull())) {

               const tbox::Pointer<xfer::RefineClasses>
               equiv_classes =
                  d_single_block_refine_alg->getEquivalenceClasses();

               TBOX_ASSERT(d_neighbor_copy_overlaps[nb][nc].size() == 0);
               createCopyOverlaps(d_neighbor_copy_overlaps[nb][nc],
                                  dst_level,
                                  src_level,
                                  ni(),
                                  nb,
                                  nc,
                                  equiv_classes);

            }
            nc++;
         }
      }
   }
}

/*
 * ************************************************************************
 *                                                                        *
 * Create and store overlaps that will be used in refineScratchData.      *
 *                                                                        *
 * ************************************************************************
 */
void MultiblockRefineSchedule::createCopyOverlaps(
   tbox::List<tbox::Pointer<hier::BoxOverlap> >& overlaps,
   const tbox::Pointer<hier::PatchLevel>& dst_patch_level,
   const tbox::Pointer<hier::PatchLevel>& src_patch_level,
   const hier::GridGeometry::Neighbor& neighbor,
   const int block_number,
   const int neighbor_counter,
   const tbox::Pointer<xfer::RefineClasses>& equiv_classes)
{
   const tbox::Dimension& dim(dst_patch_level->getDim());
   const tbox::SAMRAI_MPI& mpi(dst_patch_level->getMappedBoxLevel()->getMPI());
   const int nb = block_number;
   const int nc = neighbor_counter;
   const int num_classes = equiv_classes->getNumberOfEquivalenceClasses();
   const int rank = mpi.getRank();
   const hier::BlockId dst_block_id(block_number);
   const hier::BlockId nbr_id(neighbor.getBlockNumber());

   bool copy_all = d_fill_pattern->doesSourceLevelCommunicateToDestination(); 
   bool filling_interior_only = !d_fill_pattern->fillingCoarseFineGhosts();

   hier::BoxList unfilled_boxes(dim);
   if (!copy_all) {
      unfilled_boxes = d_neighbor_unfilled_boxes[nb][nc];
      for (hier::BoxList::Iterator ub(unfilled_boxes); ub; ub++) {
         ub().rotate(neighbor.getRotationIdentifier());
         ub().shift(neighbor.getShift() *
                    dst_patch_level->getRatioToLevelZero());
      }
   } else if (d_only_copy_source_intersections) {
      TBOX_ASSERT(d_neighbor_unfilled_boxes[nb][nc].size() == 0);

      if (!src_patch_level.isNull()) {
         tbox::Pointer<hier::GridGeometry> grid_geom(
            src_patch_level->getGridGeometry());
         hier::IntVector ratio(src_patch_level->getRatioToLevelZero());

         hier::BoxList trans_neighbor_boxes(dim);
         grid_geom->getTranslatedBlock(
            trans_neighbor_boxes,
            nb,
            nbr_id.getBlockValue());
         trans_neighbor_boxes.refine(dst_patch_level->getRatioToLevelZero());

         const hier::MappedBoxSet& global_src_mapped_boxes =
            src_patch_level->getMappedBoxLevel()->
            getGlobalizedVersion().getGlobalMappedBoxes();

         hier::BoxList neighbor_src_boxes;
         for (hier::MappedBoxSetSingleBlockIterator
              mb(global_src_mapped_boxes, nbr_id); mb.isValid(); mb++) {
            hier::Box nbr_box(mb->getBox());
            grid_geom->translateBox(nbr_box, ratio, dst_block_id, nbr_id);

            neighbor_src_boxes.appendItem(nbr_box);
         }

         if (neighbor_src_boxes.size()) {
            const hier::MappedBoxSet& dst_mapped_boxes =
               dst_patch_level->getMappedBoxLevel()->getMappedBoxes();

            for (hier::MappedBoxSetSingleBlockIterator
                 mb(dst_mapped_boxes, dst_block_id); mb.isValid(); mb++) {
               d_neighbor_unfilled_boxes[nb][nc].appendItem(mb->getBox());
            }
            d_neighbor_unfilled_boxes[nb][nc].grow(d_data_fill_gcw);
            d_neighbor_unfilled_boxes[nb][nc].intersectBoxes(neighbor_src_boxes);
            if (d_neighbor_unfilled_boxes[nb][nc].size()) {
               d_neighbor_unfilled_boxes[nb][nc].coalesceBoxes();
               unfilled_boxes = d_neighbor_unfilled_boxes[nb][nc];
            }
         }
      }
   }

   int num_overs_per_patch = (copy_all && !d_only_copy_source_intersections) ?
                             1 : unfilled_boxes.size();

   if (num_overs_per_patch && !filling_interior_only) {
      hier::MappedBoxSetSingleBlockIterator dst_local_iter(
         d_multiblock_dst_level->getMappedBoxLevel()->getMappedBoxes(),
         dst_block_id);

      if (d_finalize_ghost_patch_numbers[nb][nc].size()) {
         int pc = 0;
         for ( ; dst_local_iter.isValid(); dst_local_iter++) {

            if (d_finalize_ghost_patch_numbers[nb][nc][pc] >= 0) {

               const hier::MappedBoxId& mapped_box_id = dst_local_iter->getId();
               tbox::Pointer<hier::Patch> dst_patch(
                  d_multiblock_dst_level->getPatch(mapped_box_id));

               hier::Box dst_grow_box(dim);
               if (copy_all) {
                  dst_grow_box = hier::Box::grow(
                        dst_patch->getBox(),
                        d_data_fill_gcw);
               }

               const int num_src_patches =
                  d_finalize_ghost_num_src_patches[nb][nc][pc];
               for (int ns = 0; ns < num_src_patches; ns++) {
                  hier::LocalId local_id(d_finalize_ghost_patch_numbers[nb][nc][pc]
                                         + ns);
                  hier::MappedBoxId mbid(local_id, rank, dst_block_id);
                  tbox::Pointer<hier::Patch> src_patch =
                     d_finalize_ghost_level->getPatch(mbid);

                  hier::BoxList::Iterator itr(unfilled_boxes);
                  for (int nop = 0; nop < num_overs_per_patch; nop++, itr++) {
                     hier::Box restrict_box(dim);
                     if (!copy_all) {
                        restrict_box = *itr;
                     } else {
                        restrict_box = dst_grow_box;
                        if (d_only_copy_source_intersections) {
                           restrict_box = restrict_box * *itr;
                        }
                     }

                     if (copy_all || !itr().isEmpty()) {
                        if (!neighbor.isSingularity() ||
                            !d_enable_singularity_patches) {
                           for (int ne = 0; ne < num_classes; ne++) {

                              const RefineClasses::Data& rep_item =
                                 equiv_classes->getClassRepresentative(ne);

                              tbox::Pointer<hier::BoxOverlap> box_overlap =
                                 calculateOverlap(*dst_patch,
                                    *src_patch,
                                    rep_item,
                                    restrict_box);

                              overlaps.appendItem(box_overlap);
                           }
                        } else {

                           for (int ne = 0; ne < num_classes; ne++) {

                              const RefineClasses::Data& rep_item =
                                 equiv_classes->getClassRepresentative(ne);

                              tbox::Pointer<hier::BoxOverlap> box_overlap =
                                 calculateSingularityPatchOverlap(*src_patch,
                                    rep_item,
                                    restrict_box);

                              overlaps.appendItem(box_overlap);

                           }
                        }
                     }
                  }
               }
            }
            pc++;
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

void MultiblockRefineSchedule::fillData(
   double fill_time,
   bool do_physical_boundary_fill) const
{
   fillData(fill_time,
      do_physical_boundary_fill,
      false,
      false);
}

void MultiblockRefineSchedule::fillData(
   double fill_time,
   bool do_physical_boundary_fill,
   bool filling_coarse_scratch,
   bool filling_crse_scr_recursive) const
{
   const tbox::Dimension& dim(d_multiblock_hierarchy->getDim());

   const tbox::SAMRAI_MPI& mpi(
      d_multiblock_dst_level->getMappedBoxLevel()->getMPI());

   const int rank = mpi.getRank();

   tbox::Pointer<hier::GridGeometry> grid_geometry(
      d_multiblock_hierarchy->getGridGeometry());

   int dst_ln = d_multiblock_dst_level->getLevelNumber();

   bool is_scr_recursive = false;
   tbox::Array<hier::ComponentSelector> fill_crse_scr_selector(0);

   const int nblocks = d_multiblock_hierarchy->getGridGeometry()->getNumberBlocks();
   if (filling_coarse_scratch) {
      d_multiblock_strategy->setFillingCoarseScratch(true);
      is_scr_recursive = true;
      fill_crse_scr_selector.resizeArray(nblocks);
   }

   hier::ComponentSelector allocate_scr_vector;

   const hier::MappedBoxSet& dst_global_mapped_boxes =
      d_multiblock_dst_level->getMappedBoxLevel()->getGlobalizedVersion().getGlobalMappedBoxes();

   if (!filling_coarse_scratch) {
      allocateScratchSpace(
         allocate_scr_vector,
         d_multiblock_dst_level,
         fill_time);
   }

   bool copy_all = d_fill_pattern->doesSourceLevelCommunicateToDestination();

   if (d_local_fill_only) {
      if (!d_single_block_fill_local.isNull()) {
         if (filling_coarse_scratch) {
            d_single_block_fill_local->allocateDestinationSpace(
               *d_coarse_selector, fill_time);
         }

         /*
          * If there is only one block, then call physical boundary filling,
          * if desired, internally within RefineSchedule.  Otherwise
          * physical boundary filling will be called later on in this method.
          */
         if (nblocks == 1) {
            d_single_block_fill_local->fillData(fill_time,
                                                do_physical_boundary_fill);
         } else {
            d_single_block_fill_local->fillData(fill_time, false);
         }
      }
   } else {
      d_multiblock_coarse_schedule->fillData(fill_time,
         true,
         true,
         is_scr_recursive);

      if (filling_coarse_scratch) {
         if (dst_global_mapped_boxes.size()) {
            d_multiblock_coarse_schedule->allocateScratchSpace(
               *d_coarse_selector,
               d_multiblock_dst_level,
               fill_time);
         }   
      }

      if (!d_single_block_fill_local.isNull()) {
         d_single_block_fill_local->fillData(fill_time, false);
      }

      if (dst_global_mapped_boxes.size()) {

         bool fill_fine_gcw = d_fill_pattern->fillingCoarseFineGhosts();

         for (int nb = 0; nb < nblocks; nb++) {
            refineScratchData(d_multiblock_coarse_scratch_level,
               d_multiblock_dst_level,
               d_local_refine_overlaps[nb],
               d_crs_to_dst_connector,
               d_unfilled_boxes[nb],
               nb,
               fill_fine_gcw,
               false);

            copyScratchToDestination(d_multiblock_dst_level,
               d_unfilled_boxes[nb],
               hier::BlockId(nb),
               d_single_block_refine_alg->getEquivalenceClasses());
         }
      }

      tbox::Pointer<hier::ComponentSelector> crs_scr_vector =
         d_multiblock_coarse_schedule->d_coarse_selector;
      d_multiblock_coarse_scratch_level->
      deallocatePatchData(*crs_scr_vector);
   }

   if (nblocks > 1 && d_multiblock_strategy != NULL) {
      tbox::Array<bool> singularity_to_fill(nblocks, false);
      tbox::Array<tbox::Array<tbox::List<tbox::Pointer<hier::Patch> > > >
         singularity_patches(nblocks);
      tbox::Array<tbox::Array<hier::ComponentSelector> >
         local_selector(nblocks);
      tbox::Array<int> num_local_on_block(nblocks,0);
      bool do_neighbor_copy = false;
      for (int nb = 0; nb < nblocks; nb++) {

         hier::MappedBoxSetSingleBlockIterator dst_iter(
            dst_global_mapped_boxes, hier::BlockId(nb));

         if (dst_iter.isValid()) {

            int num_global_on_block = 0; 
            for ( ; dst_iter.isValid(); dst_iter++) {
               num_global_on_block++;
               if (dst_iter->getOwnerRank() == rank) {
                  num_local_on_block[nb]++;
               }
            }

            if (d_enable_singularity_patches) {
               singularity_patches[nb].resizeArray(num_global_on_block);
            }

            if (grid_geometry->reducedConnectivityExists(nb) &&
                d_enable_singularity_patches) {
               singularity_to_fill[nb] = true;
            }
            local_selector[nb].resizeArray(num_global_on_block);

            int nc = 0;
            for (tbox::List<hier::GridGeometry::Neighbor>::
                 Iterator ni(grid_geometry->getNeighbors(nb));
                 ni; ni++) {

               if (copy_all || d_neighbor_unfilled_boxes[nb][nc].size() > 0) {
                  do_neighbor_copy = true;
               }
               nc++;
            }      
         }
      }
      
      hier::ComponentSelector dst_vector;
      if (!(d_neighbor_single_block_refine_schedule.isNull())
          || !(d_neighbor_multiblock_coarse_schedule.isNull())) {

         if (!d_neighbor_single_block_refine_schedule.isNull()) {
            d_neighbor_single_block_refine_schedule->
            allocateDestinationSpace(dst_vector, fill_time);
         } else {
            if (!d_single_block_fill_local.isNull()) {
               d_single_block_fill_local->
               initializeDestinationVector(dst_vector);
            } else {
               initializeDestinationVector(dst_vector);
            }

            const int ncomponents = d_neighbor_ghost_level->
               getPatchDescriptor()->
               getMaxNumberRegisteredComponents();

            for (int di = 0; di < ncomponents; di++) {
               if (dst_vector.isSet(di)) {
                  if (d_neighbor_ghost_level->
                      checkAllocated(di)) {
                     dst_vector.clrFlag(di);
                  }
               }
            }

            d_neighbor_ghost_level->allocatePatchData(
                        dst_vector, fill_time);

         }

         bool do_fill_fine_gcw = d_fill_pattern->fillingCoarseFineGhosts();
         if (do_neighbor_copy && do_fill_fine_gcw) {
            if (!d_neighbor_multiblock_coarse_schedule.isNull()) {
               d_neighbor_multiblock_coarse_schedule->
               fillData(fill_time, true, true, is_scr_recursive);
            }
            if (!d_neighbor_single_block_refine_schedule.
                isNull()) {
               d_neighbor_single_block_refine_schedule->
               fillData(fill_time, false);
            }

            if (!d_neighbor_multiblock_coarse_schedule.isNull()) {

               hier::ComponentSelector scr_vector;
               d_neighbor_multiblock_coarse_schedule->
               allocateScratchSpace(
                  scr_vector,
                  d_neighbor_ghost_level,
                  fill_time);

               for (int nb = 0; nb < nblocks; nb++) {

                  hier::MappedBoxSetSingleBlockIterator dst_iter(
                     dst_global_mapped_boxes, hier::BlockId(nb));

                  if (dst_iter.isValid()) { 
                     int nc = 0;
                     for (tbox::List<hier::GridGeometry::Neighbor>::
                          Iterator ni(grid_geometry->getNeighbors(nb));
                          ni; ni++) {

                        refineScratchData(
                           d_neighbor_multiblock_coarse_level,
                           d_neighbor_ghost_level,
                           d_neighbor_refine_overlaps[nb][nc],
                           d_neighbor_crs_to_dst_connector,
                           d_neighbor_unfilled_boxes[nb][nc],
                           ni().getBlockNumber(),
                           false,
                           true);

                        copyScratchToDestination(
                           d_neighbor_ghost_level,
                           d_neighbor_unfilled_boxes[nb][nc],
                           hier::BlockId(ni().getBlockNumber()),
                           d_single_block_refine_alg->getEquivalenceClasses());

                        nc++;
                     }
                  }
               }

               d_neighbor_ghost_level->
                  deallocatePatchData(scr_vector);

               tbox::Pointer<hier::ComponentSelector> crs_scr_vector =
                  d_neighbor_multiblock_coarse_schedule->
                     d_coarse_selector;
                  d_neighbor_multiblock_coarse_level->
                  deallocatePatchData(*crs_scr_vector);

            }

            /*
             * get a src_vector to allocate on the finalize_ghost_level,
             * then copy between blocks.
             */

            hier::ComponentSelector src_vector;
            if (do_neighbor_copy && do_fill_fine_gcw) {
               if (d_neighbor_single_block_refine_schedule.isNull()) {
                  d_neighbor_multiblock_coarse_schedule->
                  initializeSourceVector(
                     src_vector);
               } else {
                  d_neighbor_single_block_refine_schedule->
                  initializeSourceVector(
                     src_vector);
               }

               d_finalize_ghost_level->
                  allocatePatchData(src_vector);
            }

            const tbox::Pointer<xfer::RefineClasses>
            equiv_classes =
               d_single_block_refine_alg->getEquivalenceClasses();
            int num_classes =
               equiv_classes->getNumberOfEquivalenceClasses();

            for (int nb = 0; nb < nblocks; nb++) {

               int nc = 0;
               for (tbox::List<hier::GridGeometry::Neighbor>::
                    Iterator ni(grid_geometry->getNeighbors(nb));
                    ni; ni++) {

                  hier::BoxList unfilled_boxes(dim);

                  if (num_local_on_block[nb] > 0) {
                     if (do_neighbor_copy && do_fill_fine_gcw) {
                        if (d_using_standard_transaction) {
                           copyBetweenBlocks(d_finalize_ghost_level,
                              d_neighbor_ghost_level,
                              (ni().getShift()) *
                              (d_multiblock_dst_level->getRatioToLevelZero()),
                              ni().getRotationIdentifier(), equiv_classes,
                              hier::BlockId(nb),
                              hier::BlockId(ni().getBlockNumber()));
                        } else {
                           fillBetweenBlocks(d_finalize_ghost_level,
                              d_neighbor_ghost_level,
                              (ni().getShift()) *
                              (d_multiblock_dst_level->getRatioToLevelZero()),
                              ni().getRotationIdentifier(), equiv_classes,
                              hier::BlockId(nb),
                              hier::BlockId(ni().getBlockNumber()));
                        }
                     }

                     if (!copy_all) {
                        unfilled_boxes = d_neighbor_unfilled_boxes[nb][nc];
                        for (hier::BoxList::Iterator ub(unfilled_boxes); ub; ub++) {
                           ub().rotate(ni().getRotationIdentifier());
                           ub().shift(ni().getShift()
                              * d_multiblock_dst_level->getRatioToLevelZero());
                        }
                     } else if (d_only_copy_source_intersections) {
                        if (d_neighbor_unfilled_boxes[nb][nc].size()) {
                           unfilled_boxes = d_neighbor_unfilled_boxes[nb][nc];
                        } 
                     }
                  }

                  if (do_neighbor_copy &&
                      num_local_on_block[nb] &&
                      do_fill_fine_gcw &&
                      d_neighbor_copy_overlaps[nb][nc].size()) {

                      int num_overs_per_patch =
                         (copy_all && !d_only_copy_source_intersections) ?
                         1 : unfilled_boxes.size();

                     tbox::ListIterator<tbox::Pointer<hier::BoxOverlap> >
                        overlap_iter(d_neighbor_copy_overlaps[nb][nc]); 

                     hier::MappedBoxSetSingleBlockIterator dst_local_iter(
                        d_multiblock_dst_level->getMappedBoxLevel()->getMappedBoxes(),
                        hier::BlockId(nb));

                     int pc = 0;
                     for ( ; dst_local_iter.isValid(); dst_local_iter++) { 

                        if (d_finalize_ghost_patch_numbers[nb][nc][pc] >= 0) {

                           const hier::MappedBoxId& mapped_box_id =
                              dst_local_iter->getId();
                           tbox::Pointer<hier::Patch> dst_patch(
                              d_multiblock_dst_level->getPatch(mapped_box_id));

                           const int num_src_patches =
                              d_finalize_ghost_num_src_patches[nb][nc][pc];
                           for (int ns = 0; ns < num_src_patches; ns++) {
                              hier::LocalId local_id(
                                 d_finalize_ghost_patch_numbers[nb][nc][pc]+ns);
                              hier::MappedBoxId mbid(local_id,
                                                     rank,
                                                     hier::BlockId(nb));

                              tbox::Pointer<hier::Patch> src_patch =
                                 d_finalize_ghost_level->getPatch(mbid);

                              hier::BoxList::Iterator itr(unfilled_boxes);
                              for (int nop = 0; nop < num_overs_per_patch; nop++, itr++) {

                                 if (copy_all || !itr().empty()) {
                                    if (!ni().isSingularity() ||
                                        !d_enable_singularity_patches) {
                                       for (int ne = 0; ne < num_classes; ne++) {
                                          if (!(overlap_iter())->isOverlapEmpty()) {
                                             for (tbox::List<int>::Iterator
                                                  l(equiv_classes->getIterator(ne));
                                                  l; l++) {

                                                const RefineClasses::Data& item =
                                                equiv_classes->getRefineItem(l());
                                                TBOX_ASSERT(item.d_class_id == ne);

                                                const int dst = item.d_dst;
                                                const int src = item.d_src;

                                                dst_patch->getPatchData(dst)->
                                                copy(*(src_patch->getPatchData(src)),
                                                   *(overlap_iter()));
                                             }
                                          }
                                          overlap_iter++;
                                       }
                                    } else {
                                       singularity_to_fill[nb] = true;
                                       tbox::Pointer<hier::Patch> patch(
                                          new hier::Patch(
                                             src_patch->getMappedBox(),
                                             src_patch->getPatchDescriptor()));
                                       tbox::Pointer<hier::PatchGeometry>
                                       patch_geom(
                                          new hier::PatchGeometry(
                                             hier::IntVector(dim, 1),
                                             hier::PatchGeometry::TwoDimBool(dim, 0),
                                             hier::PatchGeometry::TwoDimBool(dim, 0)));
                                       patch->setPatchGeometry(patch_geom);
                                       patch->setPatchLevelNumber(dst_ln);

                                       for (int ne = 0; ne < num_classes; ne++) {

                                          for (tbox::List<int>::Iterator
                                               l(equiv_classes->getIterator(ne));
                                               l; l++) {

                                             const RefineClasses::Data& item =
                                                equiv_classes->getRefineItem(l());
                                             TBOX_ASSERT(item.d_class_id == ne);

                                             const int dst = item.d_dst;
                                             const int src = item.d_src;

                                             if (!patch->checkAllocated(dst)) {
                                                patch->allocatePatchData(dst);
                                                local_selector[nb][pc].setFlag(dst);
                                             }

                                             patch->getPatchData(dst)->
                                             copy(*(src_patch->getPatchData(src)),
                                                *(overlap_iter()));
                                          }
                                          overlap_iter++; 
                                       }
   
                                       singularity_patches[nb][pc].addItem(patch);
                                    }
                                 }
                              }
                           }
                        }
                        pc++;
                     }
                  }
                  nc++;
               }
            }
            if (do_neighbor_copy && do_fill_fine_gcw) {
               d_finalize_ghost_level->
               deallocatePatchData(src_vector);
            }
            d_neighbor_ghost_level->
            deallocatePatchData(dst_vector);

            for (int nb = 0; nb < nblocks; nb++) {
               if (singularity_to_fill[nb]) {
                  fillSingularityBoundary(
                     singularity_patches[nb],
                     nb, fill_time);
               }

               if (do_physical_boundary_fill) {
                  const hier::IntVector gcw = getBoundaryFillGhostWidth();
                  hier::MappedBoxSetSingleBlockIterator dst_local_iter(
                     d_multiblock_dst_level->getMappedBoxLevel()->getMappedBoxes(),
                     hier::BlockId(nb));

                  for ( ; dst_local_iter.isValid(); dst_local_iter++) {

                     const hier::MappedBoxId& mapped_box_id =
                        dst_local_iter->getId();

                     tbox::Pointer<hier::Patch> patch(
                        d_multiblock_dst_level->getPatch(mapped_box_id)); 

                     if (patch->getPatchGeometry()->
                         intersectsPhysicalBoundary()) {
                        d_multiblock_strategy->
                        setPhysicalBoundaryConditions(*patch,
                           fill_time,
                           gcw);
                     }
                  }
               }

               for (int j = 0; j < singularity_patches[nb].size(); j++) {
                  for (tbox::List<tbox::Pointer<hier::Patch> >::Iterator
                       sp(singularity_patches[nb][j]);
                       sp; sp++) {
                     sp()->deallocatePatchData(local_selector[nb][j]);
                  }
               }
            }
         }
      }
   }

   if (!filling_coarse_scratch) {
      d_multiblock_dst_level->deallocatePatchData(allocate_scr_vector);
   }

   if (filling_coarse_scratch) {
      if (!filling_crse_scr_recursive) {
         d_multiblock_strategy->setFillingCoarseScratch(false);
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

void MultiblockRefineSchedule::fillSingularityBoundary(
   tbox::Array<tbox::List<tbox::Pointer<hier::Patch> > >& singularity_patches,
   const int block_number,
   const double fill_time) const
{
   const tbox::Dimension& dim(d_multiblock_hierarchy->getDim());

   tbox::Pointer<hier::GridGeometry> grid_geometry(
      d_multiblock_hierarchy->getGridGeometry());

   hier::IntVector ratio(d_multiblock_dst_level->getRatioToLevelZero());

   for (hier::BoxList::Iterator sb(
           grid_geometry->getSingularityBoxList(block_number));
        sb;
        sb++) {
      hier::Box singularity(sb());
      if (d_multiblock_dst_level->getLevelNumber() != 0) {
         singularity.refine(ratio);
      }

      hier::IntVector gcw(
         d_multiblock_dst_level->getPatchDescriptor()->getMaxGhostWidth(dim));

      if (d_multiblock_strategy != NULL) {
         int pc = 0;
         hier::MappedBoxSetSingleBlockIterator dst_local_iter(
            d_multiblock_dst_level->getMappedBoxLevel()->getMappedBoxes(),
            hier::BlockId(block_number));

         for ( ; dst_local_iter.isValid(); dst_local_iter++) {

            const hier::MappedBoxId& mapped_box_id =
               dst_local_iter->getId();

            tbox::Pointer<hier::Patch> patch(
               d_multiblock_dst_level->getPatch(mapped_box_id));
            tbox::Pointer<hier::PatchGeometry> pgeom =
               patch->getPatchGeometry();

            tbox::Array<hier::BoundaryBox> nboxes =
               pgeom->getNodeBoundaries();

            if (nboxes.getSize()) {
               for (int bb = 0; bb < nboxes.getSize(); bb++) {
                  hier::Box intersection = (nboxes[bb].getBox())
                     * singularity;
                  if (!(intersection.empty())) {
                     hier::Box fill_box =
                        pgeom->getBoundaryFillBox(nboxes[bb],
                           patch->getBox(),
                           gcw);
                     d_multiblock_strategy->fillSingularityBoundaryConditions(
                        *patch, singularity_patches[pc],
                        fill_time, fill_box, nboxes[bb]);
                  }
               }
            }

            if (dim == tbox::Dimension(3)) {
               tbox::Array<hier::BoundaryBox> eboxes =
                  pgeom->getEdgeBoundaries();

               if (eboxes.getSize()) {
                  for (int bb = 0; bb < eboxes.getSize(); bb++) {
                     hier::Box intersection =
                        (eboxes[bb].getBox()) * singularity;
                     if (!(intersection.empty())) {
                        hier::Box fill_box =
                           pgeom->getBoundaryFillBox(eboxes[bb],
                              patch->getBox(),
                              gcw);
                        d_multiblock_strategy->
                        fillSingularityBoundaryConditions(
                           *patch, singularity_patches[pc],
                           fill_time, fill_box, eboxes[bb]);
                     }
                  }
               }
            }
            pc++;
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

void MultiblockRefineSchedule::copyBetweenBlocks(
   tbox::Pointer<hier::PatchLevel> dst_level,
   const tbox::Pointer<hier::PatchLevel> src_level,
   const hier::IntVector& shift,
   const hier::Transformation::RotationIdentifier rotate,
   const tbox::Pointer<xfer::RefineClasses> refine_classes,
   const hier::BlockId& dst_block_id,
   const hier::BlockId& src_block_id) const
{
   TBOX_ASSERT(!dst_level.isNull());
   TBOX_ASSERT(!src_level.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS4(*d_multiblock_hierarchy,
      *dst_level,
      *src_level,
      shift);

   tbox::ConstPointer<hier::MappedBoxLevel> dst_mapped_box_level(
      dst_level->getMappedBoxLevel());

   const int num_refine_items =
      refine_classes->getNumberOfRefineItems();

   for (hier::PatchLevel::Iterator p(src_level); p; p++) {
      tbox::Pointer<hier::Patch> src_patch(*p);
      const hier::MappedBoxId& src_mapped_box_id =
         src_patch->getMappedBox().getId();

      if (src_mapped_box_id.getBlockId() == src_block_id) {
         hier::MappedBoxId dst_mapped_box_id(src_mapped_box_id.getGlobalId(),
                                             dst_block_id);

         if (dst_mapped_box_level->hasMappedBox(dst_mapped_box_id)) {
            tbox::Pointer<hier::Patch> dst_patch(
               dst_level->getPatch(dst_mapped_box_id));

            for (int nd = 0; nd < num_refine_items; nd++) {

               const RefineClasses::Data& item = refine_classes->getRefineItem(nd);
               const int src = item.d_dst; //dst for src_level
               const int dst = item.d_src; //src for dst_level

               hier::Transformation::translateAndCopyData(
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
}

/*
 * ************************************************************************
 *                                                                        *
 * Call the routines that will fill data from one block to another.       *
 *                                                                        *
 * ************************************************************************
 */

void MultiblockRefineSchedule::fillBetweenBlocks(
   tbox::Pointer<hier::PatchLevel> dst_level,
   const tbox::Pointer<hier::PatchLevel> src_level,
   const hier::IntVector& shift,
   const hier::Transformation::RotationIdentifier rotate,
   const tbox::Pointer<xfer::RefineClasses> refine_classes,
   const hier::BlockId& dst_block_id,
   const hier::BlockId& src_block_id) const

{
   TBOX_ASSERT(!dst_level.isNull());
   TBOX_ASSERT(!src_level.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS4(*d_multiblock_hierarchy,
      *dst_level,
      *src_level,
      shift);

   tbox::ConstPointer<hier::MappedBoxLevel> dst_mapped_box_level(
      dst_level->getMappedBoxLevel());

   const int num_refine_items =
      refine_classes->getNumberOfRefineItems();

   for (hier::PatchLevel::Iterator p(src_level); p; p++) {
      tbox::Pointer<hier::Patch> src_patch(*p);
      const hier::MappedBoxId& src_mapped_box_id =
         src_patch->getMappedBox().getId();
      if (src_mapped_box_id.getBlockId() == src_block_id) {
         hier::MappedBoxId dst_mapped_box_id(src_mapped_box_id.getGlobalId(),
                                             dst_block_id);

         if (dst_mapped_box_level->hasMappedBox(dst_mapped_box_id)) {

            tbox::Pointer<hier::Patch> dst_patch(
               dst_level->getPatch(dst_mapped_box_id));

            for (int nd = 0; nd < num_refine_items; nd++) {

               const RefineClasses::Data& item = refine_classes->getRefineItem(nd);
               const int src = item.d_dst; //dst for src_level
               const int dst = item.d_src; //src for dst_level

               hier::Transformation::translateAndFillData(
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
}

/*
 * ************************************************************************
 *                                                                        *
 * Initialize a component selector to contain the source data components  *
 * for this schedule.                                                     *
 *                                                                        *
 * ************************************************************************
 */

void MultiblockRefineSchedule::initializeSourceVector(
   hier::ComponentSelector& allocate_vector) const
{
   if (!d_single_block_fill_local.isNull()) {
      d_single_block_fill_local->initializeSourceVector(allocate_vector);
      return;
   }

   if (!d_multiblock_coarse_schedule.isNull()) {
      d_multiblock_coarse_schedule->initializeSourceVector(
         allocate_vector);
      return;
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

const tbox::Pointer<xfer::RefineClasses>&
MultiblockRefineSchedule::getEquivalenceClasses() const
{
   return d_single_block_refine_alg->getEquivalenceClasses();
}

/*
 * ************************************************************************
 *                                                                        *
 * Determine if the fillData operation requires source data from blocks   *
 * other than the one containing the data being filled                    *
 *                                                                        *
 * ************************************************************************
 */

bool
MultiblockRefineSchedule::needOtherSourceBlocks(
   hier::BoxList& dst_boxes,
   hier::BoxList& src_boxes,
   hier::BoxList& domain_outside_block) const
{
   const tbox::Dimension& dim(d_multiblock_hierarchy->getDim());

   const hier::IntVector constant_one_intvector(dim, 1);

   bool needed;

   hier::BoxList dst_grow_list(dst_boxes);
   dst_grow_list.grow(constant_one_intvector);

   dst_grow_list.intersectBoxes(domain_outside_block);

   if (dst_grow_list.size() == 0) {
      needed = false;
   } else {

      if (dst_boxes.getTotalSizeOfBoxes() > src_boxes.getTotalSizeOfBoxes()) {

         hier::BoxList dst_remove_list(dst_boxes);
         dst_remove_list.removeIntersections(src_boxes);

         dst_remove_list.grow(constant_one_intvector);

         dst_remove_list.intersectBoxes(domain_outside_block);

         if (dst_remove_list.size() == 0) {

            needed = false;

         } else {

            needed = true;

         }
      } else {

         hier::BoxList intersection(src_boxes);
         intersection.intersectBoxes(dst_boxes);

         if (dst_boxes.getTotalSizeOfBoxes() >
             intersection.getTotalSizeOfBoxes()) {

            hier::BoxList dst_remove_list(dst_boxes);
            dst_remove_list.removeIntersections(src_boxes);

            dst_remove_list.grow(constant_one_intvector);

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

   return needed;

}

/*
 * ************************************************************************
 *                                                                        *
 * Create a schedule that will fill data on a temporary coarse level on a *
 * neighboring block.                                                     *
 *                                                                        *
 * ************************************************************************
 */

void MultiblockRefineSchedule::createNeighborCoarseSchedule(
   const tbox::Pointer<hier::PatchLevel>& neighbor_coarse_level,
   const int next_coarser_level,
   const hier::IntVector& ratio_to_coarser,
   const tbox::Pointer<hier::PatchHierarchy>& hierarchy)
{
   NULL_USE(ratio_to_coarser);
   NULL_USE(hierarchy);

   tbox::Pointer<hier::GridGeometry> grid_geometry(
      d_multiblock_hierarchy->getGridGeometry());

   tbox::Pointer<hier::PatchLevel> coarse_hierarchy_level =
      d_multiblock_hierarchy->getPatchLevel(next_coarser_level);

   bool unfilled_boxes_exist = false;

   for (int nb = 0; nb < grid_geometry->getNumberBlocks(); nb++) {

      hier::BoxList domain_outside_block;

      tbox::List<hier::GridGeometry::Neighbor> neighbors =
         grid_geometry->getNeighbors(nb);

      for (tbox::List<hier::GridGeometry::Neighbor>::
           Iterator nei(neighbors); nei; nei++) {
         domain_outside_block.unionBoxes(hier::BoxList(nei().getTranslatedDomain()));
      }

      domain_outside_block.refine(neighbor_coarse_level->getRatioToLevelZero());

      hier::BoxList pseudo_domain(domain_outside_block);
      pseudo_domain.unionBoxes(
         hier::BoxList(neighbor_coarse_level->getPhysicalDomain(hier::BlockId(nb))));

      hier::BoxList unfilled_boxes;
      neighbor_coarse_level->getBoxes(unfilled_boxes, hier::BlockId(nb));

      findUnfilledBoxes(unfilled_boxes,
         nb,
         coarse_hierarchy_level,
         pseudo_domain,
         d_data_fill_gcw);

      if (unfilled_boxes.size()) {
         unfilled_boxes_exist = true;
         break;
      }
   }

   d_neighbor_multiblock_coarse_level =
      neighbor_coarse_level;

   grid_geometry->adjustMultiblockPatchLevelBoundaries(
      *d_neighbor_multiblock_coarse_level);

   tbox::Pointer<PatchLevelFullFillPattern> fill_pattern(
      new PatchLevelFullFillPattern());

   if (unfilled_boxes_exist) {
      d_neighbor_multiblock_coarse_schedule =
            new MultiblockRefineSchedule(
               fill_pattern,
               d_neighbor_multiblock_coarse_level,
               coarse_hierarchy_level,
               d_multiblock_hierarchy,
               d_single_block_scratch_refine_alg,
               d_transaction_factory,
               d_multiblock_strategy,
               true,
               d_enable_singularity_patches);
   } else {
      d_neighbor_multiblock_coarse_schedule =
            new MultiblockRefineSchedule(
               fill_pattern,
               d_neighbor_multiblock_coarse_level,
               coarse_hierarchy_level,
               next_coarser_level - 1,
               d_multiblock_hierarchy,
               d_single_block_scratch_refine_alg,
               d_transaction_factory,
               d_multiblock_strategy,
               true,
               d_enable_singularity_patches);

   }
}

/*
 * ************************************************************************
 *                                                                        *
 * Create a schedule to fill data on a temporary coarse level             *
 *                                                                        *
 * ************************************************************************
 */

void MultiblockRefineSchedule::createCoarseSchedule(
   const tbox::Pointer<hier::PatchLevel>& coarse_level,
   const int next_coarser_level)
{
   tbox::Pointer<hier::GridGeometry> grid_geometry(
      d_multiblock_hierarchy->getGridGeometry());

   d_multiblock_coarse_scratch_level = coarse_level;

   grid_geometry->adjustMultiblockPatchLevelBoundaries(
      *d_multiblock_coarse_scratch_level);

   tbox::Pointer<hier::PatchLevel> coarse_hierarchy_level =
      d_multiblock_hierarchy->getPatchLevel(next_coarser_level);

   bool unfilled_boxes_exist = false;

   for (int nb = 0; nb < grid_geometry->getNumberBlocks(); nb++) {

      hier::BoxList domain_outside_block;

      tbox::List<hier::GridGeometry::Neighbor> neighbors =
         grid_geometry->getNeighbors(nb);

      for (tbox::List<hier::GridGeometry::Neighbor>::
           Iterator nei(neighbors); nei; nei++) {
         domain_outside_block.unionBoxes(hier::BoxList(nei().getTranslatedDomain()));
      }

      domain_outside_block.refine(coarse_level->getRatioToLevelZero());

      hier::BoxList pseudo_domain(domain_outside_block);
      pseudo_domain.unionBoxes(
         hier::BoxList(coarse_level->getPhysicalDomain(hier::BlockId(nb))));
   
      hier::BoxList unfilled_boxes;
      coarse_level->getBoxes(unfilled_boxes, hier::BlockId(nb));

      findUnfilledBoxes(unfilled_boxes,
         nb,
         coarse_hierarchy_level,
         pseudo_domain,
         d_data_fill_gcw);

      if (unfilled_boxes.size()) {
         unfilled_boxes_exist = true;
         break;
      }
   }

   tbox::Pointer<PatchLevelFullFillPattern> fill_pattern(
      new PatchLevelFullFillPattern());

   if (!unfilled_boxes_exist) {
      d_multiblock_coarse_schedule =
         new MultiblockRefineSchedule(
            fill_pattern,
            d_multiblock_coarse_scratch_level,
            coarse_hierarchy_level,
            d_multiblock_hierarchy,
            d_single_block_scratch_refine_alg,
            d_transaction_factory,
            d_multiblock_strategy,
            true,
            d_enable_singularity_patches);
   } else {
      d_multiblock_coarse_schedule =
         new MultiblockRefineSchedule(
            fill_pattern,
            d_multiblock_coarse_scratch_level,
            coarse_hierarchy_level,
            next_coarser_level - 1,
            d_multiblock_hierarchy,
            d_single_block_scratch_refine_alg,
            d_transaction_factory,
            d_multiblock_strategy,
            true,
            d_enable_singularity_patches);
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

void MultiblockRefineSchedule::findUnfilledBoxes(
   hier::BoxList& unfilled_boxes,
   const int block_number,
   tbox::Pointer<hier::PatchLevel> coarse_hierarchy_level,
   const hier::BoxList& pseudo_domain,
   const hier::IntVector& gcw)
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(*d_multiblock_hierarchy,
      *coarse_hierarchy_level,
      gcw);

   unfilled_boxes.grow(gcw);
   unfilled_boxes.intersectBoxes(pseudo_domain);

   hier::BoxList coarse_hierarchy_boxes;
   coarse_hierarchy_level->getBoxes(coarse_hierarchy_boxes,
                                    hier::BlockId(block_number));

   unfilled_boxes.removeIntersections(coarse_hierarchy_boxes);

   tbox::Pointer<hier::GridGeometry> grid_geometry(
      d_multiblock_hierarchy->getGridGeometry());

   if (unfilled_boxes.size()) {
      for (tbox::List<hier::GridGeometry::Neighbor>::
           Iterator ni(grid_geometry->getNeighbors(block_number));
           ni; ni++) {

         int id = ni().getBlockNumber();

         hier::BoxList level_boxlist(grid_geometry->getDim());
         coarse_hierarchy_level->getBoxes(level_boxlist, hier::BlockId(id));

         if (level_boxlist.size()) {

            grid_geometry->translateBoxList(level_boxlist,
               coarse_hierarchy_level->getRatioToLevelZero(),
               hier::BlockId(block_number),
               hier::BlockId(id));

            unfilled_boxes.removeIntersections(level_boxlist);
         }

         if (unfilled_boxes.size() == 0) {
            break;
         }

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

void MultiblockRefineSchedule::constructScratchRefineAlgorithm()
{
   const tbox::Dimension& dim(d_multiblock_hierarchy->getDim());

   d_single_block_scratch_refine_alg = new xfer::RefineAlgorithm(dim);

   const tbox::Pointer<RefineClasses> refine_classes =
      d_single_block_refine_alg->getEquivalenceClasses();

   tbox::Pointer<RefineClasses> scratch_refine_classes(new RefineClasses());

   tbox::Pointer<BoxGeometryVariableFillPattern> bg_fill_pattern(
      new BoxGeometryVariableFillPattern());

   const int num_refine_items = refine_classes->getNumberOfRefineItems();

   for (int nd = 0; nd < num_refine_items; nd++) {

      RefineClasses::Data data = refine_classes->getRefineItem(nd);

      data.d_dst = data.d_scratch;
      data.d_var_fill_pattern = bg_fill_pattern;

      scratch_refine_classes->insertEquivalenceClassItem(data);

   }

   d_single_block_scratch_refine_alg->setEquivalenceClasses(
      scratch_refine_classes);

}

/*
 * ************************************************************************
 *                                                                        *
 * Copy data from scratch to destination where there is overlap between   *
 * the data.  Nothing is done if scratch and destination are the same.    *
 *                                                                        *
 * ************************************************************************
 */

void
MultiblockRefineSchedule::copyScratchToDestination(
   tbox::Pointer<hier::PatchLevel> level,
   const hier::BoxList& unfilled_boxes,
   const hier::BlockId& block_id,
   tbox::Pointer<xfer::RefineClasses> refine_classes) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*d_multiblock_hierarchy, *level);

   const tbox::Dimension& dim(d_multiblock_hierarchy->getDim());

   const int num_classes = refine_classes->getNumberOfEquivalenceClasses();
   const int num_refine_items = refine_classes->getNumberOfRefineItems();

   tbox::Pointer<hier::PatchDescriptor> descriptor =
      level->getPatchDescriptor();

   bool copy_needed = false;
   for (int nd = 0; nd < num_refine_items; nd++) {

      const xfer::RefineClasses::Data& data = refine_classes->getRefineItem(nd);

      const int dst = data.d_dst;
      const int src = data.d_scratch;

      if (dst != src) {
         copy_needed = true;
         break;
      }
   }

   if (copy_needed) {
      hier::Transformation zero_trans(hier::IntVector::getZero(dim));
      for (hier::PatchLevel::Iterator p(level); p; p++) {

         tbox::Pointer<hier::Patch> patch = *p;
         const hier::BlockId& patch_block = patch->getMappedBox().getBlockId();

         if (patch_block != block_id) {
            continue;
         }

         hier::Box patch_box(patch->getBox());

         for (int ne = 0; ne < num_classes; ne++) {

            const xfer::RefineClasses::Data& rep_item =
               refine_classes->getClassRepresentative(ne);

            const int rep_item_src_id = rep_item.d_scratch;
            const int rep_item_dst_id = rep_item.d_scratch;

            tbox::Pointer<hier::PatchDataFactory> src_pdf =
               descriptor->getPatchDataFactory(rep_item_src_id);
            tbox::Pointer<hier::PatchDataFactory> dst_pdf =
               descriptor->getPatchDataFactory(rep_item_dst_id);

            for (hier::BoxList::Iterator b(unfilled_boxes);
                 b; b++) {

               const hier::Box fill_box(b());

               const hier::Box src_mask(fill_box);

               if (!src_mask.empty()) {
                  tbox::Pointer<hier::BoxOverlap> overlap =
                     rep_item.d_var_fill_pattern->calculateOverlap(
                        *dst_pdf->getBoxGeometry(fill_box),
                        *src_pdf->getBoxGeometry(patch_box),
                        patch_box,
                        src_mask,
                        fill_box, 
                        true, zero_trans);

                  for (tbox::List<int>::Iterator
                       l(refine_classes->getIterator(ne)); l; l++) {

                     const RefineClasses::Data& item =
                        refine_classes->getRefineItem(l());
                     const int dst = item.d_dst;
                     const int src = item.d_scratch;

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

void MultiblockRefineSchedule::refineScratchData(
   tbox::Pointer<hier::PatchLevel> coarse_mb_level,
   tbox::Pointer<hier::PatchLevel> fine_level,
   const std::vector<std::vector< tbox::Pointer<hier::BoxOverlap> > >& overlaps,
   const hier::Connector& crs_to_dst,
   const hier::BoxList& unfilled_boxes,
   const int block_number,
   const bool fill_fine_gcw,
   const bool filling_neighbor) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(*d_multiblock_hierarchy,
      *coarse_mb_level,
      *fine_level);

   const tbox::SAMRAI_MPI& mpi(fine_level->getMappedBoxLevel()->getMPI());
   const tbox::Dimension& dim(d_multiblock_hierarchy->getDim());
   const hier::BlockId block_id(block_number);

   const tbox::Pointer<xfer::RefineClasses> refine_classes =
      d_single_block_refine_alg->getEquivalenceClasses();

   const int num_refine_items = refine_classes->getNumberOfRefineItems();

   const int rank = mpi.getRank();

   const hier::IntVector ratio =
      fine_level->getRatioToLevelZero() / coarse_mb_level->getRatioToLevelZero();

   hier::IntVector gcw(hier::IntVector::getZero(dim));

   if (fill_fine_gcw) {
      gcw = fine_level->getPatchDescriptor()->getMaxGhostWidth(dim);
   }

   tbox::Pointer<VariableFillPattern> bg_fill_pattern;
   if (filling_neighbor) {
      bg_fill_pattern = new BoxGeometryVariableFillPattern();
   }

   const hier::MappedBoxSet& coarse_mapped_boxes =
      coarse_mb_level->getMappedBoxLevel()->getMappedBoxes();
   TBOX_ASSERT(overlaps.size() == coarse_mapped_boxes.size()); 
   int oc = 0;
   for (hier::MappedBoxSetSingleBlockIterator ni(coarse_mapped_boxes, block_id);
        ni.isValid(); ++ni) {

      if (overlaps[oc].size()) {
         const hier::MappedBox& coarse_mapped_box = *ni;

         const hier::MappedBoxSet& dst_nabrs =
            crs_to_dst.getNeighborSet(coarse_mapped_box.getId());

         for (hier::MappedBoxSet::const_iterator dn = dst_nabrs.begin();
              dn != dst_nabrs.end(); ++dn) {

            const hier::MappedBox& dst_mapped_box = *dn;

            if (dst_mapped_box.getOwnerRank() != rank) {
               continue;
            }

            tbox::Pointer<hier::Patch> fine_patch =
               fine_level->getPatch(dst_mapped_box.getId());
            tbox::Pointer<hier::Patch> crse_patch =
               coarse_mb_level->getPatch(coarse_mapped_box.getId());

            hier::Box fine_box(fine_patch->getBox());
            fine_box.grow(gcw);

            hier::BoxList fill_boxes(unfilled_boxes);
            fill_boxes.intersectBoxes(fine_box);

            if (fill_boxes.size() > 0) {

               hier::Box refined_coarse_box(crse_patch->getBox());
               refined_coarse_box.refine(ratio);
               fill_boxes.intersectBoxes(refined_coarse_box);

               fill_boxes.coalesceBoxes();

               if (d_multiblock_strategy != NULL) {
                  d_multiblock_strategy->preprocessRefineBoxes(*fine_patch,
                     *crse_patch,
                     fill_boxes,
                     ratio);
               }

               for (int iri = 0; iri < num_refine_items; iri++) {
                  const xfer::RefineClasses::Data& ref_item =
                     refine_classes->getRefineItem(iri);

                  if (!(ref_item.d_oprefine.isNull())) {

                     tbox::Pointer<hier::BoxOverlap> refine_overlap(
                        overlaps[oc][ref_item.d_class_id]);

                     const int scratch_id = ref_item.d_scratch;

                     ref_item.d_oprefine->refine(*fine_patch, *crse_patch,
                        scratch_id, scratch_id,
                        *refine_overlap, ratio);

                  }
               }

               if (d_multiblock_strategy != NULL) {
                  d_multiblock_strategy->postprocessRefineBoxes(*fine_patch,
                     *crse_patch,
                     fill_boxes,
                     ratio);

               }
            }
         }
      }
      oc++;
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

tbox::Pointer<hier::BoxOverlap>
MultiblockRefineSchedule::calculateOverlap(
   const hier::Patch& dst_patch,
   const hier::Patch& src_patch,
   const xfer::RefineClasses::Data& item,
   const hier::Box& restrict_box) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(*d_multiblock_hierarchy, dst_patch, src_patch);

   const tbox::Dimension& dim(d_multiblock_hierarchy->getDim());

   tbox::Pointer<hier::PatchDescriptor> dst_patch_descriptor(
      dst_patch.getPatchDescriptor());
   tbox::Pointer<hier::PatchDescriptor> src_patch_descriptor(
      src_patch.getPatchDescriptor());

   const int item_dst_id = item.d_dst;
   const int item_src_id = item.d_src;

   tbox::Pointer<hier::PatchDataFactory> src_pdf(
      src_patch_descriptor->getPatchDataFactory(item_src_id));
   tbox::Pointer<hier::PatchDataFactory> dst_pdf(
      dst_patch_descriptor->getPatchDataFactory(item_dst_id));

   hier::Box dst_ghost_box(dst_patch.getBox());
   dst_ghost_box.grow(dst_pdf->getGhostCellWidth());
   hier::Box dst_fill_box(dst_ghost_box * src_patch.getBox() * restrict_box);

   bool overwrite_interior = true;

   tbox::Pointer<hier::BoxOverlap> overlap =
      item.d_var_fill_pattern->calculateOverlap(
         *dst_pdf->getBoxGeometry(dst_patch.getBox()),
         *src_pdf->getBoxGeometry(src_patch.getBox()),
         dst_patch.getBox(),
         dst_fill_box,
         dst_fill_box,
         overwrite_interior,
         hier::Transformation(hier::IntVector::getZero(dim)));

   return overlap;
}

/*
 * ************************************************************************
 *                                                                        *
 * Calculate the overlap used to fill data on a singularity patch for a   *
 * single equivalence class.                                              *
 *                                                                        *
 * ************************************************************************
 */

tbox::Pointer<hier::BoxOverlap>
MultiblockRefineSchedule::calculateSingularityPatchOverlap(
   const hier::Patch& patch,
   const xfer::RefineClasses::Data& item,
   const hier::Box& restrict_box) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*d_multiblock_hierarchy, patch);

   const tbox::Dimension& dim(d_multiblock_hierarchy->getDim());

   tbox::Pointer<hier::PatchDescriptor> patch_descriptor(
      patch.getPatchDescriptor());

   const int item_dst_id = item.d_dst;
   const int item_src_id = item.d_src;

   tbox::Pointer<hier::PatchDataFactory> src_pdf(
      patch_descriptor->getPatchDataFactory(item_src_id));
   tbox::Pointer<hier::PatchDataFactory> dst_pdf(
      patch_descriptor->getPatchDataFactory(item_dst_id));

   hier::Box dst_ghost_box(patch.getBox());
   dst_ghost_box.grow(dst_pdf->getGhostCellWidth());
   hier::Box dst_fill_box(dst_ghost_box * patch.getBox() * restrict_box);

   bool overwrite_interior = true;

   tbox::Pointer<xfer::BoxGeometryVariableFillPattern> fill_pattern(
      new xfer::BoxGeometryVariableFillPattern());

   tbox::Pointer<hier::BoxOverlap> overlap =
      fill_pattern->calculateOverlap(
         *dst_pdf->getBoxGeometry(patch.getBox()),
         *src_pdf->getBoxGeometry(patch.getBox()),
         patch.getBox(),
         dst_fill_box,
         dst_fill_box,
         overwrite_interior,
         hier::Transformation(hier::IntVector::getZero(dim)));

   return overlap;
}


/*
 * ************************************************************************
 *                                                                        *
 * Create and store overlaps that will be used in refineScratchData.      *
 *                                                                        *
 * ************************************************************************
 */
void MultiblockRefineSchedule::createRefineOverlaps(
   std::vector<std::vector< tbox::Pointer<hier::BoxOverlap> > >& overlaps,
   tbox::Pointer<hier::PatchLevel> coarse_mb_level,
   tbox::Pointer<hier::PatchLevel> fine_level,
   const hier::Connector& crs_to_dst,
   const hier::BoxList& unfilled_boxes,
   const int block_number,
   const bool fill_fine_gcw,
   const bool filling_neighbor)
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(*d_multiblock_hierarchy,
      *coarse_mb_level,
      *fine_level);

   const tbox::SAMRAI_MPI& mpi(fine_level->getMappedBoxLevel()->getMPI());
   const tbox::Dimension& dim(d_multiblock_hierarchy->getDim());
   const hier::BlockId block_id(block_number);

   tbox::Pointer<hier::PatchDescriptor> fine_patch_descriptor =
      fine_level->getPatchDescriptor();

   const tbox::Pointer<xfer::RefineClasses> refine_classes =
      d_single_block_refine_alg->getEquivalenceClasses();

   const int num_classes =
      refine_classes->getNumberOfEquivalenceClasses();

   const int rank = mpi.getRank();

   const hier::IntVector ratio =
      fine_level->getRatioToLevelZero() / coarse_mb_level->getRatioToLevelZero();

   hier::IntVector gcw(hier::IntVector::getZero(dim));

   if (fill_fine_gcw) {
      gcw = fine_level->getPatchDescriptor()->getMaxGhostWidth(dim);
   }

   tbox::Pointer<VariableFillPattern> bg_fill_pattern;
   if (filling_neighbor) {
      bg_fill_pattern = new BoxGeometryVariableFillPattern();
   }

   const hier::MappedBoxSet& coarse_mapped_boxes =
      coarse_mb_level->getMappedBoxLevel()->getMappedBoxes();
   overlaps.resize(coarse_mapped_boxes.size());
   int oc = 0;
   for (hier::MappedBoxSetSingleBlockIterator ni(coarse_mapped_boxes, block_id);
        ni.isValid(); ++ni) {
      const hier::MappedBox& coarse_mapped_box = *ni;

      const hier::MappedBoxSet& dst_nabrs =
         crs_to_dst.getNeighborSet(coarse_mapped_box.getId());

      TBOX_ASSERT(dst_nabrs.size() == 1 || dst_nabrs.size() == 0);

      for (hier::MappedBoxSet::const_iterator dn = dst_nabrs.begin();
           dn != dst_nabrs.end(); ++dn) {

         const hier::MappedBox& dst_mapped_box = *dn;

         if (dst_mapped_box.getOwnerRank() != rank) {
            continue;
         }

         tbox::Pointer<hier::Patch> fine_patch =
            fine_level->getPatch(dst_mapped_box.getId());
         tbox::Pointer<hier::Patch> crse_patch =
            coarse_mb_level->getPatch(coarse_mapped_box.getId());

         hier::Box fine_box(fine_patch->getBox());
         fine_box.grow(gcw);

         hier::BoxList fill_boxes(unfilled_boxes);
         fill_boxes.intersectBoxes(fine_box);

         if (fill_boxes.size() > 0) {

            hier::Box refined_coarse_box(crse_patch->getBox());
            refined_coarse_box.refine(ratio);
            fill_boxes.intersectBoxes(refined_coarse_box);

            fill_boxes.coalesceBoxes();

            overlaps[oc].resize(num_classes);

            for (int ne = 0; ne < num_classes; ne++) {

               const xfer::RefineClasses::Data& rep_item =
                  refine_classes->getClassRepresentative(ne);

               tbox::Pointer<hier::BoxOverlap> refine_overlap;

               if (!(rep_item.d_oprefine.isNull())) {

                  const tbox::Pointer<VariableFillPattern>& var_fill_pattern =
                     filling_neighbor ?
                     bg_fill_pattern : rep_item.d_var_fill_pattern;

                  const int scratch_id = rep_item.d_scratch;
                  tbox::Pointer<hier::PatchDataFactory> fine_pdf =
                     fine_patch_descriptor->getPatchDataFactory(scratch_id);
                  const hier::IntVector& fine_scratch_gcw =
                     fine_pdf->getGhostCellWidth();
                  hier::Box fine_data_box(fine_patch->getBox());
                  fine_data_box.grow(fine_scratch_gcw);
                  const hier::Box scratch_space(
                     fine_data_box * refined_coarse_box);

                  refine_overlap =
                     var_fill_pattern->computeFillBoxesOverlap(
                        fill_boxes,
                        fine_patch->getBox(),
                        scratch_space,
                        *fine_pdf);

               } else {
                  refine_overlap.setNull();
               }
               overlaps[oc][ne] = refine_overlap;
            }
         }
      }
      oc++;
   }
}


/*
 * ************************************************************************
 *                                                                        *
 * Allocate scratch space on a patch level.                               *
 *                                                                        *
 * ************************************************************************
 */

void MultiblockRefineSchedule::allocateScratchSpace(
   hier::ComponentSelector& allocate_vector,
   tbox::Pointer<hier::PatchLevel> level,
   double fill_time) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*d_multiblock_hierarchy, *level);

   if (!d_single_block_fill_local.isNull()) {
      d_single_block_fill_local->allocateScratchSpace(
         allocate_vector, level, fill_time);
      return;
   }

   if (!d_multiblock_coarse_schedule.isNull()) {
      d_multiblock_coarse_schedule->allocateScratchSpace(
         allocate_vector, level, fill_time);
      return;
   }
}

/*
 * ************************************************************************
 *                                                                        *
 * Get the ghost width needed for boundary filling.                       *
 *                                                                        *
 * ************************************************************************
 */

hier::IntVector
MultiblockRefineSchedule::getBoundaryFillGhostWidth() const
{
   const tbox::Dimension& dim(d_multiblock_hierarchy->getDim());

   if (!d_single_block_fill_local.isNull()) {
      return d_single_block_fill_local->getBoundaryFillGhostWidth();
   }

   if (!d_multiblock_coarse_schedule.isNull()) {
      return d_multiblock_coarse_schedule->getBoundaryFillGhostWidth();
   }

   TBOX_ASSERT(0);

   return hier::IntVector(dim, -1);
}

/*
 * ************************************************************************
 *                                                                        *
 * Initialize a ComponentSelector to store all destination components for *
 * this schedule.                                                         *
 *                                                                        *
 * ************************************************************************
 */

void MultiblockRefineSchedule::initializeDestinationVector(
   hier::ComponentSelector& dst_vector) const
{
   dst_vector.clrAllFlags();

   const tbox::Pointer<xfer::RefineClasses> refine_classes =
      d_single_block_refine_alg->getEquivalenceClasses();

   const int num_refine_items = refine_classes->getNumberOfRefineItems();

   for (int nd = 0; nd < num_refine_items; nd++) {

      dst_vector.setFlag(refine_classes->getRefineItem(nd).d_dst);

   }
}

/*
 *************************************************************************
 *
 * Reset schedule with new set of refine items.
 *
 ************************************************************************
 */

void MultiblockRefineSchedule::reset(
   const tbox::Pointer<xfer::RefineAlgorithm> refine_algorithm)
{
   TBOX_ASSERT(!refine_algorithm.isNull());

   d_single_block_refine_alg = refine_algorithm;
   constructScratchRefineAlgorithm();

   tbox::Pointer<xfer::RefineClasses> refine_classes(
      d_single_block_refine_alg->getEquivalenceClasses());

   if (!d_single_block_fill_local.isNull()) {
      d_single_block_fill_local->reset(refine_classes);
   }

   if (!d_multiblock_coarse_schedule.isNull()) {
      d_multiblock_coarse_schedule->reset(
         d_single_block_scratch_refine_alg);
   }

   if (!d_neighbor_multiblock_coarse_schedule.isNull()) {
      d_neighbor_multiblock_coarse_schedule->reset(
         d_single_block_scratch_refine_alg);
   }

   if (!d_neighbor_single_block_refine_schedule.isNull()) {
      d_neighbor_single_block_refine_schedule->reset(refine_classes);
   }

}

}
}
#endif
