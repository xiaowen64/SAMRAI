/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Multiblock binary trees of Boxes for overlap searches.
 *
 ************************************************************************/

#ifndef included_hier_MultiblockBoxTree_C
#define included_hier_MultiblockBoxTree_C

#include "SAMRAI/hier/MultiblockBoxTree.h"
#include "SAMRAI/hier/GridGeometry.h"
#include "SAMRAI/hier/BoxContainerConstIterator.h"
#include "SAMRAI/hier/BoxContainerIterator.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"
#include "SAMRAI/tbox/Statistician.h"
#include "SAMRAI/tbox/TimerManager.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace hier {

/*
 *************************************************************************
 *************************************************************************
 */

MultiblockBoxTree::MultiblockBoxTree(
   const tbox::ConstPointer<GridGeometry>& grid_geometry,
   const BoxSet& mapped_boxes,
   size_t min_number):
   d_grid_geometry(grid_geometry)
{
   NULL_USE(min_number);

   const tbox::Dimension& dim(d_grid_geometry->getDim());

   /*
    * Group Boxes by their BlockId and
    * create a tree for each BlockId.
    */

   std::map<BlockId, BoxSet> mapped_boxes_by_block;
   for (BoxSet::ConstIterator bi = mapped_boxes.begin();
        bi != mapped_boxes.end(); ++bi) {

      TBOX_ASSERT((*bi).getId().isValid());
      const BlockId& block_id = (*bi).getBlockId();
      mapped_boxes_by_block[block_id].order();
      mapped_boxes_by_block[block_id].insert(
         mapped_boxes_by_block[block_id].end(), *bi);
   }

   for (std::map<BlockId, BoxSet>::iterator blocki = mapped_boxes_by_block.begin();
        blocki != mapped_boxes_by_block.end(); ++blocki) {

      TBOX_ASSERT(blocki->first.getBlockValue() >= 0 &&
         blocki->first.getBlockValue() < grid_geometry->getNumberBlocks());

      /*
       * The following lines do this:
       * d_single_block_trees[blocki->first].generateTree(blocki->second, min_number);
       *
       * We cannot do it the concise way because it
       * requires the default constructor for BoxTree.
       */
      const std::pair<std::map<BlockId, BoxTree>::iterator, bool> insert_return_value(
         d_single_block_trees.insert(
            std::pair<BlockId, BoxTree>(
               blocki->first,
               BoxTree(dim, blocki->second))));
      TBOX_ASSERT(insert_return_value.second);
   }
}

/*
 *************************************************************************
 *************************************************************************
 */
#if 0
MultiblockBoxTree::MultiblockBoxTree(
   const tbox::ConstPointer<GridGeometry>& grid_geometry,
   const BoxSet& mapped_boxes,
   size_t min_number):
   d_grid_geometry(grid_geometry)
{
   generateTree(grid_geometry, mapped_boxes, min_number);
}
#endif
/*
 *************************************************************************
 *************************************************************************
 */

MultiblockBoxTree::MultiblockBoxTree(
   const tbox::ConstPointer<GridGeometry>& grid_geometry,
   const std::map<BlockId, BoxList>& boxes,
   size_t min_number):
   d_grid_geometry(grid_geometry)
{
   generateTree(grid_geometry, boxes, min_number);
}

/*
 *************************************************************************
 *************************************************************************
 */

MultiblockBoxTree::MultiblockBoxTree():
   d_grid_geometry(NULL)
{
}

/*
 *************************************************************************
 *************************************************************************
 */

MultiblockBoxTree::~MultiblockBoxTree()
{
}

/*
 *************************************************************************
 * Generate the tree from a given mutable vector of mapped_boxes.
 * The vector will be changed and its output state is undefined.
 *************************************************************************
 */
void MultiblockBoxTree::generateTree(
   const tbox::ConstPointer<GridGeometry>& grid_geometry,
   const BoxSet& mapped_boxes,
   size_t min_number)
{
   NULL_USE(min_number);

   d_grid_geometry = grid_geometry;

   /*
    * Group Boxes by their BlockId and create a tree for each
    * BlockId.
    */
   std::map<BlockId, BoxSet> mapped_boxes_by_block;
   for (BoxSet::ConstIterator bi = mapped_boxes.begin();
        bi != mapped_boxes.end(); ++bi) {

      const BlockId& block_id = bi->getBlockId();

      std::map<BlockId, BoxSet>::iterator iter =
         mapped_boxes_by_block.find(block_id);

      if (iter != mapped_boxes_by_block.end()) {
         iter->second.insert(iter->second.end(), *bi);
      } else {
         hier::BoxContainer boxes(*bi, true);
         mapped_boxes_by_block.insert(
            std::pair<BlockId, BoxContainer>(block_id, boxes));
      }
   }

   for (std::map<BlockId, BoxSet>::iterator blocki = mapped_boxes_by_block.begin();
        blocki != mapped_boxes_by_block.end(); ++blocki) {

      TBOX_ASSERT(blocki->first.getBlockValue() >= 0 &&
         blocki->first.getBlockValue() < grid_geometry->getNumberBlocks());

      /*
       * The following lines do this:
       * d_single_block_trees[blocki->first].generateTree(blocki->second, min_number);
       *
       * We cannot do it the concise way because it requires the
       * default constructor for BoxTree.
       */
      const std::pair<std::map<BlockId, BoxTree>::iterator, bool> insert_return_value(
         d_single_block_trees.insert(
            std::pair<BlockId, BoxTree>(
               blocki->first,
               BoxTree(grid_geometry->getDim(), blocki->second))));
      TBOX_ASSERT(insert_return_value.second);
   }
}

/*
 *************************************************************************
 * Generate the tree from a given mutable vector of mapped_boxes.
 * The vector will be changed and its output state is undefined.
 *************************************************************************
 */
void MultiblockBoxTree::generateTree(
   const tbox::ConstPointer<GridGeometry>& grid_geometry,
   const std::map<BlockId, BoxList>& boxes,
   size_t min_number)
{
   d_grid_geometry = grid_geometry;

   for (std::map<BlockId, BoxList>::const_iterator blocki = boxes.begin();
        blocki != boxes.end(); ++blocki) {

      TBOX_ASSERT(blocki->first.getBlockValue() >= 0 &&
         blocki->first.getBlockValue() < grid_geometry->getNumberBlocks());

      /*
       * The following lines do this:
       * d_single_block_trees[blocki->first].generateTree(blocki->second, min_number);
       *
       * We cannot do it the concise way because it requires the
       * default constructor for BoxTree.
       */
      const std::pair<std::map<BlockId, BoxTree>::iterator, bool> insert_return_value(
         d_single_block_trees.insert(
            std::pair<BlockId, BoxTree>(
               blocki->first,
               BoxTree(grid_geometry->getDim(), blocki->second,
                  blocki->first, min_number))));
      TBOX_ASSERT(insert_return_value.second);
   }
}

/*
 *************************************************************************
 *************************************************************************
 */

void MultiblockBoxTree::clear()
{
   d_single_block_trees.clear();
   d_grid_geometry.setNull();
}

/*
 *************************************************************************
 *************************************************************************
 */

bool MultiblockBoxTree::isInitialized() const
{
   return !d_grid_geometry.isNull();
}

/*
 *************************************************************************
 *************************************************************************
 */

const tbox::ConstPointer<GridGeometry>& MultiblockBoxTree::getGridGeometry() const
{
   return d_grid_geometry;
}

/*
 *************************************************************************
 *************************************************************************
 */

bool MultiblockBoxTree::hasBoxInBlock(
   const BlockId& block_id) const
{
   return d_single_block_trees.find(block_id) !=
          d_single_block_trees.end();
}

/*
 *************************************************************************
 *************************************************************************
 */

const BoxTree& MultiblockBoxTree::getSingleBlockBoxTree(
   const BlockId& block_id) const
{
   std::map<BlockId, BoxTree>::const_iterator mi =
      d_single_block_trees.find(block_id);

   if (mi == d_single_block_trees.end()) {
      TBOX_ERROR("hier::MultiblockBoxTree::getSingleBlockBoxTree: cannot\n"
         << "return the single-block BoxTree for block " << block_id
         << "\neither because there is no such block in the domain configuration\n"
         << "or no Boxes from that block exists in the tree.\n"
         << "You can use MultiblockBoxTree::hasBoxInBlock()\n"
         << "to determine whether the tree exists.");
   }

   return mi->second;
}

/*
 *************************************************************************
 *************************************************************************
 */
#if 0
void MultiblockBoxTree::findOverlapBoxes(
   BoxSet& overlap_mapped_boxes,
   const Box& box,
   const BlockId& block_id,
   const IntVector& refinement_ratio,
   bool include_singularity_block_neighbors) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(*d_grid_geometry, box, refinement_ratio);

   TBOX_ASSERT(block_id.getBlockValue() >= 0 &&
      block_id.getBlockValue() < d_grid_geometry->getNumberBlocks());

   /*
    * Search in the index space of block_id for overlaps.
    */

   std::map<BlockId, BoxTree>::const_iterator
   blocki(d_single_block_trees.find(block_id));

   if (blocki != d_single_block_trees.end()) {
      blocki->second.findOverlapBoxes(overlap_mapped_boxes, box);
   }

   /*
    * Search in the index spaces neighboring block_id for overlaps.
    */

   const tbox::List<GridGeometry::Neighbor>& block_neighbors(
      d_grid_geometry->getNeighbors(block_id));

   for (tbox::ListIterator<GridGeometry::Neighbor> ni(block_neighbors); ni; ni++) {

      const GridGeometry::Neighbor& neighbor(*ni);

      if (!include_singularity_block_neighbors && neighbor.isSingularity()) {
         continue;
      }

      const BlockId neighbor_block_id(neighbor.getBlockId());

      blocki = d_single_block_trees.find(neighbor_block_id);

      if (blocki == d_single_block_trees.end()) {
         continue;
      }

      Box transformed_box(box);

      d_grid_geometry->transformBox(transformed_box,
         refinement_ratio,
         neighbor_block_id,
         block_id);

      blocki->second.findOverlapBoxes(overlap_mapped_boxes, transformed_box);

   }
}
#endif
/*
 *************************************************************************
 *************************************************************************
 */
void MultiblockBoxTree::findOverlapBoxes(
   Connector& overlap_connector,
   const Box& box,
   const IntVector& refinement_ratio,
   bool include_singularity_block_neighbors) const
{
   const BlockId& block_id = box.getBlockId();

   TBOX_DIM_ASSERT_CHECK_ARGS3(*d_grid_geometry, box, refinement_ratio);

   TBOX_ASSERT(block_id.getBlockValue() >= 0 &&
               block_id.getBlockValue() < d_grid_geometry->getNumberBlocks());

   /*
    * Search in the index space of block_id for overlaps.
    */

   std::map<BlockId, BoxTree>::const_iterator
   blocki(d_single_block_trees.find(block_id));

   if (blocki != d_single_block_trees.end()) {
      blocki->second.findOverlapBoxes(overlap_connector, box);
   }

   /*
    * Search in the index spaces neighboring block_id for overlaps.
    */

   const tbox::List<GridGeometry::Neighbor>& block_neighbors(
      d_grid_geometry->getNeighbors(block_id));

   for (tbox::ListIterator<GridGeometry::Neighbor> ni(block_neighbors); ni; ni++) {

      const GridGeometry::Neighbor& neighbor(*ni);

      if (!include_singularity_block_neighbors && neighbor.isSingularity()) {
         continue;
      }

      const BlockId neighbor_block_id(neighbor.getBlockId());

      blocki = d_single_block_trees.find(neighbor_block_id);

      if (blocki == d_single_block_trees.end()) {
         continue;
      }

      Box transformed_box(box);

      d_grid_geometry->transformBox(transformed_box,
         refinement_ratio,
         neighbor_block_id,
         block_id);

      blocki->second.findOverlapBoxes(overlap_connector, transformed_box);

   }
   return;
}

/*
 *************************************************************************
 *************************************************************************
 */

void MultiblockBoxTree::findOverlapBoxes(
   std::vector<Box>& overlap_mapped_boxes,
   const Box& box,
   const BlockId& block_id,
   const IntVector& refinement_ratio,
   bool include_singularity_block_neighbors) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(*d_grid_geometry, box, refinement_ratio);

   TBOX_ASSERT(block_id.getBlockValue() >= 0 &&
      block_id.getBlockValue() < d_grid_geometry->getNumberBlocks());

   /*
    * Search in the index space of block_id for overlaps.
    */

   std::map<BlockId, BoxTree>::const_iterator blocki(d_single_block_trees.find(block_id));

   if (blocki != d_single_block_trees.end()) {
      blocki->second.findOverlapBoxes(overlap_mapped_boxes, box);
   }

   /*
    * Search in the index spaces neighboring block_id for overlaps.
    */

   const tbox::List<GridGeometry::Neighbor>& block_neighbors(
      d_grid_geometry->getNeighbors(block_id));

   for (tbox::ListIterator<GridGeometry::Neighbor> ni(block_neighbors); ni; ni++) {

      const GridGeometry::Neighbor& neighbor(*ni);

      if (!include_singularity_block_neighbors && neighbor.isSingularity()) {
         continue;
      }

      const BlockId neighbor_block_id(neighbor.getBlockId());

      blocki = d_single_block_trees.find(neighbor_block_id);

      if (blocki == d_single_block_trees.end()) {
         continue;
      }

      Box transformed_box(box);

      d_grid_geometry->transformBox(transformed_box,
         refinement_ratio,
         neighbor_block_id,
         block_id);

      blocki->second.findOverlapBoxes(overlap_mapped_boxes, transformed_box);

   }
}

/*
 *************************************************************************
 *************************************************************************
 */

void MultiblockBoxTree::findOverlapBoxes(
   std::vector<const Box *>& overlap_mapped_boxes,
   const Box& box,
   const BlockId& block_id,
   const IntVector& refinement_ratio,
   bool include_singularity_block_neighbors) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(*d_grid_geometry, box, refinement_ratio);

   TBOX_ASSERT(block_id.getBlockValue() >= 0 &&
      block_id.getBlockValue() < d_grid_geometry->getNumberBlocks());

   /*
    * Search in the index space of block_id for overlaps.
    */

   std::map<BlockId, BoxTree>::const_iterator blocki(d_single_block_trees.find(block_id));

   if (blocki != d_single_block_trees.end()) {
      blocki->second.findOverlapBoxes(overlap_mapped_boxes, box);
   }

   /*
    * Search in the index spaces neighboring block_id for overlaps.
    */

   const tbox::List<GridGeometry::Neighbor>& block_neighbors(
      d_grid_geometry->getNeighbors(block_id));

   for (tbox::ListIterator<GridGeometry::Neighbor> ni(block_neighbors); ni; ni++) {

      const GridGeometry::Neighbor& neighbor(*ni);

      if (!include_singularity_block_neighbors && neighbor.isSingularity()) {
         continue;
      }

      const BlockId neighbor_block_id(neighbor.getBlockId());

      blocki = d_single_block_trees.find(neighbor_block_id);

      if (blocki == d_single_block_trees.end()) {
         continue;
      }

      Box transformed_box(box);

      d_grid_geometry->transformBox(transformed_box,
         refinement_ratio,
         neighbor_block_id,
         block_id);

      blocki->second.findOverlapBoxes(overlap_mapped_boxes, transformed_box);

   }
}

/*
 *************************************************************************
 *************************************************************************
 */

void MultiblockBoxTree::findOverlapBoxes(
   BoxList& overlap_boxes,
   const Box& box,
   const BlockId& block_id,
   const IntVector& refinement_ratio,
   bool include_singularity_block_neighbors) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(*d_grid_geometry, box, refinement_ratio);

   TBOX_ASSERT(block_id.getBlockValue() >= 0 &&
      block_id.getBlockValue() < d_grid_geometry->getNumberBlocks());

   /*
    * Search in the index space of block_id for overlaps.
    */

   std::map<BlockId, BoxTree>::const_iterator blocki(d_single_block_trees.find(block_id));

   if (blocki != d_single_block_trees.end()) {
      blocki->second.findOverlapBoxes(overlap_boxes, box);
   }

   /*
    * Search in the index spaces neighboring block_id for overlaps.
    */

   const tbox::List<GridGeometry::Neighbor>& block_neighbors(
      d_grid_geometry->getNeighbors(block_id));

   for (tbox::ListIterator<GridGeometry::Neighbor> ni(block_neighbors); ni; ni++) {

      const GridGeometry::Neighbor& neighbor(*ni);

      if (!include_singularity_block_neighbors && neighbor.isSingularity()) {
         continue;
      }

      const BlockId neighbor_block_id(neighbor.getBlockId());

      blocki = d_single_block_trees.find(neighbor_block_id);

      if (blocki == d_single_block_trees.end()) {
         continue;
      }

      Box transformed_box(box);

      d_grid_geometry->transformBox(transformed_box,
         refinement_ratio,
         neighbor_block_id,
         block_id);

      blocki->second.findOverlapBoxes(overlap_boxes, transformed_box);

   }
}

/*
 *************************************************************************
 *************************************************************************
 */
tbox::Pointer<MultiblockBoxTree> MultiblockBoxTree::createRefinedTree(
   const IntVector& ratio) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*d_grid_geometry, ratio);
   const tbox::Dimension& dim(d_grid_geometry->getDim());
   TBOX_ASSERT(ratio >= IntVector::getOne(dim));

   MultiblockBoxTree* rval = new MultiblockBoxTree();
   rval->d_grid_geometry = d_grid_geometry;

   for (std::map<BlockId, BoxTree>::const_iterator mi = d_single_block_trees.begin();
        mi != d_single_block_trees.end(); ++mi) {

      const BlockId& block_id(mi->first);
      const BoxTree& mapped_box_tree(mi->second);

      /*
       * The next lines really just do this:
       * rval->d_single_block_trees[block_id] = *mapped_box_tree.createRefinedTree(ratio);
       *
       * We cannot do it concisely without using the default
       * constructor for BoxTree.
       */
      const std::pair<std::map<BlockId, BoxTree>::iterator, bool> insert_return_value =
         rval->d_single_block_trees.insert(
            std::pair<BlockId, BoxTree>(
               block_id,
               BoxTree(dim)));
      TBOX_ASSERT(insert_return_value.second);
      insert_return_value.first->second = *mapped_box_tree.createRefinedTree(ratio);

   }

   return tbox::Pointer<MultiblockBoxTree>(rval);
}

/*
 *************************************************************************
 *************************************************************************
 */
void MultiblockBoxTree::getBoxes(
   std::vector<Box>& mapped_boxes) const
{
   for (std::map<BlockId, BoxTree>::const_iterator mi = d_single_block_trees.begin();
        mi != d_single_block_trees.end(); ++mi) {
      mi->second.getBoxes(mapped_boxes);
   }
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
