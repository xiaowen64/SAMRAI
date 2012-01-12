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
   const GridGeometry& grid_geometry,
   const BoxContainer& boxes,
   size_t min_number):
   d_grid_geometry(&grid_geometry)
{
   NULL_USE(min_number);

   const tbox::Dimension& dim(d_grid_geometry->getDim());

   /*
    * Group Boxes by their BlockId and
    * create a tree for each BlockId.
    */

   for (BoxContainer::ConstIterator bi = boxes.begin();
        bi != boxes.end(); ++bi) {
      if (boxes.isOrdered()) {
         TBOX_ASSERT((*bi).getId().isValid());
         const BlockId& block_id = (*bi).getBlockId();
         d_single_block_trees[block_id].order();
         d_single_block_trees[block_id].insert(
            d_single_block_trees[block_id].end(), *bi);
      } else {
         const BlockId& block_id = (*bi).getBlockId();
         d_single_block_trees[block_id].pushBack(*bi); 
      }
   }

   for (std::map<BlockId, BoxContainer>::iterator blocki = d_single_block_trees.begin();
        blocki != d_single_block_trees.end(); ++blocki) {

      TBOX_ASSERT(blocki->first.getBlockValue() >= 0 &&
         blocki->first.getBlockValue() < d_grid_geometry->getNumberBlocks());

      blocki->second.makeTree();
   }
#if 0
      blocki->second.makeTree(); 

      /*
       * The following lines do this:
       * d_single_block_trees[blocki->first].generateTree(blocki->second, min_number);
       *
       * We cannot do it the concise way because it
       * requires the default constructor for BoxTree.
       */
      const std::pair<std::map<BlockId, BoxContainer>::iterator, bool> insert_return_value(
         d_single_block_trees.insert(
            std::pair<BlockId, BoxContainer>(
               blocki->first,
               blocki->second)));
      TBOX_ASSERT(insert_return_value.second);
   }
#endif
}

/*
 *************************************************************************
 *************************************************************************
 */

MultiblockBoxTree::MultiblockBoxTree(
   const GridGeometry& grid_geometry,
   const std::map<BlockId, BoxContainer>& boxes,
   size_t min_number):
   d_grid_geometry(&grid_geometry)
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
 * Generate the tree from a given container of boxes.
 *************************************************************************
 */
void MultiblockBoxTree::generateTree(
   const GridGeometry& grid_geometry,
   const BoxContainer& boxes,
   size_t min_number)
{
   NULL_USE(min_number);

   clear();

   d_grid_geometry = &grid_geometry;

   /*
    * Group Boxes by their BlockId and create a tree for each
    * BlockId.
    */
   for (BoxContainer::ConstIterator bi = boxes.begin();
        bi != boxes.end(); ++bi) {
      d_single_block_trees[bi->getBlockId()].insert(
         d_single_block_trees[bi->getBlockId()].end(),
         *bi);
   }

   for (std::map<BlockId, BoxContainer>::iterator blocki = d_single_block_trees.begin();
        blocki != d_single_block_trees.end(); ++blocki) {

      TBOX_ASSERT(blocki->first.getBlockValue() >= 0 &&
         blocki->first.getBlockValue() < d_grid_geometry->getNumberBlocks());

      blocki->second.makeTree();
   }
#if 0
      /*
       * The following lines do this:
       * d_single_block_trees[blocki->first].generateTree(blocki->second, min_number);
       *
       * We cannot do it the concise way because it requires the
       * default constructor for BoxTree.
       */
      const std::pair<std::map<BlockId, BoxContainer>::iterator, bool> insert_return_value(
         d_single_block_trees.insert(
            std::pair<BlockId, BoxContainer>(
               blocki->first,
               blocki->second)));
      TBOX_ASSERT(insert_return_value.second);
   }
#endif
}

/*
 *************************************************************************
 * Generate the tree from a given map of BoxContainers.
 *************************************************************************
 */
void MultiblockBoxTree::generateTree(
   const GridGeometry& grid_geometry,
   const std::map<BlockId, BoxContainer>& boxes,
   size_t min_number)
{
   clear();

   d_single_block_trees = boxes;
   d_grid_geometry = &grid_geometry;

   for (std::map<BlockId, BoxContainer>::iterator blocki = d_single_block_trees.begin();
        blocki != d_single_block_trees.end(); ++blocki) {

      TBOX_ASSERT(blocki->first.getBlockValue() >= 0 &&
         blocki->first.getBlockValue() < d_grid_geometry->getNumberBlocks());

      blocki->second.makeTree();
   }
#if 0

      /*
       * The following lines do this:
       * d_single_block_trees[blocki->first].generateTree(blocki->second, min_number);
       *
       * We cannot do it the concise way because it requires the
       * default constructor for BoxTree.
       */
      const std::pair<std::map<BlockId, BoxContainer>::iterator, bool> insert_return_value(
         d_single_block_trees.insert(
            std::pair<BlockId, BoxContainer>(
               blocki->first,
               blocki->second)));
      TBOX_ASSERT(insert_return_value.second);
   }
#endif
}

/*
 *************************************************************************
 * Generate the tree from a given vector of boxes.
 *************************************************************************
 */
void MultiblockBoxTree::generateNonPeriodicTree(
   const GridGeometry& grid_geometry,
   const BoxContainer& boxes,
   size_t min_number)
{
   NULL_USE(min_number);

   clear();

   d_grid_geometry = &grid_geometry;

   /*
    * Group Boxes by their BlockId and create a tree for each
    * BlockId.
    */
   for (BoxContainer::ConstIterator bi = boxes.begin();
        bi != boxes.end(); ++bi) {
      if (!bi->isPeriodicImage()) {
         d_single_block_trees[bi->getBlockId()].insert(
            d_single_block_trees[bi->getBlockId()].end(),
            *bi);
      }
   }

   for (std::map<BlockId, BoxContainer>::iterator blocki = d_single_block_trees.begin();
        blocki != d_single_block_trees.end(); ++blocki) {

      TBOX_ASSERT(blocki->first.getBlockValue() >= 0 &&
         blocki->first.getBlockValue() < d_grid_geometry->getNumberBlocks());

      blocki->second.makeTree();
   }
#if 0
      /*
       * The following lines do this:
       * d_single_block_trees[blocki->first].generateTree(blocki->second, min_number);
       *
       * We cannot do it the concise way because it requires the
       * default constructor for BoxTree.
       */
      const std::pair<std::map<BlockId, BoxContainer>::iterator, bool> insert_return_value(
         d_single_block_trees.insert(
            std::pair<BlockId, BoxContainer>(
               blocki->first,
               blocki->second)));
      TBOX_ASSERT(insert_return_value.second);
   }
#endif
}

/*
 *************************************************************************
 *************************************************************************
 */

void MultiblockBoxTree::clear()
{
   d_single_block_trees.clear();
   d_grid_geometry = NULL;
}

/*
 *************************************************************************
 *************************************************************************
 */

bool MultiblockBoxTree::isInitialized() const
{
   return d_grid_geometry != NULL;
}

/*
 *************************************************************************
 *************************************************************************
 */

const GridGeometry& MultiblockBoxTree::getGridGeometry() const
{
   return *d_grid_geometry;
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

const BoxContainer& MultiblockBoxTree::getSingleBlockBoxTree(
   const BlockId& block_id) const
{
   std::map<BlockId, BoxContainer>::const_iterator mi =
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

   std::map<BlockId, BoxContainer>::const_iterator
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
   std::vector<Box>& overlap_boxes,
   const Box& box,
   const IntVector& refinement_ratio,
   bool include_singularity_block_neighbors) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(*d_grid_geometry, box, refinement_ratio);

   const BlockId &block_id = box.getBlockId();
   TBOX_ASSERT(block_id.getBlockValue() >= 0 &&
      block_id.getBlockValue() < d_grid_geometry->getNumberBlocks());

   /*
    * Search in the index space of block_id for overlaps.
    */

   std::map<BlockId, BoxContainer>::const_iterator blocki(d_single_block_trees.find(block_id));

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

void MultiblockBoxTree::findOverlapBoxes(
   std::vector<const Box *>& overlap_boxes,
   const Box& box,
   const IntVector& refinement_ratio,
   bool include_singularity_block_neighbors) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(*d_grid_geometry, box, refinement_ratio);

   const BlockId &block_id = box.getBlockId();
   TBOX_ASSERT(block_id.getBlockValue() >= 0 &&
      block_id.getBlockValue() < d_grid_geometry->getNumberBlocks());

   /*
    * Search in the index space of block_id for overlaps.
    */

   std::map<BlockId, BoxContainer>::const_iterator blocki(d_single_block_trees.find(block_id));

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

void MultiblockBoxTree::findOverlapBoxes(
   BoxContainer& overlap_boxes,
   const Box& box,
   const IntVector& refinement_ratio,
   bool include_singularity_block_neighbors) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(*d_grid_geometry, box, refinement_ratio);

   const BlockId &block_id = box.getBlockId();
   TBOX_ASSERT(block_id.getBlockValue() >= 0 &&
      block_id.getBlockValue() < d_grid_geometry->getNumberBlocks());

   /*
    * Search in the index space of block_id for overlaps.
    */

   std::map<BlockId, BoxContainer>::const_iterator blocki(d_single_block_trees.find(block_id));

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

   for (std::map<BlockId, BoxContainer>::const_iterator mi = d_single_block_trees.begin();
        mi != d_single_block_trees.end(); ++mi) {

      const BlockId& block_id(mi->first);
//      const BoxContainer box_tree(mi->second);

      /*
       * The next lines really just do this:
       * rval->d_single_block_trees[block_id] = *box_tree.createRefinedTree(ratio);
       *
       * We cannot do it concisely without using the default
       * constructor for BoxTree.
       */
      rval->d_single_block_trees[block_id] = mi->second;
      rval->d_single_block_trees[block_id].refine(ratio);
      rval->d_single_block_trees[block_id].makeTree();
/*
      const std::pair<std::map<BlockId, BoxTree>::iterator, bool> insert_return_value =
         rval->d_single_block_trees.insert(
            std::pair<BlockId, BoxTree>(
               block_id,
               BoxTree(*box_tree.createRefinedTree(ratio))));
      TBOX_ASSERT(insert_return_value.second);
*/
   }

   return tbox::Pointer<MultiblockBoxTree>(rval);
}

/*
 *************************************************************************
 *************************************************************************
 */
void MultiblockBoxTree::getBoxes(
   std::vector<Box>& boxes) const
{
   for (std::map<BlockId, BoxContainer>::const_iterator mi = d_single_block_trees.begin();
        mi != d_single_block_trees.end(); ++mi) {
      //mi->second.getBoxes(boxes);
      for (BoxContainer::ConstIterator itr = mi->second.begin();
           itr != mi->second.end(); ++itr) {
         boxes.push_back(*itr);
      }
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
