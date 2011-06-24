/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Multiblock binary trees of MappedBoxes for overlap searches.
 *
 ************************************************************************/

#ifndef included_hier_MultiblockMappedBoxTree_C
#define included_hier_MultiblockMappedBoxTree_C

#include "SAMRAI/hier/MultiblockMappedBoxTree.h"
#include "SAMRAI/hier/GridGeometry.h"
#include "SAMRAI/hier/BoxList.h"
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

MultiblockMappedBoxTree::MultiblockMappedBoxTree(
   const tbox::ConstPointer<hier::GridGeometry> &grid_geometry,
   const MappedBoxSet& mapped_boxes,
   size_t min_number ):
   d_grid_geometry(grid_geometry)
{
   const tbox::Dimension &dim(d_grid_geometry->getDim());

   /*
    * Group MappedBoxes by their BlockId and
    * create a tree for each BlockId.
    */

   std::map<BlockId,std::vector<MappedBox> > mapped_boxes_by_block;
   for ( MappedBoxSet::const_iterator bi=mapped_boxes.begin();
         bi!=mapped_boxes.end(); ++bi ) {
      mapped_boxes_by_block[bi->getBlockId()].push_back(*bi);
   }

   for ( std::map<BlockId,std::vector<MappedBox> >::iterator blocki=mapped_boxes_by_block.begin();
         blocki!=mapped_boxes_by_block.end(); ++blocki ) {

      const BlockId &bid(blocki->first);
      std::vector<MappedBox> &mapped_boxes_for_block(blocki->second);

      TBOX_ASSERT( bid.getBlockValue() >= 0 &&
                   bid.getBlockValue() < grid_geometry->getNumberBlocks() );

      /*
       * The following lines do this:
       * d_single_block_trees[bid].generateTree(mapped_boxes_for_block, min_number);
       *
       * We cannot do it the concise way because it
       * requires the default constructor for MappedBoxTree.
       */
      const std::pair<std::map<BlockId,MappedBoxTree>::iterator,bool> insert_return_value(
         d_single_block_trees.insert(
            std::pair<BlockId,MappedBoxTree>(
               bid,
               MappedBoxTree(dim)) ) );
      TBOX_ASSERT(insert_return_value.second);
      insert_return_value.first->second.generateTree(mapped_boxes_for_block, min_number);

   }

   return;
}

/*
*************************************************************************
*************************************************************************
*/

MultiblockMappedBoxTree::MultiblockMappedBoxTree(
   const tbox::ConstPointer<hier::GridGeometry> &grid_geometry,
   const std::vector<MappedBox>& mapped_boxes,
   size_t min_number ):
   d_grid_geometry(grid_geometry)
{
   generateTree(grid_geometry, mapped_boxes, min_number);
   return;
}

/*
*************************************************************************
*************************************************************************
*/

MultiblockMappedBoxTree::MultiblockMappedBoxTree()
   : d_grid_geometry(NULL)
{
   return;
}

/*
*************************************************************************
*************************************************************************
*/

MultiblockMappedBoxTree::~MultiblockMappedBoxTree()
{
}

/*
*************************************************************************
Generate the tree from a given mutable vector of mapped_boxes.
The vector will be changed and its output state is undefined.
*************************************************************************
*/
void MultiblockMappedBoxTree::generateTree(
   const tbox::ConstPointer<hier::GridGeometry> &grid_geometry,
   const std::vector<MappedBox>& mapped_boxes,
   size_t min_number)
{
   d_grid_geometry = grid_geometry;

   /*
    * Group MappedBoxes by their BlockId and create a tree for each
    * BlockId.
    */
   std::map<BlockId,std::vector<MappedBox> > mapped_boxes_by_block;
   for ( std::vector<MappedBox>::const_iterator bi=mapped_boxes.begin();
         bi!=mapped_boxes.end(); ++bi ) {
      mapped_boxes_by_block[bi->getBlockId()].push_back(*bi);
   }

   for ( std::map<BlockId,std::vector<MappedBox> >::iterator blocki=mapped_boxes_by_block.begin();
         blocki!=mapped_boxes_by_block.end(); ++blocki ) {

      const BlockId &bid(blocki->first);
      std::vector<MappedBox> &mapped_boxes_for_block(blocki->second);

      TBOX_ASSERT( bid.getBlockValue() >= 0 &&
                   bid.getBlockValue() < grid_geometry->getNumberBlocks() );

      /*
       * The following lines do this:
       * d_single_block_trees[bid].generateTree(mapped_boxes_for_block, min_number);
       *
       * We cannot do it the concise way because it requires the
       * default constructor for MappedBoxTree.
       */
      const std::pair<std::map<BlockId,MappedBoxTree>::iterator,bool> insert_return_value(
         d_single_block_trees.insert(
            std::pair<BlockId,MappedBoxTree>(
               bid,
               MappedBoxTree(grid_geometry->getDim())) ) );
      TBOX_ASSERT(insert_return_value.second);
      insert_return_value.first->second.generateTree(mapped_boxes_for_block, min_number);

   }

   return;
}

/*
*************************************************************************
*************************************************************************
*/

void MultiblockMappedBoxTree::clear()
{
   d_single_block_trees.clear();
   d_grid_geometry.setNull();
}

/*
*************************************************************************
*************************************************************************
*/

bool MultiblockMappedBoxTree::isInitialized() const
{
   return !d_grid_geometry.isNull();
}

/*
*************************************************************************
*************************************************************************
*/

const tbox::ConstPointer<hier::GridGeometry> &MultiblockMappedBoxTree::getGridGeometry() const
{
   return d_grid_geometry;
}

/*
*************************************************************************
*************************************************************************
*/

bool MultiblockMappedBoxTree::hasMappedBoxInBlock(
   const hier::BlockId &block_id ) const
{
   return (d_single_block_trees.find(block_id) !=
           d_single_block_trees.end());
}

/*
*************************************************************************
*************************************************************************
*/

const MappedBoxTree &MultiblockMappedBoxTree::getSingleBlockMappedBoxTree(
   const hier::BlockId &block_id ) const
{
   std::map<BlockId,MappedBoxTree>::const_iterator mi =
      d_single_block_trees.find(block_id);

   if ( mi == d_single_block_trees.end() ) {
      TBOX_ERROR("MultiblockMappedBoxTree::getSingleBlockMappedBoxTree: cannot\n"
                 <<"return the single-block MappedBoxTree for block " << block_id
                 <<"\neither because there is no such block in the domain configuration\n"
                 <<"or no MappedBoxes from that block exists in the tree.\n"
                 <<"You can use MultiblockMappedBoxTree::hasMappedBoxInBlock()\n"
                 <<"to determine whether the tree exists.");
   }

   return mi->second;
}

/*
*************************************************************************
*************************************************************************
*/

void MultiblockMappedBoxTree::findOverlapMappedBoxes(
   MappedBoxSet& overlap_mapped_boxes,
   const Box& box,
   const BlockId &block_id,
   const IntVector &refinement_ratio,
   bool include_singularity_block_neighbors ) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(*d_grid_geometry, box, refinement_ratio);

   TBOX_ASSERT( block_id.getBlockValue() >= 0 &&
                block_id.getBlockValue() < d_grid_geometry->getNumberBlocks() );

   /*
    * Search in the index space of block_id for overlaps.
    */

   std::map<BlockId,MappedBoxTree>::const_iterator
      blocki(d_single_block_trees.find(block_id));

   if ( blocki != d_single_block_trees.end() ) {
      blocki->second.findOverlapMappedBoxes(overlap_mapped_boxes, box);
   }


   /*
    * Search in the index spaces neighboring block_id for overlaps.
    */

   const tbox::List<GridGeometry::Neighbor>& block_neighbors(
      d_grid_geometry->getNeighbors(block_id.getBlockValue()));

   for ( tbox::ListIterator<GridGeometry::Neighbor> ni(block_neighbors); ni; ni++ ) {

      const GridGeometry::Neighbor &neighbor(*ni);

      if ( !include_singularity_block_neighbors && neighbor.isSingularity() ) {
         continue;
      }

      const BlockId neighbor_block_id(neighbor.getBlockNumber());

      blocki = d_single_block_trees.find(neighbor_block_id);

      if ( blocki == d_single_block_trees.end() ) {
         continue;
      }

      Box transformed_box(box);

      d_grid_geometry->translateBox(transformed_box,
                                    refinement_ratio,
                                    neighbor_block_id,
                                    block_id);

      blocki->second.findOverlapMappedBoxes(overlap_mapped_boxes, transformed_box);

   }

   return;
}

/*
*************************************************************************
*************************************************************************
*/

void MultiblockMappedBoxTree::findOverlapMappedBoxes(
   std::vector<MappedBox>& overlap_mapped_boxes,
   const Box& box,
   const BlockId &block_id,
   const IntVector &refinement_ratio,
   bool include_singularity_block_neighbors ) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(*d_grid_geometry, box, refinement_ratio);

   TBOX_ASSERT( block_id.getBlockValue() >= 0 &&
                block_id.getBlockValue() < d_grid_geometry->getNumberBlocks() );

   /*
    * Search in the index space of block_id for overlaps.
    */

   std::map<BlockId,MappedBoxTree>::const_iterator blocki(d_single_block_trees.find(block_id));

   if ( blocki != d_single_block_trees.end() ) {
      blocki->second.findOverlapMappedBoxes(overlap_mapped_boxes, box);
   }


   /*
    * Search in the index spaces neighboring block_id for overlaps.
    */

   const tbox::List<GridGeometry::Neighbor>& block_neighbors(
      d_grid_geometry->getNeighbors(block_id.getBlockValue()));

   for ( tbox::ListIterator<GridGeometry::Neighbor> ni(block_neighbors); ni; ni++ ) {

      const GridGeometry::Neighbor &neighbor(*ni);

      if ( !include_singularity_block_neighbors && neighbor.isSingularity() ) {
         continue;
      }

      const BlockId neighbor_block_id(neighbor.getBlockNumber());

      blocki = d_single_block_trees.find(neighbor_block_id);

      if ( blocki == d_single_block_trees.end() ) {
         continue;
      }

      Box transformed_box(box);

      d_grid_geometry->translateBox(transformed_box,
                                    refinement_ratio,
                                    neighbor_block_id,
                                    block_id);

      blocki->second.findOverlapMappedBoxes(overlap_mapped_boxes, transformed_box);

   }

   return;
}

/*
*************************************************************************
*************************************************************************
*/

void MultiblockMappedBoxTree::findOverlapMappedBoxes(
   hier::BoxList & overlap_mapped_boxes,
   const Box& box,
   const BlockId &block_id,
   const IntVector &refinement_ratio,
   bool include_singularity_block_neighbors ) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(*d_grid_geometry, box, refinement_ratio);

   TBOX_ASSERT( block_id.getBlockValue() >= 0 &&
                block_id.getBlockValue() < d_grid_geometry->getNumberBlocks() );

   /*
    * Search in the index space of block_id for overlaps.
    */

   std::map<BlockId,MappedBoxTree>::const_iterator blocki(d_single_block_trees.find(block_id));

   if ( blocki != d_single_block_trees.end() ) {
      blocki->second.findOverlapMappedBoxes(overlap_mapped_boxes, box);
   }


   /*
    * Search in the index spaces neighboring block_id for overlaps.
    */

   const tbox::List<GridGeometry::Neighbor>& block_neighbors(
      d_grid_geometry->getNeighbors(block_id.getBlockValue()));

   for ( tbox::ListIterator<GridGeometry::Neighbor> ni(block_neighbors); ni; ni++ ) {

      const GridGeometry::Neighbor &neighbor(*ni);

      if ( !include_singularity_block_neighbors && neighbor.isSingularity() ) {
         continue;
      }

      const BlockId neighbor_block_id(neighbor.getBlockNumber());

      blocki = d_single_block_trees.find(neighbor_block_id);

      if ( blocki == d_single_block_trees.end() ) {
         continue;
      }

      Box transformed_box(box);

      d_grid_geometry->translateBox(transformed_box,
                                    refinement_ratio,
                                    neighbor_block_id,
                                    block_id);

      blocki->second.findOverlapMappedBoxes(overlap_mapped_boxes, transformed_box);

   }

   return;
}



/*
*************************************************************************
*************************************************************************
*/
tbox::Pointer<MultiblockMappedBoxTree> MultiblockMappedBoxTree::createRefinedTree(
   const IntVector& ratio) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*d_grid_geometry, ratio);
   const tbox::Dimension &dim(d_grid_geometry->getDim());
   TBOX_ASSERT(ratio >= IntVector::getOne(dim));

   MultiblockMappedBoxTree* rval = new MultiblockMappedBoxTree();
   rval->d_grid_geometry = d_grid_geometry;

   for ( std::map<BlockId,MappedBoxTree>::const_iterator mi=d_single_block_trees.begin();
         mi!=d_single_block_trees.end(); ++mi ) {

      const BlockId &block_id(mi->first);
      const MappedBoxTree &mapped_box_tree(mi->second);

      /*
       * The next lines really just do this:
       * rval->d_single_block_trees[block_id] = *mapped_box_tree.createRefinedTree(ratio);
       *
       * We cannot do it concisely without using the default
       * constructor for MappedBoxTree.
       */
      const std::pair<std::map<BlockId,MappedBoxTree>::iterator,bool> insert_return_value =
         rval->d_single_block_trees.insert(
            std::pair<BlockId,MappedBoxTree>(
               block_id,
               MappedBoxTree(dim)) );
      TBOX_ASSERT(insert_return_value.second);
      insert_return_value.first->second = *mapped_box_tree.createRefinedTree(ratio);

   }

   return tbox::Pointer<MultiblockMappedBoxTree>(rval);
}



/*
*************************************************************************
*************************************************************************
*/
void MultiblockMappedBoxTree::getMappedBoxes(
   std::vector<MappedBox>& mapped_boxes) const
{
   for ( std::map<BlockId,MappedBoxTree>::const_iterator mi=d_single_block_trees.begin();
         mi!=d_single_block_trees.end(); ++mi ) {
      mi->second.getMappedBoxes(mapped_boxes);
   }
   return;
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
