/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Multiblock binary trees of MappedBoxes for overlap searches.
 *
 ************************************************************************/

#ifndef included_hier_MultiblockMappedBoxTree
#define included_hier_MultiblockMappedBoxTree

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/MappedBoxSet.h"
#include "SAMRAI/hier/MappedBoxTree.h"
#include "SAMRAI/tbox/ConstPointer.h"
#include "SAMRAI/tbox/DescribedClass.h"

#include <vector>
#include <map>

namespace SAMRAI {
namespace hier {

class GridGeometry;

/*!
 * @brief Utility sorting MappedBoxes into tree-like form for finding
 * box overlaps.  All boxes must be specified in the transformation
 * specified by their BlockId.
 *
 * Overlap searches are done by
 * - hasOverlap()
 * - findOverlapMappedBoxes()
 *
 * Significant changes since design review:
 * - findOverlapMappedBoxes requires refinement_ratio.
 */

class MultiblockMappedBoxTree : public tbox::DescribedClass
{

public:

   /*!
    * @brief Constructs a MultiblockMappedBoxTree from set of MappedBoxes.
    *
    * @param[in] grid_geometry GridGeometry desribing the multiblock
    * environment.
    *
    * @param[in] mapped_boxes.  No empty boxes are allowed.  An assertion
    *                           failure will occur if the mapped boxes in this
    *                           input set do not all have the same BlockId.
    *
    * @param[in] min_number Split up sets of boxes while the number of
    * boxes in a subset is greater than this value.  Setting to a
    * larger value tends to make tree building faster but tree
    * searching slower, and vice versa.  @b Default: 10
    */
   explicit MultiblockMappedBoxTree(
      const tbox::ConstPointer<GridGeometry> &grid_geometry,
      const MappedBoxSet& mapped_boxes,
      size_t min_number = 10);

   /*!
    * @brief Constructs a MultiblockMappedBoxTree from vector of MappedBoxes.
    *
    * See MultiblockMappedBoxTree( const tbox::Dimension& , const MappedBoxSet& , size_t min_number );
    *
    * @param[in] grid_geometry
    *
    * @param[in] mapped_boxes.  No empty boxes are allowed.
    *
    * @param[in] min_number.  @b Default: 10
    */
   explicit MultiblockMappedBoxTree(
      const tbox::ConstPointer<GridGeometry> &grid_geometry,
      const std::vector<Box>& mapped_boxes,
      size_t min_number = 10);

   /*!
    * @brief Constructs a MultiblockMappedBoxTree from vector of MappedBoxes.
    *
    * @param[in] grid_geometry
    *
    * @param[in] boxes.  No empty boxes are allowed.
    *
    * @param[in] min_number.  @b Default: 10
    */
   explicit MultiblockMappedBoxTree(
      const tbox::ConstPointer<GridGeometry> &grid_geometry,
      const std::map<BlockId, BoxList>& boxes,
      size_t min_number = 10);

   /*!
    * @brief Default constructor constructs an uninitialized
    * MultiblockMappedBoxTree.
    */
   explicit MultiblockMappedBoxTree();

   /*!
    * @brief Destructor.
    */
   ~MultiblockMappedBoxTree();

   /*!
    * @brief Generates the tree from a MUTABLE vector of MappedBoxes.
    *
    * For efficiency reasons, mapped_boxes is changed in the process.
    * Its output state is undefined.  However, you can change
    * mapped_boxes after tree generation without invalidating the
    * tree.
    *
    * @param[in] grid_geometry
    *
    * @param[in] mapped_boxes.  No empty boxes are allowed.
    *
    * @param[in] min_number
    */
   void
   generateTree(
      const tbox::ConstPointer<GridGeometry> &grid_geometry,
      const std::vector<Box>& mapped_boxes,
      size_t min_number = 10);

   /*!
    * @brief Generates the tree from MUTABLE lists of Boxes.
    *
    * For efficiency reasons, boxes is changed in the process.
    * Its output state is undefined.  However, you can change
    * boxes after tree generation without invalidating the
    * tree.
    *
    * @param[in] grid_geometry
    *
    * @param[in] mapped_boxes.  No empty boxes are allowed.
    *
    * @param[in] min_number
    */
   void
   generateTree(
      const tbox::ConstPointer<GridGeometry> &grid_geometry,
      const std::map<BlockId, BoxList>& boxes,
      size_t min_number = 10);

   /*!
    * @brief Return whether the tree contains any MappedBoxes with the
    * given BlockId.
    *
    * If the method getSingleBlockMappedBoxTree(const BlockId&) method
    * will throw an unrecoverable error if this method returns false
    * for the given BlockId.
    */
   bool hasMappedBoxInBlock(
      const BlockId &block_id ) const;

   /*!
    * @brief Return the tree for a single block.
    *
    * If the MappedBoxes initializing the tree did not contain at
    * least one Box with the given BlockId, the corresponding
    * single-block tree does not exist, and this method throws an
    * unrecoverable error.  To check for the existance of the tree,
    * use hasMappedBoxInBlock().
    */
   const MappedBoxTree &getSingleBlockMappedBoxTree(
      const BlockId &block_id ) const;

   /*!
    * @brief Reset to uninitialized state.
    *
    * Uninitialized trees can be initialized using generateTree().
    */
   void
   clear();

   /*!
    * @brief Check whether the tree has been initialized.
    *
    * Uninitialized trees can be initialized using generateTree().
    */
   bool
   isInitialized() const;

   /*!
    * @brief Return the GridGeometry object for the multiblock
    * environment.
    */
   const tbox::ConstPointer<GridGeometry> &getGridGeometry() const;


   //@{

   //! @name Overlap checks

   /*!
    * @brief Whether the given Box has an overlap with
    * MappedBoxes in the tree.
    *
    * We also check for overlap with MappedBoxes in blocks adjacent
    * to mapped_box's block.
    *
    * @param[in] box
    *
    * @param[in] block_id Specifies the block in which box is
    * specified.
    *
    * @param[in] include_singularity_block_neighbors Whether to include
    * intersections with boxes in blocks that are neighbors of block
    * block_id across a multiblock singularity.
    */
   bool
   hasOverlap(
      const Box& box,
      const BlockId &block_id,
      bool include_singularity_block_neighbors = false) const;

   /*!
    * @brief Find all boxes that overlap the given Box.
    *
    * To avoid unneeded work, the output @b overlap_mapped_boxes
    * container is not emptied.  Overlapping MappedBoxes are simply
    * added.
    *
    * Output is sorted.
    *
    * @param[out] overlap_mapped_boxes MappedBoxes that overlap with box.
    *
    * @param[in] box
    *
    * @param[in] block_id Specifies the block in which box is
    * specified.
    *
    * @param[in] refinement_ratio Refinement ratio of box's index
    * space.
    *
    * @param[in] include_singularity_block_neighbors Whether to include
    * intersections with boxes in blocks that are neighbors of block
    * block_id across a multiblock singularity.
    */
   void
   findOverlapMappedBoxes(
      MappedBoxSet& overlap_mapped_boxes,
      const Box& box,
      const BlockId &block_id,
      const IntVector &refinement_ratio,
      bool include_singularity_block_neighbors = false) const;

   /*!
    * @brief Find all boxes that overlap the given Box.
    *
    * To avoid unneeded work, the output @b overlap_mapped_boxes
    * container is not emptied.  Overlapping MappedBoxes are simply
    * added.
    *
    * Output is unsorted.
    *
    * @param[out] overlap_mapped_boxes
    *
    * @param[in] box
    *
    * @param[in] block_id Specifies the block in which box is
    * specified.
    *
    * @param[in] refinement_ratio Refinement ratio of box's index
    * space.
    *
    * @param[in] include_singularity_block_neighbors Whether to include
    * intersections with boxes in blocks that are neighbors of block
    * block_id across a multiblock singularity.
    */
   void
   findOverlapMappedBoxes(
      std::vector<Box>& overlap_mapped_boxes,
      const Box& box,
      const BlockId &block_id,
      const IntVector &refinement_ratio,
      bool include_singularity_block_neighbors = false) const;

   /*!
    * @brief Find all boxes that overlap the given Box.
    *
    * Analogous to findOverlapMappedBoxes returning a vector of MappedBoxes
    * but avoids the copies.  If the returned overlapped mapped boxes are used
    * in a context in which the MultiblockMappedBoxTree is constant there is
    * no point in incurring the cost of copying the tree's MappedBoxes.  Just
    * return a vector of their addresses.
    *
    * Output is unsorted.
    *
    * @param[out] overlap_mapped_boxes Pointers to MappedBoxes that overlap
    * with box.
    *
    * @param[in] box
    *
    * @param[in] block_id Specifies the block in which box is
    * specified.
    *
    * @param[in] refinement_ratio Refinement ratio of box's index
    * space.
    *
    * @param[in] include_singularity_block_neighbors Whether to include
    * intersections with boxes in blocks that are neighbors of block
    * block_id across a multiblock singularity.
    */
   void
   findOverlapMappedBoxes(
      std::vector<const Box*>& overlap_mapped_boxes,
      const Box& box,
      const BlockId &block_id,
      const IntVector &refinement_ratio,
      bool include_singularity_block_neighbors = false) const;

   /*!
    * @brief Find all boxes that overlap the given Box.
    *
    * To avoid unneeded work, the output @b overlap_boxes
    * container is not emptied.  Overlapping Boxes are simply
    * added.
    *
    * Output is unsorted.
    *
    * @param[in] box
    *
    * @param[in] block_id Specifies the block in which box is
    * specified.
    *
    * @param[in] refinement_ratio Refinement ratio of box's index
    * space.
    *
    * @param[in] include_singularity_block_neighbors Whether to include
    * intersections with boxes in blocks that are neighbors of block
    * block_id across a multiblock singularity.
    */
   void
   findOverlapBoxes(
      BoxList& overlap_boxes,
      const Box& box,
      const BlockId &block_id,
      const IntVector &refinement_ratio,
      bool include_singularity_block_neighbors = false) const;

   /*!
    * @brief Get the MappedBoxes in the tree.
    *
    * @param[out] mapped_boxes
    */
   void
   getMappedBoxes(
      std::vector<Box>& mapped_boxes) const;

   /*!
    * @brief Create a similar tree with the boxes refined by a given
    * ratio.
    *
    * @param[in] ratio The boxes are refined by this ratio.
    *
    * Note that there is no coresponding version to create a coarsened
    * tree.  Coarsened trees cannot be trivially generated like
    * refined trees can.  To create a coarsened tree, you must
    * manually get the boxes, coarsen them and use them to build a new
    * tree.
    */
   tbox::Pointer<MultiblockMappedBoxTree>
   createRefinedTree(
      const IntVector& ratio) const;

   //@}

private:

   /*!
    * @brief Container of single-block MappedBoxTrees.
    *
    * For each BlockId represented in the set of MappedBoxes, there is
    * an entry in this container.
    */
   std::map<BlockId,MappedBoxTree> d_single_block_trees;

   tbox::ConstPointer<GridGeometry> d_grid_geometry;

};

}
}

#endif
