/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Multiblock binary trees of Boxes for overlap searches.
 *
 ************************************************************************/

#ifndef included_hier_MultiblockBoxTree
#define included_hier_MultiblockBoxTree

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/BoxTree.h"
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/tbox/DescribedClass.h"

#include <vector>
#include <map>

namespace SAMRAI {
namespace hier {

class GridGeometry;

/*!
 * @brief Utility sorting Boxes into tree-like form for finding
 * box overlaps.  All boxes must be specified in the transformation
 * specified by their BlockId.
 *
 * Overlap searches are done by
 * - hasOverlap()
 * - findOverlapBoxes()
 *
 * Significant changes since design review:
 * - findOverlapBoxes requires refinement_ratio.
 */

class MultiblockBoxTree:public tbox::DescribedClass
{

public:
   /*!
    * @brief Constructs a MultiblockBoxTree from set of Boxes.
    *
    * @param[in] grid_geometry GridGeometry desribing the multiblock
    * environment.
    *
    * @param[in] boxes.  No empty boxes are allowed.
    *
    * @param[in] min_number Split up sets of boxes while the number of
    * boxes in a subset is greater than this value.  Setting to a
    * larger value tends to make tree building faster but tree
    * searching slower, and vice versa.  @b Default: 10
    */
   MultiblockBoxTree(
      const GridGeometry& grid_geometry,
      const BoxContainer& boxes,
      size_t min_number = 10);

   /*!
    * @brief Constructs a MultiblockBoxTree from vector of Boxes.
    *
    * See MultiblockBoxTree( const tbox::Dimension& , const BoxContainer& , size_t min_number );
    *
    * @param[in] grid_geometry
    *
    * @param[in] boxes.  No empty boxes are allowed.
    *
    * @param[in] min_number.  @b Default: 10
    */
   MultiblockBoxTree(
      const GridGeometry& grid_geometry,
      const std::vector<Box>& boxes,
      size_t min_number = 10);

   /*!
    * @brief Constructs a MultiblockBoxTree from a collection of BoxContainers each
    * of which is associated with a specific BlockId.
    *
    * @param[in] grid_geometry
    *
    * @param[in] boxes.  No empty boxes are allowed.
    *
    * @param[in] min_number.  @b Default: 10
    */
   MultiblockBoxTree(
      const GridGeometry& grid_geometry,
      const std::map<BlockId, BoxContainer>& boxes,
      size_t min_number = 10);

   /*!
    * @brief Default constructor constructs an uninitialized
    * MultiblockBoxTree.
    */
   MultiblockBoxTree();

   /*!
    * @brief Destructor.
    */
   ~MultiblockBoxTree();

   /*!
    * @brief Generates the tree from a BoxContainer.
    *
    * @param[in] grid_geometry
    *
    * @param[in] boxes.  No empty boxes are allowed.
    *
    * @param[in] min_number
    */
   void
   generateTree(
      const GridGeometry& grid_geometry,
      const BoxContainer& boxes,
      size_t min_number = 10);

   /*!
    * @brief Generates the tree from lists of Boxes.
    *
    * @param[in] grid_geometry
    *
    * @param[in] boxes.  No empty boxes are allowed.
    *
    * @param[in] min_number
    */
   void
   generateTree(
      const GridGeometry& grid_geometry,
      const std::map<BlockId, BoxContainer>& boxes,
      size_t min_number = 10);

   /*!
    * @brief Generates the tree of non-periodic Boxes from a BoxContainer.
    *
    * @param[in] grid_geometry
    *
    * @param[in] boxes.  No empty boxes are allowed.
    *
    * @param[in] min_number
    */
   void
   generateNonPeriodicTree(
      const GridGeometry& grid_geometry,
      const BoxContainer& boxes,
      size_t min_number = 10);

   /*!
    * @brief Return whether the tree contains any Boxes with the
    * given BlockId.
    *
    * If the method getSingleBlockBoxTree(const BlockId&) method
    * will throw an unrecoverable error if this method returns false
    * for the given BlockId.
    */
   bool
   hasBoxInBlock(
      const BlockId& block_id) const;

   /*!
    * @brief Return the tree for a single block.
    *
    * If the Boxes initializing the tree did not contain at
    * least one Box with the given BlockId, the corresponding
    * single-block tree does not exist, and this method throws an
    * unrecoverable error.  To check for the existance of the tree,
    * use hasBoxInBlock().
    */
   const BoxContainer&
   getSingleBlockBoxTree(
      const BlockId& block_id) const;

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
    *
    * Do not deallocate the returned GridGeometry.
    */
   const GridGeometry&
   getGridGeometry() const;

   //@{

   //! @name Overlap checks

   /*!
    * @brief Whether the given Box has an overlap with
    * Boxes in the tree.
    *
    * We also check for overlap with Boxes in blocks adjacent
    * to box's block.
    *
    * @param[in] box
    *
    * @param[in] include_singularity_block_neighbors Whether to include
    * intersections with boxes in blocks that are neighbors of box's
    * block across a multiblock singularity.
    */
   bool
   hasOverlap(
      const Box& box,
      bool include_singularity_block_neighbors = false) const;

   /*!
    * @brief Find all boxes that overlap the given Box and insert as neighbors
    * of that box in overlap_connector.
    *
    * @param[out] overlap_connector Overlap connector with box in its base
    * BoxLevel.
    *
    * @param[in] box
    *
    * @param[in] refinement_ratio Refinement ratio of box's index
    * space.
    *
    * @param[in] include_singularity_block_neighbors Whether to include
    * intersections with boxes in blocks that are neighbors of box's
    * block across a multiblock singularity.
    */
   void
   findOverlapBoxes(
      Connector& overlap_connector,
      const Box& box,
      const IntVector& refinement_ratio,
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
    * @param[out] overlap_boxes
    *
    * @param[in] box
    *
    * @param[in] refinement_ratio Refinement ratio of box's index
    * space.
    *
    * @param[in] include_singularity_block_neighbors Whether to include
    * intersections with boxes in blocks that are neighbors of box's
    * block across a multiblock singularity.
    */
   void
   findOverlapBoxes(
      std::vector<Box>& overlap_boxes,
      const Box& box,
      const IntVector& refinement_ratio,
      bool include_singularity_block_neighbors = false) const;

   /*!
    * @brief Find all boxes that overlap the given Box.
    *
    * Analogous to findOverlapBoxes returning a vector of Boxes
    * but avoids the copies.  If the returned overlapped boxes are used
    * in a context in which the MultiblockBoxTree is constant there is
    * no point in incurring the cost of copying the tree's Boxes.  Just
    * return a vector of their addresses.
    *
    * Output is unsorted.
    *
    * @param[out] overlap_boxes Pointers to Boxes that overlap
    * with box.
    *
    * @param[in] box
    *
    * @param[in] refinement_ratio Refinement ratio of box's index
    * space.
    *
    * @param[in] include_singularity_block_neighbors Whether to include
    * intersections with boxes in blocks that are neighbors of box's
    * block across a multiblock singularity.
    */
   void
   findOverlapBoxes(
      std::vector<const Box *>& overlap_boxes,
      const Box& box,
      const IntVector& refinement_ratio,
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
    * @param[in] refinement_ratio Refinement ratio of box's index
    * space.
    *
    * @param[in] include_singularity_block_neighbors Whether to include
    * intersections with boxes in blocks that are neighbors of box's
    * block across a multiblock singularity.
    */
   void
   findOverlapBoxes(
      BoxContainer& overlap_boxes,
      const Box& box,
      const IntVector& refinement_ratio,
      bool include_singularity_block_neighbors = false) const;

   /*!
    * @brief Get the Boxes in the tree.
    *
    * @param[out] boxes
    */
   void
   getBoxes(
      std::vector<Box>& boxes) const;

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
   tbox::Pointer<MultiblockBoxTree>
   createRefinedTree(
      const IntVector& ratio) const;

   //@}

private:
   /*!
    * @brief Container of single-block BoxTrees.
    *
    * For each BlockId represented in the set of Boxes, there is
    * an entry in this container.
    */
   std::map<BlockId, BoxContainer> d_single_block_trees;

   /*
    * @brief GridGeometry object.
    *
    * Do not delete this object.
    */
   const GridGeometry *d_grid_geometry;

};

}
}

#endif
