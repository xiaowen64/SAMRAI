/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Multiblock binary trees of Boxes for overlap searches.
 *
 ************************************************************************/

namespace SAMRAI {
namespace hier {

/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
MultiblockBoxTree::clear()
{
   d_single_block_trees.clear();
}

/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
int
MultiblockBoxTree::getNumberBlocksInTree() const
{
   return static_cast<int>(d_single_block_trees.size());
}

/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
const GridGeometry*
MultiblockBoxTree::getGridGeometry() const
{
   return d_grid_geometry;
}

/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool
MultiblockBoxTree::hasBoxInBlock(
   const BlockId& block_id) const
{
   return d_single_block_trees.find(block_id) !=
          d_single_block_trees.end();
}

}
}
