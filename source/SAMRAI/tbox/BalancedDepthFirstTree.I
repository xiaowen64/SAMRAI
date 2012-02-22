/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Utility for building efficient communication tree.
 *
 ************************************************************************/

namespace SAMRAI {
namespace tbox {

SAMRAI_INLINE_KEYWORD
unsigned int BalancedDepthFirstTree::getRank() const
{
   return d_rank;
}

SAMRAI_INLINE_KEYWORD
unsigned int BalancedDepthFirstTree::getParentRank() const
{
   return d_parent;
}

SAMRAI_INLINE_KEYWORD
unsigned int BalancedDepthFirstTree::getChildRank(
   unsigned int child_number) const
{
   return d_children[child_number];
}


SAMRAI_INLINE_KEYWORD
unsigned short int BalancedDepthFirstTree::getNumberOfChildren() const
{
   return d_num_children;
}

SAMRAI_INLINE_KEYWORD
unsigned int BalancedDepthFirstTree::getInvalidRank() const
{
   return 1 << (8 * sizeof(unsigned int) - 2);
}

}
}
