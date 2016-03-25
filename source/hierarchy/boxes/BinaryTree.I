//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/boxes/BinaryTree.I $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	Utility class that provides standard binary tree functions.
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
   namespace hier {


template<int DIM> inline
int BinaryTree<DIM>::findParent(int node)
{
   return d_tree[node].parent;
}

template<int DIM> inline
int BinaryTree<DIM>::findLeftChild(int node)
{
   return d_tree[node].lft_child;
}

template<int DIM> inline
int BinaryTree<DIM>::findRightChild(int node)
{
   return d_tree[node].rgt_child;
}

template<int DIM> inline
int BinaryTree<DIM>::findRoot()
{
   return d_root;
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}
