//
// File:	BinaryTree.I
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 374 $
// Modified:	$Date: 2005-05-20 11:58:22 -0700 (Fri, 20 May 2005) $
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
