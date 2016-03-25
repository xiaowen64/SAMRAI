//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/boxes/BinaryTree.h $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Utility class that provides standard binary tree functions.
//


#ifndef included_hier_BinaryTree
#define included_hier_BinaryTree

#include "SAMRAI_config.h"

#include <iostream>

#include "Box.h"
#include "BoxArray.h"
#include "ProcessorMapping.h"
#include "tbox/Array.h"

#include "tbox/SAMRAI_MPI.h"

namespace SAMRAI {
   namespace hier {

/**
 * Class BinaryTree<DIM> is a utility class designed to support
 * more efficient implementations of the Berger-Rigoutsous clustering
 * algorithm.  
 */

template<int DIM> class BinaryTree
{
public:
   /**
    * Constructs a binary tree.  If "number_of_nodes = 0,"
    * the tree will be constructed with the same number of
    * nodes as there are MPI processes.
    */
   BinaryTree(const ProcessorMapping &mapping, 
                    const BoxArray<DIM> &boxes);

   /**
    * The destructor does nothing of public interest.
    */
   ~BinaryTree<DIM>();

   /**
    * Returns true if "box" intersects a patch that is local 
    * to the processor "node_id;" else returns false.
    * However, always returns true when node_id == findRoot().
    * (This a highly specialized call, designed for use with
    * Berger-Rigoutsous research.)
    * (You must call init() before calling this method.)
    */
   bool participates(int node_id, const Box<DIM> &box);

   /**
    * Returns the processor id of the nearest anscestor in 
    * the tree that has a local patch that intersects with box.  
    * If no such processor exists, returns the id of the root node.
    * (This a highly specialized call, designed for use with
    * Berger-Rigoutsous research.)
    * (You must call init() before calling this method.)
    */
   int findParticipatingAnscestor(const Box<DIM> &box);

   /**
    * Returns, in id_out, the processor ids of this processor's
    * descendants that will send to this processor (i.e., 
    * if "i" is listed in id_out, then when processor "i"
    * calls findParticipatingAnscestor, id_out = this processor).
    * Note that id_out may be empty on return.
    * (This a highly specialized call, designed for use with
    * Berger-Rigoutsous research.)
    * (You must call init() before calling this method.)
    */
   void findParticipatingDescendants(const Box<DIM> &box,
                                     tbox::Array<int> &id_out);

   /**
    * Performs an all-to-one sum reduction on data;
    * root ends up with the summation; a processor participates
    * in the reduction only if it has a local patch that intersects 
    * with box.
    * (This a highly specialized call, designed for use with
    * Berger-Rigoutsous research.)
    * (You must call init() before calling this method.)
    */
   void reduce(const Box<DIM> &box, int *data, int len);

   /**
    * Broadcasts "box_to_bcast" to a subset of processors.
    * The subset is composed of those processors for which
    * participatesInReduction(processor, participants_box) is true.
    * (This a highly specialized call, designed for use with
    * Berger-Rigoutsous research.)
    * (You must call init() before calling this method.)
    */
   void partialBcast(const Box<DIM> &participants_box, 
                     Box<DIM> &box_to_bcast);
   void partialBcast(const Box<DIM> &participants_box, 
                     int &value);

   /*
    * Returns a new group and communicator that contains a subset
    * of the processors that are members of the old group and 
    * communicator.  A processor from the old group will be included
    * in the new group if the processor has a local patch that
    * intersects with the box.
    */
   void buildParticipatingCommunicator(const Box<DIM> &box, 
                                       tbox::SAMRAI_MPI::comm old_comm, 
                                       tbox::SAMRAI_MPI::group &new_group,
                                       tbox::SAMRAI_MPI::comm &new_comm);


private:

   /*
    * Returns the parent of "node" (by definition, the root of the
    * tree is its own parent).
    */
   int findParent(int node);

   /*
    * Returns the left child of "node;" if node has no left child,
    * returns -1.
    */
   int findLeftChild(int node);

   /*
    * Returns the right child of "node;" if node has no right child,
    * returns -1.
    */
   int findRightChild(int node);

   /*
    * Returns the root of the tree.
    */
   int findRoot();

   struct binTreeNode {
      int id;
      int parent;
      int lft_child;
      int rgt_child;
   };

   /*
    * Recursive call; if node_id participates in the histogram
    * reduction operation, add it to the d_descendants and
    * return; else, recurse on left and right children (if
    * they exist)
    */
   void preorderComputation(int node_id, const Box<DIM> &box);


  /*
   * The copy constructor and assignment operator are not implemented.
   */
   BinaryTree(const BinaryTree<DIM>&);
   BinaryTree<DIM>& operator=(const BinaryTree<DIM>&);

   tbox::Array<binTreeNode> d_tree;
   int d_root;

   ProcessorMapping d_mapping;
   tbox::Array< Box<DIM> > d_boxes;
   tbox::Array<int> d_first;

   /*
    * work space used in findParticipatingDescendants()
    */
   tbox::Array<int> d_descendants;
   int d_descendants_len;

   int d_myid;
};

}
}

#ifndef DEBUG_NO_INLINE
#include "BinaryTree.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "BinaryTree.C"
#endif
