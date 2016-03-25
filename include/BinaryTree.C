//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-1/source/hierarchy/boxes/BinaryTree.C $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1846 $
// Modified:	$LastChangedDate: 2008-01-11 09:51:05 -0800 (Fri, 11 Jan 2008) $
// Description: Utility class that provides standard binary tree functions.
//

#ifndef included_hier_BinaryTree_C
#define included_hier_BinaryTree_C

#include "BinaryTree.h"
#include "ProcessorMapping.h"
#include "tbox/SAMRAI_MPI.h"
#include "BoxComm.h"
#include "tbox/Utilities.h"

#ifdef DEBUG_NO_INLINE
#include "BinaryTree.I"
#endif

namespace SAMRAI {
   namespace hier {


/*
*********************************************************************
*
* Constructor and Destructor.  The constructor builds the
* binary tree, which is stored in d_tree.
*
*********************************************************************
*/

template<int DIM>  BinaryTree<DIM>::BinaryTree(const ProcessorMapping &mapping, 
                                   const BoxArray<DIM> &boxes)
: d_mapping(mapping)
{
   d_root = 0;

   int number_of_nodes = tbox::SAMRAI_MPI::getNodes();
   int number_of_boxes = boxes.getNumberOfBoxes();
   d_myid = tbox::SAMRAI_MPI::getRank();

   d_descendants.resizeArray(number_of_nodes);

   /*
    * 
    */
   d_first.resizeArray(number_of_nodes+1);
   d_boxes.resizeArray(number_of_boxes);

   /*
    * Construct d_boxes so that all boxes belonging to a processor
    * are contiguous: the boxes belonging to processor j will
    * be located in d_boxes[d_first[j]] through d_boxes[d_first[j+1]]-1.
    */

   //zero the array
   int k;
   for (k=0; k<number_of_nodes+1; ++k) {
     d_first[k] = 0;
   }

   //count the number of boxes local to each processor
   for (k=0; k<number_of_boxes; ++k) {
      int owner = mapping.getProcessorAssignment(k);
      d_first[owner+1] += 1;
   }

   //prefix-sum to find where the boxes for each processor start
   for (k=1; k<number_of_nodes+1; ++k) {
      d_first[k] += d_first[k-1];
   }

   //put the boxes where they belong
   for (k=0; k<number_of_boxes; ++k) {
      int owner = mapping.getProcessorAssignment(k);
      d_boxes[ d_first[owner]  ] = boxes(k);
      d_first[owner] += 1;
   }

   //finally, shift entries in d_first to the right
   for (k=number_of_nodes; k>0; --k) {
     d_first[k] = d_first[k-1]; 
   }
   d_first[k] = 0;

   /*
    * build the binary tree
    */

   //d_tree.resizeArray(number_of_nodes);
   d_tree.resizeArray(number_of_nodes);

   for (int i=0; i<number_of_nodes; ++i) {
      d_tree[i].id = i;
      d_tree[i].parent = -1;
      d_tree[i].lft_child = -1;
      d_tree[i].rgt_child = -1;
   }

   //possibly todo: perform a random initialization
   for (int id=0; id<number_of_nodes; ++id) {
      if (id != 0) {
        d_tree[id].parent = (id-1)/2;
      }

      int lft_child = id*2+1;
      if (lft_child < number_of_nodes) {
         d_tree[id].lft_child = lft_child;
      }

      int rgt_child = lft_child + 1;
      if (rgt_child < number_of_nodes) {
         d_tree[id].rgt_child = rgt_child;
      }
   }
}

template<int DIM>  BinaryTree<DIM>::~BinaryTree()
{
}


/*
**********************************************************************
* Returns, in id_out,  the processor id of the nearest anscestor in
* the tree that has a local patch this intersects with box.  If no
* such processor exists, returns the id of the root node.
* (This a highly specialized call, designed for use with
* Berger-Rigoutsous research.)
*
**********************************************************************
*/
template<int DIM> int BinaryTree<DIM>::findParticipatingAnscestor(const Box<DIM> &box)
{
   int node = d_myid;
   do {
      node = findParent(node);
   } while (! participates(node, box) ); 
  return node;
}


template<int DIM> void BinaryTree<DIM>::findParticipatingDescendants(const Box<DIM> &box,
                                                    tbox::Array<int> &id_out)
{
   d_descendants_len = 0;

   //perform recursive pre-order computation on left child, if s/he exists
   int left_child = findLeftChild(d_myid);
   if (left_child != -1) {
      preorderComputation(left_child, box);
   }

   //perform recursive pre-order computation on right child, if s/he exists
   int right_child = findRightChild(d_myid);
   if (right_child != -1) {
      preorderComputation(right_child, box);
   }

   //copy list of participating descendants from working to output space.
   id_out.resizeArray(d_descendants_len);
   for (int j=0; j<d_descendants_len; ++j) {
      id_out[j] = d_descendants[j];
   }
}

/*
**********************************************************************
*
*
**********************************************************************
*/


template<int DIM> void BinaryTree<DIM>::preorderComputation(int node_id, const Box<DIM> &box)
{
   if (participates(node_id, box)) {
      d_descendants[d_descendants_len++] = node_id;
      return;
   }

   //perform recursive pre-order computation on left child, if s/he exists
   int left_child = findLeftChild(node_id);
   if (left_child != -1) {
      preorderComputation(left_child, box);
   }

   //perform recursive pre-order computation on right child, if s/he exists
   int right_child = findRightChild(node_id);
   if (right_child != -1) {
      preorderComputation(right_child, box);
   }
}

/*
**********************************************************************
*
*
**********************************************************************
*/


template<int DIM> bool BinaryTree<DIM>::participates(int node_id, const Box<DIM> &box)
{
   if (node_id == findRoot()) {
     return true;
   }

   for (int i=d_first[node_id]; i<d_first[node_id+1]; ++i) {
      if (d_boxes[i].intersects(box)) {
         return true;
      }
   }
   return false;
}

/*
**********************************************************************
*
* Performs an all-to-one sum reduction on data;
* root ends up with the summation; a processor participates
* in the reduction only if it has a local patch that intersects
* with box.
*
**********************************************************************
*/

template<int DIM> void BinaryTree<DIM>::reduce(const Box<DIM> &box, int *data, int len)
{

#ifdef HAVE_MPI
   if (! participates(d_myid, box)) {
      return;
   }

   tbox::Array<int> descendants;
   findParticipatingDescendants(box, descendants);

   int rcv_ct = descendants.getSize();

   //if no descendants are sending us data, skip this
   //block.  Otherwise, receive the data from each sender,
   //then sum it with out own data.
   if (rcv_ct) {
      MPI_Request *request = new MPI_Request[rcv_ct];
      MPI_Status *status = new MPI_Status[rcv_ct];
      int *buf_in = new int[rcv_ct*len];
  
      //start non-blocking receive from each descendant
      for (int j=0; j<rcv_ct; ++j) {
         int source = descendants[j];
         MPI_Irecv(buf_in+j*len, len, MPI_INT, source, 
                  0, MPI_COMM_WORLD, request+j);
      }

      //wait for all recvs to finish
      MPI_Waitall(rcv_ct, request, status);
      delete [] request;
      delete [] status;

      //sum the data with our own
      for (int j=0; j<rcv_ct; ++j) {
         int offset = len*j;
         for (int k=0; k<len; ++k) {
            data[k] += buf_in[offset+k];
         }
      }
      delete [] buf_in;
   } //if (rcv_ct) 


   //send the result (i.e., the partially reduced data) to an ancestor
   if (d_myid != findRoot()) {
      int anscestor = findParticipatingAnscestor(box);
      MPI_Send(data, len, MPI_INT, anscestor, 0, MPI_COMM_WORLD);
   }
#endif
}


/*
**********************************************************************
*
* The box_to_bcast is sent from findRoot() to all processors
* for which participates(participants_box) == true;
* The first function is for broadcasting a box; the second
* is for broadcasting a single integer.
*
**********************************************************************
*/
template<int DIM> void BinaryTree<DIM>::partialBcast(const Box<DIM> &participants_box, 
                                    Box<DIM> &box_to_bcast)
{
   if (! participates(d_myid, participants_box)) {
      return;
   }
   
   /*
    * receive the box being broadcast
    */
   if (d_myid != findRoot()) {
     int anscestor = findParticipatingAnscestor(participants_box);
     BoxComm<DIM>::recvBox(box_to_bcast, anscestor);
   }

   /*
    * send the box being broadcast to descendants in the tree
    */
   tbox::Array<int> descendants;
   findParticipatingDescendants(participants_box, descendants);
   int rcv_ct = descendants.getSize();
  
   if (rcv_ct) {
      BoxComm<DIM>::sendBox(box_to_bcast, descendants);
   }
}

template<int DIM> void BinaryTree<DIM>::partialBcast(const Box<DIM> &participants_box, 
                                    int &value)
{
#ifdef HAVE_MPI
   if (! participates(d_myid, participants_box)) {
      return;
   }
   
   /*
    * receive the value being broadcast
    */
   if (d_myid != findRoot()) {
     int anscestor = findParticipatingAnscestor(participants_box);
     MPI_Status status;
     MPI_Recv(&value, 1, MPI_INT, anscestor, 0, MPI_COMM_WORLD, &status);
   }

   /*
    * send the box being broadcast to descendants in the tree
    */
   tbox::Array<int> descendants;
   findParticipatingDescendants(participants_box, descendants);
   int send_ct = descendants.getSize();
  
   if (send_ct) {
  
      //start non-blocking send to each descendant
      MPI_Request *request = new MPI_Request[send_ct];
      MPI_Status *status = new MPI_Status[send_ct];
      for (int j=0; j<send_ct; ++j) {
         int dest = descendants[j];
         MPI_Isend(&value, 1, MPI_INT, dest, 
                   0, MPI_COMM_WORLD, request+j);
      }
      //wait for all sends to finish
      MPI_Waitall(send_ct, request, status);
      delete [] request;
      delete [] status;
   }
#endif
}

/*
**********************************************************************
*
* Returns a new group and communicator that contains a subset
* of the processors that are members of the old group and
* communicator.  A processor from the old group will be included
* in the new group if the processor has a local patch that
* intersects with the box.
*
**********************************************************************
*/

template<int DIM> void BinaryTree<DIM>::buildParticipatingCommunicator(
   const Box<DIM> &box, 
   tbox::SAMRAI_MPI::comm old_comm,
   tbox::SAMRAI_MPI::group &new_group,
   tbox::SAMRAI_MPI::comm &new_comm)
{

#ifdef HAVE_MPI
   /*
    * Get a handle to the old group of processors.  
    * Also get the number of processors in the group.
    */
   int old_np;
   MPI_Group old_group, world_group;
   MPI_Comm_size(old_comm, &old_np);
   MPI_Comm_group(old_comm, &old_group);
   MPI_Comm_group(MPI_COMM_WORLD, &world_group);

   /*
    * translate (map) the ranks of the processors in the old_group 
    * to their respective ranks in the world group.
    */
   tbox::Array<int> world_ranks_array(old_np);
   tbox::Array<int> old_ranks_array(old_np);
   int *world_ranks = world_ranks_array.getPointer();
   int *old_ranks = old_ranks_array.getPointer();
   for (int k=0; k<old_np; ++k) {
      old_ranks[k] = k;
   }
   MPI_Group_translate_ranks(old_group, old_np, old_ranks,
                             world_group, world_ranks);

   /*
    * for each processor in the old group, determine if
    * if should be in the new group.
    */
   tbox::Array<bool> in_new_group(old_np);
   int new_np = 0;
   for (int i=0; i<old_np; ++i) {
      int world_id = world_ranks[i];

      if (participates(world_id, box)) {

         in_new_group[i] = true;
         ++new_np;
      } else {
         in_new_group[i] = false;
      } 

      int root = findRoot();
      int new_root = -1;
      for (int m=0; m<old_np; ++m) {
         world_id = world_ranks[m];
         if (world_id == root) {
            new_root = m;
         }
      }

      /*
       * ensure the root processor is in the new group!
       */
      if (! in_new_group[new_root]) {
         in_new_group[new_root] = true;
         ++new_np;
      }
   }

   /*
    * form the new group and communicator
    */
   tbox::Array<int> new_ranks_array(new_np);
   int *new_ranks = new_ranks_array.getPointer();
   int idx = 0;
   for (int j=0; j<old_np; ++j) {
      if (in_new_group[j]) {
         new_ranks[idx++] = j;
      }
   }

   MPI_Group_incl(old_group, new_np, new_ranks, &new_group);
   MPI_Comm_create(old_comm, new_group, &new_comm);
#endif
}

}
}

#endif
