//
// File:        BoxTreeNode.C
// Package:     SAMRAI hierarchy
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 179 $
// Modified:    $Date: 2005-01-20 14:50:51 -0800 (Thu, 20 Jan 2005) $
// Description: Basic building block used by BoxTree class.
//

#ifndef included_hier_BoxTreeNode_C
#define included_hier_BoxTreeNode_C

#include "BoxTreeNode.h"
#include "BoxGraphUtilities.h"
#include "tbox/MPI.h"


namespace SAMRAI {
   namespace hier {

/*
 * ************************************************************************
 *                                                                       *
 * Constructs a BoxTreeNodeX that represents the physical
 * domain specified by box.
 *                                                                       *
 * ************************************************************************
 */


template<int DIM>  BoxTreeNode<DIM>::BoxTreeNode(
   tbox::Array<typename BoxTreeNode<DIM>::Triple> box_triples,
   bool have_mapping,
   int min_length,
   int dim,
   int recurse_level)
{
   d_dim = dim;
   d_rank = tbox::MPI::getRank();

   /*
    * Compute this node's domain, which is the bounding box
    * for the list of boxes.
    */
   int len = box_triples.size();
   for (int j=0; j<len; ++j) {
      d_domain += box_triples[j].box;
   }

   /*
    * If the list of boxes is small enough, we won't
    * do any recursive stuff: we'll just let the boxes
    * live here.  In this case, there is no left child,
    * no right child, and no recursive d_tree.
    */
   min_length = (min_length < 1) ? 1 : min_length;
   if (len <= min_length || dim == 0) {
      d_box_triples = box_triples;
      d_work.resizeArray( len );
      return;
   }

   /*
    * Partition the boxes into three sets:
    *    - those that belong to me
    *    - those that belong to my left child
    *    - those that belong to my right child
    */
   tbox::Array< typename BoxTreeNode<DIM>::Triple> left, right;
   int mid = (d_domain.lower(d_dim) + d_domain.upper(d_dim))/2;

   //first, compute the number of items in each list
   //(this loop would be eliminated if we had an Array::appendItem())
   int my_ct, left_ct, right_ct;
   my_ct = left_ct = right_ct = 0;
   for (int i=0; i<len; ++i) {
      if (box_triples[i].box.upper(d_dim)  <= mid) {
         ++left_ct;
      } else if (box_triples[i].box.lower(d_dim) > mid) {
         ++right_ct;
      } else {
         ++my_ct;
      }
   }

   //second, resize the arrays and insert the elements
   d_box_triples.resizeArray(my_ct); 
   left.resizeArray(left_ct);
   right.resizeArray(right_ct);
   my_ct = left_ct = right_ct = 0;

   for (int i=0; i<len; ++i) {
      if (box_triples[i].box.upper(d_dim)  <= mid) {
         left[left_ct++] = box_triples[i];
      } else if (box_triples[i].box.lower(d_dim) > mid) {
         right[right_ct++] = box_triples[i];
      } else {
         d_box_triples[my_ct++] = box_triples[i];
      }
   }

   /*
    * Recurse to build a private tree on this node;
    * this tree contains the boxes assigned to this node,
    * but they will be partitioned with referenc to
    * the next lower dimension.
    */
   if (d_dim > 0) {
      d_tree = new BoxTreeNode<DIM>(d_box_triples, have_mapping, 
                                    min_length, d_dim-1, recurse_level+1);
      d_box_triples.setNull();
   } else {
      d_work.resizeArray( d_box_triples.size() );
   }

   /*
    * Recurse to build this node's left and right children.
    */
   if (!left.isNull()) {
      d_left_child = new BoxTreeNode<DIM>(left, have_mapping, 
                                          min_length, d_dim-1, recurse_level+1);
   }

   if (!right.isNull()) {
      d_right_child = new BoxTreeNode<DIM>(right, have_mapping, 
                                           min_length, d_dim-1, recurse_level+1);
   }
}


/*
 * ************************************************************************
 *                                                                       *
 * dtor
 *                                                                       *
 * ************************************************************************
 */

template<int DIM>  BoxTreeNode<DIM>::~BoxTreeNode()
{
}

/*
 * ************************************************************************
 * 
 *  Append to "indices" the indice of the boxes owned by this
 *  node that intersect with "box."  Then recurse (1) for child
 *  nodes; (2) for the tree rooted at this node.
 * 
 * ************************************************************************
 */

template<int DIM> void BoxTreeNode<DIM>::findOverlapIndices(
   tbox::List<int> &indices,
   const Box<DIM> & box,
   bool find_local_boxes,
   int recurse_level)
{
   if (box.intersects(d_domain)) {
      if (d_tree) {
         d_tree->findOverlapIndices(indices, box, find_local_boxes);
      }

      /*
       * "d_work.isNull()" if the calling BoxTree was not passed processor
       * mapping information.  The import is that all indices int d_box_triples
       * are unique, so we don't need to check for duplicates.
       */
      else if (d_work.isNull()) {
         int len = d_box_triples.size();

         //case 1: only interested in local overlapping boxes
         if (find_local_boxes) {
            for (int i=0; i<len; ++i) {
               if (d_box_triples[i].owner == d_rank) {
                  if (box.intersects(d_box_triples[i].box)) {
                     indices.appendItem(d_box_triples[i].idx);
                  }
               }
            }
         //case 2: interested in all overlapping boxes
         } else {
            for (int i=0; i<len; ++i) {
               if (box.intersects(d_box_triples[i].box)) {
                  indices.appendItem(d_box_triples[i].idx);
               }
            }
         }
      }

      /*
       * If there is a mapping, then it's possible that same index
       * could be appended more than once, so extra care is required
       * to guard against duplicate indices.
       */
      else {

         int len = d_box_triples.size();
         int idx = 0;

         //case 1: only interested in local overlapping boxes
         if (find_local_boxes) {
            for (int i=0; i<len; ++i) {
               if (d_box_triples[i].owner == d_rank) {
                  if (box.intersects(d_box_triples[i].box)) {
                     d_work[idx++] = d_box_triples[i].idx;
                  }
               }
            }
         //case 2: interested in all overlapping boxes
         } else {
            for (int i=0; i<len; ++i) {
               if (box.intersects(d_box_triples[i].box)) {
                  d_work[idx++] = d_box_triples[i].idx;
               }
            }
         }

         for (int j=0; j<idx; ++j) {
            indices.appendItem(d_work[j]);
         }

      }

      if (d_left_child) {
         d_left_child->findOverlapIndices(indices, box, 
                                          find_local_boxes, recurse_level+1);
      }

      if (d_right_child) {
         d_right_child->findOverlapIndices(indices, box, 
                                           find_local_boxes, recurse_level+1);
      }
   }

   /*
    * check for and eliminate any duplicate indices that
    * may have crept in during the recursive calls to 
    * left or right children.  Note: this really only needs
    * to be done before exiting at the top recursion level.
    */

   int size = indices.size();
   if (size > 1) {
      tbox::Array<int> sortme(size);
      int j = 0;
      for (tbox::List<int>::Iterator l(indices); l; l++) {
         sortme[j++] = l();
      }

      int *ptr = sortme.getPointer();
      qsort((void*)ptr, size, sizeof(int),
         hier::BoxGraphUtilities<DIM>::qsortIntCompare);

      indices.clearItems();
      indices.appendItem(sortme[0]);
      for (int k=1; k<size; ++k) {
         if (sortme[k] != sortme[k-1]) {
            indices.appendItem(sortme[k]);
         }
      }
   }

}

template<int DIM> void BoxTreeNode<DIM>::findOverlapBoxes(
   BoxList<DIM> &overlap_boxes,
   const Box<DIM> & box)
{
   if (box.intersects(d_domain)) {
      if (d_tree) {
         d_tree->findOverlapBoxes(overlap_boxes, box);
      }

      else {

         int len = d_box_triples.size();
         for (int i=0; i<len; ++i) {
            if (box.intersects(d_box_triples[i].box)) {
               overlap_boxes.appendItem( d_box_triples[i].box);
            }
         }
      }

      if (d_left_child) {
         d_left_child->findOverlapBoxes(overlap_boxes, box);
      }

      if (d_right_child) {
         d_right_child->findOverlapBoxes(overlap_boxes, box);
      }
   }
}

}
}

#endif
