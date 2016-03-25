//
// File:        BoxTree.C
// Package:     SAMRAI hierarchy
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 266 $
// Modified:    $Date: 2005-03-15 16:26:07 -0800 (Tue, 15 Mar 2005) $
// Description: Utility class to reduce complexity of box calculus operations.
//

#ifndef included_hier_BoxTree_C
#define included_hier_BoxTree_C

#include "BoxTree.h"
#include "BoxGraphUtilities.h"
#include "tbox/Utilities.h"

#include <iostream>
#include <iomanip>
using namespace std;

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif


namespace SAMRAI {
   namespace hier {

static void buildTboxArrayFromList(
   tbox::Array<int> &array,
   const tbox::List<int> &list)
{
   int len = list.getNumberItems();
   array.resizeArray(len);   
   int count = 0;
   for (tbox::List<int>::Iterator j(list); j; j++) {
      array[count++] = j();
   }
}

/*
 * ************************************************************************
 *
 * ctors and dtor
 *
 * ************************************************************************
 */

template<int DIM>  BoxTree<DIM>::BoxTree(
   const BoxArray<DIM> & box_array,
   const tbox::Array< tbox::List< IntVector<DIM> > >& box_shifts,
   const ProcessorMapping & box_mapping,
   int min_length)
: d_boxes(box_array),
  d_have_mapping(true)
{
   tbox::Array<int> indices;
   BoxArray<DIM>  boxes;
   BoxGraphUtilities<DIM>::makeBoxesPlusPeriodicBoxes(
      boxes, indices, box_array, box_shifts);

   int len = boxes.getNumberOfBoxes();
   tbox::Array<typename BoxTreeNode<DIM>::Triple> box_triples(len);
   for (int i=0; i<len; ++i) {
      int idx = indices[i];
      box_triples[i].box = boxes(i);
      box_triples[i].idx = idx;
      box_triples[i].owner = box_mapping.getProcessorAssignment(idx);
   }

   d_tree = new BoxTreeNode<DIM>(box_triples, d_have_mapping, min_length);
}

template<int DIM>  BoxTree<DIM>::BoxTree(
   const BoxArray<DIM> & box_array,
   const tbox::Array< tbox::List< IntVector<DIM> > >& box_shifts,
   int min_length)
: d_boxes(box_array),
  d_have_mapping(false)
{
   tbox::Array<int> indices;
   BoxArray<DIM>  boxes;
   BoxGraphUtilities<DIM>::makeBoxesPlusPeriodicBoxes(
      boxes, indices, box_array, box_shifts);

   int len = boxes.getNumberOfBoxes();
   tbox::Array<typename BoxTreeNode<DIM>::Triple> box_triples(len);
   for (int i=0; i<len; ++i) {
      int idx = indices[i];
      box_triples[i].box = boxes(i);
      box_triples[i].idx = idx;
      box_triples[i].owner = -1;
   }

   d_tree = new BoxTreeNode<DIM>(box_triples, d_have_mapping, min_length);
}

template<int DIM>  BoxTree<DIM>::BoxTree(
   const BoxArray<DIM> & box_array,
   int min_length)
: d_boxes(box_array),
  d_have_mapping(false)
{
   int len = d_boxes.getNumberOfBoxes();
   tbox::Array<typename BoxTreeNode<DIM>::Triple> box_triples(len);
   for (int i=0; i<len; ++i) {
      box_triples[i].box = d_boxes(i);
      box_triples[i].idx = i;
      box_triples[i].owner = -1;
   }

   d_tree = new BoxTreeNode<DIM>(box_triples, d_have_mapping, min_length);
}


template<int DIM>  BoxTree<DIM>::~BoxTree()
{
}

/*
 *************************************************************************
 *
 * Routines to return an array of integer indices indicating elements
 * of the box array passed to the ctor that overlap with the specified box.  
 * Depending on which routine is called, the indices will either correspond 
 * to all boxes in the original box array or only those that are local to 
 * this processor. 
 *
 *************************************************************************
 */

template<int DIM> void BoxTree<DIM>::findOverlapIndices(
   tbox::Array<int>& indices,
   const Box<DIM>& box) const
{
   bool find_local_overlap_indices = false;
   privateFindOverlapIndices(indices, box, find_local_overlap_indices);
}

template<int DIM> void BoxTree<DIM>::findLocalOverlapIndices(
   tbox::Array<int>& indices,
   const Box<DIM>& box) const
{
   if (!d_have_mapping) {
      TBOX_ERROR("BoxTree<DIM>::findLocalOverlapIndices() error!!"
                 << "\n Processor mapping must be passed to class constructor"
                 << "to use this method. None was passed." << endl);
   }

   bool find_local_overlap_indices = true;
   privateFindOverlapIndices(indices, box, find_local_overlap_indices);
}

template<int DIM> void BoxTree<DIM>::privateFindOverlapIndices(
   tbox::Array<int>& indices,
   const Box<DIM>& box,
   bool find_local_overlaps) const
{
   indices.setNull();
   tbox::List<int> indices_list;
   d_tree->findOverlapIndices(indices_list, box, find_local_overlaps);

   buildTboxArrayFromList(indices, indices_list);

   int len = indices.getSize();
   if (len) {
      int *idx_ptr = indices.getPointer();
      qsort((void*)idx_ptr, len, sizeof(int),
            BoxGraphUtilities<DIM>::qsortIntCompare);
   }
}

/*
 *************************************************************************
 *
 * Routine that returns a list of boxes that overlap with the argument box.
 *
 *************************************************************************
 */

template<int DIM> void BoxTree<DIM>::findOverlapBoxes(
   BoxList<DIM>& overlap_boxes,
   const Box<DIM>& box) const
{
   overlap_boxes.clearItems();
   tbox::List<int> indices_list;
   d_tree->findOverlapBoxes(overlap_boxes, box);
}

/*
 * **************************************************************
 * CAUTION: the semantics of this call differ from that of
 * BoxList::removeIntersections(const BoxList<DIM> takeaway);
 * here, the list that is being modified is the list that is
 * passed as an argument; the "takeaway" list is the list that
 * was passed when the BoxTop object was constructed.
 * **************************************************************
 */

template<int DIM> void BoxTree<DIM>::removeIntersections(BoxList<DIM>& fragments) const
{
   //This will hold boxes in "list" minus the intersection with "takeaway"
   BoxArray<DIM> boxes_in(fragments);
   fragments.clearItems();

   /*
    * Compare each box in "list" against all boxes in "take_away"
    * Note: this is opposite the approach of the current implementation;
    * here, the "take_away" list is the list with which the tree was
    * instantiated.
    */
   int len = boxes_in.getNumberOfBoxes();
   for (int k=0; k<len; ++k) {
      Box<DIM> tryme = boxes_in(k);

      /*
       * Find a shorter list of boxes from "takeaway" that might intersect
       * with "tryme."  Recall that "takeaway" was the BoxList with which
       * this object was instantiated; the next call returns a list of
       * boxes that are in "takeaway," and are also nabors (i.e., overlap with)
       * the "tryme" box.
       */
      BoxList<DIM> nabors;
      findOverlapBoxes(nabors, tryme);
      int nabor_count = nabors.size();

      /*
       * if "tryme" doesn't intersect any boxes on takeaway
       * i.e., if "nabors" is empty" then keep "tryme."
       */
      if (nabor_count == 0) {
        fragments.appendItem(tryme);
      } 

      /*
       * Else, burst up tryme against the nabors.
       * Here, we call the BoxList<DIM> version of
       * removeIntersections.
       */
      else {

         //possibly todo: hand-code a removeIntersections arrays ...
         BoxList<DIM> tryme_burst(tryme);
         tryme_burst.removeIntersections(nabors);
         fragments.copyItems(tryme_burst);
      }
   }
}



}
}


#endif
