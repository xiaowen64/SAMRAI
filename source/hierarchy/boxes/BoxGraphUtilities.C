//
// File:        BoxGraphUtilities.C
// Package:     SAMRAI hierarchy
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Utility class for operations that reduce complexity of box calculus
//

#ifndef included_hier_BoxGraphUtilities_C
#define included_hier_BoxGraphUtilities_C

#include "BoxGraphUtilities.h"
#include "tbox/Utilities.h"


namespace SAMRAI {
   namespace hier {

/*
 * ************************************************************************
 * 
 * Returns an array of boxes that includes entries for each
 * box that touches a periodic boundary.
 * 
 * ************************************************************************
 */
template<int DIM> void BoxGraphUtilities<DIM>::makeBoxesPlusPeriodicBoxes(
   BoxArray<DIM>& out_boxes,
   const BoxArray<DIM>& in_boxes,
   const tbox::Array< tbox::List< IntVector<DIM> > >& shifts)
{
   /*
    * Compute the number of boxes that will be in the list;
    * this is the number of original boxes, plus all the shifts.
    */
   int p_num = in_boxes.getNumberOfBoxes();
   int b_num = countPeriodicBoxes(shifts);
   int ct = p_num + b_num;

   int test = shifts.getSize();
   if (! (test == 0 || test == p_num)) {
      TBOX_ERROR("BoxGraphUtilities<DIM>::makeBoxesPlusPeriodicBoxes() error"
                 << "The shift array must either have zero length, or be the " 
                 << "same length as in_boxes" << endl);
   }

   /*
    *  Resize the array to hold the expanded box list.
    */
   out_boxes.resizeBoxArray(ct);

   /*
    * Finally, fill in the arrays.
    */
   int idx = 0;
   for (int i = 0; i < p_num; i++) {
      /*
       * append the unshifted box
       */
      out_boxes(idx) = in_boxes(i);
      ++idx;

      /*
       * append any associated shifted boxes
       */
      if (b_num) {
         for (typename tbox::List< IntVector<DIM> >::Iterator
            sh(shifts[i]); sh; sh++) {
            out_boxes(idx) = Box<DIM>::shift(in_boxes(i), sh());
            ++idx;
         }
      }
   }
}


/*
 * ************************************************************************
 * 
 * Returns an array of boxes that includes entries for each
 * box that touches a periodic boundary.
 * 
 * ************************************************************************
 */
template<int DIM> void BoxGraphUtilities<DIM>::makeBoxesPlusPeriodicBoxes(
   BoxArray<DIM>& out_boxes,
   tbox::Array<int>& out_indices,
   const BoxArray<DIM>& in_boxes,
   const tbox::Array< tbox::List< IntVector<DIM> > >& shifts)
{
   /*
    * Compute the number of boxes that will be in the list;
    * this is the number of original boxes, plus all the shifts.
    */
   int p_num = in_boxes.getNumberOfBoxes();
   int b_num = countPeriodicBoxes(shifts);
   int ct = p_num + b_num;

   int test = shifts.getSize();
   if (! (test == 0 || test == p_num)) {
      TBOX_ERROR("BoxGraphUtilities<DIM>::makeBoxesPlusPeriodicBoxes() error"
                 << "The shift array must either have zero length, or be the "
                 << "same length as in_boxes" << endl);
   }

   /*
    *  Resize the array to hold the expanded box list.
    */
   out_boxes.resizeBoxArray(ct);
   out_indices.resizeArray(ct);

   /*
    * Finally, fill in the arrays.
    */
   int idx = 0;
   for (int i = 0; i < p_num; i++) {
      /*
       * append the unshifted box
       */
      out_boxes(idx) = in_boxes(i);
      out_indices[idx] = i;
      ++idx;

      /*
       * append any associated shifted boxes
       */
      if (b_num) {
         for (typename tbox::List< IntVector<DIM> >::Iterator
            sh(shifts[i]); sh; sh++) {
            out_boxes(idx) = Box<DIM>::shift(in_boxes(i), sh());
            out_indices[idx] = i;
            ++idx;
         }
      }
   }
}

/*
 * ************************************************************************
 * 
 *  Returns the sum of shifts[j].getNumberItems()
 * 
 * ************************************************************************
 */
template<int DIM> int BoxGraphUtilities<DIM>::countPeriodicBoxes(
   const tbox::Array< tbox::List< IntVector<DIM> > >& shifts)
{
   int count = 0;
   int n = shifts.getSize();
   for (int j=0; j<n; ++j) {
      count += shifts[j].getNumberItems();
   }
  return count;
}


/*
 * ************************************************************************
 * 
 *  for use when sorting integers using the C-library qsort
 * 
 * ************************************************************************
 */
template<int DIM> int BoxGraphUtilities<DIM>::qsortIntCompare(const void *v, const void *w)
{
   int i = *(int *)v;
   int j = *(int *)w;
   if (i > j) {
      return (1);
   }
   if (i < j) {
      return (-1);
   }
   return (0);
}


}
}

#endif
