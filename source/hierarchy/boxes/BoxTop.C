//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/hierarchy/boxes/BoxTop.C $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Utility class to reduce complexity of box calculus operations.
//

#ifndef included_hier_BoxTop_C
#define included_hier_BoxTop_C

#include "BoxTop.h"
#include "BoxList.h"
#include "BoxGraphUtilities.h"
#include "tbox/SAMRAI_MPI.h"


namespace SAMRAI {
   namespace hier {


/*
*************************************************************************
*                                                                       *
* Constructors and Destructor                                           *
*                                                                       *
*************************************************************************
*/

template<int DIM>  BoxTop<DIM>::BoxTop(
   const BoxArray<DIM>& in_boxes,
   const tbox::Array< tbox::List< IntVector<DIM> > >& shifts)
{
   BoxGraphUtilities<DIM>::makeBoxesPlusPeriodicBoxes(
      d_boxes, in_boxes, shifts);  
  setup();
}


template<int DIM>  BoxTop<DIM>::BoxTop(const BoxArray<DIM>& in_boxes)
: d_boxes(in_boxes)
{
  setup();
}

template<int DIM>  BoxTop<DIM>::~BoxTop()
{
   for (int i=0; i<DIM*2; ++i) {
      boxElt *tmp = d_sorted_lists[i];
      delete [] tmp;
   }
}


/*
************************************************************************
*                 
* Constructs an array of boxes (and/or their indices), 
* that is a subset of the array of boxes passed to the constructor.  
*                
*************************************************************************
*/

template<int DIM> void BoxTop<DIM>::findOverlappingBoxes(
   BoxArray<DIM> &overlaps,
   const Box<DIM> & box)
{
   findNaborsPrivate(box);
   tbox::Array<int> indices;
   bool build_indices = false;
   bool build_overlaps = true;
   buildShortestList(overlaps, build_overlaps, indices, build_indices, box);
}

template<int DIM> void BoxTop<DIM>::findOverlappingBoxIndices(
   tbox::Array<int> &indices,
   const Box<DIM> & box)
{
   findNaborsPrivate(box);
   BoxArray<DIM> overlaps;
   bool build_indices = true;
   bool build_overlaps = false;
   buildShortestList(overlaps, build_overlaps, indices, build_indices, box);
}

template<int DIM> void BoxTop<DIM>::findOverlappingBoxesAndIndices(
   BoxArray<DIM> &overlaps,
   tbox::Array<int> &indices,
   const Box<DIM> & box)
{
   findNaborsPrivate(box);
   bool build_indices = true;
   bool build_overlaps = true;
   buildShortestList(overlaps, build_overlaps, indices, build_indices, box);
}


/*
*************************************************************************
*                                                                       *
*  
*                                                                       *
*************************************************************************
*/

template<int DIM> void BoxTop<DIM>::removeIntersections(BoxList<DIM>& fragments)
{
   //This will hold boxes in "list" minus the intersection with "takeaway"
   BoxArray<DIM> boxes_in(fragments);
   fragments.clearItems();


   //Compare each box in "list" against all boxes in "takeaway"
   //Note: this is opposite the approach of the current implementation
   int len = boxes_in.getNumberOfBoxes();
   for (int k=0; k<len; ++k) {
      Box<DIM> tryme = boxes_in(k);

      //Find a shorter list of boxes from "takeaway" that might intersect
      //with "tryme."  Recall that "takeaway" was the list with which
      //the BoxTop was instantiated; the next call returns a list of
      //boxes that in "takeaway," and are also nabors (i.e., overlap with)
      //the "tryme" box.
      //
      //Cost of this operation is 
      //O( max{ lg(n)), nabors.getNumberItems() } ).
      BoxArray<DIM> nabors_array;
      findOverlappingBoxes(nabors_array, tryme);
      BoxList<DIM> nabors(nabors_array);

      //if "tryme" doesn't intersect any boxes on takeaway
      //(i.e., if "nabors" is empty" then keep "tryme."
      if (nabors.getNumberItems() == 0) {
        fragments.appendItem(tryme);
      } 

      //else, burst up tryme against the nabors.
      else {
         BoxList<DIM> tryme_burst(tryme);
         tryme_burst.removeIntersections(nabors);
         fragments.copyItems(tryme_burst);
      }
   }
}

/*
*************************************************************************
*                                                                       *
*  Performs the search function.                                        *
*                                                                       *
*************************************************************************
*/

template<int DIM> void BoxTop<DIM>::findNaborsPrivate(const Box<DIM>& box)
{
   int len = d_boxes.getNumberOfBoxes();
   d_shortest_list_length = len;
   d_shortest_list = d_sorted_lists[0];

   //test "box" coordinate against sorted lists to find the
   //smallest possible subset of d_boxes that might
   //intersect "box."  The resulting list is stored in d_shortest_list.
   boxElt *tmp;
   int count; 
   int offset = -1;

   for (int k=0; k<DIM; ++k) {

      //Compare box's upper coordinate against the box
      //list sorted by lower coordinates.
      ++offset;
      tmp = d_sorted_lists[offset];
      int up = box.upper(k);
      count = binSearch(tmp, len, up, ASCENDING);
      if (count < d_shortest_list_length) {
        d_shortest_list_length = count;
        d_shortest_list = tmp;
      }

      //Compare box's lower coordinate against the box
      //list sorted by upper coordinates.
      ++offset;
      tmp = d_sorted_lists[offset];
      int dn = box.lower(k);
      count = binSearch(tmp, len, dn, DESCENDING);
      if (count < d_shortest_list_length) {
        d_shortest_list_length = count;
        d_shortest_list = tmp;
      }
   }
}

/*
*************************************************************************
*                                                                       *
* Print functions, used during development and testing.                 *
* These print the arrays of sorted boxElts.                             *
*                                                                       *
*************************************************************************
*/

template<int DIM> void BoxTop<DIM>::print(std::ostream& os)
{
   os << "\n----- BoxTop::print -----\n";

   int length = d_boxes.getNumberOfBoxes();
   for (int i=0; i<DIM; ++i) {
      os << "Dimension: " << i+1 << std::endl << std::endl;

      os << "  sorted by upper:\n";
      printEltArray(d_sorted_lists[i*2+1], length, os);

      os << "  sorted by lower:\n";
      printEltArray(d_sorted_lists[i*2], length, os);
   }
}


/*
*************************************************************************
* Note: All functions from here to end-of-file are private.             *
*************************************************************************
*/


/*
*************************************************************************
*                                                                       *
* Setup the sorted arrays of boxElts.                                   *
* (Private function)                                                    *
*                                                                       *
*************************************************************************
*/

template<int DIM> void BoxTop<DIM>::setup()
{
   int length = d_boxes.getNumberOfBoxes();

   int offset = -1;
   for (int k=0; k<DIM; ++k) {
      ++offset;

      //Create and initialize array for sorting boxes by lower coordinate;
      //boxes are sorted by lower coordinate in ascending order
      d_sorted_lists[offset] = new boxElt[length];
      boxElt *tmp = d_sorted_lists[offset];
      int i;
      for (i=0; i<length; ++i) {
         tmp[i].coord = d_boxes.getBox(i).lower(k);
         tmp[i].box = &(d_boxes.getBox(i));
         tmp[i].idx = i;
      }
      //sort the arrays
      sort(tmp, length, ASCENDING);

      //Create and initialize array for sorting boxes by upper coordinate
      //boxes are sorted by upper coordinate in descending order.
      ++offset;
      d_sorted_lists[offset] = new boxElt[length];
      tmp = d_sorted_lists[offset];
      for (i=0; i<length; ++i) {
         tmp[i].coord = d_boxes.getBox(i).upper(k);
         tmp[i].box = &(d_boxes.getBox(i));
         tmp[i].idx = i;
      }

      //sort the arrays
      sort(tmp, length, DESCENDING);
   }
}

/*
*************************************************************************
*                                                                       *
* Prints an array of boxElts.  Called by print.                         *
* (Private function)                                                    *
*                                                                       *
*************************************************************************
*/

template<int DIM> void BoxTop<DIM>::printEltArray(boxElt *&data, int len, std::ostream& os)
{
   for (int j=0; j<len; ++j) {
      os << "    coord: " << data[j].coord 
         << "  box: " << *(data[j].box) << std::endl;
   }
   os << std::endl;
}


/*
*************************************************************************
*                                                                       *
*  Returns the number of elements in "tmp" that have coord fields that  *
*  are less than "edge," if ordering == ASCENDING.                      *
*  If ordering == DESCENDING, returns number of elts greater than edge. *
* (Private function)                                                    *
*                                                                       *
*************************************************************************
*/

template<int DIM> int BoxTop<DIM>::binSearch(boxElt *data, int data_len, int coord, int ordering)
{
   int idx=0;

   //Determine the number of elements in the (left hand portion) of
   //the array that are smaller than coord.
   //(The array is sorted in ASCENDING order)
   if (ordering == ASCENDING) {

      //extreme case #1 - all elts are larger than coord
      if (data[0].coord > coord) { 
         return 0;
       }

       //extreme case #2 - all elts are smaller than coord
       else if (data[data_len-1].coord < coord) {
         return data_len;
       }

      //for all other cases, use binary search
      else {
         int lo = 0;
         int hi = data_len-1;
         while (hi >= lo) {
            idx = (hi+lo)/2;
            if (data[idx].coord == coord) {
               break;
            } else if (data[idx].coord > coord) {
               hi = idx - 1;
            } else if (data[idx].coord < coord) {
              lo = idx + 1;
            }
         }

         //there may be more than one box with the same coord,
         //so must make sure we get 'em all!
         while ( (idx < data_len) && (data[idx].coord <= coord) ) {
           ++idx;
         }
      }
   }
   

   //Determine the number of elements in the (left hand portion) of
   //the array that are larger than coord.
   //(The array is sorted in DESCENDING order)
   else {  

      //extreme case #1 - all elts are smaller than coord
      if (data[0].coord < coord) {
        return 0;
      }

      //extreme case #2 - all elts are larger than coord
      else if (data[data_len-1].coord > coord) { 
         return data_len;
      }

      //for all other cases, use binary search
      else {
         int lo = 0;
         int hi = data_len-1;
         while (hi >= lo) {
            idx = (hi+lo)/2;
            if (data[idx].coord == coord) {
               break;
            } else if (data[idx].coord < coord) { 
               hi = idx - 1;
            } else if (data[idx].coord > coord) {
               lo = idx + 1;
            }
         }
         while ( (idx < data_len) && (data[idx].coord >= coord) ) {
            ++idx;
         }
      }
   }
   return idx;
}

/*
*************************************************************************
*                                                                       *
* Constructs ``BoxList d_short_list'' from ``boxElt *d_shortest_list.'' *
* On entrance, the first "d_shortest_list_length" elements of           *
* "d_shortestList" contain those boxes that may intersect the given     *
* box.  This function looks at each box, and adds it to the final list  *
* it actually intersects the given box.
*
* (Private function)                                                    *
*                                                                       *
*************************************************************************
*/

template<int DIM> void BoxTop<DIM>::buildShortestList(
   BoxArray<DIM> &overlaps,
   bool build_overlaps,
   tbox::Array<int> &indices,
   bool build_indices,
   const Box<DIM>& box)
{
   //Count the number of neighbors; the count is needed so we can
   //set the size of the arrays.
   int count = 0;
   int i;
   for (i=0; i<d_shortest_list_length; ++i) {
      bool keep = true; 
      for (int k=0; k<DIM; ++k) {
         if (   box.lower(k) > d_shortest_list[i].box->upper(k)
             || box.upper(k) < d_shortest_list[i].box->lower(k)) {
            keep = false;
            k = DIM;
         }
      }

      if (keep) {
         ++count;
      }
   }

   if (build_indices) {
      indices.resizeArray(count);
   }
   if (build_overlaps) {
     overlaps.resizeBoxArray(count);
  }

   //go through the same loop again, but this time insert the indices
   count = 0;
   for (i=0; i<d_shortest_list_length; ++i) {
      bool keep = true;
      for (int k=0; k<DIM; ++k) {
         if (   box.lower(k) > d_shortest_list[i].box->upper(k)
             || box.upper(k) < d_shortest_list[i].box->lower(k)) {
            keep = false;
            k = DIM;
         }
      }

      if (keep)  {
         if (build_indices) {
            indices[count]  = d_shortest_list[i].idx;
         }
         if (build_overlaps) {
            overlaps.getBox(count) = ( *(d_shortest_list[i].box) );
         }
         count++;
      }
   }
}

/*
*************************************************************************
*                                                                       *
* Sorts an array of boxElts.  This is a wrapper around qsort, in case   *
* we want to use a different sort function at some future time.         *
* (Private function)                                                   *
*                                                                       *
*************************************************************************
*/

template<int DIM> void BoxTop<DIM>::sort(boxElt *&array, int len, int ordering)
{
   if (ordering == ASCENDING) {
      qsort(array, len, sizeof(boxElt), boxEltCompareA);
   } else {  // ordering == DESCENDING
      qsort(array, len, sizeof(boxElt), boxEltCompareD);
   }
}


/*
*************************************************************************
*                                                                       *
* Static functions needed for qsort.                                    *
* (Private functions)                                                   *
*                                                                       *
*************************************************************************
*/

template<int DIM> int BoxTop<DIM>::boxEltCompareA(const void *v, const void *w)
{
   boxElt *a = (boxElt*) v;
   boxElt *b = (boxElt*) w;
   int i = a->coord;
   int j = b->coord;
 
   if (i > j) return 1;
   if (i < j) return -1;
   return 0;
}

template<int DIM> int BoxTop<DIM>::boxEltCompareD(const void *v, const void *w)
{
   boxElt *a = (boxElt*) v;
   boxElt *b = (boxElt*) w;
   int i = a->coord;
   int j = b->coord;
 
   if (i < j) return 1;
   if (i > j) return -1;
   return 0;
}


}
}

#endif
