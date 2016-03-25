//
// File:	Box.C
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Box representing a portion of the AMR index space
//

#ifndef included_hier_Box_C
#define included_hier_Box_C

#include "Box.h"

#include "tbox/Utilities.h"

#ifdef DEBUG_NO_INLINE
#include "Box.I"
#endif

namespace SAMRAI {
   namespace hier {


/*
*************************************************************************
*                                                                       *
* Return the dimension of the box that is the longest.                  *
*                                                                       *
*************************************************************************
*/

template<int DIM> int Box<DIM>::longestDimension() const
{
  int max = upper(0)-lower(0);
  int dim = 0;
 
  for (int i=1; i< DIM; i++)
    if ((upper(i)-lower(i)) > max) {
      max = upper(i)-lower(i);
      dim = i;
    }
  return dim;
}

/*
*************************************************************************
*									*
* Type Conversions							*
*									*
*************************************************************************
*/

template<int DIM>
tbox::DatabaseBox Box<DIM>::DatabaseBox_from_Box() const
{
   tbox::DatabaseBox new_Box;

   new_Box.setDimension(DIM);

   for (int i = 0; i < DIM; i++) {
      new_Box.lower(i) = d_lo(i);
      new_Box.upper(i) = d_hi(i);
   }

   return new_Box;
}

template<int DIM>
void Box<DIM>::set_Box_from_DatabaseBox(const tbox::DatabaseBox& box)
{
   for (int i = 0; i < box.getDimension(); i++) {
      d_lo(i) = box.lower(i);
      d_hi(i) = box.upper(i);
   }
}

/*
*************************************************************************
*									*
* Stream input/output operators: [(l0,...,ln),(u0,...,un)].		*
*									*
*************************************************************************
*/
template<int DIM>
istream& operator>>(istream& s, Box<DIM>& box)
{
   while (s.get() != '[');
   s >> box.lower();
   while (s.get() != ',') NULL_STATEMENT;
   s >> box.upper();
   while (s.get() != ']') NULL_STATEMENT;
   return(s);
}

template<int DIM>
ostream& operator<<(ostream& s, const Box<DIM>& box)
{
   if (box.empty()) {
      s << "[(),()]";
   } else {
      s << '[' << box.lower() << ',' << box.upper() << ']';
   }
   return(s);
}

template<int DIM> Box<DIM>& Box<DIM>::operator+=(const Box<DIM>& box)
{
   if (!box.empty()) {
      if (empty()) {
         *this = box;
      } else {
         d_lo.min(box.d_lo);
         d_hi.max(box.d_hi);
      }
   }
   return(*this);
}

/*
*************************************************************************
*                                                                       *
* Static member function called from coalesceWith().  It attempts to    *
* recursively coalesce intervals individual dimensions in index space.  *
* If it is possible to coalesce two intervals (defined by a proper      *
* overlap or adjacency relationship), the value true is returned.       *
* If this is impossible, false is returned.                             *
*                                                                       *
*************************************************************************
*/

template<int DIM> bool Box<DIM>::coalesceIntervals(
   const int* lo1, const int* hi1,
   const int* lo2, const int* hi2,
   const int dim)
{
   bool retval = false;
   if (dim == 1) {
      // interval 1 to the right of interval 2.
      if ( (lo1[0] <= hi2[0]+1) && (hi2[0] <= hi1[0]) ) {
         retval = true; return(retval);
      }
      // interval 1 to the left of interval 2.
      if ( (lo1[0] <= lo2[0]) && (lo2[0] <= hi1[0]+1) ) {
         retval = true; return(retval);
      }
   } else {
      for (int id = 0; id < dim; id++) {
         if ( (lo1[id] == lo2[id]) && (hi1[id] == hi2[id]) ) {
            int id2;
            int low1[DIM], high1[DIM], low2[DIM], high2[DIM];
            for (id2 = 0; id2 < id; id2++) {
               low1[id2] = lo1[id2]; high1[id2] = hi1[id2];
               low2[id2] = lo2[id2]; high2[id2] = hi2[id2];
            }
            for (id2 = id+1; id2 < dim; id2++) {
               int id1 = id2-1;
               low1[id1] = lo1[id2]; high1[id1] = hi1[id2];
               low2[id1] = lo2[id2]; high2[id1] = hi2[id2];
            }
            if ( coalesceIntervals(low1, high1, low2, high2, dim-1) ) {
               retval = true; return(retval);
            }
         }
      }
   }

   return(retval);
}

/*
*************************************************************************
*                                                                       *
* Return true if this box can be coalesced with the argument box,       *
* and set this box to the union of the boxes.  Otherwise, return false  *
* and leave this box as is.  Two boxes may be coalesced if their union  *
* is a box.  This routine attempts to coalesce the boxes along          *
* each coordinate direction using the coalesceIntervals() function.     *
*                                                                       *
*************************************************************************
*/

template <int DIM>
bool Box<DIM>::coalesceWith(const Box<DIM>& box)
{
   bool retval = false;

   if (empty() || box.empty()) {
      retval = true;
      *this += box;
   } else {
      int id;
      const int* box_lo = box.lower();
      const int* box_hi = box.upper();
      int me_lo[DIM], me_hi[DIM];
      for (id = 0; id < DIM; id++) {
         me_lo[id] = d_lo(id); me_hi[id] = d_hi(id);
      }
      if ( coalesceIntervals(box_lo, box_hi, me_lo, me_hi, DIM) ) {
         retval = true;
      } else { // test for one box containing the other...
         // test whether me contains box.
         retval = true; id = 0;
         while (retval && (id < DIM)) {
            retval = ((me_lo[id] <= box_lo[id]) && (me_hi[id] >= box_hi[id]));
            id++;
         }
         if (!retval) { // me doesn't contain box; check other way around...
            retval = true; id = 0;
            while (retval && (id < DIM)) {
               retval = (    (box_lo[id] <= me_lo[id])
                          && (box_hi[id] >= me_hi[id]) );
               id++;
            }
         }
      }
   }

   if (retval) *this += box;

   return(retval);
}

/*
*************************************************************************
*                                                                       *
* Rotates a 3-Dimensional box 45*num_rotations degrees around the given *
* and set this box to the union of the boxes.                           *
*                                                                       *
*************************************************************************
*/

template<int DIM>
void Box<DIM>::rotateAboutAxis(const int axis, const int num_rotations)
{
                                                                                
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(axis < DIM);
   assert(DIM == 3);
#endif
                                                                                
   const int a = (axis+1)% DIM;
   const int b = (axis+2)% DIM;
                                                                                
   Index<DIM> tmp_lo;
   Index<DIM> tmp_hi;
   for (int j = 0; j < num_rotations; j++) {
      tmp_lo = d_lo;
      tmp_hi = d_hi;
      d_lo(a) = tmp_lo(b);
      d_lo(b) = -tmp_hi(a)-1;
      d_hi(a) = tmp_hi(b);
      d_hi(b) = -tmp_lo(a)-1;
   }
}

/*
*************************************************************************
*                                                                       *
* Rotate a box in the manner determined by the rotation number          *
*                                                                       *
*************************************************************************
*/

template<int DIM>
void Box<DIM>::rotate(const int rotation_number)
{

   if (DIM == 2) {
      for (int j = 0; j < rotation_number; j++) {
         Index<DIM> tmp_lo(d_lo);
         Index<DIM> tmp_hi(d_hi);

         d_lo(0) = tmp_lo(1);
         d_lo(1) = -tmp_hi(0)-1;
         d_hi(0) = tmp_hi(1);
         d_hi(1) = -tmp_lo(0)-1;
      }
   }

   if (DIM == 3) {
      if (rotation_number == 0) {
         return;
      } else if (rotation_number == 1) {
         rotateAboutAxis(0,3);
         rotateAboutAxis(2,3);
      } else if (rotation_number == 2) {
         rotateAboutAxis(1,1);
         rotateAboutAxis(2,1);
      } else if (rotation_number == 3) {
         rotateAboutAxis(1,2);
         rotateAboutAxis(0,3);
      } else if (rotation_number == 4) {
         rotateAboutAxis(1,3);
      } else if (rotation_number == 5) {
         rotateAboutAxis(2,1);
      } else if (rotation_number == 6) {
         rotateAboutAxis(1,1);
      } else if (rotation_number == 7) {
         rotateAboutAxis(0,3);
      } else if (rotation_number == 8) {
         rotateAboutAxis(0,2);
         rotateAboutAxis(2,3);
      } else if (rotation_number == 9) {
         rotateAboutAxis(0,3);
         rotateAboutAxis(2,1);
      } else if (rotation_number == 10) {
         rotateAboutAxis(1,2);
      } else if (rotation_number == 11) {
         rotateAboutAxis(0,3);
         rotateAboutAxis(1,3);
      } else if (rotation_number == 12) {
         rotateAboutAxis(2,3);
      } else if (rotation_number == 13) {
         rotateAboutAxis(0,1);
      } else if (rotation_number == 14) {
         rotateAboutAxis(0,2);
         rotateAboutAxis(1,1);
      } else if (rotation_number == 15) {
         rotateAboutAxis(0,1);
         rotateAboutAxis(1,3);
      } else if (rotation_number == 16) {
         rotateAboutAxis(0,2);
         rotateAboutAxis(1,2);
      } else if (rotation_number == 17) {
         rotateAboutAxis(0,1);
         rotateAboutAxis(2,1);
      } else if (rotation_number == 18) {
         rotateAboutAxis(0,3);
         rotateAboutAxis(1,1);
      } else if (rotation_number == 19) {
         rotateAboutAxis(0,1);
         rotateAboutAxis(2,3);
      } else if (rotation_number == 20) {
         rotateAboutAxis(0,2);
      } else if (rotation_number == 21) {
         rotateAboutAxis(0,2);
         rotateAboutAxis(2,1);
      } else if (rotation_number == 22) {
         rotateAboutAxis(0,2);
         rotateAboutAxis(1,3);
      } else if (rotation_number == 23) {
         rotateAboutAxis(1,2);
         rotateAboutAxis(0,1);
      }
   }
}


}
}

#endif
