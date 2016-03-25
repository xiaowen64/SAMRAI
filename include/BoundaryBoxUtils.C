//
// File:	$URL$
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2006 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Generic utilities for boundary box calculus.
//

#ifndef included_hier_BoundaryBoxUtils_C
#define included_hier_BoundaryBoxUtils_C

#include "BoundaryBoxUtils.h"
#include "Box.h"
#include "tbox/Utilities.h"

#ifdef DEBUG_NO_INLINE
// #include "BoundaryBoxUtils.I"
#endif


namespace SAMRAI {
   namespace hier {


template<int DIM>
BoundaryBoxUtils<DIM>::BoundaryBoxUtils()
   : d_bbox(),
     d_outward(0)
{
   return;
}


template<int DIM>
BoundaryBoxUtils<DIM>::BoundaryBoxUtils(const BoundaryBox<DIM>& bbox)
   : d_bbox(bbox),
     d_outward(0)
{
   computeOutwardShift();
   return;
}


template<int DIM>
BoundaryBoxUtils<DIM>::~BoundaryBoxUtils()
{
   return;
}


template<int DIM>
void BoundaryBoxUtils<DIM>::setBoundaryBox(const BoundaryBox<DIM>& bbox)
{
   d_bbox = bbox;
   computeOutwardShift();
   return;
}


template<int DIM>
const BoundaryBox<DIM> &BoundaryBoxUtils<DIM>::getBoundaryBox() const
{
   return d_bbox;
}


template<int DIM>
const IntVector<DIM> &BoundaryBoxUtils<DIM>::getOutwardShift() const
{
   return d_outward;
}


template<int DIM>
void BoundaryBoxUtils<DIM>::computeOutwardShift()
{
   /*
    * Note that d_outward contains information that is redundant
    * with respect to the boundary box.  The values of d_outward
    * depends strictly the location of the boundary box.
    */
   const int lidx = d_bbox.getLocationIndex();

   switch ( d_bbox.getBoundaryType() ) {

   // Note: number of non-zero in d_outward is the same as boundary type.

   case 1:
   {
      int i = lidx/2;
      d_outward(i) = lidx%2 == 0 ? -1 /* lower side */ : 1 /* upper side */;
   }
   break;

   case 2:
   {
      if ( DIM == 2 ) {
         d_outward(0) = lidx%2 == 0 ? -1 : 1;
         d_outward(1) = lidx/2 == 0 ? -1 : 1;
      }
      else if ( DIM == 3 ) {
         const int dir = lidx/4;
         const int rem = lidx%4;
         if ( dir == 0 ) {
            // Nonzero in dirs 1 and 2.
            d_outward(1) = rem%2 == 0 ? -1 : 1;
            d_outward(2) = rem/2 == 0 ? -1 : 1;
         }
         else if ( dir == 1 ) {
            // Nonzero in dirs 0 and 2.
            d_outward(0) = rem/2 == 0 ? -1 : 1;
            d_outward(2) = rem%2 == 0 ? -1 : 1;
         }
         else {
            // Nonzero in dirs 0 and 1.
            d_outward(0) = rem%2 == 0 ? -1 : 1;
            d_outward(1) = rem/2 == 0 ? -1 : 1;
         }
      }
      else {
         TBOX_ERROR("BoundaryBoxUtils cannot compute\n"
                    <<"boundary direction for " << d_bbox.getBox());
      }
   }
   break;

   case 3:
   {
      if ( DIM == 3 ) {
         d_outward(0) = lidx%2 == 0     ? -1 : 1;
         d_outward(1) = (lidx%4)/2 == 0 ? -1 : 1;
         d_outward(2) = lidx/4 == 0     ? -1 : 1;
      }
      else {
         TBOX_ERROR("BoundaryBoxUtils cannot compute\n"
                    <<"boundary direction for " << d_bbox.getBox());
      }
   }
   break;

   default:
      TBOX_ERROR("BoundaryBoxUtils cannot compute\n"
                 <<"boundary direction for type " << d_bbox.getBoundaryType());
   }
   return;
}


template<int DIM>
void BoundaryBoxUtils<DIM>::stretchBoxToGhostWidth(
   Box<DIM>& box,
   const hier::IntVector<DIM> &gcw ) const
{
   TBOX_ASSERT( gcw >= hier::IntVector<DIM>(0) );

   box = d_bbox.getBox();
   for ( int d=0; d<DIM; ++d ) {
      /*
       * If gcw along dimension d is > 1, stretch it out to that width.
       * If gcw a long dimension d is 0, shrink the box down to nothing
       * in that dimension.
       */
      if ( d_outward(d) == -1 ) {
         if ( gcw(d) > 1 ) box.growLower(d, gcw(d)-1);
         else              box.lower()(d) = box.upper()(d) + 1;
      }
      else if ( d_outward(d) ==  1 ) {
         if ( gcw(d) > 1 ) box.growUpper(d, gcw(d)-1);
         else              box.upper()(d) = box.lower()(d) - 1;
      }
   }
   return;
}


template<int DIM>
void BoundaryBoxUtils<DIM>::extendBoxOutward(
   Box<DIM>& box,
   const hier::IntVector<DIM> &extension ) const
{
   box = d_bbox.getBox();
   for ( int d=0; d<DIM; ++d ) {
      if      ( d_outward(d) == -1 ) box.growLower(d, extension(d));
      else if ( d_outward(d) ==  1 ) box.growUpper(d, extension(d));
   }
   return;
}


template<int DIM>
void BoundaryBoxUtils<DIM>::shiftBoxInward(
   Box<DIM>& box,
   const hier::IntVector<DIM> &distance ) const
{
   box = d_bbox.getBox();
   for ( int d=0; d<DIM; ++d ) {
      if      ( d_outward(d) == -1 ) box.shift(d,  distance(d));
      else if ( d_outward(d) ==  1 ) box.shift(d, -distance(d));
   }
   return;
}


template<int DIM>
int BoundaryBoxUtils<DIM>::normalDir() const
{
   return d_bbox.getLocationIndex()/2;
}

}
}

#endif
