//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/boxgeometry/OutersideGeometry.C $
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	hier::Box geometry information for outerside centered objects
//

#ifndef included_pdat_OutersideGeometry_C
#define included_pdat_OutersideGeometry_C

#include "OutersideGeometry.h"
#include "BoxList.h"
#include "SideGeometry.h"
#include "SideOverlap.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#include "tbox/Utilities.h"
#endif

#ifdef DEBUG_NO_INLINE
#include "OutersideGeometry.I"
#endif
namespace SAMRAI {
    namespace pdat {

/*
*************************************************************************
*									*
* Create a side geometry object given the box and ghost cell width.	*
*									*
*************************************************************************
*/

template<int DIM>  OutersideGeometry<DIM>::OutersideGeometry(
   const hier::Box<DIM>& box,
   const hier::IntVector<DIM>& ghosts)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(ghosts.min() >= 0);
#endif
   d_box    = box;
   d_ghosts = ghosts;
}

template<int DIM>  OutersideGeometry<DIM>::~OutersideGeometry()
{
}

/*
*************************************************************************
*									*
* Attempt to calculate the intersection between two outerside centered	*
* box geometries.  The calculateOverlap() checks whether both arguments	*
* are outerside geometries; if so, it compuates the intersection.  If	*
* not, then it calls calculateOverlap() on the source object (if retry	*
* is true) to allow the source a chance to calculate the intersection.	*
* See the hier::BoxGeometry<DIM> base class for more information about the	*
* protocol.  A pointer to null is returned if the intersection canot be	*
* computed.								*
* 									*
*************************************************************************
*/

template<int DIM> 
tbox::Pointer< hier::BoxOverlap<DIM> > OutersideGeometry<DIM>::calculateOverlap(
   const hier::BoxGeometry<DIM>& dst_geometry,
   const hier::BoxGeometry<DIM>& src_geometry,
   const hier::Box<DIM>& src_mask,
   const bool overwrite_interior,
   const hier::IntVector<DIM>& src_offset,
   const bool retry) const
{
   const SideGeometry<DIM> *t_dst = 
      dynamic_cast<const SideGeometry<DIM> *>(&dst_geometry);
   const OutersideGeometry<DIM> *t_src =
      dynamic_cast<const OutersideGeometry<DIM> *>(&src_geometry);

   tbox::Pointer< hier::BoxOverlap<DIM> > over = NULL;
   if ((t_src != NULL) && (t_dst != NULL)) {
      over = doOverlap(*t_dst, *t_src, src_mask, overwrite_interior, 
		       src_offset);
   } else if (retry) {
      over = src_geometry.calculateOverlap(
         dst_geometry, src_geometry, src_mask,
         overwrite_interior, src_offset, false);
   }
   return(over);
}

/*
*************************************************************************
*									*
* Compute the overlap between a side geometry destination box and an	*
* outerside geometry source box.  The intersection algorithm is similar	*
* the side geometry algorithm except that only the borders of source	*
* are used in the intersection computation.				*
*									*
*************************************************************************
*/

template<int DIM> 
tbox::Pointer< hier::BoxOverlap<DIM> > OutersideGeometry<DIM>::doOverlap(
   const SideGeometry<DIM>& dst_geometry,
   const OutersideGeometry<DIM>& src_geometry,
   const hier::Box<DIM>& src_mask,
   const bool overwrite_interior,
   const hier::IntVector<DIM>& src_offset)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(dst_geometry.getDirectionVector() == hier::IntVector<DIM>(1));
#endif

   hier::BoxList<DIM> dst_boxes[DIM];

   // Perform a quick-and-dirty intersection to see if the boxes might overlap

   const hier::Box<DIM> src_box =
      hier::Box<DIM>::grow(src_geometry.d_box, src_geometry.d_ghosts) * src_mask;
   const hier::Box<DIM> src_shift =
      hier::Box<DIM>::shift(src_box, src_offset);
   const hier::Box<DIM> dst_ghost =
      hier::Box<DIM>::grow(dst_geometry.getBox(), dst_geometry.getGhosts());

   // Compute the intersection (if any) for each of the side directions

   const hier::Box<DIM> quick_check =
      hier::Box<DIM>::grow(src_shift, 1) * hier::Box<DIM>::grow(dst_ghost, 1);

   if (!quick_check.empty()) {

      const hier::Box<DIM> mask_shift = hier::Box<DIM>::shift(src_mask, src_offset);

      for (int d = 0; d < DIM; d++) {

         const hier::Box<DIM> msk_side =
            SideGeometry<DIM>::toSideBox(mask_shift, d);
         const hier::Box<DIM> dst_side =
            SideGeometry<DIM>::toSideBox(dst_ghost, d);
         const hier::Box<DIM> src_side =
            SideGeometry<DIM>::toSideBox(src_shift, d);

         const hier::Box<DIM> together = dst_side * src_side;

         if (!together.empty()) {

            // Add lower side intersection (if any) to the box list
            hier::Box<DIM> low_side = src_side;
            low_side.upper(d) = low_side.lower(d); //+ghosts;
            dst_boxes[d].unionBoxes(low_side * msk_side * dst_side);

            // Add upper side intersection (if any) to the box list
            hier::Box<DIM> hig_side = src_side;
            hig_side.lower(d) = hig_side.upper(d); //-ghosts;
            dst_boxes[d].unionBoxes(hig_side * msk_side * dst_side);

            // Take away the interior of over_write interior is not set
            if (!overwrite_interior) {
               dst_boxes[d].removeIntersections(
                  SideGeometry<DIM>::toSideBox(dst_geometry.getBox(), d));
            }

         }  // if (!together.empty())

      }  // loop over dim

   } // if (!quick_check.empty())

   // Create the side overlap data object using the boxes and source shift

   hier::BoxOverlap<DIM> *overlap = new SideOverlap<DIM>(dst_boxes, src_offset);
   return(tbox::Pointer< hier::BoxOverlap<DIM> >(overlap));
}

}
}
#endif
