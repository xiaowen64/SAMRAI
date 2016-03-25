//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/boxgeometry/SideGeometry.C $
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	hier::Box geometry information for side centered objects
//

#ifndef included_pdat_SideGeometry_C
#define included_pdat_SideGeometry_C

#include "SideGeometry.h"
#include "BoxList.h"
#include "SideOverlap.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#include "tbox/Utilities.h"
#endif

#ifdef DEBUG_NO_INLINE
#include "SideGeometry.I"
#endif
namespace SAMRAI {
    namespace pdat {

/*
*************************************************************************
*									*
* Create a side geometry object given the box, ghost cell width, and    *
* direction information.                                                *
*									*
*************************************************************************
*/

template<int DIM>  SideGeometry<DIM>::SideGeometry(
   const hier::Box<DIM>& box,
   const hier::IntVector<DIM>& ghosts,
   const hier::IntVector<DIM>& directions)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(ghosts.min() >= 0);
   TBOX_ASSERT(directions.min() >= 0);
#endif
   d_box    = box;
   d_ghosts = ghosts;
   d_directions = directions;
}

template<int DIM>  SideGeometry<DIM>::~SideGeometry()
{
}

/*
*************************************************************************
*									*
* Attempt to calculate the intersection between two side centered box	*
* geometries.  The calculateOverlap() checks whether both arguments are	*
* side geometries; if so, it compuates the intersection.  If not, then	*
* it calls calculateOverlap() on the source object (if retry is true)	*
* to allow the source a chance to calculate the intersection.  See the	*
* hier::BoxGeometry<DIM> base class for more information about the protocol.	*
* A pointer to null is returned if the intersection cannot be computed.	*
* 									*
*************************************************************************
*/

template<int DIM> tbox::Pointer< hier::BoxOverlap<DIM> > SideGeometry<DIM>::calculateOverlap(
   const hier::BoxGeometry<DIM>& dst_geometry,
   const hier::BoxGeometry<DIM>& src_geometry,
   const hier::Box<DIM>& src_mask,
   const bool overwrite_interior,
   const hier::IntVector<DIM>& src_offset,
   const bool retry) const
{
   const SideGeometry<DIM> *t_dst =
      dynamic_cast<const SideGeometry<DIM> *>(&dst_geometry);
   const SideGeometry<DIM> *t_src =
      dynamic_cast<const SideGeometry<DIM> *>(&src_geometry);

   tbox::Pointer< hier::BoxOverlap<DIM> > over = NULL;
   if ((t_src != NULL) && (t_dst != NULL)) {
      over = doOverlap(*t_dst, *t_src, src_mask, overwrite_interior, 
		       src_offset);
   } else if (retry) {
      over = src_geometry.calculateOverlap(dst_geometry, src_geometry,
                                           src_mask, overwrite_interior,
                                           src_offset, false);
   }
   return(over);
}

/*
*************************************************************************
*									*
* Convert an AMR-index space hier::Box into a side-index space box by a	*
* increasing the index size by one in the axis direction.		*
*									*
*************************************************************************
*/

template<int DIM> hier::Box<DIM> SideGeometry<DIM>::toSideBox(
   const hier::Box<DIM>& box,
   int side_normal)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( (side_normal >= 0) && (side_normal < DIM) );
#endif
   hier::Box<DIM> side_box;

   if (!box.empty()) {
      side_box = box;
      side_box.upper(side_normal) += 1;
   }

   return(side_box);
}

/*
*************************************************************************
*									*
* Compute the overlap between two side centered boxes.  The algorithm	*
* is fairly straight-forward.  First, we perform a quick-and-dirty	*
* intersection to see if the boxes might overlap.  If that intersection	*
* is not empty, then we need to do a better job calculating the overlap	*
* for each dimension.  Note that the AMR index space boxes must be	*
* shifted into the side centered space before we calculate the proper	*
* intersections.							*
*									*
*************************************************************************
*/

template<int DIM> tbox::Pointer< hier::BoxOverlap<DIM> > SideGeometry<DIM>::doOverlap(
   const SideGeometry<DIM>& dst_geometry,
   const SideGeometry<DIM>& src_geometry,
   const hier::Box<DIM>& src_mask,
   const bool overwrite_interior,
   const hier::IntVector<DIM>& src_offset)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(dst_geometry.getDirectionVector() 
           == src_geometry.getDirectionVector());
#endif

   hier::BoxList<DIM> dst_boxes[DIM];

   // Perform a quick-and-dirty intersection to see if the boxes might overlap

   const hier::Box<DIM> src_box =
      hier::Box<DIM>::grow(src_geometry.d_box, src_geometry.d_ghosts) * src_mask;
   const hier::Box<DIM> src_shift =
      hier::Box<DIM>::shift(src_box, src_offset);
   const hier::Box<DIM> dst_ghost =
      hier::Box<DIM>::grow(dst_geometry.d_box, dst_geometry.d_ghosts);

   // Compute the intersection (if any) for each of the side directions

   const hier::Box<DIM> quick_check =
      hier::Box<DIM>::grow(src_shift, 1) * hier::Box<DIM>::grow(dst_ghost, 1);

   if (!quick_check.empty()) {

      const hier::IntVector<DIM>& dirs = src_geometry.getDirectionVector();
      for (int d = 0; d < DIM && dirs(d); d++) {

         const hier::Box<DIM> dst_side = toSideBox(dst_ghost, d);
         const hier::Box<DIM> src_side = toSideBox(src_shift, d);
         const hier::Box<DIM> together = dst_side * src_side;

         if (!together.empty()) {

            dst_boxes[d].unionBoxes(together);
            if (!overwrite_interior) {
               const hier::Box<DIM> int_side = toSideBox(dst_geometry.d_box, d);
               dst_boxes[d].removeIntersections(together,int_side);
            } else {
               dst_boxes[d].appendItem(together);
            }

         }  // if (!together.empty())

      }  // loop over dim && dirs(d)

   }  // if (!quick_check.empty())

   // Create the side overlap data object using the boxes and source shift

   hier::BoxOverlap<DIM> *overlap = new SideOverlap<DIM>(dst_boxes, src_offset);
   return(tbox::Pointer< hier::BoxOverlap<DIM> >(overlap));
}

}
}
#endif
