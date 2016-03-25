//
// File:	EdgeGeometry.C
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	hier::Box geometry information for edge centered objects
//

#ifndef included_pdat_EdgeGeometry_C
#define included_pdat_EdgeGeometry_C

#include "EdgeGeometry.h"
#include "BoxList.h"
#include "EdgeOverlap.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

#ifdef DEBUG_NO_INLINE
#include "EdgeGeometry.I"
#endif
namespace SAMRAI {
    namespace pdat {

/*
*************************************************************************
*									*
* Create a edge geometry object given the box and ghost cell width.	*
*									*
*************************************************************************
*/

template<int DIM>  EdgeGeometry<DIM>::EdgeGeometry(
   const hier::Box<DIM>& box, const hier::IntVector<DIM>& ghosts)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(ghosts.min() >= 0);
#endif
   d_box    = box;
   d_ghosts = ghosts;
}

template<int DIM>  EdgeGeometry<DIM>::~EdgeGeometry()
{
}

/*
*************************************************************************
*									*
* Attempt to calculate the intersection between two edge centered box	*
* geometries.  The calculateOverlap() checks whether both arguments are	*
* edge geometries; if so, it compuates the intersection.  If not, then	*
* it calls calculateOverlap() on the source object (if retry is true)	*
* to allow the source a chance to calculate the intersection.  See the	*
* hier::BoxGeometry<DIM> base class for more information about the protocol.	*
* A pointer to null is returned if the intersection cannot be computed.	*
* 									*
*************************************************************************
*/

template<int DIM> tbox::Pointer< hier::BoxOverlap<DIM> > EdgeGeometry<DIM>::calculateOverlap(
   const hier::BoxGeometry<DIM>& dst_geometry,
   const hier::BoxGeometry<DIM>& src_geometry,
   const hier::Box<DIM>& src_mask,
   const bool overwrite_interior,
   const hier::IntVector<DIM>& src_offset,
   const bool retry) const
{
   const EdgeGeometry<DIM> *t_dst = 
      dynamic_cast<const EdgeGeometry<DIM> *>(&dst_geometry);
   const EdgeGeometry<DIM> *t_src =
      dynamic_cast<const EdgeGeometry<DIM> *>(&src_geometry);

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
* Convert an AMR-index space hier::Box into a edge-index space box by a	*
* cyclic shift of indices.						*
*									*
*************************************************************************
*/

template<int DIM> hier::Box<DIM> EdgeGeometry<DIM>::toEdgeBox(const hier::Box<DIM>& box, const int axis)
{
   hier::Box<DIM> edge_box;

   if (!box.empty()) {
      edge_box = box;
      for (int i = 0; i < DIM; i++) {
         if (axis != i) {
            edge_box.upper(i) += 1;
         }
      }
   }

   return(edge_box);
}

/*
*************************************************************************
*									*
* Compute the overlap between two edge centered boxes.  The algorithm	*
* is fairly straight-forward.  First, we perform a quick-and-dirty	*
* intersection to see if the boxes might overlap.  If that intersection	*
* is not empty, then we need to do a better job calculating the overlap	*
* for each dimension.  Note that the AMR index space boxes must be	*
* shifted into the edge centered space before we calculate the proper	*
* intersections.							*
*									*
*************************************************************************
*/

template<int DIM> tbox::Pointer< hier::BoxOverlap<DIM> > EdgeGeometry<DIM>::doOverlap(
   const EdgeGeometry<DIM>& dst_geometry,
   const EdgeGeometry<DIM>& src_geometry,
   const hier::Box<DIM>& src_mask,
   const bool overwrite_interior,
   const hier::IntVector<DIM>& src_offset)
{
   hier::BoxList<DIM> dst_boxes[DIM];

   // Perform a quick-and-dirty intersection to see if the boxes might overlap

   const hier::Box<DIM> src_box =
      hier::Box<DIM>::grow(src_geometry.d_box, src_geometry.d_ghosts) * src_mask;
   const hier::Box<DIM> src_shift =
      hier::Box<DIM>::shift(src_box, src_offset);
   const hier::Box<DIM> dst_ghost =
      hier::Box<DIM>::grow(dst_geometry.d_box, dst_geometry.d_ghosts);

   // Compute the intersection (if any) for each of the edge directions

   const hier::Box<DIM> quick_check =
      hier::Box<DIM>::grow(src_shift, 1) * hier::Box<DIM>::grow(dst_ghost, 1);

   if (!quick_check.empty()) {

      for (int d = 0; d < DIM; d++) {

         const hier::Box<DIM> dst_edge = toEdgeBox(dst_ghost, d);
         const hier::Box<DIM> src_edge = toEdgeBox(src_shift, d);
         const hier::Box<DIM> together = dst_edge * src_edge;

         if (!together.empty()) {

            dst_boxes[d].unionBoxes(together);
            if (!overwrite_interior) {
               const hier::Box<DIM> int_edge = toEdgeBox(dst_geometry.d_box, d);
               dst_boxes[d].removeIntersections(together,int_edge);
            } else {
               dst_boxes[d].appendItem(together);
            }

         }  // if (!together.empty())

      }  // loop over dim

   }  // if (!quick_check.empty())

   // Create the edge overlap data object using the boxes and source shift

   hier::BoxOverlap<DIM> *overlap = new EdgeOverlap<DIM>(dst_boxes, src_offset);
   return(tbox::Pointer< hier::BoxOverlap<DIM> >(overlap));
}

}
}
#endif
