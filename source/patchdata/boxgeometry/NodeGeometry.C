//
// File:	NodeGeometry.C
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	hier::Box geometry information for node centered objects
//

#ifndef included_pdat_NodeGeometry_C
#define included_pdat_NodeGeometry_C

#include "NodeGeometry.h"
#include "BoxList.h"
#include "NodeOverlap.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

#ifdef DEBUG_NO_INLINE
#include "NodeGeometry.I"
#endif
namespace SAMRAI {
    namespace pdat {

/*
*************************************************************************
*									*
* Create a node geometry object given the box and ghost cell width.	*
*									*
*************************************************************************
*/

template<int DIM>  NodeGeometry<DIM>::NodeGeometry(const hier::Box<DIM>& box,
                                       const hier::IntVector<DIM>& ghosts)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(ghosts.min() >= 0);
#endif
   d_box    = box;
   d_ghosts = ghosts;
}

template<int DIM>  NodeGeometry<DIM>::~NodeGeometry()
{
}

/*
*************************************************************************
*									*
* Attempt to calculate the intersection between two node centered box	*
* geometries.  The calculateOverlap() checks whether both arguments are	*
* node geometries; if so, it computes the intersection.  If not, then	*
* it calls calculateOverlap() on the source object (if retry is true)	*
* to allow the source a chance to calculate the intersection.  See the	*
* hier::BoxGeometry<DIM> base class for more information about the protocol.	*
* A pointer to null is returned if the intersection cannot be computed.	*
* 									*
*************************************************************************
*/

template<int DIM> tbox::Pointer< hier::BoxOverlap<DIM> > NodeGeometry<DIM>::calculateOverlap(
   const hier::BoxGeometry<DIM>& dst_geometry,
   const hier::BoxGeometry<DIM>& src_geometry,
   const hier::Box<DIM>& src_mask,
   const bool overwrite_interior,
   const hier::IntVector<DIM>& src_offset,
   const bool retry) const
{
   const NodeGeometry<DIM> *t_dst = 
      dynamic_cast<const NodeGeometry<DIM> *>(&dst_geometry);
   const NodeGeometry<DIM> *t_src =
      dynamic_cast<const NodeGeometry<DIM> *>(&src_geometry);

   tbox::Pointer< hier::BoxOverlap<DIM> > over= NULL;
   if ( (t_src != NULL) && (t_src != NULL)) {
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
* Compute the overlap between two node centered boxes.  The algorithm	*
* is fairly straight-forward.  First, the two boxes are converted into	*
* node coordinates.  Then, the boxes are intersected and, if necessary,	*
* the interior section is removed from the destination box.		*
*									*
*************************************************************************
*/

template<int DIM> tbox::Pointer< hier::BoxOverlap<DIM> > NodeGeometry<DIM>::doOverlap(
   const NodeGeometry<DIM>& dst_geometry,
   const NodeGeometry<DIM>& src_geometry,
   const hier::Box<DIM>& src_mask,
   const bool overwrite_interior,
   const hier::IntVector<DIM>& src_offset)
{
   hier::BoxList<DIM> dst_boxes;

   // Translate the source box and grow the destination box by the ghost cells

   const hier::Box<DIM> src_box =
      hier::Box<DIM>::grow(src_geometry.d_box, src_geometry.d_ghosts) * src_mask;
   const hier::Box<DIM> src_shift =
      hier::Box<DIM>::shift(src_box, src_offset);
   const hier::Box<DIM> dst_ghost =
      hier::Box<DIM>::grow(dst_geometry.d_box, dst_geometry.d_ghosts);

   // Convert the boxes into node space and compute the intersection

   const hier::Box<DIM> dst_node = NodeGeometry<DIM>::toNodeBox(dst_ghost);
   const hier::Box<DIM> src_node = NodeGeometry<DIM>::toNodeBox(src_shift);
   const hier::Box<DIM> together = dst_node * src_node;

   if (!together.empty()) {
      dst_boxes.unionBoxes(together);
      if (!overwrite_interior) {
         const hier::Box<DIM> int_node = toNodeBox(dst_geometry.d_box);
         dst_boxes.removeIntersections(together,int_node);
      } else {
         dst_boxes.appendItem(together);
      }
   }

   // Create the node overlap data object using the boxes and source shift

   hier::BoxOverlap<DIM> *overlap = new NodeOverlap<DIM>(dst_boxes, src_offset);
   return(tbox::Pointer< hier::BoxOverlap<DIM> >(overlap));
}

}
}
#endif
