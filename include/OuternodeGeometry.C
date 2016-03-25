//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/boxgeometry/OuternodeGeometry.C $
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	hier::Box geometry information for outernode centered objects
//

#ifndef included_pdat_OuternodeGeometry_C
#define included_pdat_OuternodeGeometry_C

#include "OuternodeGeometry.h"
#include "BoxList.h"
#include "NodeGeometry.h"
#include "NodeOverlap.h"
#include "tbox/Utilities.h"


#ifdef DEBUG_NO_INLINE
#include "OuternodeGeometry.I"
#endif
namespace SAMRAI {
    namespace pdat {

/*
*************************************************************************

Create a outernode geometry object given the box and ghost cell width.     

*************************************************************************
*/

template<int DIM>  OuternodeGeometry<DIM>::OuternodeGeometry(
   const hier::Box<DIM>& box,
   const hier::IntVector<DIM>& ghosts)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(ghosts.min() >= 0);
#endif
   d_box    = box;
   d_ghosts = ghosts;
}

template<int DIM>  OuternodeGeometry<DIM>::~OuternodeGeometry()
{
}

/*
*************************************************************************

Attempt to calculate the intersection between two outernode centered  
box geometries.  The calculateOverlap() checks whether both arguments 
are outernode geometries; if so, it computes the intersection.  If    
not, then it calls calculateOverlap() on the source object (if retry  
is true) to allow the source a chance to calculate the intersection.  
See the hier::BoxGeometry<DIM> base class for more information about the   
protocol.  A pointer to null is returned if the intersection canot be 
computed.                                                             

*************************************************************************
*/

template<int DIM> tbox::Pointer< hier::BoxOverlap<DIM> > 
OuternodeGeometry<DIM>::calculateOverlap(
   const hier::BoxGeometry<DIM>& dst_geometry,
   const hier::BoxGeometry<DIM>& src_geometry,
   const hier::Box<DIM>& src_mask,
   const bool overwrite_interior,
   const hier::IntVector<DIM>& src_offset,
   const bool retry) const
{
   const NodeGeometry<DIM> *t_dst_node = 
      dynamic_cast<const NodeGeometry<DIM> *>(&dst_geometry);
   const OuternodeGeometry<DIM> *t_dst_onode = 
      dynamic_cast<const OuternodeGeometry<DIM> *>(&dst_geometry);
   const NodeGeometry<DIM> *t_src_node =
      dynamic_cast<const NodeGeometry<DIM> *>(&src_geometry);
   const OuternodeGeometry<DIM> *t_src_onode =
      dynamic_cast<const OuternodeGeometry<DIM> *>(&src_geometry);

   tbox::Pointer< hier::BoxOverlap<DIM> > over = NULL;
   if ((t_src_onode != NULL) && (t_dst_node != NULL)) {
      over = doOverlap(*t_dst_node, *t_src_onode, src_mask, overwrite_interior, 
		       src_offset);
   } else if ((t_dst_onode != NULL) && (t_src_node != NULL)) {
      over = doOverlap(*t_dst_onode, *t_src_node, src_mask, overwrite_interior, 
		       src_offset);
   } else if ((t_src_onode != NULL) && (t_dst_onode != NULL)) {
      over = doOverlap(*t_dst_onode, *t_src_onode, src_mask, overwrite_interior,
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
*
* Compute the overlap between a node geometry destination box and an    
* outernode geometry source box.  The intersection algorithm is similar 
* the node geometry algorithm except that only the borders of source    
* are used in the intersection computation.                             
*
*************************************************************************
*/

template<int DIM> tbox::Pointer< hier::BoxOverlap<DIM> > 
OuternodeGeometry<DIM>::doOverlap(
   const NodeGeometry<DIM>& dst_geometry,
   const OuternodeGeometry<DIM>& src_geometry,
   const hier::Box<DIM>& src_mask,
   const bool overwrite_interior,
   const hier::IntVector<DIM>& src_offset)
{

   hier::BoxList<DIM> dst_boxes;

   // Perform a quick-and-dirty intersection to see if the boxes might overlap

   const hier::Box<DIM> src_box =
      hier::Box<DIM>::grow(src_geometry.d_box, src_geometry.d_ghosts) * src_mask;
   const hier::Box<DIM> src_box_shifted = hier::Box<DIM>::shift(src_box, src_offset);
   const hier::Box<DIM> dst_box =
      hier::Box<DIM>::grow(dst_geometry.getBox(), dst_geometry.getGhosts());

   const hier::Box<DIM> dst_node_box = NodeGeometry<DIM>::toNodeBox(dst_box);
   const hier::Box<DIM> src_node_box = NodeGeometry<DIM>::toNodeBox(src_box_shifted);

   // Compute the intersection (if any) for each of the side directions

   if ( dst_node_box.intersects(src_node_box) ) {

      const hier::Box<DIM> msk_node_box = 
         NodeGeometry<DIM>::toNodeBox( hier::Box<DIM>::shift(src_mask, src_offset) );

      for (int d = 0; d < DIM; d++) {

	 hier::Box<DIM> trimmed_src_node_box = src_node_box;
	 for ( int dh=d+1; dh< DIM; ++dh ) {
	    /*
	      For dimensions higher than d, narrow the box down to avoid
	      representing edge and corner nodes multiple times.
            */
	    ++trimmed_src_node_box.lower(dh);
	    --trimmed_src_node_box.upper(dh);
	 }

	 // Add lower side intersection (if any) to the box list
	 hier::Box<DIM> low_node_box = trimmed_src_node_box;
	 low_node_box.upper(d) = low_node_box.lower(d);
	 dst_boxes.unionBoxes(low_node_box * msk_node_box * dst_node_box);

	 // Add upper side intersection (if any) to the box list
	 hier::Box<DIM> hig_node_box = trimmed_src_node_box;
	 hig_node_box.lower(d) = hig_node_box.upper(d);
	 dst_boxes.unionBoxes(hig_node_box * msk_node_box * dst_node_box);

	 // Take away the interior if over_write interior is not set

	 if (!overwrite_interior) {
	    dst_boxes.removeIntersections(
	       NodeGeometry<DIM>::toNodeBox(dst_geometry.getBox()));
	 }

      }  // loop over dim

   }  // src and dst boxes intersect

   // Create the outernode overlap data object using the boxes and source shift

   hier::BoxOverlap<DIM> *overlap = new NodeOverlap<DIM>(dst_boxes, src_offset);
   return(tbox::Pointer< hier::BoxOverlap<DIM> >(overlap));
}

/*
*************************************************************************
*
* Compute the overlap between an outernode geometry destination box and a    
* node geometry source box.  The intersection algorithm is similar 
* the node geometry algorithm except that only the borders of the dest    
* are used in the intersection computation.                             
*
*************************************************************************
*/

template<int DIM> tbox::Pointer< hier::BoxOverlap<DIM> > 
OuternodeGeometry<DIM>::doOverlap(
   const OuternodeGeometry<DIM>& dst_geometry,
   const NodeGeometry<DIM>& src_geometry,
   const hier::Box<DIM>& src_mask,
   const bool overwrite_interior,
   const hier::IntVector<DIM>& src_offset)
{

   hier::BoxList<DIM> src_boxes;

   // Perform a quick-and-dirty intersection to see if the boxes might overlap

   const hier::Box<DIM> src_box =
      hier::Box<DIM>::grow(src_geometry.getBox(), src_geometry.getGhosts()) * src_mask;
   const hier::Box<DIM> src_box_shifted = hier::Box<DIM>::shift(src_box, src_offset);
   const hier::Box<DIM> dst_box =
      hier::Box<DIM>::grow(dst_geometry.d_box, dst_geometry.d_ghosts);

   const hier::Box<DIM> dst_node_box = NodeGeometry<DIM>::toNodeBox(dst_box);
   const hier::Box<DIM> src_node_box = NodeGeometry<DIM>::toNodeBox(src_box_shifted);

   // Compute the intersection (if any) for each of the side directions

   if ( dst_node_box.intersects(src_node_box) ) { 

      const hier::Box<DIM> msk_node_box = 
         NodeGeometry<DIM>::toNodeBox( hier::Box<DIM>::shift(src_mask, src_offset) );

      for (int d = 0; d < DIM; d++) {

	 hier::Box<DIM> trimmed_dst_node_box = dst_node_box;
	 for ( int dh=d+1; dh < DIM; ++dh ) {
	    /*
	      For dimensions higher than d, narrow the box down to avoid
	      representing edge and corner nodes multiple times.
            */
	    ++trimmed_dst_node_box.lower(dh);
	    --trimmed_dst_node_box.upper(dh);
	 }

	 // Add lower side intersection (if any) to the box list
	 hier::Box<DIM> low_node_box = trimmed_dst_node_box;
	 low_node_box.upper(d) = low_node_box.lower(d);
	 src_boxes.unionBoxes(low_node_box * msk_node_box * src_node_box);

	 // Add upper side intersection (if any) to the box list
	 hier::Box<DIM> hig_node_box = trimmed_dst_node_box;
	 hig_node_box.lower(d) = hig_node_box.upper(d);
	 src_boxes.unionBoxes(hig_node_box * msk_node_box * src_node_box);

	 // Take away the interior of over_write interior is not set

	 if (!overwrite_interior) {
	    src_boxes.removeIntersections(
	       NodeGeometry<DIM>::toNodeBox(dst_geometry.getBox()));
	 }

      }  // loop over dim

   }  // src and dst boxes intersect

   // Create the side overlap data object using the boxes and source shift

   hier::BoxOverlap<DIM> *overlap = new NodeOverlap<DIM>(src_boxes, src_offset);
   return(tbox::Pointer< hier::BoxOverlap<DIM> >(overlap));
}


/*
*************************************************************************
*
* Compute the overlap between an outernode geometry destination box and an    
* outernode geometry source box.  The intersection algorithm is similar 
* the node geometry algorithm except that only the borders of source    
* are used in the intersection computation.                             
*
*************************************************************************
*/

template<int DIM> tbox::Pointer< hier::BoxOverlap<DIM> > 
OuternodeGeometry<DIM>::doOverlap(
   const OuternodeGeometry<DIM>& dst_geometry,
   const OuternodeGeometry<DIM>& src_geometry,
   const hier::Box<DIM>& src_mask,
   const bool overwrite_interior,
   const hier::IntVector<DIM>& src_offset)
{

   hier::BoxList<DIM> dst_boxes;

   // Perform a quick-and-dirty intersection to see if the boxes might overlap

   const hier::Box<DIM> src_box =
      hier::Box<DIM>::grow(src_geometry.d_box, src_geometry.d_ghosts) * src_mask;
   const hier::Box<DIM> src_box_shifted = hier::Box<DIM>::shift(src_box, src_offset);
   const hier::Box<DIM> dst_box =
      hier::Box<DIM>::grow(dst_geometry.getBox(), dst_geometry.getGhosts());

   const hier::Box<DIM> dst_node_box = NodeGeometry<DIM>::toNodeBox(dst_box);
   const hier::Box<DIM> src_node_box = NodeGeometry<DIM>::toNodeBox(src_box_shifted);

   // Compute the intersection (if any) for each of the side directions

   if ( dst_node_box.intersects(src_node_box) ) { 

      const hier::Box<DIM> msk_node_box =
         NodeGeometry<DIM>::toNodeBox( hier::Box<DIM>::shift(src_mask, src_offset) );

      int dst_d, src_d;

      for ( dst_d=0; dst_d<DIM; ++dst_d ) {

	 hier::Box<DIM> trimmed_dst_node_box = dst_node_box;
	 for ( int dh=dst_d+1; dh<DIM; ++dh ) {
	    ++trimmed_dst_node_box.lower(dh);
	    --trimmed_dst_node_box.upper(dh);
	 }

	 hier::Box<DIM> lo_dst_node_box = trimmed_dst_node_box;
	 lo_dst_node_box.upper(dst_d) = lo_dst_node_box.lower(dst_d);

	 hier::Box<DIM> hi_dst_node_box = trimmed_dst_node_box;
	 hi_dst_node_box.lower(dst_d) = hi_dst_node_box.upper(dst_d);


	 for (src_d = 0; src_d < DIM; ++src_d) {

	    hier::Box<DIM> trimmed_src_node_box = src_node_box;
	    for ( int dh=src_d+1; dh<DIM; ++dh ) {
	       ++trimmed_src_node_box.lower(dh);
	       --trimmed_src_node_box.upper(dh);
	    }

	    hier::Box<DIM> lo_src_node_box = trimmed_src_node_box;
	    lo_src_node_box.upper(src_d) = lo_src_node_box.lower(src_d);

	    hier::Box<DIM> hi_src_node_box = trimmed_src_node_box;
	    hi_src_node_box.lower(src_d) = hi_src_node_box.upper(src_d);

	    dst_boxes.unionBoxes(lo_src_node_box*msk_node_box*lo_dst_node_box);
	    dst_boxes.unionBoxes(hi_src_node_box*msk_node_box*lo_dst_node_box);
	    dst_boxes.unionBoxes(lo_src_node_box*msk_node_box*hi_dst_node_box);
	    dst_boxes.unionBoxes(hi_src_node_box*msk_node_box*hi_dst_node_box);

	    // Take away the interior of over_write interior is not set

	    if (!overwrite_interior) {
	       dst_boxes.removeIntersections(
	       NodeGeometry<DIM>::toNodeBox(dst_geometry.d_box));
	    }

	 }  // loop over src dim

      }  // loop over dst dim

   }  // if src and dst boxes intersect

   // Create the side overlap data object using the boxes and source shift

   hier::BoxOverlap<DIM> *overlap = new NodeOverlap<DIM>(dst_boxes, src_offset);
   return(tbox::Pointer< hier::BoxOverlap<DIM> >(overlap));

}

}
}
#endif
