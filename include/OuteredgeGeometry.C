//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/boxgeometry/OuteredgeGeometry.C $
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Release:	$Name$
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Box geometry information for outeredge centered objects
//

#ifndef included_pdat_OuteredgeGeometry_C
#define included_pdat_OuteredgeGeometry_C

#include "OuteredgeGeometry.h"

#include "BoxList.h"
#include "EdgeGeometry.h"
#include "EdgeOverlap.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#include "tbox/Utilities.h"
#endif

#ifdef DEBUG_NO_INLINE
#include "OuteredgeGeometry.I"
#endif

namespace SAMRAI {
    namespace pdat {

/*
*************************************************************************
*									*
* Create an outeredge geometry object given box and ghost cell width.	*
*									*
*************************************************************************
*/

template<int DIM> OuteredgeGeometry<DIM>::OuteredgeGeometry(
   const hier::Box<DIM>& box, 
   const hier::IntVector<DIM>& ghosts)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(ghosts.min() >= 0);
#endif
   d_box    = box;
   d_ghosts = ghosts;
}

template<int DIM> OuteredgeGeometry<DIM>::~OuteredgeGeometry()
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
* hier::BoxGeometry base class for more information about the protocol.	*
* A pointer to null is returned if the intersection cannot be computed.	*
* 									*
*************************************************************************
*/

template<int DIM> tbox::Pointer<hier::BoxOverlap<DIM> > 
OuteredgeGeometry<DIM>::calculateOverlap(
   const hier::BoxGeometry<DIM>& dst_geometry,
   const hier::BoxGeometry<DIM>& src_geometry,
   const hier::Box<DIM>& src_mask,
   const bool overwrite_interior,
   const hier::IntVector<DIM>& src_offset,
   const bool retry) const
{
   const pdat::EdgeGeometry<DIM> *t_dst_edge = 
      dynamic_cast<const pdat::EdgeGeometry<DIM> *>(&dst_geometry);
   const OuteredgeGeometry<DIM> *t_dst_oedge =
      dynamic_cast<const OuteredgeGeometry<DIM> *>(&dst_geometry);
   const OuteredgeGeometry<DIM> *t_src =
      dynamic_cast<const OuteredgeGeometry<DIM> *>(&src_geometry);

   tbox::Pointer<hier::BoxOverlap<DIM> > over = NULL;

   if ((t_src != NULL) && (t_dst_edge != NULL)) {
      over = doOverlap(*t_dst_edge, *t_src, src_mask, overwrite_interior, 
		       src_offset);
   } else if ((t_src != NULL) && (t_dst_oedge != NULL)) {
      over = doOverlap(*t_dst_oedge, *t_src, src_mask, overwrite_interior, 
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
* Compute the overlap between an edge and an outeredge centered boxes.  *
* The algorithm is similar to the standard edge intersection algorithm  *
* except we operate only on the boundaries of the source box.           *
*									*
*************************************************************************
*/
template<int DIM> tbox::Pointer<hier::BoxOverlap<DIM> > 
OuteredgeGeometry<DIM>::doOverlap(
   const pdat::EdgeGeometry<DIM>& dst_geometry,
   const OuteredgeGeometry<DIM>& src_geometry,
   const hier::Box<DIM>& src_mask,
   const bool overwrite_interior,
   const hier::IntVector<DIM>& src_offset)
{
   hier::BoxList<DIM> dst_boxes[DIM];

   // Perform a quick-and-dirty intersection to see if the boxes might overlap

   const hier::Box<DIM> src_box =
      hier::Box<DIM>::grow(src_geometry.d_box, src_geometry.d_ghosts) * src_mask;
   const hier::Box<DIM> src_box_shifted = hier::Box<DIM>::shift(src_box, src_offset);
   const hier::Box<DIM> dst_box =
      hier::Box<DIM>::grow(dst_geometry.getBox(), dst_geometry.getGhosts());

   // Compute the intersection (if any) for each of the edge directions

   bool quick_boxes_intersect = 
      ( hier::Box<DIM>::grow(src_box_shifted, 1) ).intersects(
                                                   hier::Box<DIM>::grow(dst_box, 1) );
   if (quick_boxes_intersect) {

      for (int axis = 0; axis < DIM; ++axis) {

         const hier::Box<DIM> dst_edge_box = 
            pdat::EdgeGeometry<DIM>::toEdgeBox(dst_box, axis);
         const hier::Box<DIM> src_edge_box = 
            pdat::EdgeGeometry<DIM>::toEdgeBox(src_box_shifted, axis);

         bool boxes_intersect = dst_edge_box.intersects(src_edge_box);

         if (boxes_intersect) {

            for (int face_normal = 0; face_normal < DIM; ++face_normal) {

               if ( face_normal != axis ) {
 
                  for (int side = 0; side < 2; ++side) {
                     hier::Box<DIM> outeredge_src_box = toOuteredgeBox(src_box_shifted, 
                                                                       axis, 
                                                                       face_normal,
                                                                       side);
                     dst_boxes[axis].unionBoxes(outeredge_src_box * dst_edge_box);
                  }

               }  // data is not defined when face_normal == axis

            }  // iterate over face normal directions
  
            if (!overwrite_interior) {
               const hier::Box<DIM> interior_edges = 
                  pdat::EdgeGeometry<DIM>::toEdgeBox(dst_geometry.getBox(), 
                                                     axis);
               dst_boxes[axis].removeIntersections(interior_edges);
            }

         }  // if source and destination edge boxes overlap in axis direction

      }  // iterate over axis directions

   }  // if quick check passes

   // Create the edge overlap data object using the boxes and source shift
   hier::BoxOverlap<DIM> *overlap = 
      new pdat::EdgeOverlap<DIM>(dst_boxes, src_offset);
   return(tbox::Pointer<hier::BoxOverlap<DIM> >(overlap));
}


/*
*************************************************************************
*									*
* Compute the overlap between two outeredge centered boxes.             *
* The algorithm is similar to the standard edge intersection algorithm  *
* except we operate only on the boundaries of the source box.           *
*									*
*************************************************************************
*/

template<int DIM> tbox::Pointer<hier::BoxOverlap<DIM> > 
OuteredgeGeometry<DIM>::doOverlap(
   const OuteredgeGeometry<DIM>& dst_geometry,
   const OuteredgeGeometry<DIM>& src_geometry,
   const hier::Box<DIM>& src_mask,
   const bool overwrite_interior,
   const hier::IntVector<DIM>& src_offset)
{

   hier::BoxList<DIM> dst_boxes[DIM];

   // Perform a quick-and-dirty intersection to see if the boxes might overlap

   const hier::Box<DIM> src_box =
      hier::Box<DIM>::grow(src_geometry.d_box, src_geometry.d_ghosts) * src_mask;
   const hier::Box<DIM> src_box_shifted = hier::Box<DIM>::shift(src_box, src_offset);
   const hier::Box<DIM> dst_box =
      hier::Box<DIM>::grow(dst_geometry.getBox(), dst_geometry.getGhosts());

   // Compute the intersection (if any) for each of the edge directions

   bool quick_boxes_intersect =
      ( hier::Box<DIM>::grow(src_box_shifted, 1) ).intersects(
                                                   hier::Box<DIM>::grow(dst_box, 1) );
   if (quick_boxes_intersect) {

      for (int axis = 0; axis < DIM; ++axis) {

         const hier::Box<DIM> dst_edge_box = 
            pdat::EdgeGeometry<DIM>::toEdgeBox(dst_box, axis);
         const hier::Box<DIM> src_edge_box = 
            pdat::EdgeGeometry<DIM>::toEdgeBox(src_box_shifted, axis);

         bool boxes_intersect = dst_edge_box.intersects(src_edge_box);

         if (boxes_intersect) {

            for (int src_face_normal = 0; src_face_normal < DIM; ++src_face_normal) {

               if ( src_face_normal != axis ) {

                  hier::Box<DIM> outeredge_src_box_lo = toOuteredgeBox(src_box_shifted,
                                                                       axis,
                                                                       src_face_normal,
                                                                       0);
                  hier::Box<DIM> outeredge_src_box_up = toOuteredgeBox(src_box_shifted,
                                                                       axis,
                                                                       src_face_normal,
                                                                       1);

                  for (int dst_face_normal = 0; dst_face_normal < DIM; ++dst_face_normal) {

                     if ( dst_face_normal != axis ) {
 
                        hier::Box<DIM> outeredge_dst_box_lo = toOuteredgeBox(dst_box,
                                                                             axis,
                                                                             dst_face_normal,
                                                                             0);
                        hier::Box<DIM> outeredge_dst_box_up = toOuteredgeBox(dst_box,
                                                                             axis,
                                                                             dst_face_normal,
                                                                             1);

                        dst_boxes[axis].unionBoxes(outeredge_src_box_lo * outeredge_dst_box_lo);
                        dst_boxes[axis].unionBoxes(outeredge_src_box_lo * outeredge_dst_box_up);
                        dst_boxes[axis].unionBoxes(outeredge_src_box_up * outeredge_dst_box_lo);
                        dst_boxes[axis].unionBoxes(outeredge_src_box_up * outeredge_dst_box_up);

                     }  // dst data undefined when dst_face_normal == axis

                  }  // iterate over dst face normal directions

               }  // src data undefined when src_face_normal == axis

            }  // iterate over src face normal directions

         }  // if source and destination edge boxes overlap in axis direction

         if (!overwrite_interior) {
            const hier::Box<DIM> interior_edges = 
               pdat::EdgeGeometry<DIM>::toEdgeBox(dst_geometry.getBox(), 
                                                  axis);
            dst_boxes[axis].removeIntersections(interior_edges);
         }

      }  // iterate over axis directions

   }  // if quick check passes

   // Create the edge overlap data object using the boxes and source shift
   hier::BoxOverlap<DIM> *overlap = 
      new pdat::EdgeOverlap<DIM>(dst_boxes, src_offset);
   return(tbox::Pointer<hier::BoxOverlap<DIM> >(overlap));
}


/*
*************************************************************************
*                                                                       *
* Convert an AMR-index space hier::Box into a edge-index space box      *
* for an outeredge region.                                              *
*                                                                       *
*************************************************************************
*/

template<int DIM>
hier::Box<DIM> OuteredgeGeometry<DIM>::toOuteredgeBox(
   const hier::Box<DIM>& box,
   int axis,
   int face_normal, 
   int side)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(0 <= axis && axis < DIM); 
   TBOX_ASSERT(0 <= face_normal && face_normal < DIM); 
   TBOX_ASSERT(face_normal != axis);
   TBOX_ASSERT(side == 0 || side == 1);
#endif

   hier::Box<DIM> oedge_box;

   /*
    * If data is defined (i.e., face_normal != axis), then
    *    1) Make an edge box for the given axis.
    *    2) Trim box as needed to avoid redundant edge indices
    *       for different face normal directions. 
    *    3) Restrict box to lower or upper face for given
    *       face normal direction.
    */

   if ( (face_normal != axis) && !box.empty() ) {

      oedge_box = EdgeGeometry<DIM>::toEdgeBox(box, axis);

      for (int d = 0; d < DIM; ++d) {

         if ( d != axis ) {  // do not trim in axis direction

            for (int dh = d+1; dh < DIM; ++dh) { // trim in higher dimensions 

               if ( dh != axis && dh != face_normal ) {  
                  // do not trim in axis or face_normal direction 

                  ++oedge_box.lower(dh);
                  --oedge_box.upper(dh);

               }

            }
 
         }

      }

      if (side == 0 ) {  // lower side in face normal direction
         oedge_box.upper(face_normal) = oedge_box.lower(face_normal);
      } else {  // side == 1; upper side in face normal direction
         oedge_box.lower(face_normal) = oedge_box.upper(face_normal);
      }

   }

   return(oedge_box);
}


}
}
#endif

