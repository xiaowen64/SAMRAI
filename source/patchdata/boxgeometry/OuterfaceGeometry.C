//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/boxgeometry/OuterfaceGeometry.C $
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	hier::Box geometry information for outerface centered objects
//

#ifndef included_pdat_OuterfaceGeometry_C
#define included_pdat_OuterfaceGeometry_C

#include "OuterfaceGeometry.h"
#include "BoxList.h"
#include "FaceGeometry.h"
#include "FaceOverlap.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#include "tbox/Utilities.h"
#endif

#ifdef DEBUG_NO_INLINE
#include "OuterfaceGeometry.I"
#endif
namespace SAMRAI {
    namespace pdat {

/*
*************************************************************************
*									*
* Create a face geometry object given the box and ghost cell width.	*
*									*
*************************************************************************
*/

template<int DIM>  OuterfaceGeometry<DIM>::OuterfaceGeometry(
   const hier::Box<DIM>& box, const hier::IntVector<DIM>& ghosts)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(ghosts.min() >= 0);
#endif
   d_box    = box;
   d_ghosts = ghosts;
}

template<int DIM>  OuterfaceGeometry<DIM>::~OuterfaceGeometry()
{
}

/*
*************************************************************************
*									*
* Attempt to calculate the intersection between two outerface centered	*
* box geometries.  The calculateOverlap() checks whether both arguments	*
* are outerface geometries; if so, it compuates the intersection.  If	*
* not, then it calls calculateOverlap() on the source object (if retry	*
* is true) to allow the source a chance to calculate the intersection.	*
* See the hier::BoxGeometry<DIM> base class for more information about  *
* the protocol.  A pointer to null is returned if the intersection      *
* cannot be computed.                                                   *
* 									*
*************************************************************************
*/

template<int DIM> 
tbox::Pointer< hier::BoxOverlap<DIM> > OuterfaceGeometry<DIM>::calculateOverlap(
   const hier::BoxGeometry<DIM>& dst_geometry,
   const hier::BoxGeometry<DIM>& src_geometry,
   const hier::Box<DIM>& src_mask,
   const bool overwrite_interior,
   const hier::IntVector<DIM>& src_offset,
   const bool retry) const
{
   const FaceGeometry<DIM> *t_dst = 
      dynamic_cast<const FaceGeometry<DIM> *>(&dst_geometry);
   const OuterfaceGeometry<DIM> *t_src =
      dynamic_cast<const OuterfaceGeometry<DIM> *>(&src_geometry);

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
* Compute the overlap between a face geometry destination box and an	*
* outerface geometry source box.  The intersection algorithm is similar	*
* the face geometry algorithm except that only the borders of source	*
* are used in the intersection computation.				*
*									*
*************************************************************************
*/

template<int DIM> 
tbox::Pointer< hier::BoxOverlap<DIM> > OuterfaceGeometry<DIM>::doOverlap(
   const FaceGeometry<DIM>& dst_geometry,
   const OuterfaceGeometry<DIM>& src_geometry,
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
      hier::Box<DIM>::grow(dst_geometry.getBox(), dst_geometry.getGhosts());

   // Compute the intersection (if any) for each of the face directions

   const hier::Box<DIM> quick_check =
      hier::Box<DIM>::grow(src_shift, 1) * hier::Box<DIM>::grow(dst_ghost, 1);

   if (!quick_check.empty()) {

      const hier::Box<DIM> mask_shift = hier::Box<DIM>::shift(src_mask, src_offset);

      for (int d = 0; d < DIM; d++) {

         const hier::Box<DIM> msk_face =
            FaceGeometry<DIM>::toFaceBox(mask_shift, d);
         const hier::Box<DIM> dst_face =
            FaceGeometry<DIM>::toFaceBox(dst_ghost, d);
         const hier::Box<DIM> src_face =
            FaceGeometry<DIM>::toFaceBox(src_shift, d);

         const hier::Box<DIM> together = dst_face * src_face;

         if (!together.empty()) {

            // Add lower face intersection (if any) to the box list
            hier::Box<DIM> low_face = src_face;
            low_face.upper(0) = low_face.lower(0);  //+ghosts;
            dst_boxes[d].unionBoxes(low_face * msk_face * dst_face);

            // Add upper face intersection (if any) to the box list
            hier::Box<DIM> hig_face = src_face;
            hig_face.lower(0) = hig_face.upper(0);  //-ghosts;
            dst_boxes[d].unionBoxes(hig_face * msk_face * dst_face);

            // Take away the interior of over_write interior is not set
            if (!overwrite_interior) {
               dst_boxes[d].removeIntersections(
                  FaceGeometry<DIM>::toFaceBox(dst_geometry.getBox(), d));
            }

         }  // if (!together.empty())

      }  // loop over dim

   } // if (!quick_check.empty())


   // Create the face overlap data object using the boxes and source shift

   hier::BoxOverlap<DIM> *overlap = new FaceOverlap<DIM>(dst_boxes, src_offset);
   return(tbox::Pointer< hier::BoxOverlap<DIM> >(overlap));
}

}
}
#endif
