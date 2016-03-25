//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/boxgeometry/FaceGeometry.C $
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	hier::Box geometry information for face centered objects
//

#ifndef included_pdat_FaceGeometry_C
#define included_pdat_FaceGeometry_C

#include "FaceGeometry.h"
#include "BoxList.h"
#include "FaceOverlap.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#include "tbox/Utilities.h"
#endif

#ifdef DEBUG_NO_INLINE
#include "FaceGeometry.I"
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

template<int DIM>  FaceGeometry<DIM>::FaceGeometry(
   const hier::Box<DIM>& box, const hier::IntVector<DIM>& ghosts)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(ghosts.min() >= 0);
#endif
   d_box    = box;
   d_ghosts = ghosts;
}

template<int DIM>  FaceGeometry<DIM>::~FaceGeometry()
{
}

/*
*************************************************************************
*									*
* Attempt to calculate the intersection between two face centered box	*
* geometries.  The calculateOverlap() checks whether both arguments are	*
* face geometries; if so, it compuates the intersection.  If not, then	*
* it calls calculateOverlap() on the source object (if retry is true)	*
* to allow the source a chance to calculate the intersection.  See the	*
* hier::BoxGeometry<DIM> base class for more information about the protocol.	*
* A pointer to null is returned if the intersection cannot be computed.	*
* 									*
*************************************************************************
*/

template<int DIM> 
tbox::Pointer< hier::BoxOverlap<DIM> > FaceGeometry<DIM>::calculateOverlap(
   const hier::BoxGeometry<DIM>& dst_geometry,
   const hier::BoxGeometry<DIM>& src_geometry,
   const hier::Box<DIM>& src_mask,
   const bool overwrite_interior,
   const hier::IntVector<DIM>& src_offset,
   const bool retry) const
{
   const FaceGeometry<DIM> *t_dst = 
      dynamic_cast<const FaceGeometry<DIM> *>(&dst_geometry);
   const FaceGeometry<DIM> *t_src =
      dynamic_cast<const FaceGeometry<DIM> *>(&src_geometry);

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
* Convert an AMR-index space hier::Box into a face-index space box by a	*
* cyclic shift of indices.						*
*									*
*************************************************************************
*/

template<int DIM> hier::Box<DIM> 
FaceGeometry<DIM>::toFaceBox(
   const hier::Box<DIM>& box, 
   int face_normal)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( (face_normal >= 0) && (face_normal < DIM) );
#endif

   hier::Box<DIM> face_box;

   if (!box.empty()) {
      const int x = face_normal;
      face_box.lower(0) = box.lower(x);
      face_box.upper(0) = box.upper(x)+1;
      for (int i = 1; i < DIM; i++) {
         const int y = (face_normal + i) % DIM;
         face_box.lower(i) = box.lower(y);
         face_box.upper(i) = box.upper(y);
      }
   }

   return(face_box);
}

/*
*************************************************************************
*									*
* Compute the overlap between two face centered boxes.  The algorithm	*
* is fairly straight-forward.  First, we perform a quick-and-dirty	*
* intersection to see if the boxes might overlap.  If that intersection	*
* is not empty, then we need to do a better job calculating the overlap	*
* for each dimension.  Note that the AMR index space boxes must be	*
* shifted into the face centered space before we calculate the proper	*
* intersections.							*
*									*
*************************************************************************
*/

template<int DIM> 
tbox::Pointer< hier::BoxOverlap<DIM> > FaceGeometry<DIM>::doOverlap(
   const FaceGeometry<DIM>& dst_geometry,
   const FaceGeometry<DIM>& src_geometry,
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

   // Compute the intersection (if any) for each of the face directions

   const hier::Box<DIM> quick_check =
      hier::Box<DIM>::grow(src_shift, 1) * hier::Box<DIM>::grow(dst_ghost, 1);

   if (!quick_check.empty()) {
      for (int d = 0; d < DIM; d++) {
         const hier::Box<DIM> dst_face = toFaceBox(dst_ghost, d);
         const hier::Box<DIM> src_face = toFaceBox(src_shift, d);
         const hier::Box<DIM> together = dst_face * src_face;
         if (!together.empty()) {
            dst_boxes[d].unionBoxes(together);
            if (!overwrite_interior) {
               const hier::Box<DIM> int_face = toFaceBox(dst_geometry.d_box, d);
               dst_boxes[d].removeIntersections(together,int_face);
            } else {
               dst_boxes[d].appendItem(together);
            }
         }  // if (!together.empty())
      }  // loop over dim
   }  // !quick_check.empty()

   // Create the face overlap data object using the boxes and source shift

   hier::BoxOverlap<DIM> *overlap = new FaceOverlap<DIM>(dst_boxes, src_offset);
   return(tbox::Pointer< hier::BoxOverlap<DIM> >(overlap));
}

}
}
#endif
