//
// File:	OuteredgeGeometry.C
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Release:	$Name$
// Revision:	$Revision: 692 $
// Modified:	$Date: 2005-10-28 13:37:46 -0700 (Fri, 28 Oct 2005) $
// Description:	Box geometry information for outeredge centered objects
//

#ifndef included_pdat_OuteredgeGeometry_C
#define included_pdat_OuteredgeGeometry_C

#include "OuteredgeGeometry.h"

#include "BoxList.h"
#include "EdgeGeometry.h"
#include "EdgeOverlap.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

#ifdef DEBUG_NO_INLINE
#include "OuteredgeGeometry.I"
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

template<int DIM> OuteredgeGeometry<DIM>::OuteredgeGeometry(
   const hier::Box<DIM>& box, const hier::IntVector<DIM>& ghosts)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(ghosts.min() >= 0);
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
   const hier::Box<DIM> src_shift =
      hier::Box<DIM>::shift(src_box, src_offset);
   const hier::Box<DIM> dst_ghost =
      hier::Box<DIM>::grow(dst_geometry.getBox(), dst_geometry.getGhosts());

   // Compute the intersection (if any) for each of the edge directions

   const hier::Box<DIM> quick_check =
      hier::Box<DIM>::grow(src_shift, 1) * 
      hier::Box<DIM>::grow(dst_ghost, 1);

   if (!quick_check.empty()) {

      // a = axis
      for (int a = 0; a < DIM; a++) {

         const hier::Box<DIM> dst_box = 
            pdat::EdgeGeometry<DIM>::toEdgeBox(dst_ghost, a);
         const hier::Box<DIM> src_box = 
            pdat::EdgeGeometry<DIM>::toEdgeBox(src_shift, a);
         const hier::Box<DIM> together = dst_box * src_box;

         if (!together.empty()) {

            // f = face normal
            for (int f = 0; f < DIM; f++) {

               hier::Box<DIM> boxlo = src_box;
               boxlo.lower(f) = src_box.lower(f);
               boxlo.upper(f) = src_box.lower(f);
               hier::Box<DIM> boxup = src_box;
               boxup.lower(f) = src_box.upper(f);
               boxup.upper(f) = src_box.upper(f);
               
               // Data NULL when f = a
               if (f == a) {

                  boxlo.setEmpty();
                  boxup.setEmpty();

               } else {

                  trimBoxes(boxlo, boxup, a, f);
                  
                  // Add upper/lower side intersection (if any) to the box list
                  dst_boxes[a].unionBoxes(boxlo * dst_box);
                  dst_boxes[a].unionBoxes(boxup * dst_box);
               }
            } // loop over f
  
            if (!overwrite_interior) {
               const hier::Box<DIM> int_edge = 
                  pdat::EdgeGeometry<DIM>::toEdgeBox(dst_geometry.getBox(), 
                                                     a);
               dst_boxes[a].removeIntersections(int_edge);
            }

         }  // if (!together.empty())

      }  // loop over a

   }  // if (!quick_check.empty())


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

   const hier::Box<DIM> src_mask_box =
      hier::Box<DIM>::grow(src_geometry.d_box, src_geometry.d_ghosts) * src_mask;
   const hier::Box<DIM> src_shift =
      hier::Box<DIM>::shift(src_mask_box, src_offset);
   const hier::Box<DIM> dst_ghost =
      hier::Box<DIM>::grow(dst_geometry.getBox(), dst_geometry.getGhosts());

   // Compute the intersection (if any) for each of the edge directions

   const hier::Box<DIM> quick_check =
      hier::Box<DIM>::grow(src_shift, 1) * hier::Box<DIM>::grow(dst_ghost, 1);

   if (!quick_check.empty()) {

      // a = axis
      for (int a = 0; a < DIM; a++) {

         const hier::Box<DIM> dst_box = 
            pdat::EdgeGeometry<DIM>::toEdgeBox(dst_ghost, a);
         const hier::Box<DIM> src_box = 
            pdat::EdgeGeometry<DIM>::toEdgeBox(src_shift, a);
         const hier::Box<DIM> together = dst_box * src_box;

         if (!together.empty()) {

            // sf = source face normal
            for (int sf = 0; sf < DIM; sf++) {

               hier::Box<DIM> src_boxlo = src_box;
               src_boxlo.upper(sf) = src_box.lower(sf);
               hier::Box<DIM> src_boxup = src_box;
               src_boxup.lower(sf) = src_box.upper(sf);

               if (sf == a) {
                  src_boxlo.setEmpty();
                  src_boxup.setEmpty();
               } else {
                  trimBoxes(src_boxlo, src_boxup, a, sf);                  
               }

               // df = destination face normal
               for (int df = 0; df < DIM; df++) {
                  
                  hier::Box<DIM> dst_boxlo = dst_box;
                  dst_boxlo.upper(df) = dst_box.lower(df);
                  hier::Box<DIM> dst_boxup = dst_box;
                  dst_boxup.lower(df) = dst_box.upper(df);
                  
                  if (df == a) {
                     dst_boxlo.setEmpty();
                     dst_boxup.setEmpty();                     
                  } else {                     
                     trimBoxes(dst_boxlo, dst_boxup, a, df);
                  }
               
                  // Determine overlap
                  dst_boxes[a].unionBoxes(src_boxlo*dst_boxlo);
                  dst_boxes[a].unionBoxes(src_boxup*dst_boxlo);
                  dst_boxes[a].unionBoxes(src_boxlo*dst_boxup);
                  dst_boxes[a].unionBoxes(src_boxup*dst_boxup);

               } // loop over df
            } // loop over sf
         } // !together.empty()

         if (!overwrite_interior) {
            const hier::Box<DIM> int_edge = 
               pdat::EdgeGeometry<DIM>::toEdgeBox(dst_geometry.getBox(), a);
            dst_boxes[a].removeIntersections(int_edge);
         }

      } // loop over a
   } // !quick_check.isEmpty()
   
   // Create the edge overlap data object using the boxes and source shift
   hier::BoxOverlap<DIM> *overlap = 
      new pdat::EdgeOverlap<DIM>(dst_boxes, src_offset);
   return(tbox::Pointer<hier::BoxOverlap<DIM> >(overlap));
}

/*
*************************************************************************
*									*
* Trim the databoxes for particular axis,face_nrml directions to avoid  *
* duplicating data on common edges.					*
*									*
*************************************************************************
*/

template<int DIM> void 
OuteredgeGeometry<DIM>::trimBoxes(hier::Box<DIM>& boxlo, 
                                  hier::Box<DIM>& boxup,
                                  const int axis,
                                  const int face_nrml)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(axis >= 0 && axis < DIM);
   assert(face_nrml >= 0 && face_nrml < DIM);
#endif
   /*
    * In 1D and 2D, no trimming is necessary.  Simply set empty
    * boxes when face_nrml = axis.
    */
   if (DIM < 3) {
      
      if (face_nrml == axis) {
         boxlo.setEmpty();
         boxup.setEmpty();
      }
      
   } else {
      
      /*
       * For 3D boxes, trim boxlo, boxup boxes in the largest non-null 
       * dimension. i.e.
       *
       * if axis = 0, face_nrml = 0 - data is NULL (boxlo = boxhi = empty)
       * if axis = 0, face_nrml = 1 - trim in Z
       * if axis = 0, face_nrml = 2 - no trimming
       *
       * if axis = 1, face_nrml = 0 - trim in Z
       * if axis = 1, face_nrml = 1 - data is NULL (boxlo = boxhi = empty)
       * if axis = 1, face_nrml = 2 - no trimming
       *
       * if axis = 2, face_nrml = 0 - trim in Y
       * if axis = 2, face_nrml = 1 - no trimming
       * if axis = 2, face_nrml = 2 - data is NULL (boxlo = boxhi = empty)
       */

      int trim_dir;
   
      if (axis == 0) {
         if (face_nrml == 0) {
            boxlo.setEmpty();
            boxup.setEmpty();
         } else if (face_nrml == 1) {
            trim_dir = 2;
            boxlo.lower(trim_dir) += 1;
            boxlo.upper(trim_dir) -= 1;
            boxup.lower(trim_dir) += 1;
            boxup.upper(trim_dir) -= 1;
         } 
      }
   
      if (axis == 1) {
         if (face_nrml == 0) {
            trim_dir = 2;
            boxlo.lower(trim_dir) += 1;
            boxlo.upper(trim_dir) -= 1;
            boxup.lower(trim_dir) += 1;
            boxup.upper(trim_dir) -= 1;
         } else if (face_nrml == 1) {
            boxlo.setEmpty();
            boxup.setEmpty();
         } 
      }

      if (axis == 2) {
         if (face_nrml == 0) {
            trim_dir = 1;
            boxlo.lower(trim_dir) += 1;
            boxlo.upper(trim_dir) -= 1;
            boxup.lower(trim_dir) += 1;
            boxup.upper(trim_dir) -= 1;
         } else if (face_nrml == 2) {
            boxlo.setEmpty();
            boxup.setEmpty();
         } 
      }
   }   

}

}
}
#endif

