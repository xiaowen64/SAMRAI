/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_OuterfaceGeometry_C
#define included_pdat_OuterfaceGeometry_C

#include "SAMRAI/pdat/OuterfaceGeometry.h"
#include "SAMRAI/pdat/FaceGeometry.h"
#include "SAMRAI/pdat/FaceOverlap.h"
#include "SAMRAI/hier/BoxContainerConstIterator.h"
#include "SAMRAI/tbox/Utilities.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/pdat/OuterfaceGeometry.I"
#endif
namespace SAMRAI {
namespace pdat {

/*
 *************************************************************************
 *
 * Create a face geometry object given the box and ghost cell width.
 *
 *************************************************************************
 */

OuterfaceGeometry::OuterfaceGeometry(
   const hier::Box& box,
   const hier::IntVector& ghosts):
   d_box(box),
   d_ghosts(ghosts)
{
   TBOX_ASSERT(ghosts.min() >= 0);
}

OuterfaceGeometry::~OuterfaceGeometry()
{
}

/*
 *************************************************************************
 *
 * Attempt to calculate the intersection between two outerface centered
 * box geometries.  The calculateOverlap() checks whether both arguments
 * are outerface geometries; if so, it compuates the intersection.  If
 * not, then it calls calculateOverlap() on the source object (if retry
 * is true) to allow the source a chance to calculate the intersection.
 * See the hier::BoxGeometry base class for more information about
 * the protocol.  A pointer to null is returned if the intersection
 * cannot be computed.
 *
 *************************************************************************
 */

tbox::Pointer<hier::BoxOverlap> OuterfaceGeometry::calculateOverlap(
   const hier::BoxGeometry& dst_geometry,
   const hier::BoxGeometry& src_geometry,
   const hier::Box& src_mask,
   const hier::Box& fill_box,
   const bool overwrite_interior,
   const hier::Transformation& transformation,
   const bool retry,
   const hier::BoxContainer& dst_restrict_boxes) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(d_box, src_mask);

   const FaceGeometry* t_dst =
      dynamic_cast<const FaceGeometry *>(&dst_geometry);
   const OuterfaceGeometry* t_src =
      dynamic_cast<const OuterfaceGeometry *>(&src_geometry);

   tbox::Pointer<hier::BoxOverlap> over(NULL);
   if ((t_src != NULL) && (t_dst != NULL)) {
      over = doOverlap(*t_dst, *t_src, src_mask, fill_box, overwrite_interior,
            transformation, dst_restrict_boxes);
   } else if (retry) {
      over = src_geometry.calculateOverlap(
            dst_geometry, src_geometry, src_mask, fill_box,
            overwrite_interior, transformation, false, dst_restrict_boxes);
   }
   return over;
}

/*
 *************************************************************************
 *
 * Compute the overlap between a face geometry destination box and an
 * outerface geometry source box.  The intersection algorithm is similar
 * the face geometry algorithm except that only the borders of source
 * are used in the intersection computation.
 *
 *************************************************************************
 */

tbox::Pointer<hier::BoxOverlap> OuterfaceGeometry::doOverlap(
   const FaceGeometry& dst_geometry,
   const OuterfaceGeometry& src_geometry,
   const hier::Box& src_mask,
   const hier::Box& fill_box,
   const bool overwrite_interior,
   const hier::Transformation& transformation,
   const hier::BoxContainer& dst_restrict_boxes)
{
   const tbox::Dimension& dim(src_mask.getDim());

   tbox::Array<hier::BoxContainer> dst_boxes(dim.getValue());
   const hier::IntVector& src_offset = transformation.getOffset();

   // Perform a quick-and-dirty intersection to see if the boxes might overlap

   const hier::Box src_box(
      hier::Box::grow(src_geometry.d_box, src_geometry.d_ghosts) * src_mask);
   hier::Box src_shift(src_box);
   transformation.transform(src_shift);
   const hier::Box dst_ghost(
      hier::Box::grow(dst_geometry.getBox(), dst_geometry.getGhosts()));

   // Compute the intersection (if any) for each of the face directions

   const hier::IntVector one_vector(dim, 1);

   const hier::Box quick_check(
      hier::Box::grow(src_shift, one_vector) * hier::Box::grow(dst_ghost,
         one_vector));

   if (!quick_check.empty()) {

      const hier::Box mask_shift(hier::Box::shift(src_mask, src_offset));

      for (int d = 0; d < dim.getValue(); d++) {

         const hier::Box msk_face(
            FaceGeometry::toFaceBox(mask_shift, d));
         const hier::Box dst_face(
            FaceGeometry::toFaceBox(dst_ghost, d));
         const hier::Box src_face(
            FaceGeometry::toFaceBox(src_shift, d));
         const hier::Box fill_face(
            FaceGeometry::toFaceBox(fill_box, d));

         const hier::Box together(dst_face * src_face * fill_face);

         if (!together.empty()) {

            // Add lower face intersection (if any) to the box list
            hier::Box low_face(src_face);
            low_face.upper(0) = low_face.lower(0);  //+ghosts;

            hier::Box low_overlap(low_face * msk_face * dst_face);
            if (!low_overlap.empty()) { 
               dst_boxes[d].pushBack(low_overlap);
            }

            // Add upper face intersection (if any) to the box list
            hier::Box hig_face(src_face);
            hig_face.lower(0) = hig_face.upper(0);  //-ghosts;

            hier::Box hig_overlap(hig_face * msk_face * dst_face);
            if (!hig_overlap.empty()) {
               dst_boxes[d].pushBack(hig_overlap);
            }

            // Take away the interior of over_write interior is not set
            if (!overwrite_interior) {
               dst_boxes[d].removeIntersections(
                  FaceGeometry::toFaceBox(dst_geometry.getBox(), d));
            }

         }  // if (!together.empty())

         if (dst_restrict_boxes.size() && dst_boxes[d].size()) {
            hier::BoxContainer face_restrict_boxes;
            for (hier::BoxContainer::ConstIterator b(dst_restrict_boxes);
                 b != dst_restrict_boxes.end(); ++b) {
               face_restrict_boxes.pushBack(FaceGeometry::toFaceBox(b(), d));
            }
            dst_boxes[d].intersectBoxes(face_restrict_boxes);
         }

      }  // loop over dim

   } // if (!quick_check.empty())

   // Create the face overlap data object using the boxes and source shift

   hier::BoxOverlap* overlap = new FaceOverlap(dst_boxes, transformation);
   return tbox::Pointer<hier::BoxOverlap>(overlap);
}

/*
 *************************************************************************
 *
 * Set up a FaceOverlap oject using the given boxes and offset
 *
 *************************************************************************
 */
tbox::Pointer<hier::BoxOverlap>
OuterfaceGeometry::setUpOverlap(
   const hier::BoxContainer& boxes,
   const hier::Transformation& transformation) const
{
   const tbox::Dimension& dim(transformation.getOffset().getDim());
   tbox::Array<hier::BoxContainer> dst_boxes(dim.getValue());

   for (hier::BoxContainer::ConstIterator b(boxes); b != boxes.end(); ++b) {
      for (int d = 0; d < dim.getValue(); d++) {
         hier::Box face_box(FaceGeometry::toFaceBox(b(), d));
         dst_boxes[d].pushBack(face_box);
      }
   }

   // Create the face overlap data object using the boxes and source shift
   hier::BoxOverlap* overlap = new FaceOverlap(dst_boxes, transformation);
   return tbox::Pointer<hier::BoxOverlap>(overlap);

}

}
}
#endif
