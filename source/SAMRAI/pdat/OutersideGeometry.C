/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_OutersideGeometry_C
#define included_pdat_OutersideGeometry_C

#include "SAMRAI/pdat/OutersideGeometry.h"
#include "SAMRAI/pdat/SideGeometry.h"
#include "SAMRAI/pdat/SideOverlap.h"
#include "SAMRAI/hier/BoxList.h"
#include "SAMRAI/tbox/Utilities.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/pdat/OutersideGeometry.I"
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

OutersideGeometry::OutersideGeometry(
   const hier::Box& box,
   const hier::IntVector& ghosts):
   d_box(box),
   d_ghosts(ghosts)

{
   TBOX_DIM_ASSERT_CHECK_ARGS2(box, ghosts);
   TBOX_ASSERT(ghosts.min() >= 0);

}

OutersideGeometry::~OutersideGeometry()
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
 * See the hier::BoxGeometry base class for more information about the	*
 * protocol.  A pointer to null is returned if the intersection canot be	*
 * computed.								*
 *                                                                      *
 *************************************************************************
 */

tbox::Pointer<hier::BoxOverlap> OutersideGeometry::calculateOverlap(
   const hier::BoxGeometry& dst_geometry,
   const hier::BoxGeometry& src_geometry,
   const hier::Box& src_mask,
   const hier::Box& fill_box,
   const bool overwrite_interior,
   const hier::Transformation& transformation,
   const bool retry,
   const hier::BoxList& dst_restrict_boxes) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(d_box, src_mask);

   const SideGeometry* t_dst =
      dynamic_cast<const SideGeometry *>(&dst_geometry);
   const OutersideGeometry* t_src =
      dynamic_cast<const OutersideGeometry *>(&src_geometry);

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
 *									*
 * Compute the overlap between a side geometry destination box and an	*
 * outerside geometry source box.  The intersection algorithm is similar	*
 * the side geometry algorithm except that only the borders of source	*
 * are used in the intersection computation.				*
 *									*
 *************************************************************************
 */

tbox::Pointer<hier::BoxOverlap> OutersideGeometry::doOverlap(
   const SideGeometry& dst_geometry,
   const OutersideGeometry& src_geometry,
   const hier::Box& src_mask,
   const hier::Box& fill_box,
   const bool overwrite_interior,
   const hier::Transformation& transformation,
   const hier::BoxList& dst_restrict_boxes)
{
   const hier::IntVector& src_offset = transformation.getOffset();
   TBOX_DIM_ASSERT_CHECK_ARGS2(src_mask, src_offset);

   const tbox::Dimension& dim(src_mask.getDim());

   TBOX_ASSERT(dst_geometry.getDirectionVector() == hier::IntVector::getOne(dim));

   tbox::Array<hier::BoxList> dst_boxes(dim.getValue());

   // Perform a quick-and-dirty intersection to see if the boxes might overlap

   const hier::Box src_box(
      hier::Box::grow(src_geometry.d_box, src_geometry.d_ghosts) * src_mask);
   hier::Box src_shift(src_box);
   transformation.transform(src_shift);
   const hier::Box dst_ghost(
      hier::Box::grow(dst_geometry.getBox(), dst_geometry.getGhosts()));

   // Compute the intersection (if any) for each of the side directions

   const hier::IntVector one_vector(dim, 1);

   const hier::Box quick_check(
      hier::Box::grow(src_shift, one_vector) * hier::Box::grow(dst_ghost,
         one_vector));

   if (!quick_check.empty()) {

      const hier::Box mask_shift = hier::Box::shift(src_mask, src_offset);

      for (int d = 0; d < dim.getValue(); d++) {

         const hier::Box msk_side(
            SideGeometry::toSideBox(mask_shift, d));
         const hier::Box dst_side(
            SideGeometry::toSideBox(dst_ghost, d));
         const hier::Box src_side(
            SideGeometry::toSideBox(src_shift, d));
         const hier::Box fill_side(
            SideGeometry::toSideBox(fill_box, d));

         const hier::Box together(dst_side * src_side);

         if (!together.empty()) {

            // Add lower side intersection (if any) to the box list
            hier::Box low_side(src_side);
            low_side.upper(d) = low_side.lower(d); //+ghosts;
            dst_boxes[d].unionBoxes(low_side * msk_side * dst_side * fill_side);

            // Add upper side intersection (if any) to the box list
            hier::Box hig_side(src_side);
            hig_side.lower(d) = hig_side.upper(d); //-ghosts;
            dst_boxes[d].unionBoxes(hig_side * msk_side * dst_side * fill_side);

            // Take away the interior of over_write interior is not set
            if (!overwrite_interior) {
               dst_boxes[d].removeIntersections(
                  SideGeometry::toSideBox(dst_geometry.getBox(), d));
            }

         }  // if (!together.empty())

         if (dst_restrict_boxes.size() && dst_boxes[d].size()) {
            hier::BoxList side_restrict_boxes;
            for (hier::BoxList::Iterator b(dst_restrict_boxes); b; b++) {
               side_restrict_boxes.appendItem(SideGeometry::toSideBox(b(), d));
            }
            dst_boxes[d].intersectBoxes(side_restrict_boxes);
         }

      }  // loop over dim

   } // if (!quick_check.empty())

   // Create the side overlap data object using the boxes and source shift

   hier::BoxOverlap* overlap = new SideOverlap(dst_boxes, transformation);
   return tbox::Pointer<hier::BoxOverlap>(overlap);
}

/*
 *************************************************************************
 *                                                                       *
 * Set up a SideOverlap oject using the given boxes and offset           *
 *                                                                       *
 *************************************************************************
 */
tbox::Pointer<hier::BoxOverlap>
OutersideGeometry::setUpOverlap(
   const hier::BoxList& boxes,
   const hier::Transformation& transformation) const
{
   const tbox::Dimension& dim(transformation.getOffset().getDim());
   tbox::Array<hier::BoxList> dst_boxes(dim.getValue());

   for (hier::BoxList::Iterator b(boxes); b; b++) {
      for (int d = 0; d < dim.getValue(); d++) {
         hier::Box side_box(SideGeometry::toSideBox(b(), d));
         dst_boxes[d].appendItem(side_box);
      }
   }

   // Create the side overlap data object using the boxes and source shift
   hier::BoxOverlap* overlap = new SideOverlap(dst_boxes, transformation);
   return tbox::Pointer<hier::BoxOverlap>(overlap);

}

}
}
#endif
