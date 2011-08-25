/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Abstract fill pattern class to provide interface for stencils
 *
 ************************************************************************/

#ifndef included_xfer_BoxGeometryVariableFillPattern_C
#define included_xfer_BoxGeometryVariableFillPattern_C

#include "SAMRAI/xfer/BoxGeometryVariableFillPattern.h"

#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace xfer {

const std::string BoxGeometryVariableFillPattern::s_name_id =
   "BOX_GEOMETRY_FILL_PATTERN";

/*
 *************************************************************************
 *                                                                       *
 * Default contructor only sets the string name identifier               *
 *                                                                       *
 *************************************************************************
 */

BoxGeometryVariableFillPattern::BoxGeometryVariableFillPattern()
{
}

/*
 *************************************************************************
 *									*
 * Destructor                                                            *
 *                                                                      *
 *************************************************************************
 */

BoxGeometryVariableFillPattern::~BoxGeometryVariableFillPattern()
{
}

/*
 *************************************************************************
 *                                                                       *
 * Calculate the overlap using the implemented calculateOverlap() method *
 * for the destination geometry.                                         *
 *                                                                       *
 *************************************************************************
 */

tbox::Pointer<hier::BoxOverlap>
BoxGeometryVariableFillPattern::calculateOverlap(
   const hier::BoxGeometry& dst_geometry,
   const hier::BoxGeometry& src_geometry,
   const hier::Box& dst_patch_box,
   const hier::Box& src_mask,
   const hier::Box& fill_box,
   const bool overwrite_interior,
   const hier::Transformation& transformation) const
{
   NULL_USE(dst_patch_box);
   TBOX_DIM_ASSERT_CHECK_ARGS2(dst_patch_box, src_mask);
   return dst_geometry.calculateOverlap(src_geometry, src_mask, fill_box,
      overwrite_interior, transformation);
}

/*
 *************************************************************************
 *                                                                       *
 * Return the string name identifier.                                    *
 *                                                                       *
 *************************************************************************
 */
const std::string& BoxGeometryVariableFillPattern::getPatternName() const
{
   return s_name_id;
}

/*
 *************************************************************************
 *                                                                       *
 * getStencilWidth() throws an error if called.  Only overridding        *
 * versions of this method in concrete subclasses should be called.      *
 *                                                                       *
 *************************************************************************
 */
const hier::IntVector& BoxGeometryVariableFillPattern::getStencilWidth()
{
   TBOX_ERROR(
      "BoxGeometryVariableFillPattern::getStencilWidth() should not be\n"
      << "called.  This pattern creates overlaps based on\n"
      << "the BoxGeometry objects and is not restricted to a\n"
      << "specific stencil.\n");

   /*
    * Dummy return value that will never get reached.
    */
   return hier::IntVector::getZero(tbox::Dimension(1));
}

/*
 *************************************************************************
 *                                                                       *
 * Compute BoxOverlap that specifies data to be filled by refinement     *
 * operator.                                                             *
 *                                                                       *
 *************************************************************************
 */
tbox::Pointer<hier::BoxOverlap>
BoxGeometryVariableFillPattern::computeFillBoxesOverlap(
   const hier::BoxList& fill_boxes,
   const hier::Box& patch_box,
   const hier::Box& data_box,
   const hier::PatchDataFactory& pdf) const
{
   /*
    * For this (default) case, the overlap is simply the intersection of
    * fill_boxes and data_box.
    */
   hier::Transformation transformation(
      hier::IntVector::getZero(patch_box.getDim()));

   hier::BoxList overlap_boxes(fill_boxes);
   overlap_boxes.intersectBoxes(data_box);

   return pdf.getBoxGeometry(patch_box)->setUpOverlap(overlap_boxes,
      transformation);
}

}
}
#endif
