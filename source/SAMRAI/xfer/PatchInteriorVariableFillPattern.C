/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Fill pattern class that fills patch interiors only
 *
 ************************************************************************/

#ifndef included_xfer_PatchInteriorVariableFillPattern_C
#define included_xfer_PatchInteriorVariableFillPattern_C

#include "SAMRAI/xfer/PatchInteriorVariableFillPattern.h"

#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace xfer {

const std::string PatchInteriorVariableFillPattern::s_name_id =
   "PATCH_INTERIOR_FILL_PATTERN";

/*
 *************************************************************************
 *                                                                       *
 * Default contructor only sets the string name identifier               *
 *                                                                       *
 *************************************************************************
 */

PatchInteriorVariableFillPattern::PatchInteriorVariableFillPattern(
   const tbox::Dimension& dim):
   d_dim(dim)
{
}

/*
 *************************************************************************
 *									*
 * Destructor                                                            *
 *                                                                      *
 *************************************************************************
 */

PatchInteriorVariableFillPattern::~PatchInteriorVariableFillPattern()
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
PatchInteriorVariableFillPattern::calculateOverlap(
   const hier::BoxGeometry& dst_geometry,
   const hier::BoxGeometry& src_geometry,
   const hier::Box& dst_patch_box,
   const hier::Box& src_mask,
   const hier::Box& fill_box,
   const bool overwrite_interior,
   const hier::Transformation& transformation) const
{
   NULL_USE(dst_patch_box);
   NULL_USE(overwrite_interior);
   TBOX_DIM_ASSERT_CHECK_ARGS2(dst_patch_box, src_mask);
   hier::BoxList dst_restrict_boxes(dst_patch_box);
   return dst_geometry.calculateOverlap(src_geometry, src_mask, fill_box,
      true, transformation,
      dst_restrict_boxes);
}

/*
 *************************************************************************
 *                                                                       *
 * Return the string name identifier.                                    *
 *                                                                       *
 *************************************************************************
 */
const std::string& PatchInteriorVariableFillPattern::getPatternName() const
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
const hier::IntVector& PatchInteriorVariableFillPattern::getStencilWidth()
{
   return hier::IntVector::getZero(d_dim);
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
PatchInteriorVariableFillPattern::computeFillBoxesOverlap(
   const hier::BoxList& fill_boxes,
   const hier::Box& patch_box,
   const hier::Box& data_box,
   const hier::PatchDataFactory& pdf) const
{
   /*
    * For this case, the overlap is simply the intersection of
    * fill_boxes, data_box, and patch_box.
    */
   hier::Transformation transformation(
      hier::IntVector::getZero(patch_box.getDim()));

   hier::BoxList overlap_boxes(fill_boxes);
   overlap_boxes.intersectBoxes(data_box * patch_box);

   return pdf.getBoxGeometry(patch_box)->setUpOverlap(overlap_boxes,
      transformation);
}

}
}
#endif
