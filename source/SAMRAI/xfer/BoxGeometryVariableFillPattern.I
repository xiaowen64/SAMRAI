/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Abstract fill pattern class to provide interface for stencils
 *
 ************************************************************************/

namespace SAMRAI {
namespace xfer {

/*
 *************************************************************************
 *
 * Calculate the overlap using the implemented calculateOverlap() method
 * for the destination geometry.
 *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
boost::shared_ptr<hier::BoxOverlap>
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
 *
 * Return the string name identifier.
 *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
const std::string&
BoxGeometryVariableFillPattern::getPatternName() const
{
   return s_name_id;
}

}
}
