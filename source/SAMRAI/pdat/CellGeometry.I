/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

namespace SAMRAI {
namespace pdat {

SAMRAI_INLINE_KEYWORD
const hier::Box&
CellGeometry::getBox() const
{
   return d_box;
}

SAMRAI_INLINE_KEYWORD
const hier::IntVector&
CellGeometry::getGhosts() const
{
   return d_ghosts;
}

SAMRAI_INLINE_KEYWORD
hier::Box
CellGeometry::toCellBox(
   const hier::Box& box)
{
   return box;
}

SAMRAI_INLINE_KEYWORD
boost::shared_ptr<hier::BoxOverlap>
CellGeometry::doOverlap(
   const CellGeometry& dst_geometry,
   const CellGeometry& src_geometry,
   const hier::Box& src_mask,
   const hier::Box& fill_box,
   const bool overwrite_interior,
   const hier::Transformation& transformation,
   const hier::BoxContainer& dst_restrict_boxes)
{
   hier::BoxContainer dst_boxes;
   dst_geometry.computeDestinationBoxes(dst_boxes,
      src_geometry,
      src_mask,
      fill_box,
      overwrite_interior,
      transformation,
      dst_restrict_boxes);

   // Create the cell overlap data object using the boxes and source shift

   return boost::make_shared<CellOverlap>(dst_boxes, transformation);
}

}
}
