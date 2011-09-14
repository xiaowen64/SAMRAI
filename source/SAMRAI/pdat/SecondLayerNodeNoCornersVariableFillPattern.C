/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Fill pattern class to provide interface for stencils
 *
 ************************************************************************/

#ifndef included_pdat_SecondLayerNodeNoCornersVariableFillPattern_C
#define included_pdat_SecondLayerNodeNoCornersVariableFillPattern_C

#include "SAMRAI/pdat/SecondLayerNodeNoCornersVariableFillPattern.h"

#include "SAMRAI/hier/BoxContainerIterator.h"
#include "SAMRAI/pdat/NodeGeometry.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace pdat {

const std::string SecondLayerNodeNoCornersVariableFillPattern::s_name_id =
   "SECOND_LAYER_NODE_NO_CORNERS_FILL_PATTERN";

/*
 *************************************************************************
 *                                                                       *
 * Constructor                                                           *
 *                                                                       *
 *************************************************************************
 */

SecondLayerNodeNoCornersVariableFillPattern::
SecondLayerNodeNoCornersVariableFillPattern(
   const tbox::Dimension& dim):
   d_dim(dim)
{
}

/*
 *************************************************************************
 *                                                                       *
 * Destructor                                                            *
 *                                                                       *
 *************************************************************************
 */

SecondLayerNodeNoCornersVariableFillPattern::~
SecondLayerNodeNoCornersVariableFillPattern()
{
}

/*
 *************************************************************************
 *                                                                       *
 * Calculate the overlap according to the desired pattern                *
 *                                                                       *
 *************************************************************************
 */

tbox::Pointer<hier::BoxOverlap>
SecondLayerNodeNoCornersVariableFillPattern::calculateOverlap(
   const hier::BoxGeometry& dst_geometry,
   const hier::BoxGeometry& src_geometry,
   const hier::Box& dst_patch_box,
   const hier::Box& src_mask,
   const hier::Box& fill_box,
   const bool overwrite_interior,
   const hier::Transformation& transformation) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(dst_patch_box, src_mask);
   NULL_USE(overwrite_interior);

   hier::BoxList dst_boxes(dst_patch_box.getDim());

   hier::Box dst_node_box(pdat::NodeGeometry::toNodeBox(dst_patch_box));
   hier::Box src_node_mask(pdat::NodeGeometry::toNodeBox(src_mask));

   bool corner_overlap = ((dst_node_box * src_node_mask).size() == 1)
      ? true : false;

   if (!corner_overlap) {
      hier::BoxList stencil_boxes(dst_patch_box.getDim());
      computeStencilBoxes(stencil_boxes, dst_patch_box);

      const NodeGeometry* t_dst =
         dynamic_cast<const NodeGeometry *>(&dst_geometry);
      const NodeGeometry* t_src =
         dynamic_cast<const NodeGeometry *>(&src_geometry);

      TBOX_ASSERT(t_dst);
      TBOX_ASSERT(t_src);

      t_dst->computeDestinationBoxes(dst_boxes, *t_src, src_mask, fill_box,
         false, transformation);

      dst_boxes.intersectBoxes(stencil_boxes);
   }

   hier::BoxOverlap* overlap = new NodeOverlap(dst_boxes, transformation);

   return tbox::Pointer<hier::BoxOverlap>(overlap);

}

/*
 *************************************************************************
 *                                                                       *
 * Return the stencil width (1)                                          *
 *                                                                       *
 *************************************************************************
 */

const hier::IntVector& SecondLayerNodeNoCornersVariableFillPattern::
getStencilWidth()
{
   return hier::IntVector::getOne(d_dim);
}

/*
 *************************************************************************
 *                                                                       *
 * Return the string name identifier                                     *
 *                                                                       *
 *************************************************************************
 */

const std::string& SecondLayerNodeNoCornersVariableFillPattern::getPatternName()
const
{
   return s_name_id;
}

/*
 *************************************************************************
 *                                                                       *
 * Compute the boxes for the stencil around a given patch box            *
 *                                                                       *
 *************************************************************************
 */

void SecondLayerNodeNoCornersVariableFillPattern::computeStencilBoxes(
   hier::BoxList& stencil_boxes,
   const hier::Box& dst_box) const
{
   TBOX_ASSERT(stencil_boxes.size() == 0);

   const tbox::Dimension& dim = dst_box.getDim();
   hier::Box dst_node_box(pdat::NodeGeometry::toNodeBox(dst_box));

   for (unsigned short i = 0; i < dim.getValue(); i++) {
      hier::Box low_box(dst_node_box);
      low_box.lower(i) = dst_node_box.lower(i) - 1;
      low_box.upper(i) = low_box.lower(i);
      stencil_boxes.pushFront(low_box);

      hier::Box high_box(dst_node_box);
      high_box.lower(i) = dst_node_box.upper(i) + 1;
      high_box.upper(i) = high_box.lower(i);
      stencil_boxes.pushFront(high_box);
   }
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
SecondLayerNodeNoCornersVariableFillPattern::computeFillBoxesOverlap(
   const hier::BoxList& fill_boxes,
   const hier::Box& patch_box,
   const hier::Box& data_box,
   const hier::PatchDataFactory& pdf) const
{
   NULL_USE(pdf);
   const tbox::Dimension& dim = patch_box.getDim();

   hier::BoxList stencil_boxes(dim);
   computeStencilBoxes(stencil_boxes, patch_box);

   hier::BoxList overlap_boxes(fill_boxes);

   /*
    * This is the equivalent of converting every box in overlap_boxes
    * to a node centering, which must be done before intersecting with
    * stencil_boxes, which is node-centered.
    */
   for (hier::BoxList::Iterator b(overlap_boxes); b; b++) {
      b().growUpper(hier::IntVector::getOne(dim));
   }

   overlap_boxes.intersectBoxes(pdat::NodeGeometry::toNodeBox(data_box));

   overlap_boxes.intersectBoxes(stencil_boxes);

   overlap_boxes.coalesce();

   hier::BoxOverlap* overlap =
      new pdat::NodeOverlap(
         overlap_boxes,
         hier::Transformation(hier::IntVector::getZero(dim)));
   return tbox::Pointer<hier::BoxOverlap>(overlap);
}

}
}
#endif
