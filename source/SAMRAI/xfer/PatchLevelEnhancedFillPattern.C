/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Abstract fill pattern class to provide interface for stencils
 *
 ************************************************************************/

#ifndef included_xfer_PatchLevelEnhancedFillPattern_C
#define included_xfer_PatchLevelEnhancedFillPattern_C

#include "SAMRAI/xfer/PatchLevelEnhancedFillPattern.h"
#include "SAMRAI/hier/BoxContainerIterator.h"
#include "SAMRAI/hier/RealBoxConstIterator.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/tbox/MathUtilities.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/xfer/PatchLevelEnhancedFillPattern.I"
#endif

namespace SAMRAI {
namespace xfer {

/*
 *************************************************************************
 *
 * Default constructor
 *
 *************************************************************************
 */

PatchLevelEnhancedFillPattern::PatchLevelEnhancedFillPattern():
   d_max_fill_boxes(0)
{
}

/*
 *************************************************************************
 *
 * Destructor
 *
 *************************************************************************
 */

PatchLevelEnhancedFillPattern::~PatchLevelEnhancedFillPattern()
{
}

/*
 *************************************************************************
 *
 * computeFillBoxesAndNeighborhoodSets
 *
 *************************************************************************
 */
void PatchLevelEnhancedFillPattern::computeFillBoxesAndNeighborhoodSets(
   hier::BoxLevel& fill_mapped_boxes,
   hier::Connector& dst_to_fill,
   const hier::BoxLevel& dst_mapped_box_level,
   const hier::Connector& dst_to_dst,
   const hier::Connector& dst_to_src,
   const hier::Connector& src_to_dst,
   const hier::IntVector& fill_ghost_width)
{
   NULL_USE(dst_to_dst);
   NULL_USE(dst_to_src);
   NULL_USE(src_to_dst);
   TBOX_DIM_ASSERT_CHECK_ARGS2(dst_mapped_box_level, fill_ghost_width);

   tbox::ConstPointer<hier::GridGeometry> grid_geometry(
      dst_mapped_box_level.getGridGeometry());

   const hier::BoxSet& dst_mapped_boxes =
      dst_mapped_box_level.getBoxes();

   hier::LocalId last_id = dst_mapped_box_level.getLastLocalId();
   for (hier::RealBoxConstIterator ni(dst_mapped_boxes);
        ni.isValid(); ++ni) {
      const hier::Box& dst_mapped_box = *ni;
      hier::BoxList fill_boxes(dst_mapped_box);
      fill_boxes.front().grow(fill_ghost_width);

      const tbox::List<hier::GridGeometry::Neighbor>& neighbors =
         grid_geometry->getNeighbors(dst_mapped_box.getBlockId());

      hier::BoxList constructed_fill_boxes(dst_mapped_box.getDim());

      for (tbox::List<hier::GridGeometry::Neighbor>::Iterator ni(neighbors);
           ni; ni++) {

         if (ni().isSingularity()) {

            hier::BoxList encon_boxes(ni().getTransformedDomain());
            encon_boxes.refine(dst_mapped_box_level.getRefinementRatio());
            encon_boxes.intersectBoxes(fill_boxes);
            encon_boxes.removeIntersections(constructed_fill_boxes);

            if (encon_boxes.size()) {

               dst_to_fill.makeEmptyLocalNeighborhood(dst_mapped_box.getId());
               for (hier::BoxList::Iterator ei(encon_boxes);
                    ei != encon_boxes.end(); ei++) {

                  hier::Box fill_mapped_box(
                     *ei,
                     ++last_id,
                     dst_mapped_box.getOwnerRank(),
                     dst_mapped_box.getBlockId());

                  fill_mapped_boxes.addBoxWithoutUpdate(fill_mapped_box);

                  dst_to_fill.insertLocalNeighbor(fill_mapped_box,
                     dst_mapped_box.getId());

                  constructed_fill_boxes.pushBack(*ei);
               }
            }
         }
      }

      d_max_fill_boxes = tbox::MathUtilities<int>::Max(
            d_max_fill_boxes,
            constructed_fill_boxes.size());
   }
   fill_mapped_boxes.finalize();
}

}
}
#endif
