/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Abstract fill pattern class to provide interface for stencils 
 *
 ************************************************************************/

#ifndef included_xfer_PatchLevelEnhancedFillPattern_C
#define included_xfer_PatchLevelEnhancedFillPattern_C

#include "SAMRAI/xfer/PatchLevelEnhancedFillPattern.h"
#include "SAMRAI/hier/RealMappedBoxConstIterator.h"
#include "SAMRAI/hier/MappedBox.h"
#include "SAMRAI/tbox/MathUtilities.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/xfer/PatchLevelEnhancedFillPattern.I"
#endif

namespace SAMRAI {
namespace xfer {

/*
 *************************************************************************
 *                                                                       *
 * Default constructor                                                   *
 *                                                                       *
 *************************************************************************
 */

PatchLevelEnhancedFillPattern::PatchLevelEnhancedFillPattern():
   d_max_fill_boxes(0)
{
}

/*
 *************************************************************************
 *									*
 * Destructor                                                            *
 *									*
 *************************************************************************
 */

PatchLevelEnhancedFillPattern::~PatchLevelEnhancedFillPattern()
{
}

/*
 *************************************************************************
 *                                                                       *
 * computeFillMappedBoxesAndNeighborhoodSets                             *
 *                                                                       *
 *************************************************************************
 */
void PatchLevelEnhancedFillPattern::computeFillMappedBoxesAndNeighborhoodSets(
   hier::MappedBoxSet& fill_mapped_boxes,
   hier::NeighborhoodSet& dst_to_fill_edges,
   const hier::MappedBoxLevel& dst_mapped_box_level,
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

   const hier::MappedBoxSet& dst_mapped_boxes =
      dst_mapped_box_level.getMappedBoxes();

   hier::LocalId last_id = dst_mapped_box_level.getLastLocalId();
   for (hier::RealMappedBoxConstIterator ni(dst_mapped_boxes);
        ni.isValid(); ++ni) {
      const hier::MappedBox& dst_mapped_box = *ni;
      hier::BoxList fill_boxes(dst_mapped_box.getBox());
      fill_boxes.getFirstItem().grow(fill_ghost_width);

      const tbox::List<hier::GridGeometry::Neighbor>& neighbors =
         grid_geometry->getNeighbors(dst_mapped_box.getBlockId().getBlockValue());

      hier::BoxList constructed_fill_boxes(dst_mapped_box.getDim());

      for (tbox::List<hier::GridGeometry::Neighbor>::Iterator ni(neighbors);
           ni; ni++) {

         if (ni().isSingularity()) {

            hier::BoxList encon_boxes(ni().getTranslatedDomain());
            encon_boxes.refine(dst_mapped_box_level.getRefinementRatio());
            encon_boxes.intersectBoxes(fill_boxes);
            encon_boxes.removeIntersections(constructed_fill_boxes);

            if (encon_boxes.size()) {

               NeighborSet& fill_nabrs =
                  dst_to_fill_edges[dst_mapped_box.getId()];

               for (hier::BoxList::Iterator ei(encon_boxes); ei; ei++) {

                  hier::MappedBox fill_mapped_box(
                     *ei,
                     ++last_id,
                     dst_mapped_box.getOwnerRank(),
                     dst_mapped_box.getBlockId());

                  fill_mapped_boxes.insert(fill_mapped_boxes.end(),
                                           fill_mapped_box);


                  fill_nabrs.insert(fill_mapped_box);

                  constructed_fill_boxes.appendItem(*ei);
               }
            }
         }
      }

      d_max_fill_boxes = tbox::MathUtilities<int>::Max(
                            d_max_fill_boxes,
                            constructed_fill_boxes.size());
   }
}

}
}
#endif
