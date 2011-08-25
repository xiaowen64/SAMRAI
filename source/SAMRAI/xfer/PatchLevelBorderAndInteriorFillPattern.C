/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Abstract fill pattern class to provide interface for stencils
 *
 ************************************************************************/

#ifndef included_xfer_PatchLevelBorderAndInteriorFillPattern_C
#define included_xfer_PatchLevelBorderAndInteriorFillPattern_C

#include "SAMRAI/xfer/PatchLevelBorderAndInteriorFillPattern.h"
#include "SAMRAI/hier/BoxList.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/tbox/MathUtilities.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/xfer/PatchLevelBorderAndInteriorFillPattern.I"
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

PatchLevelBorderAndInteriorFillPattern::PatchLevelBorderAndInteriorFillPattern():
   d_max_fill_boxes(0)
{
}

/*
 *************************************************************************
 *									*
 * Destructor                                                           *
 *									*
 *************************************************************************
 */

PatchLevelBorderAndInteriorFillPattern::~PatchLevelBorderAndInteriorFillPattern()
{
}

/*
 *************************************************************************
 *                                                                       *
 * computeFillMappedBoxesAndNeighborhoodSets                             *
 *                                                                       *
 *************************************************************************
 */
void
PatchLevelBorderAndInteriorFillPattern::computeFillMappedBoxesAndNeighborhoodSets(
   hier::BoxSet& fill_mapped_boxes,
   hier::NeighborhoodSet& dst_to_fill_edges,
   const hier::MappedBoxLevel& dst_mapped_box_level,
   const hier::Connector& dst_to_dst,
   const hier::Connector& dst_to_src,
   const hier::Connector& src_to_dst,
   const hier::IntVector& fill_ghost_width)
{
   NULL_USE(dst_to_src);
   NULL_USE(src_to_dst);
   TBOX_DIM_ASSERT_CHECK_ARGS2(dst_mapped_box_level, fill_ghost_width);

   const hier::BoxSet& dst_mapped_boxes =
      dst_mapped_box_level.getMappedBoxes();

   /*
    * Grow each patch box and remove the level from it, except the
    * patch box itself.  (Do not fill ghost cells that are normally
    * filled by same mapped_box_level.  Do fill ghost cells that are
    * normally filled by coarser mapped_box_level.)
    */
   hier::LocalId last_id = dst_mapped_box_level.getLastLocalId();
   for (hier::BoxSet::const_iterator ni = dst_mapped_boxes.begin();
        ni != dst_mapped_boxes.end(); ++ni) {

      const hier::BoxId& gid(ni->getId());
      const hier::Box& dst_mapped_box =
         *dst_mapped_box_level.getMappedBox(gid);
      hier::BoxList fill_boxes(dst_mapped_box);
      fill_boxes.getFirstItem().grow(fill_ghost_width);
      const NeighborSet& nabrs =
         dst_to_dst.getNeighborSet(dst_mapped_box.getId());

      for (NeighborSet::const_iterator na = nabrs.begin();
           na != nabrs.end(); ++na) {
         if (!ni->isSpatiallyEqual(*na)) {
            if (dst_mapped_box.getBlockId() == na->getBlockId()) {
               fill_boxes.removeIntersections(*na);
            } else {

               tbox::ConstPointer<hier::GridGeometry> grid_geometry(
                  dst_mapped_box_level.getGridGeometry());

               const hier::BlockId& dst_block_id = dst_mapped_box.getBlockId();
               const hier::BlockId& nbr_block_id = na->getBlockId();

               TBOX_ASSERT(grid_geometry->areNeighbors(dst_block_id,
                     nbr_block_id));

               hier::Transformation::RotationIdentifier rotation =
                  grid_geometry->getRotationIdentifier(dst_block_id,
                     nbr_block_id);
               hier::IntVector offset(
                  grid_geometry->getOffset(dst_block_id, nbr_block_id));

               offset *= (dst_mapped_box_level.getRefinementRatio());

               hier::Transformation transformation(rotation, offset);

               hier::Box nbr_box(*na);
               transformation.transform(nbr_box);

               fill_boxes.removeIntersections(nbr_box);

            }
         }
      }

      if (!fill_boxes.isEmpty()) {
         d_max_fill_boxes = tbox::MathUtilities<int>::Max(d_max_fill_boxes,
               fill_boxes.size());

         NeighborSet& fill_nabrs = dst_to_fill_edges[gid];
         for (hier::BoxList::Iterator li(fill_boxes); li; li++) {
            hier::Box fill_mapped_box(*li,
                                      ++last_id,
                                      dst_mapped_box.getOwnerRank(),
                                      dst_mapped_box.getBlockId());
            fill_mapped_boxes.insert(fill_mapped_boxes.end(), fill_mapped_box);
            fill_nabrs.insert(fill_mapped_box);
         }
      }
   }
}

}
}
#endif
