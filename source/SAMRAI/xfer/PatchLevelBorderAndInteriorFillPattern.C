/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Abstract fill pattern class to provide interface for stencils 
 *
 ************************************************************************/

#ifndef included_xfer_PatchLevelBorderAndInteriorFillPattern_C
#define included_xfer_PatchLevelBorderAndInteriorFillPattern_C

#include "SAMRAI/xfer/PatchLevelBorderAndInteriorFillPattern.h"
#include "SAMRAI/hier/BoxList.h"
#include "SAMRAI/hier/MappedBox.h"
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
   hier::MappedBoxSet& fill_mapped_boxes,
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

   const hier::MappedBoxSet& dst_mapped_boxes =
      dst_mapped_box_level.getMappedBoxes();

   /*
    * Grow each patch box and remove the level from it, except the
    * patch box itself.  (Do not fill ghost cells that are normally
    * filled by same mapped_box_level.  Do fill ghost cells that are
    * normally filled by coarser mapped_box_level.)
    */
   hier::LocalId last_id = dst_mapped_box_level.getLastLocalId();
   for (hier::MappedBoxSet::const_iterator ni = dst_mapped_boxes.begin();
        ni != dst_mapped_boxes.end(); ++ni) {

      const hier::MappedBoxId& gid(ni->getId());
      const hier::MappedBox& dst_mapped_box =
         *dst_mapped_box_level.getMappedBox(gid);
      hier::BoxList fill_boxes(dst_mapped_box.getBox());
      fill_boxes.getFirstItem().grow(fill_ghost_width);
      const NeighborSet& nabrs =
         dst_to_dst.getNeighborSet(dst_mapped_box.getId());

      for (NeighborSet::const_iterator na = nabrs.begin();
           na != nabrs.end(); ++na) {
         if (*ni != *na && dst_mapped_box.getBlockId() == na->getBlockId()) {
            fill_boxes.removeIntersections(na->getBox());
         }
      }

      if (!fill_boxes.isEmpty()) {
         d_max_fill_boxes = tbox::MathUtilities<int>::Max(d_max_fill_boxes,
               fill_boxes.size());

         NeighborSet& fill_nabrs = dst_to_fill_edges[gid];
         for (hier::BoxList::Iterator li(fill_boxes); li; li++) {
            hier::MappedBox fill_mapped_box(*li,
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
