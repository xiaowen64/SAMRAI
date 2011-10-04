/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Abstract fill pattern class to provide interface for stencils
 *
 ************************************************************************/

#ifndef included_xfer_PatchLevelInteriorFillPattern_C
#define included_xfer_PatchLevelInteriorFillPattern_C

#include "SAMRAI/xfer/PatchLevelInteriorFillPattern.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainerOrderedConstIterator.h"
#include "SAMRAI/tbox/MathUtilities.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/xfer/PatchLevelInteriorFillPattern.I"
#endif

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
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

PatchLevelInteriorFillPattern::PatchLevelInteriorFillPattern():
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

PatchLevelInteriorFillPattern::~PatchLevelInteriorFillPattern()
{
}

/*
 *************************************************************************
 *
 * computeFillBoxesAndNeighborhoodSets
 *
 *************************************************************************
 */

void PatchLevelInteriorFillPattern::computeFillBoxesAndNeighborhoodSets(
   hier::BoxSet& fill_mapped_boxes,
   hier::NeighborhoodSet& dst_to_fill_edges,
   const hier::BoxLevel& dst_mapped_box_level,
   const hier::Connector& dst_to_dst,
   const hier::Connector& dst_to_src,
   const hier::Connector& src_to_dst,
   const hier::IntVector& fill_ghost_width)
{
   NULL_USE(dst_to_dst);
   NULL_USE(dst_to_src);
   NULL_USE(src_to_dst);
   NULL_USE(fill_ghost_width);
   TBOX_DIM_ASSERT_CHECK_ARGS2(dst_mapped_box_level, fill_ghost_width);

   const hier::BoxSet& dst_mapped_boxes =
      dst_mapped_box_level.getBoxes();

   /*
    * Fill just the interior.  Disregard gcw.
    */
   for (hier::BoxSet::OrderedConstIterator ni = dst_mapped_boxes.orderedBegin();
        ni != dst_mapped_boxes.orderedEnd(); ++ni) {
      const hier::BoxId& gid = ni->getId();
      const hier::Box& dst_mapped_box =
         *dst_mapped_box_level.getBox(gid);
      fill_mapped_boxes.insert(fill_mapped_boxes.orderedEnd(), dst_mapped_box);
      dst_to_fill_edges.insertNeighbor(gid, dst_mapped_box);
   }
}

/*
 *************************************************************************
 *
 * computeDestinationFillBoxesOnSourceProc
 *
 *************************************************************************
 */

void PatchLevelInteriorFillPattern::computeDestinationFillBoxesOnSourceProc(
   FillSet& dst_fill_boxes_on_src_proc,
   const hier::BoxLevel& dst_mapped_box_level,
   const hier::Connector& src_to_dst,
   const hier::IntVector& fill_ghost_width)
{
   NULL_USE(fill_ghost_width);
   TBOX_DIM_ASSERT_CHECK_ARGS2(dst_mapped_box_level, fill_ghost_width);

   const tbox::Dimension& dim(fill_ghost_width.getDim());
   /*
    * src_to_dst initialized only when there is a src mapped_box_level.
    * Without the src mapped_box_level, we do not need to compute
    * dst_fill_boxes_on_src_proc.
    *
    * For PatchLevelInteriorFillPattern, the src owner can compute fill
    * boxes for all its dst neighbors using local data.  This info is
    * stored in dst_fill_boxes_on_src_proc.
    */
   NeighborSet tmp_nabrs(dim), all_dst_nabrs(dim);
   src_to_dst.getNeighborhoodSets().getNeighbors(tmp_nabrs);
   tmp_nabrs.unshiftPeriodicImageBoxes(
      all_dst_nabrs,
      dst_mapped_box_level.getRefinementRatio());
   tmp_nabrs.clear();
   for (NeighborSet::OrderedConstIterator na = all_dst_nabrs.orderedBegin();
        na != all_dst_nabrs.orderedEnd(); ++na) {
      hier::BoxSet& fill_boxes =
         dst_fill_boxes_on_src_proc.getNeighborSet(na->getId(), dim);
      fill_boxes.insert(*na);
      d_max_fill_boxes = tbox::MathUtilities<int>::Max(d_max_fill_boxes,
            static_cast<int>(fill_boxes.size()));
   }
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif

#endif
