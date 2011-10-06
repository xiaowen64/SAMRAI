/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Abstract fill pattern class to provide interface for stencils
 *
 ************************************************************************/

#ifndef included_xfer_PatchLevelFullFillPattern_C
#define included_xfer_PatchLevelFullFillPattern_C

#include "SAMRAI/xfer/PatchLevelFullFillPattern.h"
#include "SAMRAI/hier/RealBoxConstIterator.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/tbox/MathUtilities.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/xfer/PatchLevelFullFillPattern.I"
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

PatchLevelFullFillPattern::PatchLevelFullFillPattern():
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

PatchLevelFullFillPattern::~PatchLevelFullFillPattern()
{
}

/*
 *************************************************************************
 *
 * computeFillBoxesAndNeighborhoodSets
 *
 *************************************************************************
 */

void PatchLevelFullFillPattern::computeFillBoxesAndNeighborhoodSets(
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

   const hier::BoxSet& dst_mapped_boxes =
      dst_mapped_box_level.getBoxes();

   for (hier::RealBoxConstIterator ni(dst_mapped_boxes);
        ni.isValid(); ++ni) {
      const hier::Box& dst_mapped_box = *ni;
      hier::Box fill_mapped_box(dst_mapped_box);
      fill_mapped_box.grow(fill_ghost_width);
      fill_mapped_boxes.addBoxWithoutUpdate(fill_mapped_box);
      dst_to_fill.insertLocalNeighbor(fill_mapped_box, dst_mapped_box.getId());
      TBOX_ASSERT(dst_to_fill.numLocalNeighbors(dst_mapped_box.getId()) == 1);
   }
   fill_mapped_boxes.finalize();
}

/*
 *************************************************************************
 *
 * computeDestinationFillBoxesOnSourceProc
 *
 *************************************************************************
 */

void PatchLevelFullFillPattern::computeDestinationFillBoxesOnSourceProc(
   FillSet& dst_fill_boxes_on_src_proc,
   const hier::BoxLevel& dst_mapped_box_level,
   const hier::Connector& src_to_dst,
   const hier::IntVector& fill_ghost_width)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(dst_mapped_box_level, fill_ghost_width);

   const tbox::Dimension& dim(fill_ghost_width.getDim());

   /*
    * src_to_dst initialized only when there is a src mapped_box_level.
    * Without the src mapped_box_level, we do not need to compute
    * dst_fill_boxes_on_src_proc.
    *
    * For PatchLevelFullFillPattern, the src owner can compute fill boxes
    * for all its dst neighbors using local data.  This info is
    * stored in dst_fill_boxes_on_src_proc.
    */
   hier::BoxSet tmp_nabrs(dim), all_dst_nabrs(dim);
   src_to_dst.getLocalNeighbors(tmp_nabrs);
   tmp_nabrs.unshiftPeriodicImageBoxes(
      all_dst_nabrs,
      dst_mapped_box_level.getRefinementRatio());
   tmp_nabrs.clear();
   for (hier::BoxSet::OrderedConstIterator na = all_dst_nabrs.orderedBegin();
        na != all_dst_nabrs.orderedEnd(); ++na) {
      hier::BoxSet& fill_boxes = dst_fill_boxes_on_src_proc.getNeighborSet(na->getId(), dim);
      hier::Box fill_box(*na);
      fill_box.grow(fill_ghost_width);
      fill_boxes.insert(fill_box);
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
