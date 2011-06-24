/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Abstract fill pattern class to provide interface for stencils 
 *
 ************************************************************************/

#ifndef included_xfer_PatchLevelFullFillPattern_C
#define included_xfer_PatchLevelFullFillPattern_C

#include "SAMRAI/xfer/PatchLevelFullFillPattern.h"
#include "SAMRAI/hier/RealMappedBoxConstIterator.h"
#include "SAMRAI/hier/MappedBox.h"
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
 *                                                                       *
 * Default constructor                                                   *
 *                                                                       *
 *************************************************************************
 */

PatchLevelFullFillPattern::PatchLevelFullFillPattern():
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

PatchLevelFullFillPattern::~PatchLevelFullFillPattern()
{
}

/*
 *************************************************************************
 *                                                                       *
 * computeFillMappedBoxesAndNeighborhoodSets                                        *
 *                                                                       *
 *************************************************************************
 */

void PatchLevelFullFillPattern::computeFillMappedBoxesAndNeighborhoodSets(
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

   const hier::MappedBoxSet& dst_mapped_boxes =
      dst_mapped_box_level.getMappedBoxes();

   for (hier::RealMappedBoxConstIterator ni(dst_mapped_boxes);
        ni.isValid(); ++ni) {
      const hier::MappedBox& dst_mapped_box = *ni;
      hier::MappedBox fill_mapped_box(dst_mapped_box);
      fill_mapped_box.getBox().grow(fill_ghost_width);
      fill_mapped_boxes.insert(fill_mapped_boxes.end(), fill_mapped_box);
      dst_to_fill_edges[dst_mapped_box.getId()].insert(fill_mapped_box);
      TBOX_ASSERT(dst_to_fill_edges[dst_mapped_box.getId()].size() == 1);
   }
}

/*
 *************************************************************************
 *                                                                       *
 * computeDestinationFillBoxesOnSourceProc                               *
 *                                                                       *
 *************************************************************************
 */

void PatchLevelFullFillPattern::computeDestinationFillBoxesOnSourceProc(
   FillSet& dst_fill_boxes_on_src_proc,
   const hier::MappedBoxLevel& dst_mapped_box_level,
   const hier::Connector& src_to_dst,
   const hier::IntVector& fill_ghost_width)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(dst_mapped_box_level, fill_ghost_width);

   /*
    * src_to_dst initialized only when there is a src mapped_box_level.
    * Without the src mapped_box_level, we do not need to compute
    * dst_fill_boxes_on_src_proc.
    *
    * For PatchLevelFullFillPattern, the src owner can compute fill boxes
    * for all its dst neighbors using local data.  This info is
    * stored in dst_fill_boxes_on_src_proc.
    */
   NeighborSet tmp_nabrs, all_dst_nabrs;
   src_to_dst.getNeighborhoodSets().getNeighbors(tmp_nabrs);
   tmp_nabrs.unshiftPeriodicImageMappedBoxes(
      all_dst_nabrs,
      dst_mapped_box_level.getRefinementRatio());
   tmp_nabrs.clear();
   for (NeighborSet::const_iterator na = all_dst_nabrs.begin();
        na != all_dst_nabrs.end(); ++na) {
      MappedBoxVector& fill_boxes = dst_fill_boxes_on_src_proc[na->getGlobalId()];
      hier::MappedBox fill_box(*na);
      fill_box.getBox().grow(fill_ghost_width);
      fill_boxes.push_back(fill_box);
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
