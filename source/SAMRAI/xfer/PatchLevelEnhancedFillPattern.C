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
   NULL_USE(dst_mapped_box_level);
   NULL_USE(dst_to_dst);
   NULL_USE(dst_to_src);
   NULL_USE(src_to_dst);
   NULL_USE(fill_ghost_width);
   NULL_USE(fill_mapped_boxes);
   NULL_USE(dst_to_fill_edges);

   TBOX_ERROR("computeFillMappedBoxesAndNeighborhoodSets() should not be called for PatchLevelEnhancedFillPattern.");

}

}
}
#endif
