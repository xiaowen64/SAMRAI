/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Abstract fill pattern class to provide interface for stencils
 *
 ************************************************************************/

#ifndef included_xfer_PatchLevelFillPattern_C
#define included_xfer_PatchLevelFillPattern_C

#include "SAMRAI/xfer/PatchLevelFillPattern.h"

namespace SAMRAI {
namespace xfer {

/*
 *************************************************************************
 *
 * Default constructor
 *
 *************************************************************************
 */

PatchLevelFillPattern::PatchLevelFillPattern()
{
}

/*
 *************************************************************************
 *
 * Destructor
 *
 *************************************************************************
 */

PatchLevelFillPattern::~PatchLevelFillPattern()
{
}

/*
 *************************************************************************
 * Default computeDestinationFillBoxesOnSourceProc() is a no-op.
 * A concrete implementation should only be required if
 * needsToCommunicateDestinationFillBoxes() returns false.
 *************************************************************************
 */

void PatchLevelFillPattern::computeDestinationFillBoxesOnSourceProc(
   FillSet& dst_fill_boxes_on_src_proc,
   const hier::BoxLevel& dst_mapped_box_level,
   const hier::Connector& src_to_dst,
   const hier::IntVector& fill_ghost_width)
{
   NULL_USE(dst_mapped_box_level);
   NULL_USE(src_to_dst);
   NULL_USE(fill_ghost_width);
   NULL_USE(dst_fill_boxes_on_src_proc);
   if (!needsToCommunicateDestinationFillBoxes()) {
      TBOX_ERROR(
         "The concrete PatchLevelFillPattern::computeDestinationFillBoxesOnSourceProc:\n"
         << "must be implemented whenever the concrete\n"
         << "method needsToCommunicateDestinationFillBoxes() returns false.");
   }
}

}
}
#endif
