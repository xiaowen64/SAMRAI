/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Strategy interface to user routines for refining AMR data. 
 *
 ************************************************************************/

#ifndef included_xfer_MultiblockRefinePatchStrategy_C
#define included_xfer_MultiblockRefinePatchStrategy_C

#include "SAMRAI/xfer/MultiblockRefinePatchStrategy.h"

namespace SAMRAI {
namespace xfer {

/*
 *************************************************************************
 *									*
 * The default constructor and virtual destructor do nothing             *
 * particularly interesting.		                                *
 *									*
 *************************************************************************
 */

MultiblockRefinePatchStrategy::MultiblockRefinePatchStrategy(
   const tbox::Dimension& dim):
   xfer::RefinePatchStrategy(dim)
{
   d_filling_coarse_scratch = false;
}

MultiblockRefinePatchStrategy::~MultiblockRefinePatchStrategy()
{
}

void MultiblockRefinePatchStrategy::setPhysicalBoundaryConditions(
   hier::Patch& patch,
   const double fill_time,
   const hier::IntVector& ghost_width_to_fill)

{
   TBOX_DIM_ASSERT_CHECK_ARGS3(*this, patch, ghost_width_to_fill);

   NULL_USE(patch);
   NULL_USE(fill_time);
   NULL_USE(ghost_width_to_fill);
}

}
}
#endif
