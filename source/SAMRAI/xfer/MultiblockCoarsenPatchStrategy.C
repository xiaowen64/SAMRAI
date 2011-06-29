/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Strategy interface to user routines for coarsening AMR data. 
 *
 ************************************************************************/

#ifndef included_xfer_MultiblockCoarsenPatchStrategy_C
#define included_xfer_MultiblockCoarsenPatchStrategy_C

#include "SAMRAI/xfer/MultiblockCoarsenPatchStrategy.h"

namespace SAMRAI {
namespace xfer {

  const int MultiblockCoarsenPatchStrategy::MULTIBLOCK_UNDEFINED_BLOCK_NUMBER = -1;

MultiblockCoarsenPatchStrategy::MultiblockCoarsenPatchStrategy(
   const tbox::Dimension& dim):
   xfer::CoarsenPatchStrategy(dim)
{
   d_block_number = MULTIBLOCK_UNDEFINED_BLOCK_NUMBER;
}

MultiblockCoarsenPatchStrategy::~MultiblockCoarsenPatchStrategy()
{
}

}
}
#endif
