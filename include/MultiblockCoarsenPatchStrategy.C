//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/transfer/multiblock/MultiblockCoarsenPatchStrategy.C $
// Package:	SAMRAI data transfer
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	Strategy interface to user routines for coarsening AMR data.
//
 
#ifndef included_xfer_MultiblockCoarsenPatchStrategy_C
#define included_xfer_MultiblockCoarsenPatchStrategy_C

#include "MultiblockCoarsenPatchStrategy.h"

namespace SAMRAI {
    namespace xfer {

template<int DIM>  MultiblockCoarsenPatchStrategy<DIM>::MultiblockCoarsenPatchStrategy()
{
   d_block_number = MULTIBLOCK_UNDEFINED_BLOCK_NUMBER;
}

template<int DIM>  MultiblockCoarsenPatchStrategy<DIM>::~MultiblockCoarsenPatchStrategy()
{
}
 
}
}
#endif
