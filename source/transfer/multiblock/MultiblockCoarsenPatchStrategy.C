//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/transfer/multiblock/MultiblockCoarsenPatchStrategy.C $
// Package:	SAMRAI data transfer
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
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
