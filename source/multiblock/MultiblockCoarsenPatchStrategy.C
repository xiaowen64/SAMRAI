//
// File:	MultiblockCoarsenPatchStrategy.C
// Package:	SAMRAI data transfer
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 713 $
// Modified:	$Date: 2005-11-08 13:48:27 -0800 (Tue, 08 Nov 2005) $
// Description:	Strategy interface to user routines for coarsening AMR data.
//
 
#ifndef included_mblk_MultiblockCoarsenPatchStrategy_C
#define included_mblk_MultiblockCoarsenPatchStrategy_C

#include "MultiblockCoarsenPatchStrategy.h"

namespace SAMRAI {
    namespace mblk {

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
