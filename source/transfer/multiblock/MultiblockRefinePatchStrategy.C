//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/transfer/multiblock/MultiblockRefinePatchStrategy.C $
// Package:	SAMRAI multiblock package
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Strategy interface to user routines for refining AMR data.
//

#ifndef included_xfer_MultiblockRefinePatchStrategy_C
#define included_xfer_MultiblockRefinePatchStrategy_C

#include "MultiblockRefinePatchStrategy.h"

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

template<int DIM>  MultiblockRefinePatchStrategy<DIM>::MultiblockRefinePatchStrategy()
{
   d_filling_coarse_scratch = false;
   d_block_number = MULTIBLOCK_UNDEFINED_BLOCK_NUMBER;
}

template<int DIM>  MultiblockRefinePatchStrategy<DIM>::~MultiblockRefinePatchStrategy()
{
}

template<int DIM>
void MultiblockRefinePatchStrategy<DIM>::setPhysicalBoundaryConditions(
   hier::Patch<DIM>& patch,
   const double fill_time,
   const hier::IntVector<DIM>& ghost_width_to_fill)
                                                                                
{
   NULL_USE(patch);
   NULL_USE(fill_time);
   NULL_USE(ghost_width_to_fill);
}

}
}
#endif
