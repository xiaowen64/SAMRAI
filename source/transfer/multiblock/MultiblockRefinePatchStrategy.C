//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/transfer/multiblock/MultiblockRefinePatchStrategy.C $
// Package:	SAMRAI multiblock package
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
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
