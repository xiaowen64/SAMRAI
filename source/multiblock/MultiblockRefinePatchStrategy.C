//
// File:	MultiblockRefinePatchStrategy.C
// Package:	SAMRAI multiblock package
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 653 $
// Modified:	$Date: 2005-10-06 16:28:14 -0700 (Thu, 06 Oct 2005) $
// Description:	Strategy interface to user routines for refining AMR data.
//

#ifndef included_mblk_MultiblockRefinePatchStrategy_C
#define included_mblk_MultiblockRefinePatchStrategy_C

#include "MultiblockRefinePatchStrategy.h"

namespace SAMRAI {
    namespace mblk {

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
