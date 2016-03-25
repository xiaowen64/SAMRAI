//
// File:        MultiblockRefineAlgorithm.C
// Package:     SAMRAI multiblock package
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 653 $
// Modified:    $Date: 2005-10-06 16:28:14 -0700 (Thu, 06 Oct 2005) $
// Description: Base class for geometry management on patches
//

#ifndef included_mblk_MultiblockRefineAlgorithm_C
#define included_mblk_MultiblockRefineAlgorithm_C

#include "MultiblockRefineAlgorithm.h"

namespace SAMRAI {
   namespace mblk {

/*
*************************************************************************
*                                                                       *
* Constructor and destructor for multiblock refine algorithm.  The      *
* constructor simply initializes the refine algorithm data member.      *
*                                                                       *
*************************************************************************
*/

template<int DIM>
MultiblockRefineAlgorithm<DIM>::MultiblockRefineAlgorithm(
   tbox::Pointer< xfer::RefineAlgorithm<DIM> > refine_alg,
   tbox::Pointer< MultiblockPatchHierarchy<DIM> > multiblock)
{
   d_single_block_refine_alg = refine_alg;
   d_multiblock_hierarchy = multiblock;
}

template<int DIM>
MultiblockRefineAlgorithm<DIM>::~MultiblockRefineAlgorithm()
{
}

/*
*************************************************************************
*                                                                       *
* Create a communication schedule that will copy data from the          *
* interiors of the specified level into the ghost cells and             *
* interiors of the same level.                                          *
*                                                                       *
*************************************************************************
*/

template<int DIM>
tbox::Pointer< MultiblockRefineSchedule<DIM> > 
MultiblockRefineAlgorithm<DIM>::createSchedule(
   tbox::Pointer< MultiblockPatchLevel<DIM> > level,
   MultiblockRefinePatchStrategy<DIM>* refine_strategy) const
{
   return(new MultiblockRefineSchedule<DIM>(level,
                                             level,
                                             d_multiblock_hierarchy,
                                             d_single_block_refine_alg,
                                                refine_strategy));

}

/*
*************************************************************************
*                                                                       *
* Create a communication schedule that will copy data from the          *
* interiors of the source level into the ghost cell and interiors       *
* of the destination level.                                             *
*                                                                       *
*************************************************************************
*/

template<int DIM>
tbox::Pointer< MultiblockRefineSchedule<DIM> >
MultiblockRefineAlgorithm<DIM>::createSchedule(
   tbox::Pointer< MultiblockPatchLevel<DIM> > dst_level,
   tbox::Pointer< MultiblockPatchLevel<DIM> > src_level,
   MultiblockRefinePatchStrategy<DIM>* refine_strategy) const
{
   return(new MultiblockRefineSchedule<DIM>(dst_level,
                                             src_level,
                                             d_multiblock_hierarchy,
                                             d_single_block_refine_alg,
                                                refine_strategy));

}

/*
*************************************************************************
*                                                                       *
* Create a communication schedule that copies data from the interiors   *
* of the same level and coarser levels into the interior and boundary   *
* cells of the given level.                                             *
*                                                                       *
*************************************************************************
*/

template<int DIM>
tbox::Pointer< MultiblockRefineSchedule<DIM> > 
MultiblockRefineAlgorithm<DIM>::createSchedule(
   tbox::Pointer< MultiblockPatchLevel<DIM> > level,
   const int next_coarser_level,
   tbox::Pointer< MultiblockPatchHierarchy<DIM> > multiblock,
   MultiblockRefinePatchStrategy<DIM>* refine_strategy) const
{
   return(new MultiblockRefineSchedule<DIM>(level,
                                             level,
                                             next_coarser_level,
                                             multiblock,
                                             d_single_block_refine_alg,
                                                refine_strategy));
}

/*
*************************************************************************
*                                                                       *
* Create a communication schedule that copies data from the interiors   *
* of the old level and coarser levels into the ghost cells and interior *
* cells of the given new level.                                         *
*                                                                       *
*************************************************************************
*/

template<int DIM>
tbox::Pointer< MultiblockRefineSchedule<DIM> >
MultiblockRefineAlgorithm<DIM>::createSchedule(
   tbox::Pointer< MultiblockPatchLevel<DIM> > dst_level,
   tbox::Pointer< MultiblockPatchLevel<DIM> > src_level,
   const int next_coarser_level,
   tbox::Pointer< MultiblockPatchHierarchy<DIM> > multiblock,
   MultiblockRefinePatchStrategy<DIM>* refine_strategy) const
{
   return(new MultiblockRefineSchedule<DIM>(dst_level,
                                             src_level,
                                             next_coarser_level,
                                             multiblock,
                                             d_single_block_refine_alg,
                                                refine_strategy));
}

/*
*************************************************************************
*                                                                       *
* Register a refinement operation with the algorithm                    *
*                                                                       *
*************************************************************************
*/

template<int DIM>
void MultiblockRefineAlgorithm<DIM>::registerRefine(
   const int dst,
   const int src,
   const int scratch,
   tbox::Pointer< xfer::RefineOperator<DIM> > oprefine)
{
   d_single_block_refine_alg->registerRefine(dst, src, scratch, oprefine);
}

template<int DIM>
void MultiblockRefineAlgorithm<DIM>::registerRefine(
   const int dst,
   const int src,
   const int src_told,
   const int src_tnew,
   const int scratch,
   tbox::Pointer< xfer::RefineOperator<DIM> > oprefine,
   tbox::Pointer< xfer::TimeInterpolateOperator<DIM> > optime)
{
   d_single_block_refine_alg->registerRefine(dst, src, src_told, src_tnew,
                                scratch, oprefine, optime);
}


}
}
#endif
