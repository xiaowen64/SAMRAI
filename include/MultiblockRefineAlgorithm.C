//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/transfer/multiblock/MultiblockRefineAlgorithm.C $
// Package:     SAMRAI multiblock package
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1917 $
// Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description: Base class for geometry management on patches
//

#ifndef included_xfer_MultiblockRefineAlgorithm_C
#define included_xfer_MultiblockRefineAlgorithm_C

#include "MultiblockRefineAlgorithm.h"
#include "StandardRefineTransactionFactory.h"

namespace SAMRAI {
   namespace xfer {

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
   tbox::Pointer< hier::MultiblockPatchHierarchy<DIM> > multiblock)
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
   tbox::Pointer< hier::MultiblockPatchLevel<DIM> > level,
   MultiblockRefinePatchStrategy<DIM>* refine_strategy,
   tbox::Pointer< xfer::RefineTransactionFactory<DIM> > transaction_factory)
   const 
{

   tbox::Pointer< xfer::RefineTransactionFactory<DIM> > trans_factory =
      transaction_factory;

   if (trans_factory.isNull()) {
      trans_factory = new xfer::StandardRefineTransactionFactory<DIM>;
   }

   return(new MultiblockRefineSchedule<DIM>("DEFAULT_FILL",
                                            level,
                                            level,
                                            d_multiblock_hierarchy,
                                            d_single_block_refine_alg,
                                            trans_factory,
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
   tbox::Pointer< hier::MultiblockPatchLevel<DIM> > dst_level,
   tbox::Pointer< hier::MultiblockPatchLevel<DIM> > src_level,
   MultiblockRefinePatchStrategy<DIM>* refine_strategy,
   tbox::Pointer< xfer::RefineTransactionFactory<DIM> > transaction_factory)
   const
{

   tbox::Pointer< xfer::RefineTransactionFactory<DIM> > trans_factory =
      transaction_factory;

   if (trans_factory.isNull()) {
      trans_factory = new xfer::StandardRefineTransactionFactory<DIM>;
   }

   return(new MultiblockRefineSchedule<DIM>("DEFAULT_FILL",
                                            dst_level,
                                            src_level,
                                            d_multiblock_hierarchy,
                                            d_single_block_refine_alg,
                                            trans_factory,
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
   tbox::Pointer< hier::MultiblockPatchLevel<DIM> > level,
   const int next_coarser_level,
   tbox::Pointer< hier::MultiblockPatchHierarchy<DIM> > multiblock,
   MultiblockRefinePatchStrategy<DIM>* refine_strategy,                         
   tbox::Pointer< xfer::RefineTransactionFactory<DIM> > transaction_factory)
   const
{

   tbox::Pointer< xfer::RefineTransactionFactory<DIM> > trans_factory =
      transaction_factory;

   if (trans_factory.isNull()) {
      trans_factory = new xfer::StandardRefineTransactionFactory<DIM>;
   }

   return(new MultiblockRefineSchedule<DIM>("DEFAULT_FILL",
                                            level,
                                            level,
                                            next_coarser_level,
                                            multiblock,
                                            d_single_block_refine_alg,
                                            trans_factory,
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
   tbox::Pointer< hier::MultiblockPatchLevel<DIM> > dst_level,
   tbox::Pointer< hier::MultiblockPatchLevel<DIM> > src_level,
   const int next_coarser_level,
   tbox::Pointer< hier::MultiblockPatchHierarchy<DIM> > multiblock,
   MultiblockRefinePatchStrategy<DIM>* refine_strategy,                         
   tbox::Pointer< xfer::RefineTransactionFactory<DIM> > transaction_factory)
   const
{

   tbox::Pointer< xfer::RefineTransactionFactory<DIM> > trans_factory =
      transaction_factory;

   if (trans_factory.isNull()) {
      trans_factory = new xfer::StandardRefineTransactionFactory<DIM>;
   }

   return(new MultiblockRefineSchedule<DIM>("DEFAULT_FILL",
                                            dst_level,
                                            src_level,
                                            next_coarser_level,
                                            multiblock,
                                            d_single_block_refine_alg,
                                            trans_factory,
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
