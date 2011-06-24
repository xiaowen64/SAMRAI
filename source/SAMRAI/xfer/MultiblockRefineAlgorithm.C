/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Base class for geometry management on patches 
 *
 ************************************************************************/

#ifndef included_xfer_MultiblockRefineAlgorithm_C
#define included_xfer_MultiblockRefineAlgorithm_C

#include "SAMRAI/xfer/MultiblockRefineAlgorithm.h"

#include "SAMRAI/xfer/PatchLevelFullFillPattern.h"
#include "SAMRAI/xfer/StandardRefineTransactionFactory.h"

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

MultiblockRefineAlgorithm::MultiblockRefineAlgorithm(
   tbox::Pointer<hier::PatchHierarchy> hierarchy,
   const tbox::Dimension& dim):
   d_single_block_refine_alg(new xfer::RefineAlgorithm(dim)),
   d_multiblock_hierarchy(hierarchy),
   d_enable_singularity_patches(true)
{
}

MultiblockRefineAlgorithm::~MultiblockRefineAlgorithm()
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

tbox::Pointer<MultiblockRefineSchedule>
MultiblockRefineAlgorithm::createSchedule(
   tbox::Pointer<hier::PatchLevel> level,
   MultiblockRefinePatchStrategy* refine_strategy,
   tbox::Pointer<xfer::RefineTransactionFactory> transaction_factory)
const
{

   tbox::Pointer<xfer::RefineTransactionFactory> trans_factory =
      transaction_factory;

   if (trans_factory.isNull()) {
      trans_factory = new xfer::StandardRefineTransactionFactory;
   }

   tbox::Pointer<PatchLevelFullFillPattern> fill_pattern(
      new PatchLevelFullFillPattern());

   return tbox::Pointer<MultiblockRefineSchedule>(new MultiblockRefineSchedule(
                                                     fill_pattern,
                                                     level,
                                                     level,
                                                     d_multiblock_hierarchy,
                                                     d_single_block_refine_alg,
                                                     trans_factory,
                                                     refine_strategy,
                                                     false, 
                                                     d_enable_singularity_patches));
}

tbox::Pointer<MultiblockRefineSchedule>
MultiblockRefineAlgorithm::createSchedule(
   tbox::Pointer<PatchLevelFillPattern> fill_pattern,
   tbox::Pointer<hier::PatchLevel> level,
   MultiblockRefinePatchStrategy* refine_strategy,
   tbox::Pointer<xfer::RefineTransactionFactory> transaction_factory)
const
{

   tbox::Pointer<xfer::RefineTransactionFactory> trans_factory =
      transaction_factory;

   if (trans_factory.isNull()) {
      trans_factory = new xfer::StandardRefineTransactionFactory;
   }

   return tbox::Pointer<MultiblockRefineSchedule>(new MultiblockRefineSchedule(
                                                     fill_pattern,
                                                     level,
                                                     level,
                                                     d_multiblock_hierarchy,
                                                     d_single_block_refine_alg,
                                                     trans_factory,
                                                     refine_strategy,
                                                     false,
                                                     d_enable_singularity_patches));
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

tbox::Pointer<MultiblockRefineSchedule>
MultiblockRefineAlgorithm::createSchedule(
   tbox::Pointer<hier::PatchLevel> dst_level,
   tbox::Pointer<hier::PatchLevel> src_level,
   MultiblockRefinePatchStrategy* refine_strategy,
   tbox::Pointer<xfer::RefineTransactionFactory> transaction_factory)
const
{
   TBOX_ASSERT(!dst_level.isNull());
   TBOX_ASSERT(!src_level.isNull());
#ifdef DEBUG_CHECK_DIM_ASSERTIONS
   TBOX_DIM_ASSERT_CHECK_ARGS2(*dst_level, *src_level);
   if (refine_strategy) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*dst_level, *refine_strategy);
   }
#endif

   tbox::Pointer<xfer::RefineTransactionFactory> trans_factory =
      transaction_factory;

   if (trans_factory.isNull()) {
      trans_factory = new xfer::StandardRefineTransactionFactory;
   }

   tbox::Pointer<PatchLevelFullFillPattern> fill_pattern(
      new PatchLevelFullFillPattern());

   return tbox::Pointer<MultiblockRefineSchedule>(new MultiblockRefineSchedule(
                                                     fill_pattern,
                                                     dst_level,
                                                     src_level,
                                                     d_multiblock_hierarchy,
                                                     d_single_block_refine_alg,
                                                     trans_factory,
                                                     refine_strategy,
                                                     false,
                                                     d_enable_singularity_patches));
}

tbox::Pointer<MultiblockRefineSchedule>
MultiblockRefineAlgorithm::createSchedule(
   tbox::Pointer<PatchLevelFillPattern> fill_pattern,
   tbox::Pointer<hier::PatchLevel> dst_level,
   tbox::Pointer<hier::PatchLevel> src_level,
   MultiblockRefinePatchStrategy* refine_strategy,
   tbox::Pointer<xfer::RefineTransactionFactory> transaction_factory)
const
{
   TBOX_ASSERT(!dst_level.isNull());
   TBOX_ASSERT(!src_level.isNull());
#ifdef DEBUG_CHECK_DIM_ASSERTIONS
   TBOX_DIM_ASSERT_CHECK_ARGS2(*dst_level, *src_level);
   if (refine_strategy) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*dst_level, *refine_strategy);
   }
#endif

   tbox::Pointer<xfer::RefineTransactionFactory> trans_factory =
      transaction_factory;

   if (trans_factory.isNull()) {
      trans_factory = new xfer::StandardRefineTransactionFactory;
   }

   return tbox::Pointer<MultiblockRefineSchedule>(new MultiblockRefineSchedule(
                                                     fill_pattern,
                                                     dst_level,
                                                     src_level,
                                                     d_multiblock_hierarchy,
                                                     d_single_block_refine_alg,
                                                     trans_factory,
                                                     refine_strategy,
                                                     false,
                                                     d_enable_singularity_patches));
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

tbox::Pointer<MultiblockRefineSchedule>
MultiblockRefineAlgorithm::createSchedule(
   tbox::Pointer<hier::PatchLevel> level,
   const int next_coarser_level,
   tbox::Pointer<hier::PatchHierarchy> hierarchy,
   MultiblockRefinePatchStrategy* refine_strategy,
   tbox::Pointer<xfer::RefineTransactionFactory> transaction_factory)
const
{
#ifdef DEBUG_CHECK_DIM_ASSERTIONS
   TBOX_DIM_ASSERT_CHECK_ARGS2(*level, *hierarchy);
   if (refine_strategy) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*level, *refine_strategy);
   }
#endif

   tbox::Pointer<xfer::RefineTransactionFactory> trans_factory =
      transaction_factory;

   if (trans_factory.isNull()) {
      trans_factory = new xfer::StandardRefineTransactionFactory;
   }

   tbox::Pointer<PatchLevelFullFillPattern> fill_pattern(
      new PatchLevelFullFillPattern());

   return tbox::Pointer<MultiblockRefineSchedule>(new MultiblockRefineSchedule(
                                                     fill_pattern,
                                                     level,
                                                     level,
                                                     next_coarser_level,
                                                     hierarchy,
                                                     d_single_block_refine_alg,
                                                     trans_factory,
                                                     refine_strategy,
                                                     false,
                                                     d_enable_singularity_patches));
}

tbox::Pointer<MultiblockRefineSchedule>
MultiblockRefineAlgorithm::createSchedule(
   tbox::Pointer<PatchLevelFillPattern> fill_pattern,
   tbox::Pointer<hier::PatchLevel> level,
   const int next_coarser_level,
   tbox::Pointer<hier::PatchHierarchy> hierarchy,
   MultiblockRefinePatchStrategy* refine_strategy,
   tbox::Pointer<xfer::RefineTransactionFactory> transaction_factory)
const
{
#ifdef DEBUG_CHECK_DIM_ASSERTIONS
   TBOX_DIM_ASSERT_CHECK_ARGS2(*level, *hierarchy);
   if (refine_strategy) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*level, *refine_strategy);
   }
#endif

   tbox::Pointer<xfer::RefineTransactionFactory> trans_factory =
      transaction_factory;

   if (trans_factory.isNull()) {
      trans_factory = new xfer::StandardRefineTransactionFactory;
   }

   return tbox::Pointer<MultiblockRefineSchedule>(new MultiblockRefineSchedule(
                                                     fill_pattern,
                                                     level,
                                                     level,
                                                     next_coarser_level,
                                                     hierarchy,
                                                     d_single_block_refine_alg,
                                                     trans_factory,
                                                     refine_strategy,
                                                     false, 
                                                     d_enable_singularity_patches));
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

tbox::Pointer<MultiblockRefineSchedule>
MultiblockRefineAlgorithm::createSchedule(
   tbox::Pointer<hier::PatchLevel> dst_level,
   tbox::Pointer<hier::PatchLevel> src_level,
   const int next_coarser_level,
   tbox::Pointer<hier::PatchHierarchy> hierarchy,
   MultiblockRefinePatchStrategy* refine_strategy,
   tbox::Pointer<xfer::RefineTransactionFactory> transaction_factory)
const
{
#ifdef DEBUG_CHECK_DIM_ASSERTIONS
   if (!src_level.isNull()) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*dst_level, *src_level);
   }
   if (refine_strategy) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*dst_level, *refine_strategy);
   }
#endif

   tbox::Pointer<xfer::RefineTransactionFactory> trans_factory =
      transaction_factory;

   if (trans_factory.isNull()) {
      trans_factory = new xfer::StandardRefineTransactionFactory;
   }

   tbox::Pointer<PatchLevelFullFillPattern> fill_pattern(
      new PatchLevelFullFillPattern());

   return tbox::Pointer<MultiblockRefineSchedule>(new MultiblockRefineSchedule(
                                                     fill_pattern,
                                                     dst_level,
                                                     src_level,
                                                     next_coarser_level,
                                                     hierarchy,
                                                     d_single_block_refine_alg,
                                                     trans_factory,
                                                     refine_strategy,
                                                     false,
                                                     d_enable_singularity_patches));
}

tbox::Pointer<MultiblockRefineSchedule>
MultiblockRefineAlgorithm::createSchedule(
   tbox::Pointer<PatchLevelFillPattern> fill_pattern,
   tbox::Pointer<hier::PatchLevel> dst_level,
   tbox::Pointer<hier::PatchLevel> src_level,
   const int next_coarser_level,
   tbox::Pointer<hier::PatchHierarchy> hierarchy,
   MultiblockRefinePatchStrategy* refine_strategy,
   tbox::Pointer<xfer::RefineTransactionFactory> transaction_factory)
const
{
#ifdef DEBUG_CHECK_DIM_ASSERTIONS
   if (!src_level.isNull()) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*dst_level, *src_level);
   }
   if (!hierarchy.isNull()) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*dst_level, *hierarchy);
   }
   if (refine_strategy) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*dst_level, *refine_strategy);
   }
#endif

   tbox::Pointer<xfer::RefineTransactionFactory> trans_factory =
      transaction_factory;

   if (trans_factory.isNull()) {
      trans_factory = new xfer::StandardRefineTransactionFactory;
   }

   return tbox::Pointer<MultiblockRefineSchedule>(new MultiblockRefineSchedule(
                                                     fill_pattern,
                                                     dst_level,
                                                     src_level,
                                                     next_coarser_level,
                                                     hierarchy,
                                                     d_single_block_refine_alg,
                                                     trans_factory,
                                                     refine_strategy,
                                                     false,
                                                     d_enable_singularity_patches));
}

/*
 *************************************************************************
 *                                                                       *
 * Register a refinement operation with the algorithm                    *
 *                                                                       *
 *************************************************************************
 */

void MultiblockRefineAlgorithm::registerRefine(
   const int dst,
   const int src,
   const int scratch,
   tbox::Pointer<hier::RefineOperator> oprefine,
   tbox::Pointer<VariableFillPattern> var_fill_pattern)
{
   d_single_block_refine_alg->registerRefine(dst, src, scratch, oprefine,
                                             var_fill_pattern);
}

void MultiblockRefineAlgorithm::registerRefine(
   const int dst,
   const int src,
   const int src_told,
   const int src_tnew,
   const int scratch,
   tbox::Pointer<hier::RefineOperator> oprefine,
   tbox::Pointer<hier::TimeInterpolateOperator> optime,
   tbox::Pointer<VariableFillPattern> var_fill_pattern)
{
   d_single_block_refine_alg->registerRefine(dst, src, src_told, src_tnew,
      scratch, oprefine, optime, var_fill_pattern);
}

void MultiblockRefineAlgorithm::resetSchedule(
   tbox::Pointer<xfer::MultiblockRefineSchedule> schedule) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!schedule.isNull());
#endif
   if (d_single_block_refine_alg->getEquivalenceClasses()->
       checkConsistency(schedule->getEquivalenceClasses())) {
      schedule->reset(d_single_block_refine_alg);
   } else {
      TBOX_ERROR("MultiblockRefineAlgorithm::resetSchedule error..."
         << "\n Items in schedule passed to resetSchedule are"
         << "\n inconsistent with those in the MultiblockRefineAlgorithm."
         << std::endl);
   }
}

}
}
#endif
