/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Coarsening algorithm for data transfer between AMR levels 
 *
 ************************************************************************/

#ifndef included_xfer_MultiblockCoarsenAlgorithm_C
#define included_xfer_MultiblockCoarsenAlgorithm_C

#include "SAMRAI/xfer/MultiblockCoarsenAlgorithm.h"

#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/xfer/StandardCoarsenTransactionFactory.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace xfer {

/*
 *************************************************************************
 *                                                                       *
 * The constructor creates a new xfer::CoarsenClasses object             *
 * and caches a boolean indiating whether to copy data to the            *
 * destination space on the coarse level before coarsening.              *
 *                                                                       *
 *************************************************************************
 */

MultiblockCoarsenAlgorithm::MultiblockCoarsenAlgorithm(
   tbox::Pointer<hier::PatchHierarchy> hierarchy,
   bool fill_coarse_data):
   d_coarsen_classes(new xfer::CoarsenClasses(d_fill_coarse_data)),
   d_hierarchy(hierarchy),
   d_fill_coarse_data(fill_coarse_data)
{
}

/*
 *************************************************************************
 *									*
 * The destructor implicitly deallocates the list data.                  *
 *									*
 *************************************************************************
 */

MultiblockCoarsenAlgorithm::~MultiblockCoarsenAlgorithm()
{
}

/*
 *************************************************************************
 *                                                                       *
 * Register a coarsening operation with the coarsening algorithm.        *
 *                                                                       *
 *************************************************************************
 */

void MultiblockCoarsenAlgorithm::registerCoarsen(
   const int dst,
   const int src,
   const tbox::Pointer<hier::CoarsenOperator> opcoarsen,
   const hier::IntVector& gcw_to_coarsen,
   tbox::Pointer<VariableFillPattern> var_fill_pattern)
{
   const tbox::Dimension& dim(d_hierarchy->getDim());

   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(dim, gcw_to_coarsen);

   xfer::CoarsenClasses::Data data(dim);

   data.d_dst = dst;
   data.d_src = src;
   data.d_gcw_to_coarsen = gcw_to_coarsen;
   data.d_opcoarsen = opcoarsen;

   hier::VariableDatabase* var_db = hier::VariableDatabase::getDatabase();

   tbox::Pointer<hier::Variable> var;
   if (!var_db->mapIndexToVariable(dst, var)) {
      TBOX_ERROR("MultiblockCoarsenAlgorithm::registerCoarsen error..."
         << "\nNo variable associated with dst patch data index." << std::endl);
   }

   data.d_fine_bdry_reps_var = var->fineBoundaryRepresentsVariable();
   if (!(var_fill_pattern.isNull())) {
      data.d_var_fill_pattern = var_fill_pattern;
   } else {
      data.d_var_fill_pattern = new BoxGeometryVariableFillPattern();
   }

   d_coarsen_classes->insertEquivalenceClassItem(data);
}

void MultiblockCoarsenAlgorithm::registerCoarsen(
   const int dst,
   const int src,
   const tbox::Pointer<hier::CoarsenOperator> opcoarsen,
   tbox::Pointer<VariableFillPattern> var_fill_pattern)
{
   const tbox::Dimension& dim(d_hierarchy->getDim());

   registerCoarsen(dst,
      src,
      opcoarsen,
      hier::IntVector::getZero(dim),
      var_fill_pattern);
}

/*
 *************************************************************************
 *                                                                       *
 * Create a communication schedule that will coarsen data from fine      *
 * patch level to the coarse patch level.                                *
 *                                                                       *
 *************************************************************************
 */

tbox::Pointer<MultiblockCoarsenSchedule>
MultiblockCoarsenAlgorithm::createSchedule(
   tbox::Pointer<hier::PatchLevel> crse_level,
   tbox::Pointer<hier::PatchLevel> fine_level,
   MultiblockCoarsenPatchStrategy* patch_strategy,
   MultiblockRefinePatchStrategy* refine_strategy,
   tbox::Pointer<xfer::CoarsenTransactionFactory> transaction_factory) const
{
   tbox::Pointer<xfer::CoarsenTransactionFactory> trans_factory =
      transaction_factory;

   if (trans_factory.isNull()) {
      trans_factory = new StandardCoarsenTransactionFactory();
   }

   return tbox::Pointer<MultiblockCoarsenSchedule>(
             new MultiblockCoarsenSchedule(crse_level,
                fine_level,
                d_coarsen_classes,
                d_hierarchy,
                trans_factory,
                patch_strategy,
                refine_strategy,
                d_fill_coarse_data));
}

}
}
#endif
