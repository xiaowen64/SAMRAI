/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Coarsening algorithm for data transfer between AMR levels
 *
 ************************************************************************/

#ifndef included_xfer_CoarsenAlgorithm_C
#define included_xfer_CoarsenAlgorithm_C

#include "SAMRAI/xfer/CoarsenAlgorithm.h"

#include "SAMRAI/xfer/StandardCoarsenTransactionFactory.h"
#include "SAMRAI/hier/PatchDataFactory.h"
#include "SAMRAI/hier/PatchDescriptor.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace xfer {

/*
 *************************************************************************
 *
 * The constructor creates a new CoarsenClasses object
 * and caches a boolean indiating whether to copy data to the
 * destination space on the coarse level before coarsening.
 *
 *************************************************************************
 */

CoarsenAlgorithm::CoarsenAlgorithm(
   const tbox::Dimension& dim,
   bool fill_coarse_data):
   d_dim(dim),
   d_coarsen_classes(new xfer::CoarsenClasses(d_fill_coarse_data)),
   d_fill_coarse_data(fill_coarse_data),
   d_schedule_created(false)
{
}

/*
 *************************************************************************
 *
 * The destructor implicitly deallocates the list data.
 *
 *************************************************************************
 */

CoarsenAlgorithm::~CoarsenAlgorithm()
{
}

/*
 *************************************************************************
 *
 * Register a coarsening operation with the coarsening algorithm.
 *
 *************************************************************************
 */

void CoarsenAlgorithm::registerCoarsen(
   const int dst,
   const int src,
   const tbox::Pointer<hier::CoarsenOperator> opcoarsen,
   const hier::IntVector& gcw_to_coarsen,
   tbox::Pointer<VariableFillPattern> var_fill_pattern)
{
#ifdef DEBUG_CHECK_DIM_ASSERTIONS
   if (!opcoarsen.isNull()) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *opcoarsen);
   }
#endif

   if (d_schedule_created) {
      TBOX_ERROR(
         "CoarsenAlgorithm::registerCoarsen error..."
         << "\nCannot call registerCoarsen with this coarsen algorithm"
         << "\nobject since it has already been used to create a coarsen schedule."
         << std::endl);
   }

   xfer::CoarsenClasses::Data data(d_dim);

   data.d_dst = dst;
   data.d_src = src;
   data.d_fine_bdry_reps_var = hier::VariableDatabase::getDatabase()->
      getPatchDescriptor()->getPatchDataFactory(dst)->
      fineBoundaryRepresentsVariable();
   data.d_gcw_to_coarsen = gcw_to_coarsen;
   data.d_opcoarsen = opcoarsen;
   data.d_tag = -1;
   if (!(var_fill_pattern.isNull())) {
      data.d_var_fill_pattern = var_fill_pattern;
   } else {
      data.d_var_fill_pattern = new BoxGeometryVariableFillPattern();
   }

   d_coarsen_classes->insertEquivalenceClassItem(data);
}

void CoarsenAlgorithm::registerCoarsen(
   const int dst,
   const int src,
   const tbox::Pointer<hier::CoarsenOperator> opcoarsen,
   tbox::Pointer<VariableFillPattern> var_fill_pattern)
{
   registerCoarsen(dst, src, opcoarsen,
      hier::IntVector::getZero(d_dim), var_fill_pattern);
}

/*
 *************************************************************************
 *
 * Create a communication schedule that will coarsen data from fine
 * patch level to the coarse patch level.
 *
 *************************************************************************
 */

tbox::Pointer<xfer::CoarsenSchedule>
CoarsenAlgorithm::createSchedule(
   tbox::Pointer<hier::PatchLevel> crse_level,
   tbox::Pointer<hier::PatchLevel> fine_level,
   xfer::CoarsenPatchStrategy* patch_strategy,
   tbox::Pointer<xfer::CoarsenTransactionFactory> transaction_factory)
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS2(d_dim, *crse_level, *fine_level);

   d_schedule_created = true;

   tbox::Pointer<xfer::CoarsenTransactionFactory> trans_factory =
      transaction_factory;

   if (trans_factory.isNull()) {
      trans_factory = new xfer::StandardCoarsenTransactionFactory();
   }

   return tbox::Pointer<xfer::CoarsenSchedule>(new xfer::CoarsenSchedule(
                                                  crse_level,
                                                  fine_level,
                                                  d_coarsen_classes,
                                                  trans_factory,
                                                  patch_strategy,
                                                  d_fill_coarse_data));
}

/*
 **************************************************************************
 *
 * Reconfigure coarsen schedule to perform operations in this algorithm.
 *
 **************************************************************************
 */

bool CoarsenAlgorithm::checkConsistency(
   tbox::Pointer<xfer::CoarsenSchedule> schedule) const
{
   TBOX_ASSERT(!schedule.isNull());

   return d_coarsen_classes->
          classesMatch(schedule->getEquivalenceClasses());
}

void CoarsenAlgorithm::resetSchedule(
   tbox::Pointer<xfer::CoarsenSchedule> schedule) const
{

   TBOX_ASSERT(!schedule.isNull());

   if (d_coarsen_classes->classesMatch(schedule->getEquivalenceClasses())) {
      schedule->reset(d_coarsen_classes);
   } else {
      TBOX_ERROR("CoarsenAlgorithm::resetSchedule error..."
         << "\n CoarsenClasses object passed to reset routine"
         << "\n does not match that owned by existing schedule."
         << std::endl);
   }
}

/*
 *************************************************************************
 *
 * Print coarsen algorithm data to the specified output stream.
 *
 *************************************************************************
 */

void CoarsenAlgorithm::printClassData(
   std::ostream& stream) const
{
   stream << "CoarsenAlgorithm::printClassData()" << std::endl;
   stream << "----------------------------------------" << std::endl;
   stream << "d_fill_coarse_data = " << d_fill_coarse_data << std::endl;

   d_coarsen_classes->printClassData(stream);
}

/*
 *************************************************************************
 *
 * Return the dimension
 *
 *************************************************************************
 */

const tbox::Dimension& CoarsenAlgorithm::getDim() const
{
   return d_dim;
}

}
}
#endif
