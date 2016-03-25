//
// File:	CoarsenAlgorithm.C
// Package:	SAMRAI data transfer
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 651 $
// Modified:	$Date: 2005-10-05 14:54:35 -0700 (Wed, 05 Oct 2005) $
// Description:	Coarsening algorithm for data transfer between AMR levels
//
 
#ifndef included_xfer_CoarsenAlgorithm_C
#define included_xfer_CoarsenAlgorithm_C

#include "CoarsenAlgorithm.h"

#include "PatchDataFactory.h"
#include "PatchDescriptor.h"
#include "VariableDatabase.h"
#include "tbox/Utilities.h"
#include "StandardCoarsenTransactionFactory.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

namespace SAMRAI {
    namespace xfer {

/*
*************************************************************************
*                                                                       *
* The constructor creates a new CoarsenClasses<DIM> object             *
* and caches a boolean indiating whether to copy data to the            *
* destination space on the coarse level before coarsening.              *
*                                                                       *
*************************************************************************
*/

template<int DIM>  CoarsenAlgorithm<DIM>::CoarsenAlgorithm(bool fill_coarse_data)
{
   d_fill_coarse_data = fill_coarse_data;
   d_coarsen_classes = new xfer::CoarsenClasses<DIM>(d_fill_coarse_data);
   d_schedule_created = false;
}

/*
*************************************************************************
*									*
* The destructor implicitly deallocates the list data.                  *
*									*
*************************************************************************
*/

template<int DIM>  CoarsenAlgorithm<DIM>::~CoarsenAlgorithm()
{
}
 
/*
*************************************************************************
*									*
* Register a coarsening operation with the coarsening algorithm.        *
*									*
*************************************************************************
*/

template<int DIM> void CoarsenAlgorithm<DIM>::registerCoarsen(
   const int dst,
   const int src,
   const tbox::Pointer< xfer::CoarsenOperator<DIM> > opcoarsen,
   const hier::IntVector<DIM>& gcw_to_coarsen)
{
   if (d_schedule_created) {
      TBOX_ERROR("CoarsenAlgorithm<DIM>::registerCoarsen error..."
                 << "\nCannot call registerCoarsen with this coarsen algorithm"
                 << "\nobject since it has already been used to create a coarsen schedule."
                 << endl);
   }

   typename xfer::CoarsenClasses<DIM>::Data data;

   data.d_dst                = dst;
   data.d_src                = src;
   data.d_fine_bdry_reps_var = hier::VariableDatabase<DIM>::getDatabase()->
                                  getPatchDescriptor()->getPatchDataFactory(dst)->
                                     fineBoundaryRepresentsVariable();
   data.d_gcw_to_coarsen     = gcw_to_coarsen;
   data.d_opcoarsen          = opcoarsen;
   data.d_tag                = -1;

   d_coarsen_classes->insertEquivalenceClassItem(data);
}

/*
*************************************************************************
*									*
* Create a communication schedule that will coarsen data from fine      *
* patch level to the coarse patch level.                                *
*									*
*************************************************************************
*/

template<int DIM> tbox::Pointer< xfer::CoarsenSchedule<DIM> > 
CoarsenAlgorithm<DIM>::createSchedule(
   tbox::Pointer< hier::PatchLevel<DIM> > crse_level,
   tbox::Pointer< hier::PatchLevel<DIM> > fine_level,
   xfer::CoarsenPatchStrategy<DIM>* patch_strategy,
   tbox::Pointer< xfer::CoarsenTransactionFactory<DIM> > transaction_factory)
{
   d_schedule_created = true;

   tbox::Pointer< xfer::CoarsenTransactionFactory<DIM> > trans_factory = 
      transaction_factory;

   if (trans_factory.isNull()) {
      trans_factory = new xfer::StandardCoarsenTransactionFactory<DIM>;
   }
 
   return(new xfer::CoarsenSchedule<DIM>(crse_level, 
                                         fine_level, 
                                         d_coarsen_classes,
                                         trans_factory,
                                         patch_strategy,
                                         d_fill_coarse_data));
}

/*
**************************************************************************
*                                                                        *
* Reconfigure coarsen schedule to perform operations in this algorithm.  *
*                                                                        *
**************************************************************************
*/

template<int DIM> bool CoarsenAlgorithm<DIM>::checkConsistency(
   tbox::Pointer< xfer::CoarsenSchedule<DIM> > schedule) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!schedule.isNull());
#endif
   return( d_coarsen_classes->
           checkConsistency(schedule->getEquivalenceClasses()) );
}

template<int DIM> void CoarsenAlgorithm<DIM>::resetSchedule(
   tbox::Pointer< xfer::CoarsenSchedule<DIM> > schedule) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!schedule.isNull());
#endif
   if (d_coarsen_classes->checkConsistency(schedule->getEquivalenceClasses())) {
      schedule->reset(d_coarsen_classes);
   } else {
      TBOX_ERROR("CoarsenAlgorithm<DIM>::resetSchedule error..."
                 << "\n CoarsenClasses<DIM> object passed to reset routine"
                 << "\n inconsistent with that owned by existing schedule."
                 << endl);
   }
}

/*
*************************************************************************
*									*
* Print coarsen algorithm data to the specified output stream.		*
*									*
*************************************************************************
*/

template<int DIM> void CoarsenAlgorithm<DIM>::printClassData(ostream& stream) const
{
   stream << "CoarsenAlgorithm<DIM>::printClassData()" << endl;
   stream << "----------------------------------------" << endl;
   stream << "d_fill_coarse_data = " << d_fill_coarse_data << endl;

   d_coarsen_classes->printClassData(stream);
}

}
}
#endif
