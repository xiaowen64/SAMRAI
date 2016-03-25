//
// File:	LocallyActiveDataCoarsenAlgorithm.C
// Package:	SAMRAI data transfer
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 601 $
// Modified:	$Date: 2005-09-06 11:23:15 -0700 (Tue, 06 Sep 2005) $
// Description:	Coarsening algorithm for locally-active data transfer between AMR levels
//

#ifndef included_xfer_LocallyActiveDataCoarsenAlgorithm_C
#define included_xfer_LocallyActiveDataCoarsenAlgorithm_C
 
#include "LocallyActiveDataCoarsenAlgorithm.h"

#include "PatchDataFactory.h"
#include "PatchDescriptor.h"
#include "LocallyActiveVariableDatabase.h"
#include "tbox/Utilities.h"
#include "StandardLocallyActiveDataCoarsenTransactionFactory.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

namespace SAMRAI {
    namespace xfer {

/*
*************************************************************************
*                                                                       *
* The constructor creates a new CoarsenClasses object                   *
* and caches a boolean indiating whether to copy data to the            *
* destination space on the coarse level before coarsening.              *
*                                                                       *
*************************************************************************
*/

template<int DIM>
LocallyActiveDataCoarsenAlgorithm<DIM>::LocallyActiveDataCoarsenAlgorithm(
   bool fill_coarse_data)
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

template<int DIM>
LocallyActiveDataCoarsenAlgorithm<DIM>::~LocallyActiveDataCoarsenAlgorithm()
{
}
 
/*
*************************************************************************
*									*
* Register a coarsening operation with the coarsening algorithm.        *
*									*
*************************************************************************
*/

template<int DIM>
void LocallyActiveDataCoarsenAlgorithm<DIM>::registerCoarsen(
   const int dst,
   const int src,
   const tbox::Pointer< xfer::CoarsenOperator<DIM> > opcoarsen,
   const hier::IntVector<DIM>& gcw_to_coarsen)
{
   if (d_schedule_created) {
      TBOX_ERROR("LocallyActiveDataCoarsenAlgorithm<DIM>::registerCoarsen error..."
                 << "\nCannot call registerCoarsen with this coarsen algorithm"
                 << "\nobject since it has already been used to create a coarsen schedule."
                 << endl);
   }

   typename xfer::CoarsenClasses<DIM>::Data data;

   data.d_dst                = dst;
   data.d_src                = src;
   data.d_fine_bdry_reps_var = hier::LocallyActiveVariableDatabase<DIM>::getDatabase()->
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

template<int DIM>
tbox::Pointer< xfer::LocallyActiveDataCoarsenSchedule<DIM> > 
LocallyActiveDataCoarsenAlgorithm<DIM>::createSchedule(
   tbox::Pointer< hier::PatchLevel<DIM> > crse_level,
   tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > crse_level_mgr,
   tbox::Pointer< hier::PatchLevel<DIM> > fine_level,
   tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > fine_level_mgr,
   xfer::LocallyActiveDataCoarsenPatchStrategy<DIM>* patch_strategy,
   tbox::Pointer< xfer::LocallyActiveDataCoarsenTransactionFactory<DIM> > 
      transaction_factory)
{
   d_schedule_created = true;

   tbox::Pointer< xfer::LocallyActiveDataCoarsenTransactionFactory<DIM> > trans_factory =
      transaction_factory;

   if (trans_factory.isNull()) {
      trans_factory = new xfer::StandardLocallyActiveDataCoarsenTransactionFactory<DIM>;
   }

   return(new xfer::LocallyActiveDataCoarsenSchedule<DIM>(crse_level, 
                                                          crse_level_mgr,
                                                          fine_level,
                                                          fine_level_mgr,
                                                          d_coarsen_classes,
                                                          trans_factory, 
                                                          patch_strategy,
                                                          d_fill_coarse_data));
}

/*
*************************************************************************
*									*
* Print coarsen algorithm data to the specified output stream.		*
*									*
*************************************************************************
*/

template<int DIM>
void LocallyActiveDataCoarsenAlgorithm<DIM>::printClassData(ostream& stream) const
{
   stream << "LocallyActiveDataCoarsenAlgorithm<DIM>::printClassData()" << endl;
   stream << "----------------------------------------" << endl;
   stream << "d_fill_coarse_data = " << d_fill_coarse_data << endl;

   d_coarsen_classes->printClassData(stream);
}

}
}
#endif
