//
// File:	LocallyActiveDataRefineAlgorithm.C
// Package:	SAMRAI data transfer
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 697 $
// Modified:	$Date: 2005-11-03 12:27:48 -0800 (Thu, 03 Nov 2005) $
// Description:	Refine algorithm for locally-active data transfer between AMR levels
//

#ifndef included_xfer_LocallyActiveDataRefineAlgorithm_C
#define included_xfer_LocallyActiveDataRefineAlgorithm_C

#include "LocallyActiveDataRefineAlgorithm.h"

#include "PatchDataFactory.h"
#include "PatchDescriptor.h"
#include "LocallyActiveVariableDatabase.h"
#include "tbox/Utilities.h"
#include "StandardLocallyActiveDataRefineTransactionFactory.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

namespace SAMRAI {
    namespace xfer {

/*
*************************************************************************
*                                                                       *
* Default constructor creates a new RefineClasses<DIM> object.          *
*                                                                       *
*************************************************************************
*/

template<int DIM>
LocallyActiveDataRefineAlgorithm<DIM>::LocallyActiveDataRefineAlgorithm()
{
   d_refine_classes = new xfer::RefineClasses<DIM>();
   d_schedule_created = false;
}

/*
*************************************************************************
*									*
* The destructor implicitly deletes the list storage associated with	*
* the refine algorithm.							*
*									*
*************************************************************************
*/

template<int DIM>
LocallyActiveDataRefineAlgorithm<DIM>::~LocallyActiveDataRefineAlgorithm()
{
}
 
/*
*************************************************************************
*									*
* Register a refine operation that will not require time interpolation.	*
*									*
*************************************************************************
*/

template<int DIM>
void LocallyActiveDataRefineAlgorithm<DIM>::registerRefine(
   const int dst,
   const int src,
   const int scratch,
   tbox::Pointer< xfer::RefineOperator<DIM> > oprefine)
{
   if (d_schedule_created) {
      TBOX_ERROR("LocallyActiveDataRefineAlgorithm<DIM>::registerRefine error..."
                 << "\nCannot call registerRefine with this refine algorithm"
                 << "\nobject since it has already been used to create a refine schedule."
                 << endl);
   }

   typename xfer::RefineClasses<DIM>::Data data;

   data.d_dst               = dst;
   data.d_src               = src;
   data.d_src_told          = -1;
   data.d_src_tnew          = -1;
   data.d_scratch           = scratch;
   data.d_fine_bdry_reps_var = hier::LocallyActiveVariableDatabase<DIM>::getDatabase()->
                                  getPatchDescriptor()->getPatchDataFactory(dst)->
                                     fineBoundaryRepresentsVariable();
   data.d_time_interpolate  = false;
   data.d_oprefine          = oprefine;
   data.d_optime            = NULL;
   data.d_tag               = -1;

   d_refine_classes->insertEquivalenceClassItem(data); 
}

/*
*************************************************************************
*                                                                       *
* Register a refine operation that will require time interpolation.     *
*                                                                       *
*************************************************************************
*/

template<int DIM> void LocallyActiveDataRefineAlgorithm<DIM>::registerRefine(
   const int dst,
   const int src,
   const int src_told,
   const int src_tnew,
   const int scratch,
   tbox::Pointer< xfer::RefineOperator<DIM> > oprefine,
   tbox::Pointer< xfer::TimeInterpolateOperator<DIM> > optime)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!optime.isNull());
#endif

   if (d_schedule_created) {
      TBOX_ERROR("LocallyActiveDataRefineAlgorithm<DIM>::registerRefine error..."
                 << "\nCannot call registerRefine with this refine algorithm"
                 << "\nobject since it has already been used to create a refine schedule."
                 << endl);
   }

   typename xfer::RefineClasses<DIM>::Data data;

   data.d_dst                = dst;
   data.d_src                = src;
   data.d_src_told           = src_told;
   data.d_src_tnew           = src_tnew;
   data.d_scratch            = scratch;
   data.d_fine_bdry_reps_var = hier::VariableDatabase<DIM>::getDatabase()->
                                  getPatchDescriptor()->getPatchDataFactory(dst)->
                                     fineBoundaryRepresentsVariable();
   data.d_time_interpolate   = true;
   data.d_oprefine           = oprefine;
   data.d_optime             = optime;
   data.d_tag                = optime;

   d_refine_classes->insertEquivalenceClassItem(data);
}

/*
*************************************************************************
*                                                                       *
* Create a communication schedule that will move data from the          *
* interiors of the give level into the ghost cells and                  *
* interiors of the same level.                                          *
*                                                                       *
*************************************************************************
*/

template<int DIM>
tbox::Pointer< xfer::LocallyActiveDataRefineSchedule<DIM> >
LocallyActiveDataRefineAlgorithm<DIM>::createSchedule(
   tbox::Pointer< hier::PatchLevel<DIM> > level,
   tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > level_mgr,
   xfer::LocallyActiveDataRefinePatchStrategy<DIM>* patch_strategy,
   tbox::Pointer< xfer::LocallyActiveDataRefineTransactionFactory<DIM> > transaction_factory)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!level.isNull());
   assert(!level_mgr.isNull());
#endif

   d_schedule_created = true;

   tbox::Pointer< xfer::LocallyActiveDataRefineTransactionFactory<DIM> > trans_factory =
      transaction_factory;

   if (trans_factory.isNull()) {
      trans_factory = new xfer::StandardLocallyActiveDataRefineTransactionFactory<DIM>;
   }

   return(new xfer::LocallyActiveDataRefineSchedule<DIM>(level,
                                                         level_mgr,
                                                         level,
                                                         level_mgr,
                                                         d_refine_classes,
                                                         trans_factory,
                                                         patch_strategy));
}

/*
*************************************************************************
*                                                                       *
* Create a communication schedule that will move data from the          *
* interiors of the source level into the ghost cell and interiors       *
* of the destination level.                                             *
*                                                                       *
*************************************************************************
*/

template<int DIM>
tbox::Pointer< xfer::LocallyActiveDataRefineSchedule<DIM> >
LocallyActiveDataRefineAlgorithm<DIM>::createSchedule(
   tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
   tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > dst_level_mgr,
   tbox::Pointer< hier::PatchLevel<DIM> > src_level,
   tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > src_level_mgr,
   xfer::LocallyActiveDataRefinePatchStrategy<DIM>* patch_strategy,
   bool use_time_interpolation,
   tbox::Pointer< xfer::LocallyActiveDataRefineTransactionFactory<DIM> > transaction_factory)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!dst_level.isNull());
   assert(!dst_level_mgr.isNull());
   assert(!src_level.isNull());
   assert(!src_level_mgr.isNull());
#endif

   d_schedule_created = true;

   tbox::Pointer< xfer::LocallyActiveDataRefineTransactionFactory<DIM> > trans_factory =
      transaction_factory;

   if (trans_factory.isNull()) {
      trans_factory = new xfer::StandardLocallyActiveDataRefineTransactionFactory<DIM>;
   }

   return(new xfer::LocallyActiveDataRefineSchedule<DIM>(dst_level,
                                                         dst_level_mgr,
                                                         src_level,
                                                         src_level_mgr,
                                                         d_refine_classes,
                                                         trans_factory,
                                                         patch_strategy,
                                                         use_time_interpolation));

}

/*
*************************************************************************
*									*
* Create a communication schedule that moves data from the interiors	*
* of the level and coarser levels in the hierarchy into the interior    *
* and boundary cells of the given level.				*
*									*
*************************************************************************
*/

template<int DIM>
tbox::Pointer< xfer::LocallyActiveDataRefineSchedule<DIM> >
LocallyActiveDataRefineAlgorithm<DIM>::createSchedule(
   tbox::Pointer< hier::PatchLevel<DIM> > level,
   tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > level_mgr,
   const int next_coarser_level,
   tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
   xfer::LocallyActiveDataRefinePatchStrategy<DIM>* patch_strategy,
   bool use_time_interpolation,
   tbox::Pointer< xfer::LocallyActiveDataRefineTransactionFactory<DIM> > transaction_factory)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!level.isNull());
   assert(!level_mgr.isNull());
   assert((next_coarser_level == -1) || !hierarchy.isNull());
#endif

   d_schedule_created = true;

   tbox::Pointer< xfer::LocallyActiveDataRefineTransactionFactory<DIM> > trans_factory =
      transaction_factory;

   if (trans_factory.isNull()) {
      trans_factory = new xfer::StandardLocallyActiveDataRefineTransactionFactory<DIM>;
   }

   return(new xfer::LocallyActiveDataRefineSchedule<DIM>(level,
                                                         level_mgr,
                                                         level,
                                                         level_mgr,
                                                         next_coarser_level,
                                                         hierarchy,
                                                         d_refine_classes,
                                                         trans_factory,
                                                         patch_strategy,
                                                         use_time_interpolation));
}

/*
*************************************************************************
*                                                                       *
* Create a communication schedule that moves data from the interiors    *
* of the source level and coarser levels in the hierarchy into the      *
* ghost cells and interior cells of the destination level.              *
*                                                                       *
*************************************************************************
*/

template<int DIM>
tbox::Pointer< xfer::LocallyActiveDataRefineSchedule<DIM> >
LocallyActiveDataRefineAlgorithm<DIM>::createSchedule(
   tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
   tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > dst_level_mgr,
   tbox::Pointer< hier::PatchLevel<DIM> > src_level,
   tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > src_level_mgr,
   const int next_coarser_level,
   tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
   xfer::LocallyActiveDataRefinePatchStrategy<DIM>* patch_strategy,
   bool use_time_interpolation,
   tbox::Pointer< xfer::LocallyActiveDataRefineTransactionFactory<DIM> > transaction_factory)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!dst_level.isNull());
   assert(!dst_level_mgr.isNull());
   if (!src_level.isNull()) assert(!src_level_mgr.isNull());
   assert((next_coarser_level == -1) || !hierarchy.isNull());
#endif

   d_schedule_created = true;

   tbox::Pointer< xfer::LocallyActiveDataRefineTransactionFactory<DIM> > trans_factory =
      transaction_factory;

   if (trans_factory.isNull()) {
      trans_factory = new xfer::StandardLocallyActiveDataRefineTransactionFactory<DIM>;
   }

   return(new xfer::LocallyActiveDataRefineSchedule<DIM>(dst_level,
                                                         dst_level_mgr,
                                                         src_level,
                                                         src_level_mgr, 
                                                         next_coarser_level,
                                                         hierarchy,
                                                         d_refine_classes,
                                                         trans_factory,
                                                         patch_strategy,
                                                         use_time_interpolation));
}

/*
*************************************************************************
*                                                                       *
* Print refine algorithm data to the specified output stream.           *
*                                                                       *
*************************************************************************
*/

template<int DIM>
void LocallyActiveDataRefineAlgorithm<DIM>::printClassData(ostream& stream) const
{
   stream << "LocallyActiveDataRefineAlgorithm<DIM>::printClassData()" << endl;
   stream << "----------------------------------------" << endl;
   d_refine_classes->printClassData(stream);
}

}
}
#endif
