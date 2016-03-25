//
// File:	RefineAlgorithm.C
// Package:	SAMRAI data transfer
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 697 $
// Modified:	$Date: 2005-11-03 12:27:48 -0800 (Thu, 03 Nov 2005) $
// Description:	Refine algorithm for data transfer between AMR levels
//

#ifndef included_xfer_RefineAlgorithm_C
#define included_xfer_RefineAlgorithm_C

#include "RefineAlgorithm.h"

#include "PatchDataFactory.h"
#include "PatchDescriptor.h"
#include "VariableDatabase.h"
#include "tbox/Utilities.h"
#include "StandardRefineTransactionFactory.h"

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

template<int DIM>  RefineAlgorithm<DIM>::RefineAlgorithm()
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

template<int DIM>  RefineAlgorithm<DIM>::~RefineAlgorithm()
{
}
 
/*
*************************************************************************
*									*
* Register a refine operation that will not require time interpolation.	*
*									*
*************************************************************************
*/

template<int DIM> void RefineAlgorithm<DIM>::registerRefine(
   const int dst,
   const int src,
   const int scratch,
   tbox::Pointer< xfer::RefineOperator<DIM> > oprefine)
{
   if (d_schedule_created) {
      TBOX_ERROR("RefineAlgorithm<DIM>::registerRefine error..."
                 << "\nCannot call registerRefine with this refine algorithm"
                 << "\nobject since it has already been used to create a refine schedule."
                 << endl);
   }

   typename xfer::RefineClasses<DIM>::Data data;

   data.d_dst                = dst;
   data.d_src                = src;
   data.d_src_told           = -1;
   data.d_src_tnew           = -1;
   data.d_scratch            = scratch;
   data.d_fine_bdry_reps_var = hier::VariableDatabase<DIM>::getDatabase()->
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
*									*
* Register a refine operation that will require time interpolation.	*
*									*
*************************************************************************
*/

template<int DIM> void RefineAlgorithm<DIM>::registerRefine(
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
      TBOX_ERROR("RefineAlgorithm<DIM>::registerRefine error..."
                 << "\nCannot call registerRefine with this RefineAlgorithm"
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
*									*
* Create a communication schedule that will move data from the          *
* interiors of the given level into the ghost cells and                 *
* interiors of the same level.                                          *
*									*
*************************************************************************
*/

template<int DIM> tbox::Pointer< xfer::RefineSchedule<DIM> >
RefineAlgorithm<DIM>::createSchedule(
   tbox::Pointer< hier::PatchLevel<DIM> > level,
   xfer::RefinePatchStrategy<DIM>* patch_strategy,
   tbox::Pointer< xfer::RefineTransactionFactory<DIM> > transaction_factory)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!level.isNull());
#endif

   d_schedule_created = true;

   tbox::Pointer< xfer::RefineTransactionFactory<DIM> > trans_factory =
      transaction_factory;

   if (trans_factory.isNull()) {
      trans_factory = new xfer::StandardRefineTransactionFactory<DIM>;
   }

   return(new xfer::RefineSchedule<DIM>(level,
                                        level, 
                                        d_refine_classes, 
                                        trans_factory,
                                        patch_strategy));
}

/*
*************************************************************************
*									*
* Create a communication schedule that will move data from the          *
* interiors of the source level into the ghost cell and interiors       *
* of the destination level.						*
*									*
*************************************************************************
*/

template<int DIM> tbox::Pointer< xfer::RefineSchedule<DIM> >
RefineAlgorithm<DIM>::createSchedule(
   tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
   tbox::Pointer< hier::PatchLevel<DIM> > src_level,
   xfer::RefinePatchStrategy<DIM>* patch_strategy,
   bool use_time_interpolation,
   tbox::Pointer< xfer::RefineTransactionFactory<DIM> > transaction_factory)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!dst_level.isNull());
   assert(!src_level.isNull());
#endif

   d_schedule_created = true;

   tbox::Pointer< xfer::RefineTransactionFactory<DIM> > trans_factory =
      transaction_factory;

   if (trans_factory.isNull()) {
      trans_factory = new xfer::StandardRefineTransactionFactory<DIM>;
   }

   return(new xfer::RefineSchedule<DIM>(dst_level,
                                        src_level,
                                        d_refine_classes,
                                        trans_factory,
                                        patch_strategy,
                                        use_time_interpolation));
}

/*
*************************************************************************
*                                                                       *
* Create a communication schedule that moves data from the interiors    *
* of the level and coarser levels in the hierarchy into the interior    *
* and boundary cells of the given level.                                *
*                                                                       *
*************************************************************************
*/

template<int DIM> tbox::Pointer< xfer::RefineSchedule<DIM> >
RefineAlgorithm<DIM>::createSchedule(
   tbox::Pointer< hier::PatchLevel<DIM> > level,
   const int next_coarser_level,
   tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
   xfer::RefinePatchStrategy<DIM>* patch_strategy,
   bool use_time_interpolation,
   tbox::Pointer< xfer::RefineTransactionFactory<DIM> > transaction_factory)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!level.isNull());
   assert((next_coarser_level < 0) || !hierarchy.isNull());
#endif

   d_schedule_created = true;

   tbox::Pointer< xfer::RefineTransactionFactory<DIM> > trans_factory =
      transaction_factory;

   if (trans_factory.isNull()) {
      trans_factory = new xfer::StandardRefineTransactionFactory<DIM>;
   }

   return(new xfer::RefineSchedule<DIM>(level,
                                        level,
                                        next_coarser_level,
                                        hierarchy,
                                        d_refine_classes,
                                        trans_factory,
                                        patch_strategy,
                                        use_time_interpolation));
}

/*
*************************************************************************
*									*
* Create a communication schedule that moves data from the interiors	*
* of the source level and coarser levels in the hierarchy into the      *
* ghost cells and interior cells of the destination level.              *
*									*
*************************************************************************
*/

template<int DIM> tbox::Pointer< xfer::RefineSchedule<DIM> >
RefineAlgorithm<DIM>::createSchedule(
   tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
   tbox::Pointer< hier::PatchLevel<DIM> > src_level,
   const int next_coarser_level,
   tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
   xfer::RefinePatchStrategy<DIM>* patch_strategy,
   bool use_time_interpolation,
   tbox::Pointer< xfer::RefineTransactionFactory<DIM> > transaction_factory)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!dst_level.isNull());
   assert((next_coarser_level < 0) || !hierarchy.isNull());
#endif

   d_schedule_created = true;

   tbox::Pointer< xfer::RefineTransactionFactory<DIM> > trans_factory =
      transaction_factory;

   if (trans_factory.isNull()) {
      trans_factory = new xfer::StandardRefineTransactionFactory<DIM>;
   }

   return(new xfer::RefineSchedule<DIM>(dst_level,
                                        src_level,
                                        next_coarser_level,
                                        hierarchy,
                                        d_refine_classes,
                                        trans_factory,
                                        patch_strategy,
                                        use_time_interpolation)); 
}

/*
**************************************************************************
*                                                                        *
* Reconfigure refine schedule to perform operations in this algorithm.   *
*                                                                        *
**************************************************************************
*/

template<int DIM> bool RefineAlgorithm<DIM>::checkConsistency(
   tbox::Pointer< xfer::RefineSchedule<DIM> > schedule) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!schedule.isNull());
#endif
   return( d_refine_classes->
           checkConsistency(schedule->getEquivalenceClasses()) );
} 

template<int DIM> void RefineAlgorithm<DIM>::resetSchedule(
   tbox::Pointer< xfer::RefineSchedule<DIM> > schedule) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!schedule.isNull());
#endif
   if ( d_refine_classes->
        checkConsistency(schedule->getEquivalenceClasses()) ) {
      schedule->reset(d_refine_classes);
   } else {
      TBOX_ERROR("RefineAlgorithm<DIM>::resetSchedule error..."
                 << "\n Items in xfer::RefineClasses<DIM> object passed to reset" 
                 << "\n routine are inconsistent with those in existing schedule."
                 << endl);
   }
}

template<int DIM>
const tbox::Pointer< xfer::RefineClasses<DIM> >&
RefineAlgorithm<DIM>::getEquivalenceClasses() const
{
   return(d_refine_classes);
}
                                                                               
template<int DIM>
void RefineAlgorithm<DIM>::setEquivalenceClasses(
   const tbox::Pointer< xfer::RefineClasses<DIM> > refine_classes)
{
   d_refine_classes.setNull();
   d_refine_classes = refine_classes;
}

/*
*************************************************************************
*                                                                       *
* Print refine algorithm data to the specified output stream.           *
*                                                                       *
*************************************************************************
*/

template<int DIM> void RefineAlgorithm<DIM>::printClassData(ostream& stream) const
{
   stream << "RefineAlgorithm<DIM>::printClassData()" << endl;
   stream << "----------------------------------------" << endl;
   d_refine_classes->printClassData(stream);
}

}
}
#endif
