//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/transfer/datamovers/locally_active/LocallyActiveDataRefineTransactionFactory.C $
// Package:	SAMRAI transfer
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Interface for factory objects that create transactions for
//              locally-active data refine schedules.
//

#ifndef included_xfer_LocallyActiveDataRefineTransactionFactory_C
#define included_xfer_LocallyActiveDataRefineTransactionFactory_C

#include "LocallyActiveDataRefineTransactionFactory.h"

namespace SAMRAI {
    namespace xfer {

/*
*************************************************************************
*                                                                       *
* Default constructor and destructor.                                   * 
*                                                                       *
*************************************************************************
*/

template<int DIM> 
LocallyActiveDataRefineTransactionFactory<DIM>::LocallyActiveDataRefineTransactionFactory()
{
}

template<int DIM> 
LocallyActiveDataRefineTransactionFactory<DIM>::~LocallyActiveDataRefineTransactionFactory()
{
}

/*
*************************************************************************
*                                                                       *
* Default no-op implementations of optional virtual functions.          *
*                                                                       *
*************************************************************************
*/

template<int DIM> void 
LocallyActiveDataRefineTransactionFactory<DIM>::setTransactionTime(
   double fill_time)
{
   (void) fill_time;
}

template<int DIM> void 
LocallyActiveDataRefineTransactionFactory<DIM>::preprocessScratchSpace(
   tbox::Pointer< hier::PatchLevel<DIM> > level,
   double fill_time,
   const hier::LocallyActiveDataPatchLevelManager<DIM>& allocate_mgr) const
{
   (void) level;
   (void) fill_time;
   (void) allocate_mgr;
}

}
}
#endif
