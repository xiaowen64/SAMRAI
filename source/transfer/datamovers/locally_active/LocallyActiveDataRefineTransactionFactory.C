//
// File:	LocallyActiveDataRefineTransactionFactory.C
// Package:	SAMRAI transfer
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 579 $
// Modified:	$Date: 2005-08-22 15:01:30 -0700 (Mon, 22 Aug 2005) $
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
