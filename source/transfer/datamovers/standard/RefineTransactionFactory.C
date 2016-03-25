//
// File:	RefineTransactionFactory.C
// Package:	SAMRAI transfer
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 651 $
// Modified:	$Date: 2005-10-05 14:54:35 -0700 (Wed, 05 Oct 2005) $
// Description:	Interface for factory objects that create transactions for
//              refine schedules.
//

#ifndef included_xfer_RefineTransactionFactory_C
#define included_xfer_RefineTransactionFactory_C

#include "RefineTransactionFactory.h"

namespace SAMRAI {
    namespace xfer {

/*
*************************************************************************
*                                                                       *
* Default constructor and destructor.                                   * 
*                                                                       *
*************************************************************************
*/

template<int DIM> RefineTransactionFactory<DIM>::RefineTransactionFactory()
{
}

template<int DIM> RefineTransactionFactory<DIM>::~RefineTransactionFactory()
{
}

/*
*************************************************************************
*                                                                       *
* Default no-op implementations of optional virtual functions.          *
*                                                                       *
*************************************************************************
*/

template<int DIM> void RefineTransactionFactory<DIM>::setTransactionTime(
   double fill_time)
{
   (void) fill_time;
}

template<int DIM> void RefineTransactionFactory<DIM>::preprocessScratchSpace(
   tbox::Pointer< hier::PatchLevel<DIM> > level,
   double fill_time,
   const hier::ComponentSelector& preprocess_vector) const
{
   (void) level;
   (void) fill_time;
   (void) preprocess_vector;
}

}
}
#endif
