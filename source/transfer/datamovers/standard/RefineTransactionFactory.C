//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-1/source/transfer/datamovers/standard/RefineTransactionFactory.C $
// Package:	SAMRAI transfer
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1846 $
// Modified:	$LastChangedDate: 2008-01-11 09:51:05 -0800 (Fri, 11 Jan 2008) $
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
   NULL_USE(fill_time);
}

template<int DIM> void RefineTransactionFactory<DIM>::preprocessScratchSpace(
   tbox::Pointer< hier::PatchLevel<DIM> > level,
   double fill_time,
   const hier::ComponentSelector& preprocess_vector) const
{
   NULL_USE(level);
   NULL_USE(fill_time);
   NULL_USE(preprocess_vector);
}

}
}
#endif
