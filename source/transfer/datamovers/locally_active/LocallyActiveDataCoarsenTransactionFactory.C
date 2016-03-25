//
// File:	LocallyActiveDataCoarsenTransactionFactory.C
// Package:	SAMRAI transfer
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 575 $
// Modified:	$Date: 2005-08-19 17:54:08 -0700 (Fri, 19 Aug 2005) $
// Description:	Interface for factory objects that create transactions for
//              locally-active data coarsen schedules.
//

#ifndef included_xfer_LocallyActiveDataCoarsenTransactionFactory_C
#define included_xfer_LocallyActiveDataCoarsenTransactionFactory_C

#include "LocallyActiveDataCoarsenTransactionFactory.h"

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
LocallyActiveDataCoarsenTransactionFactory<DIM>::LocallyActiveDataCoarsenTransactionFactory()
{
}

template<int DIM> 
LocallyActiveDataCoarsenTransactionFactory<DIM>::~LocallyActiveDataCoarsenTransactionFactory()
{
}

}
}
#endif
