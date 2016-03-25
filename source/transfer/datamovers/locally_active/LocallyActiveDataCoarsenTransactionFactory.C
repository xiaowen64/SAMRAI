//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/transfer/datamovers/locally_active/LocallyActiveDataCoarsenTransactionFactory.C $
// Package:	SAMRAI transfer
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
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
