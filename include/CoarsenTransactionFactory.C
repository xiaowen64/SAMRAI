//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/transfer/datamovers/standard/CoarsenTransactionFactory.C $
// Package:	SAMRAI transfer
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	Interface for factory objects that create transactions for
//              coarsen schedules.
//

#ifndef included_xfer_CoarsenTransactionFactory_C
#define included_xfer_CoarsenTransactionFactory_C

#include "CoarsenTransactionFactory.h"

namespace SAMRAI {
    namespace xfer {

/*
*************************************************************************
*                                                                       *
* Default constructor and destructor.                                   * 
*                                                                       *
*************************************************************************
*/

template<int DIM> CoarsenTransactionFactory<DIM>::CoarsenTransactionFactory()
{
}

template<int DIM> CoarsenTransactionFactory<DIM>::~CoarsenTransactionFactory()
{
}

}
}
#endif
