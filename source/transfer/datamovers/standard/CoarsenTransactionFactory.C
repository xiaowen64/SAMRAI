//
// File:	CoarsenTransactionFactory.C
// Package:	SAMRAI transfer
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 651 $
// Modified:	$Date: 2005-10-05 14:54:35 -0700 (Wed, 05 Oct 2005) $
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
