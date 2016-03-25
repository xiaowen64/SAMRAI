//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/transfer/operators/RefineOperator.C $
// Package:	SAMRAI transfer
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description: Abstract base class for spatial refinement operators.
//

#ifndef included_xfer_RefineOperator_C
#define included_xfer_RefineOperator_C

#include "RefineOperator.h"

namespace SAMRAI {
    namespace xfer {

template<int DIM>  RefineOperator<DIM>::RefineOperator()
{
}

template<int DIM>  RefineOperator<DIM>::~RefineOperator()
{
}

}
}
#endif
