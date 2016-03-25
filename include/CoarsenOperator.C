//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/transfer/operators/CoarsenOperator.C $
// Package:	SAMRAI transfer 
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Abstract base class for spatial coarsening operators.
//

#ifndef included_xfer_CoarsenOperator_C
#define included_xfer_CoarsenOperator_C

#include "CoarsenOperator.h"

namespace SAMRAI {
    namespace xfer {

template<int DIM>  CoarsenOperator<DIM>::CoarsenOperator()
{
}

template<int DIM>  CoarsenOperator<DIM>::~CoarsenOperator()
{
}

}
}
#endif
