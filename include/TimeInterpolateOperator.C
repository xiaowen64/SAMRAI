//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/transfer/operators/TimeInterpolateOperator.C $
// Package:	SAMRAI transfer package
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description: Abstract base class for time interpolation operators.
//

#ifndef included_xfer_TimeInterpolateOperator_C
#define included_xfer_TimeInterpolateOperator_C

#include "TimeInterpolateOperator.h"

namespace SAMRAI {
    namespace xfer {

template<int DIM>  TimeInterpolateOperator<DIM>::TimeInterpolateOperator()
{
}

template<int DIM>  TimeInterpolateOperator<DIM>::~TimeInterpolateOperator()
{
}

}
}
#endif
