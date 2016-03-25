//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/transfer/operators/TimeInterpolateOperator.C $
// Package:	SAMRAI transfer package
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
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
