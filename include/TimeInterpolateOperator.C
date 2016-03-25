//
// File:	TimeInterpolateOperator.C
// Package:	SAMRAI transfer package
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
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
