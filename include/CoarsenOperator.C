//
// File:	CoarsenOperator.C
// Package:	SAMRAI transfer 
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
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
