//
// File:	LocallyActiveDataCoarsenPatchStrategy.C.sed
// Package:	SAMRAI data transfer
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 306 $
// Modified:	$Date: 2005-04-26 13:49:29 -0700 (Tue, 26 Apr 2005) $
// Description:	Strategy interface to user routines for coarsening locally-active AMR data.
//

#ifndef included_xfer_LocallyActiveCoarsenPatchStrategy_C
#define included_xfer_LocallyActiveCoarsenPatchStrategy_C
 
#include "LocallyActiveDataCoarsenPatchStrategy.h"

namespace SAMRAI {
    namespace xfer {

template<int DIM>
LocallyActiveDataCoarsenPatchStrategy<DIM>::LocallyActiveDataCoarsenPatchStrategy()
{
}

template<int DIM>
LocallyActiveDataCoarsenPatchStrategy<DIM>::~LocallyActiveDataCoarsenPatchStrategy()
{
}

}
}
#endif
