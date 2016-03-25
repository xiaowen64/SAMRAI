//
// File:	CoarsenPatchStrategy.C
// Package:	SAMRAI data transfer
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 651 $
// Modified:	$Date: 2005-10-05 14:54:35 -0700 (Wed, 05 Oct 2005) $
// Description:	Strategy interface to user routines for coarsening AMR data.
//
 
#ifndef included_xfer_CoarsenPatchStrategy_C
#define included_xfer_CoarsenPatchStrategy_C

#include "CoarsenPatchStrategy.h"

namespace SAMRAI {
    namespace xfer {

template<int DIM>  CoarsenPatchStrategy<DIM>::CoarsenPatchStrategy()
{
}

template<int DIM>  CoarsenPatchStrategy<DIM>::~CoarsenPatchStrategy()
{
}
 
}
}
#endif
