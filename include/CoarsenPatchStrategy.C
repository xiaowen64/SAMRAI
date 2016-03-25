//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/transfer/datamovers/standard/CoarsenPatchStrategy.C $
// Package:	SAMRAI data transfer
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
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
