//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/transfer/datamovers/standard/CoarsenPatchStrategy.C $
// Package:	SAMRAI data transfer
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
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
