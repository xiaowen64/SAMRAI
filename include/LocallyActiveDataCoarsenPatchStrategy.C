//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/transfer/datamovers/locally_active/LocallyActiveDataCoarsenPatchStrategy.C $
// Package:	SAMRAI data transfer
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
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
