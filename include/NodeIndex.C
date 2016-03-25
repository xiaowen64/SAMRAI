//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/node/NodeIndex.C $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	hier::Index for node centered patch data types
//

#ifndef included_pdat_NodeIndex_C
#define included_pdat_NodeIndex_C

#include "NodeIndex.h"

#ifdef DEBUG_NO_INLINE
#include "NodeIndex.I"
#endif
namespace SAMRAI {
    namespace pdat {

template<int DIM> hier::IntVector<DIM> NodeIndex<DIM>::s_offsets[2 << DIM];
template<int DIM> bool NodeIndex<DIM>::s_offsets_are_set = false;

}
}
#endif
