//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/node/NodeIndex.C $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
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
