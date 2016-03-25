//
// File:	NodeIndex.C
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
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
