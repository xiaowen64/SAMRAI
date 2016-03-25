//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/variables/PatchDataFactory.I $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2139 $
// Modified:	$LastChangedDate: 2008-04-17 17:35:40 -0700 (Thu, 17 Apr 2008) $
// Description:	Factory abstract base class for creating patch data objects
//

#include "tbox/Utilities.h"

#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace hier {

template<int DIM> inline
PatchDataFactory<DIM>::PatchDataFactory(const IntVector<DIM>& ghosts)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(ghosts.min() >= 0);
#endif
   d_ghosts = ghosts;
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
