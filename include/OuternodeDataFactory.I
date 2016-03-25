//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/outernode/OuternodeDataFactory.I $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1973 $
// Modified:	$LastChangedDate: 2008-02-11 16:39:15 -0800 (Mon, 11 Feb 2008) $
// Description: Factory class for creating outernode data objects
//



#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
    namespace pdat {

template <int DIM, class TYPE> inline
int OuternodeDataFactory<DIM,TYPE>::getDefaultDepth() const
{
   return(d_depth);
}

template <int DIM, class TYPE> inline
void OuternodeDataFactory<DIM,TYPE>::setDefaultDepth(const int depth)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(depth > 0);
#endif
   d_depth = depth;
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}

