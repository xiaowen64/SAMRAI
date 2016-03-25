//
// File:        BoxGraph.I
// Package:     SAMRAI hierarchy
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Utility class to determines topographical box relationships
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
   namespace hier {


template<int DIM> inline
const tbox::Array<int>& BoxGraph<DIM>::getSrcOverlapIndices(int dst_index)
{
   return d_adj[dst_index];
}

template<int DIM> inline
const tbox::Array<int>& BoxGraph<DIM>::getLocalSrcOverlapIndices(int dst_index)
{
   return d_adj_local[dst_index];
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif


}
}

