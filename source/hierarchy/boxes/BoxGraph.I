//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/boxes/BoxGraph.I $
// Package:     SAMRAI hierarchy
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1917 $
// Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
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

