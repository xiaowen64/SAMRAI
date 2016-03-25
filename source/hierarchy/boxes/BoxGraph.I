//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/hierarchy/boxes/BoxGraph.I $
// Package:     SAMRAI hierarchy
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
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

