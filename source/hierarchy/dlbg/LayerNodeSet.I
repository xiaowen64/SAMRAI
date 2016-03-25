/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/hierarchy/dlbg/LayerNodeSet.I $
 * Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 1704 $
 * Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
 * Description: Set of layer nodes in a distributed box graph.
 */

#include "LayerNodeSet.h"


#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
namespace hier {


template<int DIM> inline
typename LayerNodeSet<DIM>::ParallelState
LayerNodeSet<DIM>::getParallelState() const
{
   return d_parallel_state;
}


#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}
