/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/dlbg/LayerNodeSet.I $
 * Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 1917 $
 * Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
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
