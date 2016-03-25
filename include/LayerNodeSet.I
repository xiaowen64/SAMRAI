/*
 * File:        $RCSfile$
 * Copyright:   (c) 1997-2005 The Regents of the University of California
 * Revision:    $Revision: 346 $
 * Modified:    $Date: 2005-05-09 12:43:12 -0700 (Mon, 09 May 2005) $
 * Description: Set of layer nodes in a distributed box graph.
 */

#include "LayerNodeSet.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

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
