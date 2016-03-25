/*
 * File:        $RCSfile$
 * Copyright:   (c) 1997-2005 The Regents of the University of California
 * Revision:    $Revision: 346 $
 * Modified:    $Date: 2005-05-09 12:43:12 -0700 (Mon, 09 May 2005) $
 * Description: Box graph representing hierarchy.
 */

#ifndef included_hier_LayerHierarchy_C
#define included_hier_LayerHierarchy_C

#include "LayerHierarchy.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

#ifdef DEBUG_NO_INLINE
#include "LayerHierarchy.I"
#endif

namespace SAMRAI {
namespace hier {


template<int DIM>
LayerHierarchy<DIM>::LayerHierarchy()
{
   return;
}


template<int DIM>
LayerHierarchy<DIM>::~LayerHierarchy()
{
   d_num_levels = 0;
   d_node_sets.setNull();
   d_peer_edges.setNull();
   d_fine_edges.setNull();
   d_coarse_edges.setNull();
   return;
}


}
}
#endif
