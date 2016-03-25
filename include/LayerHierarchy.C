/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/dlbg/LayerHierarchy.C $
 * Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 1917 $
 * Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
 * Description: Box graph representing hierarchy.
 */

#ifndef included_hier_LayerHierarchy_C
#define included_hier_LayerHierarchy_C

#include "LayerHierarchy.h"


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
