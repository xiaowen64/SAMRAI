/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/hierarchy/dlbg/LayerNode.C $
 * Copyright:   (c) 1997-2003 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 1704 $
 * Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
 * Description: Node in the distribued box graph.
 */

#ifndef included_hier_LayerNode_C
#define included_hier_LayerNode_C


#include "LayerNode.h"

#include <iostream>
#include <iomanip>

#ifdef DEBUG_NO_INLINE
#include "LayerNode.I"
#endif

namespace SAMRAI {
   namespace hier {

template<int DIM>
LayerNode<DIM>::~LayerNode()
{
   return;
}



template<int DIM>
std::ostream &operator<<( std::ostream &co, const LayerNode<DIM> &r )
{
   co << r.d_owner_rank << '#' << r.d_local_index << ':' << r.getBox();
   return co;
}




}
}
#endif
