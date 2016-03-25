/*
 * File:        $RCSfile$
 * Copyright:   (c) 1997-2003 The Regents of the University of California
 * Revision:    $Revision: 346 $
 * Modified:    $Date: 2005-05-09 12:43:12 -0700 (Mon, 09 May 2005) $
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

using namespace std;


template<int DIM>
LayerNode<DIM>::~LayerNode()
{
   return;
}



template<int DIM>
ostream &operator<<( ostream &co, const LayerNode<DIM> &r )
{
   co << r.d_owner_rank << '#' << r.d_local_index << ':' << r.getBox();
   return co;
}




}
}
#endif
