/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/hierarchy/dlbg/LayerEdgeSet.I $
 * Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 1704 $
 * Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
 * Description: Set of edges in distributed box graph.
 */

#include "LayerEdgeSet.h"


#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
namespace hier {


template<int DIM> inline
typename LayerEdgeSet<DIM>::ParallelState
LayerEdgeSet<DIM>::getParallelState() const
{
   return d_parallel_state;
}


template<int DIM> inline
const hier::IntVector<DIM> &LayerEdgeSet<DIM>::getMaxGhostCellWidth() const
{
   return d_gcw;
}


template<int DIM> inline
const hier::IntVector<DIM> &LayerEdgeSet<DIM>::getBaseRefinementRatio() const
{
   return d_base->getRefinementRatio();
}

template<int DIM> inline
const hier::IntVector<DIM> &LayerEdgeSet<DIM>::getHeadRefinementRatio() const
{
   return d_head_ratio;
}


template<int DIM> inline
const typename LayerEdgeSet<DIM>::NodeContainer
&LayerEdgeSet<DIM>::getNodeContainer(const int rank) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if ( d_base == NULL ) {
      TBOX_ERROR("Cannot get node container.\n"
                 <<"The node layer has not been set.");
   }
#endif
   return d_base->getNodeContainer(rank);
}


#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}
