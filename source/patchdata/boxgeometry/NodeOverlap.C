//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/boxgeometry/NodeOverlap.C $
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	hier::Box intersection information for node centered objects
//

#ifndef included_pdat_NodeOverlap_C
#define included_pdat_NodeOverlap_C

#include "NodeOverlap.h"

#ifdef DEBUG_NO_INLINE
#include "NodeOverlap.I"
#endif
namespace SAMRAI {
    namespace pdat {

template<int DIM>  NodeOverlap<DIM>::NodeOverlap(
   const hier::BoxList<DIM>& boxes, const hier::IntVector<DIM>& src_offset)
{
   d_dst_boxes        = boxes;
   d_is_overlap_empty = d_dst_boxes.isEmpty();
   d_offset           = src_offset;
}

template<int DIM>  NodeOverlap<DIM>::~NodeOverlap()
{
}

template<int DIM> bool NodeOverlap<DIM>::isOverlapEmpty() const
{
   return(d_is_overlap_empty);
}

}
}
#endif
