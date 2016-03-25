//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/boxgeometry/NodeOverlap.C $
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
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
