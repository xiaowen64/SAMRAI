//
// File:	NodeOverlap.C
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
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
