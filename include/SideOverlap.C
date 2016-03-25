//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/boxgeometry/SideOverlap.C $
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	hier::Box intersection information for side centered objects
//

#ifndef included_pdat_SideOverlap_C
#define included_pdat_SideOverlap_C

#include "SideOverlap.h"

#ifdef DEBUG_NO_INLINE
#include "SideOverlap.I"
#endif
namespace SAMRAI {
    namespace pdat {

template<int DIM>  SideOverlap<DIM>::SideOverlap(
   const hier::BoxList<DIM> boxes[DIM], const hier::IntVector<DIM>& src_offset)
{
   d_is_overlap_empty = true;
   for (int d = 0; d < DIM; d++) {
      d_dst_boxes[d] = boxes[d];
      if (!d_dst_boxes[d].isEmpty()) d_is_overlap_empty = false;
   }
   d_offset = src_offset;
}

template<int DIM>  SideOverlap<DIM>::~SideOverlap()
{
}

template<int DIM> bool SideOverlap<DIM>::isOverlapEmpty() const
{
   return(d_is_overlap_empty);
}

}
}
#endif
