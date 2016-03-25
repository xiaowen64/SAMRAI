//
// File:	SideOverlap.C
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
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
