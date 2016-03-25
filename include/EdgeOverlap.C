//
// File:	EdgeOverlap.C
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	hier::Box intersection information for edge centered objects
//

#ifndef included_pdat_EdgeOverlap_C
#define included_pdat_EdgeOverlap_C

#include "EdgeOverlap.h"

#ifdef DEBUG_NO_INLINE
#include "EdgeOverlap.I"
#endif
namespace SAMRAI {
    namespace pdat {


template<int DIM>  EdgeOverlap<DIM>::EdgeOverlap(
   const hier::BoxList<DIM> boxes[DIM], const hier::IntVector<DIM>& src_offset)
{
   d_is_overlap_empty = true;
   for (int d = 0; d < DIM; d++) {
      d_dst_boxes[d] = boxes[d];
      if (!d_dst_boxes[d].isEmpty()) d_is_overlap_empty = false;
   }
   d_offset = src_offset;
}

template<int DIM>  EdgeOverlap<DIM>::~EdgeOverlap()
{
}

template<int DIM> bool EdgeOverlap<DIM>::isOverlapEmpty() const
{
   return(d_is_overlap_empty);
}

}
}
#endif
