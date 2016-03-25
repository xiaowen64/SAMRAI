//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/boxgeometry/FaceOverlap.C $
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	hier::Box intersection information for face centered objects
//

#ifndef included_pdat_FaceOverlap_C
#define included_pdat_FaceOverlap_C

#include "FaceOverlap.h"

#ifdef DEBUG_NO_INLINE
#include "FaceOverlap.I"
#endif
namespace SAMRAI {
    namespace pdat {

template<int DIM>  FaceOverlap<DIM>::FaceOverlap(
   const hier::BoxList<DIM> boxes[DIM], const hier::IntVector<DIM>& src_offset)
{
   d_is_overlap_empty = true;
   for (int d = 0; d < DIM; d++) {
      d_dst_boxes[d] = boxes[d];
      if (!d_dst_boxes[d].isEmpty()) d_is_overlap_empty = false;
   }
   d_offset = src_offset;
}

template<int DIM>  FaceOverlap<DIM>::~FaceOverlap()
{
}

template<int DIM> bool FaceOverlap<DIM>::isOverlapEmpty() const
{
   return(d_is_overlap_empty);
}

}
}
#endif
