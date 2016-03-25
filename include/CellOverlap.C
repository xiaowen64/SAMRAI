//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/boxgeometry/CellOverlap.C $
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	hier::Box intersection information for cell centered objects
//

#ifndef included_pdat_CellOverlap_C
#define included_pdat_CellOverlap_C

#include "CellOverlap.h"

#ifdef DEBUG_NO_INLINE
#include "CellOverlap.I"
#endif
namespace SAMRAI {
    namespace pdat {

template<int DIM>  CellOverlap<DIM>::CellOverlap(
   const hier::BoxList<DIM>& boxes, const hier::IntVector<DIM>& src_offset)
{
   d_dst_boxes        = boxes;
   d_is_overlap_empty = d_dst_boxes.isEmpty();
   d_offset           = src_offset;
}

template<int DIM>  CellOverlap<DIM>::~CellOverlap()
{
}

template<int DIM> bool CellOverlap<DIM>::isOverlapEmpty() const
{
   return(d_is_overlap_empty);
}

template<int DIM> void CellOverlap<DIM>::print(std::ostream& os) const
{
   for (typename hier::BoxList<DIM>::Iterator b(d_dst_boxes); b; b++) { 
      const hier::Box<DIM>& box = b();
      os << "      box: " << box << std::endl;
   }
}

}
}
#endif
