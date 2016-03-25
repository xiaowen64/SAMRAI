//
// File:	CellOverlap.C
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
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

template<int DIM> void CellOverlap<DIM>::print(ostream& os) const
{
   for (typename hier::BoxList<DIM>::Iterator b(d_dst_boxes); b; b++) { 
      const hier::Box<DIM>& box = b();
      os << "      box: " << box << endl;
   }
}

}
}
#endif
