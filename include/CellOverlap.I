//
// File:	CellOverlap.I
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	hier::Box intersection information for cell centered objects
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace pdat {

template<int DIM> inline
const hier::BoxList<DIM>&
CellOverlap<DIM>::getDestinationBoxList() const
{
   return(d_dst_boxes);
}

template<int DIM> inline
const hier::IntVector<DIM>& CellOverlap<DIM>::getSourceOffset() const
{
   return(d_offset);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
