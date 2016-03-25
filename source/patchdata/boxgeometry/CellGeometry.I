//
// File:	CellGeometry.I
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	hier::Box geometry information for cell centered objects
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace pdat {

template<int DIM> inline
const hier::Box<DIM>& CellGeometry<DIM>::getBox() const
{
   return(d_box);
}

template<int DIM> inline
const hier::IntVector<DIM>& CellGeometry<DIM>::getGhosts() const
{
   return(d_ghosts);
}

template<int DIM> inline
hier::Box<DIM> CellGeometry<DIM>::toCellBox(const hier::Box<DIM>& box)
{
   return(box);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
