//
// File:	SideGeometry.I
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	hier::Box geometry information for side centered objects
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace pdat {

template<int DIM> inline
const hier::Box<DIM>& SideGeometry<DIM>::getBox() const
{
   return(d_box);
}

template<int DIM> inline
const hier::IntVector<DIM>& SideGeometry<DIM>::getGhosts() const
{
   return(d_ghosts);
}

template<int DIM> inline
const hier::IntVector<DIM>& SideGeometry<DIM>::getDirectionVector() const
{
   return(d_directions);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
