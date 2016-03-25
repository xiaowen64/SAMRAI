//
// File:	NodeGeometry.I
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	hier::Box geometry information for node centered objects
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace pdat {

template<int DIM> inline
const hier::Box<DIM>& NodeGeometry<DIM>::getBox() const
{
   return(d_box);
}

template<int DIM> inline
const hier::IntVector<DIM>& NodeGeometry<DIM>::getGhosts() const
{
   return(d_ghosts);
}

template<int DIM> inline
hier::Box<DIM> NodeGeometry<DIM>::toNodeBox(const hier::Box<DIM>& box)
{
   return(box.empty() ? box : hier::Box<DIM>(box.lower(), box.upper()+1));
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
