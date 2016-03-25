//
// File:	OutersideGeometry.I
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	hier::Box geometry information for outerside centered objects
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace pdat {

template<int DIM> inline
const hier::Box<DIM>& OutersideGeometry<DIM>::getBox() const
{
   return(d_box);
}

template<int DIM> inline
const hier::IntVector<DIM>& OutersideGeometry<DIM>::getGhosts() const
{
   return(d_ghosts);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
