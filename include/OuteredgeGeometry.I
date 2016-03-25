//
// File:	OuteredgeGeometry.I
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Release:	$Name$
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Box geometry information for outer edge centered objects
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
    namespace pdat {

template<int DIM> inline
const hier::Box<DIM>& OuteredgeGeometry<DIM>::getBox() const
{
   return(d_box);
}

template<int DIM> inline
const hier::IntVector<DIM>& OuteredgeGeometry<DIM>::getGhosts() const
{
   return(d_ghosts);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}

