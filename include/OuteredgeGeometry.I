//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/boxgeometry/OuteredgeGeometry.I $
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1973 $
// Modified:	$LastChangedDate: 2008-02-11 16:39:15 -0800 (Mon, 11 Feb 2008) $
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

