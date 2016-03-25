//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/boxgeometry/OutersideGeometry.I $
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
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
