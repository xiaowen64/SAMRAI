//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/boxgeometry/NodeGeometry.I $
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
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
