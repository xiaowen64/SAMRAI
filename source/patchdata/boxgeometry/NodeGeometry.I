//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/boxgeometry/NodeGeometry.I $
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
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
