//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/variables/BoxGeometry.I $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	Box geometry description for overlap computations
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace hier {

template<int DIM> inline
BoxGeometry<DIM>::BoxGeometry()
{
}

template<int DIM> inline
tbox::Pointer< BoxOverlap<DIM> >
BoxGeometry<DIM>::calculateOverlap(
   const BoxGeometry<DIM>& src_geometry,
   const Box<DIM>& src_mask,
   const bool overwrite_interior,
   const IntVector<DIM>& src_offset) const
{
   return(this->calculateOverlap(
      *this, src_geometry, src_mask, overwrite_interior, src_offset, true));
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
