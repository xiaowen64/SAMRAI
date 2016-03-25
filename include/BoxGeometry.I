//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/hierarchy/variables/BoxGeometry.I $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
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
