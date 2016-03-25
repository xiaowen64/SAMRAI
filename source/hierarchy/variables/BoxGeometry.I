//
// File:	BoxGeometry.I
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 301 $
// Modified:	$Date: 2005-04-25 10:32:08 -0700 (Mon, 25 Apr 2005) $
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
