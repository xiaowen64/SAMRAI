//
// File:	EdgeOverlap.I
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	hier::Box intersection information for edge centered objects
//

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace pdat {

template<int DIM> inline
const hier::BoxList<DIM>&
EdgeOverlap<DIM>::getDestinationBoxList(const int axis) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((axis >= 0) && (axis < DIM));
#endif
   return(d_dst_boxes[axis]);
}

template<int DIM> inline
const hier::IntVector<DIM>& EdgeOverlap<DIM>::getSourceOffset() const
{
   return(d_offset);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
