//
// File:	BoxOverlap.I
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 301 $
// Modified:	$Date: 2005-04-25 10:32:08 -0700 (Mon, 25 Apr 2005) $
// Description:	Base class that describes intersections between AMR boxes
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace hier {

template<int DIM> inline
BoxOverlap<DIM>::BoxOverlap()
{
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
