//
// File:	BoxOverlap.C
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 301 $
// Modified:	$Date: 2005-04-25 10:32:08 -0700 (Mon, 25 Apr 2005) $
// Description:	Base class that describes intersections between AMR boxes
//

#ifndef included_hier_BoxOverlap_C
#define included_hier_BoxOverlap_C

#include "BoxOverlap.h"

#ifdef DEBUG_NO_INLINE
#include "BoxOverlap.I"
#endif
namespace SAMRAI {
    namespace hier {

template<int DIM>  BoxOverlap<DIM>::~BoxOverlap()
{
}

template<int DIM> void BoxOverlap<DIM>::print(ostream& os) const
{
   os << "print() method not implemented for this overlap type" << endl;
}

}
}
#endif
