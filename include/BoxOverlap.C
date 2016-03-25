//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/hierarchy/variables/BoxOverlap.C $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
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

template<int DIM> void BoxOverlap<DIM>::print(std::ostream& os) const
{
   os << "print() method not implemented for this overlap type" << std::endl;
}

}
}
#endif
