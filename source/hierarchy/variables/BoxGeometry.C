//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/variables/BoxGeometry.C $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	Box geometry description for overlap computations
//

#ifndef included_hier_BoxGeometry_C
#define included_hier_BoxGeometry_C

#include "BoxGeometry.h"

#ifdef DEBUG_NO_INLINE
#include "BoxGeometry.I"
#endif

namespace SAMRAI {
   namespace hier {


template<int DIM>  BoxGeometry<DIM>::~BoxGeometry()
{
}

}
}
#endif
