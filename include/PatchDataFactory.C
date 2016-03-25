//
// File:	PatchDataFactory.C
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 301 $
// Modified:	$Date: 2005-04-25 10:32:08 -0700 (Mon, 25 Apr 2005) $
// Description:	Factory abstract base class for creating patch data objects
//

#ifndef included_hier_PatchDataFactory_C
#define included_hier_PatchDataFactory_C

#include "PatchDataFactory.h"

#ifdef DEBUG_NO_INLINE
#include "PatchDataFactory.I"
#endif
namespace SAMRAI {
    namespace hier {

template<int DIM>  PatchDataFactory<DIM>::~PatchDataFactory()
{
}

}
}
#endif
