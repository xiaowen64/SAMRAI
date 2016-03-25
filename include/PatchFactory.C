//
// File:	PatchFactory.C
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Abstract factory class for creating patch classes
//

#ifndef included_hier_PatchFactory_C
#define included_hier_PatchFactory_C

#include "PatchFactory.h"

#ifdef DEBUG_NO_INLINE
#include "PatchFactory.I"
#endif
namespace SAMRAI {
    namespace hier {

template<int DIM>  PatchFactory<DIM>::~PatchFactory()
{
}

template<int DIM> tbox::Pointer< Patch<DIM> > PatchFactory<DIM>::allocate(
   const Box<DIM>& box,
   tbox::Pointer< PatchDescriptor<DIM> > descriptor) const
{
   return(tbox::Pointer< Patch<DIM> >(new Patch<DIM>(box, descriptor)));
}

}
}
#endif
