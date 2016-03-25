//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/patches/PatchFactory.C $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
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
