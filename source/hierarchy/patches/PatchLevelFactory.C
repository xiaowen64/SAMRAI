//
// File:	PatchLevelFactory.C
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Abstract factory class for creating patch level objects
//

#ifndef included_hier_PatchLevelFactory_C
#define included_hier_PatchLevelFactory_C

#include "PatchLevelFactory.h"

#ifdef DEBUG_NO_INLINE
#include "PatchLevelFactory.I"
#endif
namespace SAMRAI {
    namespace hier {

template<int DIM>  PatchLevelFactory<DIM>::~PatchLevelFactory()
{
}

template<int DIM> tbox::Pointer< PatchLevel<DIM> > PatchLevelFactory<DIM>::allocate(
   const BoxArray<DIM>& boxes,
   const ProcessorMapping& mapping,
   const IntVector<DIM>& ratio_to_level_zero,
   const tbox::Pointer< GridGeometry<DIM> > grid_geometry,
   const tbox::Pointer< PatchDescriptor<DIM> > descriptor,
   tbox::Pointer< PatchFactory<DIM> > factory) const
{
   PatchLevel<DIM> *pl =
      new PatchLevel<DIM>(boxes,
                           mapping,
                           ratio_to_level_zero,
                           grid_geometry,
                           descriptor,
                           factory);
   return(tbox::Pointer< PatchLevel<DIM> >(pl));
}

template<int DIM> tbox::Pointer< PatchLevel<DIM> > PatchLevelFactory<DIM>::allocate(
   tbox::Pointer<tbox::Database> database,
   const tbox::Pointer< GridGeometry<DIM> > grid_geometry,
   const tbox::Pointer< PatchDescriptor<DIM> > descriptor,
   tbox::Pointer< PatchFactory<DIM> > factory,
   const ComponentSelector component_selector) const
{
   PatchLevel<DIM> *pl =
      new PatchLevel<DIM>(database,
                           grid_geometry,
                           descriptor,
                           factory,
                           component_selector);
   return(tbox::Pointer< PatchLevel<DIM> >(pl));
}

}
}
#endif
