//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/patches/PatchLevelFactory.C $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
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
   tbox::Pointer< PatchFactory<DIM> > factory,
   const bool defer_boundary_box_creation) const
{
   PatchLevel<DIM> *pl =
      new PatchLevel<DIM>(boxes,
                          mapping,
                          ratio_to_level_zero,
                          grid_geometry,
                          descriptor,
                          factory,
                          defer_boundary_box_creation);
   return(tbox::Pointer< PatchLevel<DIM> >(pl));
}

template<int DIM> tbox::Pointer< PatchLevel<DIM> > PatchLevelFactory<DIM>::allocate(
   tbox::Pointer<tbox::Database> database,
   const tbox::Pointer< GridGeometry<DIM> > grid_geometry,
   const tbox::Pointer< PatchDescriptor<DIM> > descriptor,
   tbox::Pointer< PatchFactory<DIM> > factory,
   const ComponentSelector component_selector,
   const bool defer_boundary_box_creation) const
{
   PatchLevel<DIM> *pl =
      new PatchLevel<DIM>(database,
                          grid_geometry,
                          descriptor,
                          factory,
                          component_selector,
                          defer_boundary_box_creation);
   return(tbox::Pointer< PatchLevel<DIM> >(pl));
}

}
}
#endif
