//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/patches/PatchHierarchy.I $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2043 $
// Modified:	$LastChangedDate: 2008-03-12 09:14:32 -0700 (Wed, 12 Mar 2008) $
// Description:	An AMR hierarchy of patch levels
//


#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace hier {

template<int DIM> inline
tbox::Pointer< BasePatchLevel<DIM> >
PatchHierarchy<DIM>::getPatchLevel(const int l) const 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((l >= 0) && (l < d_number_levels));
#endif
   return(d_patch_levels[l]);
}

template<int DIM> inline
tbox::Pointer< PatchDescriptor<DIM> >
PatchHierarchy<DIM>::getPatchDescriptor() const
{
   return(d_patch_descriptor);
}

template<int DIM> inline
bool PatchHierarchy<DIM>::levelExists(const int l) const
{
   return((l < d_number_levels) && !(d_patch_levels[l].isNull()));
}

template<int DIM> inline
bool PatchHierarchy<DIM>::finerLevelExists(const int l) const
{
   return((l+1 < d_number_levels) && !(d_patch_levels[l+1].isNull()));
}

template<int DIM> inline
int PatchHierarchy<DIM>::getNumberOfLevels() const
{
   return(d_number_levels);
}

template<int DIM> inline
int PatchHierarchy<DIM>::getFinestLevelNumber() const
{
   return(d_number_levels-1);
}

template<int DIM> inline
void PatchHierarchy<DIM>::setPatchFactory(
   tbox::Pointer< PatchFactory<DIM> > factory)
{
   d_patch_factory = factory;
}

template<int DIM> inline
void PatchHierarchy<DIM>::setPatchLevelFactory(
   tbox::Pointer< PatchLevelFactory<DIM> > factory)
{
   d_patch_level_factory = factory;
}

template<int DIM> inline
tbox::Pointer< GridGeometry<DIM> >
PatchHierarchy<DIM>::getGridGeometry() const 
{
   return(d_grid_geometry);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}
