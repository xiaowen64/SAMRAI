//
// File:	PatchLevel.I
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	A collection of patches at one level of the AMR hierarchy
//

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
    namespace hier {

template<int DIM> inline
int PatchLevel<DIM>::getLevelNumber() const
{
   return(d_level_number);
}

template<int DIM> inline
int PatchLevel<DIM>::getNextCoarserHierarchyLevelNumber() const
{
   return(d_next_coarser_level_number);
}

template<int DIM> inline
bool PatchLevel<DIM>::inHierarchy() const
{
   return(d_in_hierarchy);
}

template<int DIM> inline
tbox::Pointer< Patch<DIM> > PatchLevel<DIM>::getPatch(const int p) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((p >= 0) && (p < d_number_patches));
#endif
   return(d_patches[p]);
}

template<int DIM> inline
int PatchLevel<DIM>::getNumberOfPatches() const
{
   return(d_number_patches);
}

template<int DIM> inline
tbox::Pointer< PatchDescriptor<DIM> > PatchLevel<DIM>::getPatchDescriptor() const
{
   return(d_descriptor);
}

template<int DIM> inline
tbox::Pointer< PatchFactory<DIM> > PatchLevel<DIM>::getPatchFactory() const
{
   return(d_factory);
}

template<int DIM> inline
tbox::Pointer< GridGeometry<DIM> > PatchLevel<DIM>::getGridGeometry() const
{
   return(d_geometry);
}

template<int DIM> inline
const BoxArray<DIM>& PatchLevel<DIM>::getPhysicalDomain() const
{
   return(d_physical_domain);
}

template<int DIM> inline
const IntVector<DIM>& PatchLevel<DIM>::getRatio() const
{
   return(d_ratio_to_level_zero);
}

template<int DIM> inline
const IntVector<DIM>& PatchLevel<DIM>::getRatioToCoarserLevel() const
{
   return(d_ratio_to_coarser_level);
}

template<int DIM> inline
void PatchLevel<DIM>::setRatioToCoarserLevel(const IntVector<DIM>& ratio)
{
   d_ratio_to_coarser_level = ratio;
}
 
template<int DIM> inline
const BoxArray<DIM>& PatchLevel<DIM>::getBoxes() const
{
   return(d_boxes);
}

template<int DIM> inline
const tbox::Array< tbox::List< IntVector<DIM> > >& PatchLevel<DIM>::getShiftsForLevel() const
{
  return d_shifts;
}


template<int DIM> inline
const ProcessorMapping& PatchLevel<DIM>::getProcessorMapping() const
{
   return(d_mapping);
}

template<int DIM> inline
int PatchLevel<DIM>::getMappingForPatch(const int p) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((p >= 0) && (p < d_number_patches));
#endif
   return(d_mapping.getProcessorAssignment(p));
}

template<int DIM> inline
const Box<DIM>& PatchLevel<DIM>::getBoxForPatch(const int p) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((p >= 0) && (p < d_number_patches));
#endif
   return(d_boxes.getBox(p));
}

template<int DIM> inline
const tbox::List< IntVector<DIM> >&
PatchLevel<DIM>::getShiftsForPatch(const int p) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((p >= 0) && (p < d_number_patches));
#endif
   return(d_shifts[p]);
}

template<int DIM> inline 
bool PatchLevel<DIM>::patchTouchesRegularBoundary(const int p) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((p >= 0) && (p < d_number_patches));
#endif
   return(d_patch_touches_regular_boundary[p]);
}

template<int DIM> inline 
bool PatchLevel<DIM>::patchTouchesPeriodicBoundary(const int p) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((p >= 0) && (p < d_number_patches));
#endif
   return(d_patch_touches_periodic_boundary[p]);
}

template<int DIM> inline
PatchLevelIterator<DIM>::PatchLevelIterator()
{
   d_patch          = 1;
   d_number_patches = 0;
   d_local_box_indices = NULL;
}

template<int DIM> inline
PatchLevelIterator<DIM>::PatchLevelIterator(
   const PatchLevelIterator<DIM>& iter)
{
   d_patch          = iter.d_patch;
   d_number_patches = iter.d_number_patches;
   d_local_box_indices = iter.d_local_box_indices;
}

template<int DIM> inline
PatchLevelIterator<DIM>&
PatchLevelIterator<DIM>::operator=(const PatchLevelIterator<DIM>& iter)
{
   d_patch          = iter.d_patch;
   d_number_patches = iter.d_number_patches;
   d_local_box_indices = iter.d_local_box_indices;
   return(*this);
}

template<int DIM> inline
PatchLevelIterator<DIM>::~PatchLevelIterator()
{
}

template<int DIM> inline
int PatchLevelIterator<DIM>::operator*() const
{
   return((*d_local_box_indices)[d_patch]);
}

template<int DIM> inline
int PatchLevelIterator<DIM>::operator()() const
{
   return((*d_local_box_indices)[d_patch]);
}

template<int DIM> inline
PatchLevelIterator<DIM>::operator bool() const
{
   return(d_patch < d_number_patches);
}

#ifndef LACKS_BOOL_VOID_RESOLUTION 
template<int DIM> inline
PatchLevelIterator<DIM>::operator const void*() const
{
   return(d_patch < d_number_patches ? this : NULL);
}
#endif

template<int DIM> inline
bool PatchLevelIterator<DIM>::operator!() const
{
   return(d_patch >= d_number_patches);
}

template<int DIM> inline
void PatchLevelIterator<DIM>::operator++(int)
{
   if (d_patch<d_number_patches) {
      ++d_patch;
   }
}

template<int DIM> inline
bool PatchLevelIterator<DIM>::operator==(
   const PatchLevelIterator<DIM>& iter) const
{
   return(d_patch == iter.d_patch);
}

template<int DIM> inline
bool PatchLevelIterator<DIM>::operator!=(
   const PatchLevelIterator<DIM>& iter) const
{
   return(d_patch != iter.d_patch);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}
