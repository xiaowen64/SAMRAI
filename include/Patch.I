//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/patches/Patch.I $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Patch container class for patch data objects
//

#include "PatchData.h"
#include "VariableDatabase.h"

#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace hier {

template<int DIM> inline
const Box<DIM>& Patch<DIM>::getBox() const
{
   return(d_box);
}

template<int DIM> inline
int Patch<DIM>::getPatchNumber() const
{
   return(d_patch_number);
}

template<int DIM> inline
void Patch<DIM>::setPatchNumber(const int p)
{
   d_patch_number = p;
}

template<int DIM> inline
void Patch<DIM>::setPatchLevelNumber(const int level_number)
{
   d_patch_level_number = level_number;
}

template<int DIM> inline
int Patch<DIM>::getPatchLevelNumber() const
{
   return(d_patch_level_number);
}

template<int DIM> inline
void Patch<DIM>::setPatchInHierarchy(bool in_hierarchy)
{
   d_patch_in_hierarchy = in_hierarchy;
}

template<int DIM> inline
bool Patch<DIM>::inHierarchy() const
{
   return(d_patch_in_hierarchy);
}

template<int DIM> inline
tbox::Pointer< PatchDescriptor<DIM> > Patch<DIM>::getPatchDescriptor() const
{
   return(d_descriptor);
}

template<int DIM> inline
tbox::Pointer< PatchData<DIM> > Patch<DIM>::getPatchData(const int id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((id >= 0) && (id < d_patch_data.getSize()));
#endif
   return(d_patch_data[id]);
}

template<int DIM> inline
tbox::Pointer< PatchData<DIM> > Patch<DIM>::getPatchData(
   const tbox::Pointer< Variable<DIM> > variable,
   const tbox::Pointer<VariableContext> context) const
{
   int id = VariableDatabase<DIM>::getDatabase()->
            mapVariableAndContextToIndex(variable, context);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((id >= 0) && (id < d_patch_data.getSize()));
#endif
   return(d_patch_data[id]);
}

template<int DIM> inline
void Patch<DIM>::setPatchData(
   const int id, 
   tbox::Pointer< PatchData<DIM> > data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((id >= 0) && (id < d_patch_data.getSize()));
#endif
   d_patch_data[id] = data;
}

template<int DIM> inline
bool Patch<DIM>::checkAllocated(const int id) const
{
   return((id < d_patch_data.getSize()) && !(d_patch_data[id].isNull()));
}

template<int DIM> inline
void Patch<DIM>::setPatchGeometry(tbox::Pointer< PatchGeometry<DIM> > geometry)
{
   d_patch_geometry = geometry;
}

template<int DIM> inline
tbox::Pointer< PatchGeometry<DIM> > Patch<DIM>::getPatchGeometry() const
{
   return(d_patch_geometry);
}

template<int DIM> inline
void Patch<DIM>::setTime(const double timestamp, const int id)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((id >= 0) && (id < d_patch_data.getSize()));
   TBOX_ASSERT(!d_patch_data[id].isNull());
#endif
   d_patch_data[id]->setTime(timestamp);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
