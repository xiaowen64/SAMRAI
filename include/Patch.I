//
// File:	Patch.I
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Patch container class for patch data objects
//

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif
#ifndef included_hier_PatchData
#include "PatchData.h"
#endif
#ifndef included_hier_VariableDatabase
#include "VariableDatabase.h"
#endif

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
   assert((id >= 0) && (id < d_patch_data.getSize()));
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
   assert((id >= 0) && (id < d_patch_data.getSize()));
#endif
   return(d_patch_data[id]);
}

template<int DIM> inline
void Patch<DIM>::setPatchData(
   const int id, 
   tbox::Pointer< PatchData<DIM> > data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((id >= 0) && (id < d_patch_data.getSize()));
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
   assert((id >= 0) && (id < d_patch_data.getSize()));
   assert(!d_patch_data[id].isNull());
#endif
   d_patch_data[id]->setTime(timestamp);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
