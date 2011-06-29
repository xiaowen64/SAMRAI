/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Patch container class for patch data objects 
 *
 ************************************************************************/

#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/hier/VariableDatabase.h"

namespace SAMRAI {
namespace hier {

SAMRAI_INLINE_KEYWORD
const MappedBox& Patch::getMappedBox() const
{
   return d_mapped_box;
}

SAMRAI_INLINE_KEYWORD
const GlobalId& Patch::getGlobalId() const
{
   return d_mapped_box.getGlobalId();
}

SAMRAI_INLINE_KEYWORD
const Box& Patch::getBox() const
{
   return d_mapped_box.getBox();
}

SAMRAI_INLINE_KEYWORD
const LocalId &Patch::getLocalId() const
{
   return d_mapped_box.getLocalId();
}

SAMRAI_INLINE_KEYWORD
void Patch::setPatchLevelNumber(
   const int level_number)
{
   d_patch_level_number = level_number;
}

SAMRAI_INLINE_KEYWORD
int Patch::getPatchLevelNumber() const
{
   return d_patch_level_number;
}

SAMRAI_INLINE_KEYWORD
void Patch::setPatchInHierarchy(
   bool in_hierarchy)
{
   d_patch_in_hierarchy = in_hierarchy;
}

SAMRAI_INLINE_KEYWORD
bool Patch::inHierarchy() const
{
   return d_patch_in_hierarchy;
}

SAMRAI_INLINE_KEYWORD
tbox::Pointer<PatchDescriptor> Patch::getPatchDescriptor() const
{
   return d_descriptor;
}

SAMRAI_INLINE_KEYWORD
tbox::Pointer<PatchData> Patch::getPatchData(
   const int id) const
{
   TBOX_ASSERT((id >= 0) && (id < d_patch_data.getSize()));

   return d_patch_data[id];
}

SAMRAI_INLINE_KEYWORD
tbox::Pointer<PatchData> Patch::getPatchData(
   const tbox::Pointer<Variable> variable,
   const tbox::Pointer<VariableContext> context) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *variable);

   int id = VariableDatabase::getDatabase()->
      mapVariableAndContextToIndex(variable, context);

   TBOX_ASSERT((id >= 0) && (id < d_patch_data.getSize()));

   return d_patch_data[id];
}

SAMRAI_INLINE_KEYWORD
void Patch::setPatchData(
   const int id,
   tbox::Pointer<PatchData> data)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *data);
   TBOX_ASSERT((id >= 0) && (id < d_patch_data.getSize()));

   d_patch_data[id] = data;
}

SAMRAI_INLINE_KEYWORD
bool Patch::checkAllocated(
   const int id) const
{
   return (id < d_patch_data.getSize()) && !(d_patch_data[id].isNull());
}

SAMRAI_INLINE_KEYWORD
void Patch::setPatchGeometry(
   tbox::Pointer<PatchGeometry> geometry)
{
   d_patch_geometry = geometry;
}

SAMRAI_INLINE_KEYWORD
tbox::Pointer<PatchGeometry> Patch::getPatchGeometry() const
{
   return d_patch_geometry;
}

SAMRAI_INLINE_KEYWORD
void Patch::setTime(
   const double timestamp,
   const int id)
{
   TBOX_ASSERT((id >= 0) && (id < d_patch_data.getSize()));
   TBOX_ASSERT(!d_patch_data[id].isNull());

   d_patch_data[id]->setTime(timestamp);
}

SAMRAI_INLINE_KEYWORD
const tbox::Dimension& Patch::getDim() const
{
   return d_mapped_box.getDim();
}

}
}
