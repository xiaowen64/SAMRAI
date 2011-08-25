/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Registry of MappedBoxLevelHandles incident from a common MappedBoxLevel.
 *
 ************************************************************************/
#ifndef included_hier_MappedBoxLevelHandle_C
#define included_hier_MappedBoxLevelHandle_C

#include "SAMRAI/hier/MappedBoxLevelHandle.h"
#include "SAMRAI/hier/MappedBoxLevel.h"

namespace SAMRAI {
namespace hier {

/*
 ************************************************************************
 ************************************************************************
 */
MappedBoxLevelHandle::MappedBoxLevelHandle(
   const MappedBoxLevel* mapped_box_level):
   d_mapped_box_level(mapped_box_level)
{
}

/*
 ************************************************************************
 ************************************************************************
 */
MappedBoxLevelHandle::~MappedBoxLevelHandle()
{
   detachMyMappedBoxLevel();
}

/*
 ************************************************************************
 ************************************************************************
 */
bool MappedBoxLevelHandle::isAttached() const
{
   return d_mapped_box_level != NULL;
}

/*
 ************************************************************************
 ************************************************************************
 */
const MappedBoxLevel& MappedBoxLevelHandle::getMappedBoxLevel() const
{
   if (d_mapped_box_level == NULL) {
      TBOX_ERROR(
         "MappedBoxLevelHandle::getMappedBoxLevel Attempted to access a MappedBoxLevel\n"
         << "that has been detached from its handle.  Detachment happens\n"
         << "when the MappedBoxLevel changes in a way that can invalidate\n"
         << "Connector data.  Therefore, Connectors should not attempt\n"
         << "to access the MappedBoxLevel using a detatched handle.");
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   // Sanity check: the MappedBoxLevel's handle should be this handle.
   if (d_mapped_box_level->getMappedBoxLevelHandle().getPointer() != this) {
      TBOX_ERROR("Library error in MappedBoxLevelHandle::getMappedBoxLevel");
   }
#endif
   return *d_mapped_box_level;
}

/*
 ************************************************************************
 * To be called by the object's MappedBoxLevel when the
 * MappedBoxLevel changes in a way that can invalidate the
 * Connector data.
 ************************************************************************
 */
void MappedBoxLevelHandle::detachMyMappedBoxLevel()
{
   d_mapped_box_level = NULL;
}

}
}
#endif
