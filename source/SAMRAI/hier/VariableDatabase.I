/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Singleton database class for managing variables and contexts.
 *
 ************************************************************************/

namespace SAMRAI {
namespace hier {

/*
 *************************************************************************
 *
 * Inline accessory functions to retrieve data members.
 *
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
boost::shared_ptr<PatchDescriptor>
VariableDatabase::getPatchDescriptor() const
{
   return d_patch_descriptor;
}

SAMRAI_INLINE_KEYWORD
int VariableDatabase::getNumberOfRegisteredPatchDataIndices() const
{
   return d_num_registered_patch_data_ids;
}

SAMRAI_INLINE_KEYWORD
int VariableDatabase::getNumberOfRegisteredVariableContexts() const
{
   // currently, we do not allow removal of variable contexts
   // so this suffices
   return d_max_context_id + 1;
}

/*
 *************************************************************************
 *
 * Inline accessory functions to manage patch data ids for restart.
 *
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
ComponentSelector VariableDatabase::getPatchDataRestartTable() const
{
   return d_patchdata_restart_table;
}

SAMRAI_INLINE_KEYWORD
bool VariableDatabase::isPatchDataRegisteredForRestart(
   int index) const
{
   return d_patchdata_restart_table.isSet(index);
}

SAMRAI_INLINE_KEYWORD
void VariableDatabase::registerPatchDataForRestart(
   int index)
{
   d_patchdata_restart_table.setFlag(index);
}

SAMRAI_INLINE_KEYWORD
void VariableDatabase::unregisterPatchDataForRestart(
   int index)
{
   d_patchdata_restart_table.clrFlag(index);
}

}
}
