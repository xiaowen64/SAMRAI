//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/hierarchy/variables/VariableDatabase.I $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Singleton database class for managing variables and contexts.
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace hier {

/*
*************************************************************************
*                                                                       *
* Inline accessory functions to retrieve data members.                  *
*                                                                       *
*************************************************************************
*/

template<int DIM> inline 
int VariableDatabase<DIM>::idUndefined() const
{
   return(-1);
}

template<int DIM> inline
tbox::Pointer< hier::PatchDescriptor<DIM> >
VariableDatabase<DIM>::getPatchDescriptor() const
{
   return(d_patch_descriptor);
}

template<int DIM> inline
int VariableDatabase<DIM>::getNumberRegisteredPatchDataIndices() const
{
   return(d_num_registered_patch_data_ids);
}

template<int DIM> inline
int VariableDatabase<DIM>::getNumberRegisteredVariableContexts() const
{
   // currently, we do not allow removal of variable contexts
   // so this suffices
   return(d_max_context_id+1);
}

/*
*************************************************************************
*                                                                       *
* Inline accessory functions to manage patch data ids for restart.      *
*                                                                       *
*************************************************************************
*/

template<int DIM> inline
ComponentSelector VariableDatabase<DIM>::getPatchDataRestartTable() const
{
   return(d_patchdata_restart_table);
}

template<int DIM> inline
bool VariableDatabase<DIM>::isPatchDataRegisteredForRestart(int index) const
{
   return(d_patchdata_restart_table.isSet(index));
}

template<int DIM> inline
void VariableDatabase<DIM>::registerPatchDataForRestart(int index)
{
   d_patchdata_restart_table.setFlag(index);
}

template<int DIM> inline
void VariableDatabase<DIM>::unregisterPatchDataForRestart(int index)
{
   d_patchdata_restart_table.clrFlag(index);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
