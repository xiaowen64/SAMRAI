//
// File:	VariableDatabase.I
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 407 $
// Modified:	$Date: 2005-06-01 09:49:54 -0700 (Wed, 01 Jun 2005) $
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

/*
*************************************************************************
*                                                                       *
* Inline accessory functions to manage patch data ids for restart.      *
*                                                                       *
*************************************************************************
*/

template<int DIM> inline
ComponentSelector VariableDatabase<DIM>::getPatchDataRestartTable()
{
   return(d_patchdata_restart_table);
}

template<int DIM> inline
bool VariableDatabase<DIM>::isPatchDataRegisteredForRestart(int index)
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
