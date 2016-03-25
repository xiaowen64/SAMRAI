//
// File:	LocallyActiveVariableDatabase.I
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 595 $
// Modified:	$Date: 2005-08-30 20:41:20 -0700 (Tue, 30 Aug 2005) $
// Description:	Singleton database for variables defined on subset of hierarchy patches.
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
    namespace hier {

/*
*************************************************************************
*                                                                       *
* Inline accessory functions to search and return database contents.    *
*                                                                       *
*************************************************************************
*/
template<int DIM> inline
tbox::Pointer< PatchDescriptor<DIM> >
LocallyActiveVariableDatabase<DIM>::getPatchDescriptor() const
{
   return( VariableDatabase<DIM>::getPatchDescriptor() ); 
}

template<int DIM> inline
tbox::Pointer<hier::VariableContext> 
LocallyActiveVariableDatabase<DIM>::getSharedContext() const
{
   return(d_locally_active_context);
}

template<int DIM> inline 
bool LocallyActiveVariableDatabase<DIM>::checkVariableExists(
   const tbox::Pointer< hier::Variable<DIM> > variable) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!variable.isNull());
#endif
   return( (d_max_locally_active_variable_id >=
            variable->getInstanceIdentifier()) &&
           !d_locally_active_variables[variable->
            getInstanceIdentifier()].isNull() );
}

template<int DIM> inline
bool LocallyActiveVariableDatabase<DIM>::checkVariablePatchDataIndex(
   const tbox::Pointer< hier::Variable<DIM> > variable,
   int data_id)
{
   return( hier::VariableDatabase<DIM>::checkVariablePatchDataIndex(variable, data_id) ); 
}

template<int DIM> inline
int
LocallyActiveVariableDatabase<DIM>::mapVariableToIndex(
   const tbox::Pointer< hier::Variable<DIM> > variable) const
{
   int data_id = VariableDatabase<DIM>::idUndefined();
   if (checkVariableExists(variable)) {
      data_id = 
      d_locally_active_variable_descriptor_indices[
         variable->getInstanceIdentifier()];
   }
   return(data_id);
}

template<int DIM> inline
bool
LocallyActiveVariableDatabase<DIM>::mapIndexToVariable(
   const int index,
   tbox::Pointer< hier::Variable<DIM> >& variable) const
{
   return( hier::VariableDatabase<DIM>::mapIndexToVariable(index, variable) );
}

template<int DIM> inline
bool 
LocallyActiveVariableDatabase<DIM>::validLevel(
   const tbox::Pointer< hier::PatchLevel<DIM> > level) const
{
   return(!level.isNull() && 
          level->getLevelNumber() >= 0 &&
          level->inHierarchy());
}

template<int DIM> inline
const hier::LocallyActiveDataPatchLevelManager<DIM>*
LocallyActiveVariableDatabase<DIM>::getLevelManager(
   const hier::PatchLevel<DIM>& pl) const
{
   const int ln = pl.getLevelNumber();
   if ( ln >= 0 &&
        d_patch_level_active_data_manager.size() > ln &&
        !d_patch_level_active_data_manager[ln].isNull() &&
        d_patch_level_active_data_manager[ln]->checkLevel(pl) ) {
      return(d_patch_level_active_data_manager[ln]);
   } else {
      return( (hier::LocallyActiveDataPatchLevelManager<DIM>*)NULL );
   }
}

template<int DIM> inline
const hier::LocallyActiveDataPatchLevelManager<DIM>*
LocallyActiveVariableDatabase<DIM>::getLevelManager(
   const hier::PatchLevel<DIM>* pl) const
{
   const int ln = pl->getLevelNumber();
   if ( ln >= 0 &&
        d_patch_level_active_data_manager.size() > ln &&
        !d_patch_level_active_data_manager[ln].isNull() &&
        d_patch_level_active_data_manager[ln]->checkLevel(pl) ) {
      return(d_patch_level_active_data_manager[ln]);
   } else {
      return( (hier::LocallyActiveDataPatchLevelManager<DIM>*)NULL );
   }
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}
