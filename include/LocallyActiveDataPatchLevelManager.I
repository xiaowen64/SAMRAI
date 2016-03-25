//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/hierarchy/variables/LocallyActiveDataPatchLevelManager.I $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Class for managing locally-active data on a single patch level.
//


#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
    namespace hier {

/*
*************************************************************************
*                                                                       *
* Functions to check which level object was initialized with.           *
*                                                                       *
*************************************************************************
*/

template<int DIM> inline
bool LocallyActiveDataPatchLevelManager<DIM>::checkLevel(
   const hier::PatchLevel<DIM>& level) const
{
   return(&level == d_patch_level.getPointer());
}

template<int DIM> inline
bool LocallyActiveDataPatchLevelManager<DIM>::checkLevel(
   const tbox::Pointer< hier::PatchLevel<DIM> > level) const
{
   return(level.getPointer() == d_patch_level.getPointer());
}

template<int DIM> inline
bool LocallyActiveDataPatchLevelManager<DIM>::checkLevel(
   const hier::PatchLevel<DIM>* level) const
{
   return(level == d_patch_level.getPointer());
}

template<int DIM> inline
tbox::Pointer< hier::PatchLevel<DIM> > 
LocallyActiveDataPatchLevelManager<DIM>::getPatchLevel() const
{
   return(d_patch_level);
}

/*
*************************************************************************
*                                                                       *
* Accessory functions to access active patch data information.          *
*                                                                       *
*************************************************************************
*/

template<int DIM> inline
bool 
LocallyActiveDataPatchLevelManager<DIM>::getPatchDataActive(
   const PatchDataId& patch_data_id,
   const PatchNumber& patch_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_patch_level.isNull());
   TBOX_ASSERT((patch_number.pn >= 0) && (patch_number.pn < d_number_patches));
#endif

   bool ret_val = false;

   if (!d_patch_level.isNull()) {
      ret_val = d_active_patch_data[ patch_number.pn ]->isSet( patch_data_id.pd );
   }

   return(ret_val);

}

template<int DIM> inline
const hier::ComponentSelector&
LocallyActiveDataPatchLevelManager<DIM>::getAllPatchDataActive(
   const PatchNumber& patch_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_patch_level.isNull());
   TBOX_ASSERT((patch_number.pn >= 0) && (patch_number.pn < d_number_patches));
#endif
   if (!d_patch_level.isNull()) {
      return(*d_active_patch_data[ patch_number.pn ]);
   } else {
      return(*(new hier::ComponentSelector()));
   }
}

/*
*************************************************************************
*                                                                       *
* Accessory functions to set active patch data information.             *
*                                                                       *
*************************************************************************
*/

template<int DIM> inline
void 
LocallyActiveDataPatchLevelManager<DIM>::setPatchDataActive(
   const PatchDataId& patch_data_id,
   const PatchNumber& patch_number)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_patch_level.isNull());
   TBOX_ASSERT((patch_number.pn >= 0) && (patch_number.pn < d_number_patches));
#endif
   if (!d_patch_level.isNull()) {
      d_active_patch_data[ patch_number.pn ]->setFlag( patch_data_id.pd );
   }
}

template<int DIM> inline
void 
LocallyActiveDataPatchLevelManager<DIM>::setPatchDataActive(
   const hier::ComponentSelector& active_indices,
   const PatchNumber& patch_number)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_patch_level.isNull());
   TBOX_ASSERT((patch_number.pn >= 0) && (patch_number.pn < d_number_patches));
#endif
   if (!d_patch_level.isNull()) {
      tbox::Pointer<hier::ComponentSelector> active_data_copy =
         new hier::ComponentSelector(active_indices);
      d_active_patch_data[ patch_number.pn ] = active_data_copy;
   }
}

template<int DIM> inline 
void
LocallyActiveDataPatchLevelManager<DIM>::setAllPatchDataActive(
   const PatchNumber& patch_number)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_patch_level.isNull());
   TBOX_ASSERT((patch_number.pn >= 0) && (patch_number.pn < d_number_patches));
#endif
   if (!d_patch_level.isNull()) {
      d_active_patch_data[ patch_number.pn ]->setAllFlags();
   }
}

template<int DIM> inline 
void
LocallyActiveDataPatchLevelManager<DIM>::setPatchDataActive(
   const hier::ComponentSelector& active_indices)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_patch_level.isNull());
#endif
   if (!d_patch_level.isNull()) {
      tbox::Pointer<hier::ComponentSelector> active_data_copy = 
         new hier::ComponentSelector(active_indices);
      for (int ip = 0; ip < d_number_patches; ip++) {
         d_active_patch_data[ip] = active_data_copy;
      }
   }
}

/*
*************************************************************************
*                                                                       *
* Accessory functions to set inactive patch data information.           *
*                                                                       *
*************************************************************************
*/

template<int DIM> inline
void
LocallyActiveDataPatchLevelManager<DIM>::setPatchDataInactive(
   const PatchDataId& patch_data_id,
   const PatchNumber& patch_number)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_patch_level.isNull());
   TBOX_ASSERT((patch_number.pn >= 0) && (patch_number.pn < d_number_patches));
#endif
   if (!d_patch_level.isNull()) {
      d_active_patch_data[ patch_number.pn ]->clrFlag( patch_data_id.pd );
   }
}

template<int DIM> inline
void
LocallyActiveDataPatchLevelManager<DIM>::setAllPatchDataInactive(
   const PatchNumber& patch_number)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_patch_level.isNull());
   TBOX_ASSERT((patch_number.pn >= 0) && (patch_number.pn < d_number_patches));
#endif
   if (!d_patch_level.isNull()) {
      d_active_patch_data[ patch_number.pn ]->clrAllFlags();
   }
}

template<int DIM> inline
void 
LocallyActiveDataPatchLevelManager<DIM>::setAllPatchDataInactive()
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_patch_level.isNull());
#endif
   if (!d_patch_level.isNull()) {
      for (int ip = 0; ip < d_number_patches; ip++) {
         d_active_patch_data[ip]->clrAllFlags();
      }
   }
}

template<int DIM> inline
void 
LocallyActiveDataPatchLevelManager<DIM>::clearAllActiveDataInfo()
{
   d_patch_level.setNull();
   d_number_patches = -1;
   d_active_patch_data.resizeArray(0);
}

/*
*************************************************************************
*                                                                       *
* Accessory functions to allocate/deallocate active patch data.         *
*                                                                       *
*************************************************************************
*/

template<int DIM> inline
void 
LocallyActiveDataPatchLevelManager<DIM>::allocateAllPatchData(
   double timestamp,
   tbox::Pointer<tbox::Arena> pool) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_patch_level.isNull());
#endif
   if (!d_patch_level.isNull()) {
      for (typename hier::PatchLevel<DIM>::Iterator ip(d_patch_level); ip; ip++) {
         d_patch_level->getPatch(ip())->
         allocatePatchData(*(d_active_patch_data[ip()]), timestamp, pool);
      }
   }
}

template<int DIM> inline
void
LocallyActiveDataPatchLevelManager<DIM>::allocateAllPatchData(
   const PatchNumber& patch_number,
   double timestamp,
   tbox::Pointer<tbox::Arena> pool) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_patch_level.isNull());
   TBOX_ASSERT((patch_number.pn >= 0) && (patch_number.pn < d_number_patches));
#endif
   if (!d_patch_level.isNull()) {
      d_patch_level->getPatch( patch_number.pn )->
         allocatePatchData(*(d_active_patch_data[ patch_number.pn ]), 
                           timestamp, pool);
   }
}

template<int DIM> inline
void 
LocallyActiveDataPatchLevelManager<DIM>::deallocateAllPatchData() const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_patch_level.isNull());
#endif
   if (!d_patch_level.isNull()) {
      for (typename hier::PatchLevel<DIM>::Iterator ip(d_patch_level); ip; ip++) {
         d_patch_level->getPatch(ip())->
           deallocatePatchData(*(d_active_patch_data[ip()]));
      }
   }
}

template<int DIM> inline
void
LocallyActiveDataPatchLevelManager<DIM>::deallocateAllPatchData(
   const PatchNumber& patch_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_patch_level.isNull());
   TBOX_ASSERT((patch_number.pn >= 0) && (patch_number.pn < d_number_patches));
#endif
   if (!d_patch_level.isNull()) {
      d_patch_level->getPatch( patch_number.pn )->
         deallocatePatchData(*(d_active_patch_data[ patch_number.pn ]));
   }
}

/*
*************************************************************************
*                                                                       *
* Functions to create patch level iterator for locally-active data.     *
*                                                                       *
*************************************************************************
*/

template<int DIM> inline
typename LocallyActiveDataPatchLevelManager<DIM>::Iterator
LocallyActiveDataPatchLevelManager<DIM>::getIterator(
   const PatchDataId& patch_data_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_patch_level.isNull());
#endif
   return(
      (LocallyActiveDataPatchLevelIterator<DIM>(patch_data_id,
                                                d_patch_level,
                                                &d_active_patch_data)) );
}

template<int DIM> inline
typename LocallyActiveDataPatchLevelManager<DIM>::Iterator
LocallyActiveDataPatchLevelManager<DIM>::getIterator(
   const tbox::Pointer< hier::Variable<DIM> > variable) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_patch_level.isNull());
#endif
   return(
      (LocallyActiveDataPatchLevelIterator<DIM>(variable,
                                                d_patch_level)) );
}

/*
*************************************************************************
*                                                                       *
* Implementation of hier::LocallyActiveDataPatchLevelIterator<DIM>.     *
*                                                                       *
*************************************************************************
*/

template<int DIM> inline
LocallyActiveDataPatchLevelIterator<DIM>::
LocallyActiveDataPatchLevelIterator()
{
   d_patch             = 1;
   d_data_index        = -1;
   d_active_patch_info = (tbox::Array< tbox::Pointer<hier::ComponentSelector> >*)NULL;
   d_number_patches    = 0;
   d_mapping           = (hier::ProcessorMapping*)NULL;
}

template<int DIM> inline
LocallyActiveDataPatchLevelIterator<DIM>::
LocallyActiveDataPatchLevelIterator(
   const LocallyActiveDataPatchLevelIterator<DIM>& iter)
{
   d_patch             = iter.d_patch;
   d_data_index        = iter.d_data_index;
   d_active_patch_info = iter.d_active_patch_info;
   d_number_patches    = iter.d_number_patches;
   d_mapping           = iter.d_mapping;
}

template<int DIM> inline
LocallyActiveDataPatchLevelIterator<DIM>&
LocallyActiveDataPatchLevelIterator<DIM>::operator=(
   const LocallyActiveDataPatchLevelIterator<DIM>& iter)
{
   d_patch             = iter.d_patch;
   d_data_index        = iter.d_data_index;
   d_active_patch_info = iter.d_active_patch_info;
   d_number_patches    = iter.d_number_patches;
   d_mapping           = iter.d_mapping;
   return(*this);
}

template<int DIM> inline
int LocallyActiveDataPatchLevelIterator<DIM>::operator*() const
{
   return(d_patch);
}

template<int DIM> inline
int LocallyActiveDataPatchLevelIterator<DIM>::operator()() const
{
   return(d_patch);
}

template<int DIM> inline
LocallyActiveDataPatchLevelIterator<DIM>::operator bool() const
{
   return(d_patch < d_number_patches);
}

#ifndef LACKS_BOOL_VOID_RESOLUTION
template<int DIM> inline
LocallyActiveDataPatchLevelIterator<DIM>::operator const void*() const
{
   return(d_patch < d_number_patches ? this : NULL);
}
#endif

template<int DIM> inline
bool LocallyActiveDataPatchLevelIterator<DIM>::operator!() const
{
   return(d_patch >= d_number_patches);
}

template<int DIM> inline
bool LocallyActiveDataPatchLevelIterator<DIM>::operator==(
   const LocallyActiveDataPatchLevelIterator<DIM>& iter) const
{
   return(d_patch == iter.d_patch);
}

template<int DIM> inline
bool LocallyActiveDataPatchLevelIterator<DIM>::operator!=(
   const LocallyActiveDataPatchLevelIterator<DIM>& iter) const
{
   return(d_patch != iter.d_patch);
}

template<int DIM> inline
bool 
LocallyActiveDataPatchLevelIterator<DIM>::notActivePatch(int patch_number) const
{
   return( (d_patch<d_number_patches) &&
           !( d_mapping->isMappingLocal(d_patch) &&
              ((*d_active_patch_info)[d_patch]->isSet(d_data_index)) ) );
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}
