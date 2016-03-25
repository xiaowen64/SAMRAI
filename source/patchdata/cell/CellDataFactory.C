//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/cell/CellDataFactory.C $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description: Factory class for creating cell data objects
//

#ifndef included_pdat_CellDataFactory_C
#define included_pdat_CellDataFactory_C

#include "CellDataFactory.h"
#include "tbox/ArenaManager.h"
#include "tbox/Utilities.h"
#include "CellData.h"
#include "CellGeometry.h"
#include "Patch.h"

#ifdef DEBUG_NO_INLINE
#include "CellDataFactory.I"
#endif
namespace SAMRAI {
    namespace pdat {

/*
*************************************************************************
*									*
* The constructor simply caches the default ghost cell width and depth.	*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
CellDataFactory<DIM,TYPE>::CellDataFactory(
   int depth,
   const hier::IntVector<DIM>& ghosts)
:  hier::PatchDataFactory<DIM>(ghosts),
   d_depth(depth),
   d_mb_trans(NULL)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(depth > 0);
   TBOX_ASSERT(ghosts.min() >= 0);
#endif
   d_mb_trans = NULL;
}

template<int DIM, class TYPE>
CellDataFactory<DIM,TYPE>::~CellDataFactory()
{
   if (d_mb_trans) {
      delete d_mb_trans;
   }
}

/*
*************************************************************************
*									*
* Clone the factory and copy the default parameters to the new factory.	*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
tbox::Pointer< hier::PatchDataFactory<DIM> >
CellDataFactory<DIM,TYPE>::cloneFactory(const hier::IntVector<DIM>& ghosts)
{
   return(new CellDataFactory<DIM,TYPE>(d_depth, ghosts));
}

/*
*************************************************************************
*									*
* Allocate the concrete cell data classes.  If no arena is specified,	*
* then the standard memory arena is used.				*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
tbox::Pointer< hier::PatchData<DIM> >
CellDataFactory<DIM,TYPE>::allocate(const hier::Box<DIM>& box,
                                      tbox::Pointer<tbox::Arena> pool) const
{
   if (pool.isNull()) {
      pool = tbox::ArenaManager::getManager()->getStandardAllocator();
   }

   hier::PatchData<DIM> *patchdata =
      new (pool) CellData<DIM,TYPE>(box, this -> d_depth, this -> d_ghosts, pool);
   return(tbox::Pointer< hier::PatchData<DIM> >(patchdata, pool));
}

template<int DIM, class TYPE>
tbox::Pointer< hier::PatchData<DIM> >
CellDataFactory<DIM,TYPE>::allocate(const hier::Patch<DIM>& patch,
                                    tbox::Pointer<tbox::Arena> pool) const
{
   return (allocate(patch.getBox(), pool));
}

/*
*************************************************************************
*									*
* Return the box geometry type for cell data objects.			*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
tbox::Pointer< hier::BoxGeometry<DIM> >
CellDataFactory<DIM,TYPE>::getBoxGeometry(const hier::Box<DIM>& box) const
{
   hier::BoxGeometry<DIM> *boxgeometry = new CellGeometry<DIM>(box, this -> d_ghosts);
   return(tbox::Pointer< hier::BoxGeometry<DIM> >(boxgeometry));
}

/*
*************************************************************************
*									*
* Calculate the amount of memory needed to allocate the data object.	*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
size_t CellDataFactory<DIM,TYPE>::getSizeOfMemory(const hier::Box<DIM>& box) const
{
   const size_t obj =
      tbox::Arena::align(sizeof(CellData<DIM,TYPE>));
   const size_t data =
      CellData<DIM,TYPE>::getSizeOfData(box, d_depth, this -> d_ghosts);
   return(obj+data);
}

/*
*************************************************************************
*									*
* Determine whether this is a valid copy operation to/from CellData     *
* between the supplied datatype.                                        *
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
bool CellDataFactory<DIM,TYPE>::validCopyTo(
   const tbox::Pointer<hier::PatchDataFactory<DIM> >& dst_pdf) const
{

   bool valid_copy = false;

   /*
    * Only valid option is CellData.
    */
   tbox::Pointer< CellDataFactory<DIM,TYPE> > cdf = dst_pdf;
   if (!cdf.isNull()) {
      valid_copy = true;
   }
   return(valid_copy);
}


}
}
#endif
