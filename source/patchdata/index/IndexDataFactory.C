//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/index/IndexDataFactory.C $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Release:	0.1
// Revision:	$LastChangedRevision: 1776 $
// Modified:	$LastChangedDate: 2007-12-13 16:40:01 -0800 (Thu, 13 Dec 2007) $
// Description: hier::Patch data factory for irregularly indexed patch data
//

#ifndef included_pdat_IndexDataFactory_C
#define included_pdat_IndexDataFactory_C

#include "IndexDataFactory.h"
#include "tbox/Arena.h"
#include "tbox/ArenaManager.h"
#include "tbox/Utilities.h"
#include "Box.h"
#include "CellGeometry.h"
#include "IndexData.h"

namespace SAMRAI {
    namespace pdat {

/*
*************************************************************************
*									*
* The constructor simply caches the default ghost cell width.		*
*									*
*************************************************************************
*/

template<int DIM, class TYPE> IndexDataFactory<DIM,TYPE>::IndexDataFactory(const hier::IntVector<DIM>& ghosts) 
   :   hier::PatchDataFactory<DIM>(ghosts)
{
}

template<int DIM, class TYPE> IndexDataFactory<DIM,TYPE>::~IndexDataFactory()
{
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
IndexDataFactory<DIM,TYPE>::cloneFactory(const hier::IntVector<DIM>& ghosts)
{
   return(new IndexDataFactory<DIM,TYPE>(ghosts));
}

/*
*************************************************************************
*									*
* Allocate the concrete irregular data class.  If no arena is given,	*
* then use the standard memory arena.					*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
tbox::Pointer< hier::PatchData<DIM> >
IndexDataFactory<DIM,TYPE>::allocate(
   const hier::Box<DIM>& box, tbox::Pointer<tbox::Arena> pool) const
{
   if (pool.isNull()) {
      pool = tbox::ArenaManager::getManager()->getStandardAllocator();
   }
   hier::PatchData<DIM> *pd = new (pool) IndexData<DIM,TYPE>(box, this -> d_ghosts);
   return(tbox::Pointer< hier::PatchData<DIM> >(pd, pool));
}

template<int DIM, class TYPE>
tbox::Pointer< hier::PatchData<DIM> >
IndexDataFactory<DIM,TYPE>::allocate(const hier::Patch<DIM>& patch,
                                     tbox::Pointer<tbox::Arena> pool) const
{
   return (allocate(patch.getBox(), pool));
}


/*
*************************************************************************
*									*
* Return the box geometry type for index data objects.			*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
tbox::Pointer< hier::BoxGeometry<DIM> >
IndexDataFactory<DIM,TYPE>::getBoxGeometry(const hier::Box<DIM>& box) const
{
   hier::BoxGeometry<DIM> *boxgeometry = new CellGeometry<DIM>(box, this -> d_ghosts);
   return(tbox::Pointer< hier::BoxGeometry<DIM> >(boxgeometry));
}

/*
*************************************************************************
*									*
* Calculate the amount of memory needed to allocate the object.		*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
size_t IndexDataFactory<DIM,TYPE>::getSizeOfMemory(const hier::Box<DIM>& box) const
{
   NULL_USE(box);
   return(tbox::Arena::align(sizeof(IndexData<DIM,TYPE>)));
}

/*
*************************************************************************
*									*
* Determine whether this is a valid copy operation to/from IndexData    *
* between the supplied datatype.                                        *
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
bool IndexDataFactory<DIM,TYPE>::validCopyTo(
   const tbox::Pointer<hier::PatchDataFactory<DIM> >& dst_pdf) const
{

   bool valid_copy = false;

   /*
    * Valid option is another IndexData object of the same dimension 
    * and type.
    */
   if (!valid_copy) {
      tbox::Pointer< IndexDataFactory<DIM,TYPE> > idf = dst_pdf;
      if (!idf.isNull()) {
         valid_copy = true;
      }
   }

   return(valid_copy);
}

}
}
#endif
