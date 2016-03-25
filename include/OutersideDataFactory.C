//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/outerside/OutersideDataFactory.C $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1776 $
// Modified:	$LastChangedDate: 2007-12-13 16:40:01 -0800 (Thu, 13 Dec 2007) $
// Description: Factory class for creating outerside data objects
//

#ifndef included_pdat_OutersideDataFactory_C
#define included_pdat_OutersideDataFactory_C

#include "OutersideDataFactory.h"
#include "tbox/Arena.h"
#include "tbox/ArenaManager.h"
#include "tbox/Utilities.h"
#include "Box.h"
#include "OutersideData.h"
#include "OutersideGeometry.h"
#include "Patch.h"
#include "SideDataFactory.h"


#ifdef DEBUG_NO_INLINE
#include "OutersideDataFactory.I"
#endif
namespace SAMRAI {
    namespace pdat {

/*
*************************************************************************
*									*
* The constructor simply caches the default depth of the patch data.	*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
OutersideDataFactory<DIM,TYPE>::OutersideDataFactory(int depth)
   :  hier::PatchDataFactory<DIM>(hier::IntVector<DIM>(0)),
      d_depth(depth),
      d_no_ghosts(hier::IntVector<DIM>(0))
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(depth > 0);
#endif
}

template<int DIM, class TYPE>
OutersideDataFactory<DIM,TYPE>::~OutersideDataFactory()
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
OutersideDataFactory<DIM,TYPE>::cloneFactory(const hier::IntVector<DIM>& ghosts)
{
   return(new OutersideDataFactory<DIM,TYPE>(d_depth));
}

/*
*************************************************************************
*									*
* Allocate the concrete outerside data classes.  If no arena is		*
* specified, then the standard memory arena is used.			*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
tbox::Pointer< hier::PatchData<DIM> > OutersideDataFactory<DIM,TYPE>::allocate(
   const hier::Box<DIM>& box,
   tbox::Pointer<tbox::Arena> pool) const
{
   if (pool.isNull()) {
      pool = tbox::ArenaManager::getManager()->getStandardAllocator();
   }

   hier::PatchData<DIM> *patchdata =
      new (pool) OutersideData<DIM,TYPE>(box, d_depth, pool);
   return(tbox::Pointer< hier::PatchData<DIM> >(patchdata, pool));
}

template<int DIM, class TYPE>
tbox::Pointer< hier::PatchData<DIM> >
OutersideDataFactory<DIM,TYPE>::allocate(const hier::Patch<DIM>& patch,
                                         tbox::Pointer<tbox::Arena> pool) const
{
   return (allocate(patch.getBox(), pool));
}

/*
*************************************************************************
*									*
* Return the box geometry type for outerside data objects.		*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
tbox::Pointer< hier::BoxGeometry<DIM> >
OutersideDataFactory<DIM,TYPE>::getBoxGeometry(const hier::Box<DIM>& box) const
{
   hier::BoxGeometry<DIM> *boxgeometry = new OutersideGeometry<DIM>(box, 0);
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
size_t OutersideDataFactory<DIM,TYPE>::getSizeOfMemory(
   const hier::Box<DIM>& box) const
{
   const size_t obj = tbox::Arena::align(sizeof(OutersideData<DIM,TYPE>));
   const size_t data = OutersideData<DIM,TYPE>::getSizeOfData(box, d_depth);
   return(obj+data);
}

/*
*************************************************************************
*									*
* Determine whether this is a valid copy operation to/from NodeData     *
* between the supplied datatype.                                        *
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
bool OutersideDataFactory<DIM,TYPE>::validCopyTo(
   const tbox::Pointer<hier::PatchDataFactory<DIM> >& dst_pdf) const
{

   bool valid_copy = false;

   /*
    * Valid options are SideData and OutersideData.
    */
   if (!valid_copy) {
      tbox::Pointer< SideDataFactory<DIM,TYPE> > sdf = dst_pdf;
      if (!sdf.isNull()) {
         valid_copy = true;
      }
   }

   if (!valid_copy) {
      tbox::Pointer< OutersideDataFactory<DIM,TYPE> > osdf = dst_pdf;
      if (!osdf.isNull()) {
         valid_copy = true;
      }
   }

   return(valid_copy);
}

}
}
#endif
