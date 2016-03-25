//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/outerface/OuterfaceDataFactory.C $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description: Factory class for creating outerface data objects
//

#ifndef included_pdat_OuterfaceDataFactory_C
#define included_pdat_OuterfaceDataFactory_C

#include "OuterfaceDataFactory.h"
#include "tbox/Arena.h"
#include "tbox/ArenaManager.h"
#include "tbox/Utilities.h"
#include "Box.h"
#include "OuterfaceData.h"
#include "OuterfaceGeometry.h"
#include "Patch.h"
#include "FaceDataFactory.h"


#ifdef DEBUG_NO_INLINE
#include "OuterfaceDataFactory.I"
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
OuterfaceDataFactory<DIM,TYPE>::OuterfaceDataFactory(int depth)
:  hier::PatchDataFactory<DIM>(hier::IntVector<DIM>(0)),
   d_depth(depth),
   d_no_ghosts(hier::IntVector<DIM>(0))
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(depth > 0);
#endif


}

template<int DIM, class TYPE>
OuterfaceDataFactory<DIM,TYPE>::~OuterfaceDataFactory()
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
OuterfaceDataFactory<DIM,TYPE>::cloneFactory(const hier::IntVector<DIM>& ghosts)
{
   NULL_USE(ghosts);
   return(new OuterfaceDataFactory<DIM,TYPE>(d_depth));
}

/*
*************************************************************************
*									*
* Allocate the concrete outerface data classes.  If no arena is		*
* specified, then the standard memory arena is used.			*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
tbox::Pointer< hier::PatchData<DIM> > OuterfaceDataFactory<DIM,TYPE>::allocate(
   const hier::Box<DIM>& box,
   tbox::Pointer<tbox::Arena> pool) const
{
   if (pool.isNull()) {
      pool = tbox::ArenaManager::getManager()->getStandardAllocator();
   }

   hier::PatchData<DIM> *patchdata =
      new (pool) OuterfaceData<DIM,TYPE>(box, d_depth, pool);
   return(tbox::Pointer< hier::PatchData<DIM> >(patchdata, pool));
}

template<int DIM, class TYPE>
tbox::Pointer< hier::PatchData<DIM> >
OuterfaceDataFactory<DIM,TYPE>::allocate(const hier::Patch<DIM>& patch,
                                         tbox::Pointer<tbox::Arena> pool) const
{
   return (allocate(patch.getBox(), pool));
}


/*
*************************************************************************
*									*
* Return the box geometry type for outerface data objects.		*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
tbox::Pointer< hier::BoxGeometry<DIM> >
OuterfaceDataFactory<DIM,TYPE>::getBoxGeometry(const hier::Box<DIM>& box) const
{
   hier::BoxGeometry<DIM> *boxgeometry = new OuterfaceGeometry<DIM>(box, 0);
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
size_t OuterfaceDataFactory<DIM,TYPE>::getSizeOfMemory(
   const hier::Box<DIM>& box) const
{
   const size_t obj = tbox::Arena::align(sizeof(OuterfaceData<DIM,TYPE>));
   const size_t data = OuterfaceData<DIM,TYPE>::getSizeOfData(box, d_depth);
   return(obj+data);
}

/*
*************************************************************************
*									*
* Determine whether this is a valid copy operation to/from OuterfaceData*
* between the supplied datatype.                                        *
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
bool OuterfaceDataFactory<DIM,TYPE>::validCopyTo(
   const tbox::Pointer<hier::PatchDataFactory<DIM> >& dst_pdf) const
{

   bool valid_copy = false;

   /*
    * Valid options are FaceData and OuterfaceData.
    */
   if (!valid_copy) {
      tbox::Pointer< FaceDataFactory<DIM,TYPE> > fdf = dst_pdf;
      if (!fdf.isNull()) {
         valid_copy = true;
      }
   }

   if (!valid_copy) {
      tbox::Pointer< OuterfaceDataFactory<DIM,TYPE> > ofdf = dst_pdf;
      if (!ofdf.isNull()) {
         valid_copy = true;
      }
   }

   return(valid_copy);
}

}
}
#endif
