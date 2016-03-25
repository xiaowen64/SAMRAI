//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/face/FaceDataFactory.C $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description: Factory class for creating face data objects
//

#ifndef included_pdat_FaceDataFactory_C
#define included_pdat_FaceDataFactory_C

#include "FaceDataFactory.h"
#include "tbox/Arena.h"
#include "tbox/ArenaManager.h"
#include "tbox/Utilities.h"
#include "Box.h"
#include "FaceData.h"
#include "FaceGeometry.h"
#include "OuterfaceDataFactory.h"
#include "Patch.h"


#ifdef DEBUG_NO_INLINE
#include "FaceDataFactory.I"
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
FaceDataFactory<DIM,TYPE>::FaceDataFactory(
   int depth,
   const hier::IntVector<DIM>& ghosts,
   bool fine_boundary_represents_var)
:  hier::PatchDataFactory<DIM>(ghosts),
   d_depth(depth),
   d_fine_boundary_represents_var(fine_boundary_represents_var),
   d_mb_trans(NULL)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(depth > 0);
   TBOX_ASSERT(ghosts.min() >= 0);
#endif
}

template<int DIM, class TYPE>
FaceDataFactory<DIM,TYPE>::~FaceDataFactory()
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
FaceDataFactory<DIM,TYPE>::cloneFactory(const hier::IntVector<DIM>& ghosts)
{
   return(new FaceDataFactory<DIM,TYPE>(d_depth, 
                                          ghosts,
                                          d_fine_boundary_represents_var));
}

/*
*************************************************************************
*									*
* Allocate the concrete face data classes.  If no arena is specified,	*
* then the standard memory arena is used.				*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
tbox::Pointer< hier::PatchData<DIM> >
FaceDataFactory<DIM,TYPE>::allocate(
   const hier::Box<DIM>& box, tbox::Pointer<tbox::Arena> pool) const
{
   if (pool.isNull()) {
      pool = tbox::ArenaManager::getManager()->getStandardAllocator();
   }

   hier::PatchData<DIM> *patchdata =
      new (pool) FaceData<DIM,TYPE>(box, d_depth, this -> d_ghosts, pool);
   return(tbox::Pointer< hier::PatchData<DIM> >(patchdata, pool));
}

template<int DIM, class TYPE>
tbox::Pointer< hier::PatchData<DIM> >
FaceDataFactory<DIM,TYPE>::allocate(const hier::Patch<DIM>& patch,
                                    tbox::Pointer<tbox::Arena> pool) const
{
   return (allocate(patch.getBox(), pool));
}

/*
*************************************************************************
*									*
* Return the box geometry type for face data objects.			*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
tbox::Pointer< hier::BoxGeometry<DIM> >
FaceDataFactory<DIM,TYPE>::getBoxGeometry(const hier::Box<DIM>& box) const
{
   hier::BoxGeometry<DIM> *boxgeometry = new FaceGeometry<DIM>(box, this -> d_ghosts);
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
size_t FaceDataFactory<DIM,TYPE>::getSizeOfMemory(const hier::Box<DIM>& box) const
{
   const size_t obj =
      tbox::Arena::align(sizeof(FaceData<DIM,TYPE>));
   const size_t data =
      FaceData<DIM,TYPE>::getSizeOfData(box, d_depth, this -> d_ghosts);
   return(obj+data);
}

/*
*************************************************************************
*									*
* Determine whether this is a valid copy operation to/from FaceData     *
* between the supplied datatype.                                        *
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
bool FaceDataFactory<DIM,TYPE>::validCopyTo(
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
