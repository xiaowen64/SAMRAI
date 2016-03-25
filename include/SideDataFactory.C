//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/side/SideDataFactory.C $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1776 $
// Modified:	$LastChangedDate: 2007-12-13 16:40:01 -0800 (Thu, 13 Dec 2007) $
// Description: Factory class for creating side data objects
//

#ifndef included_pdat_SideDataFactory_C
#define included_pdat_SideDataFactory_C

#include "SideDataFactory.h"
#include "tbox/Arena.h"
#include "tbox/ArenaManager.h"
#include "Box.h"
#include "tbox/Utilities.h"
#include "SideData.h"
#include "SideGeometry.h"
#include "OutersideDataFactory.h"
#include "Patch.h"


#ifdef DEBUG_NO_INLINE
#include "SideDataFactory.I"
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
SideDataFactory<DIM,TYPE>::SideDataFactory(
   const int depth,
   const hier::IntVector<DIM>& ghosts,
   bool fine_boundary_represents_var,
   const hier::IntVector<DIM>& directions)
:  hier::PatchDataFactory<DIM>(ghosts),
   d_depth(depth),
   d_fine_boundary_represents_var(fine_boundary_represents_var),
   d_directions(directions),
   d_mb_trans(NULL)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(depth > 0);
   TBOX_ASSERT(ghosts.min() >= 0);
   TBOX_ASSERT(directions.min() >= 0);
#endif
}

template<int DIM, class TYPE>
SideDataFactory<DIM,TYPE>::~SideDataFactory()
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
SideDataFactory<DIM,TYPE>::cloneFactory(const hier::IntVector<DIM>& ghosts)
{
   return(new SideDataFactory<DIM,TYPE>(d_depth, 
                                          ghosts, 
                                          d_fine_boundary_represents_var,
                                          d_directions));
}

/*
*************************************************************************
*									*
* Allocate the concrete side data classes.  If no arena is specified,	*
* then the standard memory arena is used.				*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
tbox::Pointer< hier::PatchData<DIM> >
SideDataFactory<DIM,TYPE>::allocate(
   const hier::Box<DIM>& box, tbox::Pointer<tbox::Arena> pool) const
{
   if (pool.isNull()) {
      pool = tbox::ArenaManager::getManager()->getStandardAllocator();
   }

   hier::PatchData<DIM> *patchdata =
      new (pool) SideData<DIM,TYPE>(box,
                                      d_depth,
                                      this -> d_ghosts,
                                      d_directions,
                                      pool);
   return(tbox::Pointer< hier::PatchData<DIM> >(patchdata, pool));
}

template<int DIM, class TYPE>
tbox::Pointer< hier::PatchData<DIM> >
SideDataFactory<DIM,TYPE>::allocate(const hier::Patch<DIM>& patch,
                                    tbox::Pointer<tbox::Arena> pool) const
{
   return (allocate(patch.getBox(), pool));
}

/*
*************************************************************************
*									*
* Return the box geometry type for side data objects.			*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
tbox::Pointer< hier::BoxGeometry<DIM> >
SideDataFactory<DIM,TYPE>::getBoxGeometry(const hier::Box<DIM>& box) const
{
   hier::BoxGeometry<DIM> *boxgeometry = new SideGeometry<DIM>(box,
                                                           this -> d_ghosts,
                                                           d_directions);
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
size_t SideDataFactory<DIM,TYPE>::getSizeOfMemory(const hier::Box<DIM>& box) const
{
   const size_t obj =
      tbox::Arena::align(sizeof(SideData<DIM,TYPE>));
   const size_t data =
      SideData<DIM,TYPE>::getSizeOfData(box, d_depth, this -> d_ghosts, d_directions);
   return(obj+data);
}

/*
*************************************************************************
*									*
* Determine whether this is a valid copy operation to/from SideData     *
* between the supplied datatype.                                        *
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
bool SideDataFactory<DIM,TYPE>::validCopyTo(
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
