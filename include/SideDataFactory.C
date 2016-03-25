//
// File:	SideDataFactory.C
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
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

#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

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
:  hier::PatchDataFactory<DIM>(),
   d_depth(depth),
   d_ghosts(ghosts),
   d_fine_boundary_represents_var(fine_boundary_represents_var),
   d_directions(directions)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(depth > 0);
   assert(ghosts.min() >= 0);
   assert(directions.min() >= 0);
#endif
}

template<int DIM, class TYPE>
SideDataFactory<DIM,TYPE>::~SideDataFactory()
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
SideDataFactory<DIM,TYPE>::cloneFactory()
{
   return(new SideDataFactory<DIM,TYPE>(d_depth, 
                                          d_ghosts, 
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
                                      d_ghosts,
                                      d_directions,
                                      pool);
   return(tbox::Pointer< hier::PatchData<DIM> >(patchdata, pool));
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
                                                           d_ghosts,
                                                           d_directions);
   return(tbox::Pointer< hier::BoxGeometry<DIM> >(boxgeometry));
}

/*
*************************************************************************
*									*
* Get and set the default ghost cell widths for the side data objects	*
* created with this factory.						*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
const hier::IntVector<DIM>&
SideDataFactory<DIM,TYPE>::getDefaultGhostCellWidth() const
{
   return(d_ghosts);
}

template<int DIM, class TYPE>
void SideDataFactory<DIM,TYPE>::setDefaultGhostCellWidth(
   const hier::IntVector<DIM>& ghosts)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(ghosts.min() >= 0);
#endif
   d_ghosts = ghosts;
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
      SideData<DIM,TYPE>::getSizeOfData(box, d_depth, d_ghosts, d_directions);
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
