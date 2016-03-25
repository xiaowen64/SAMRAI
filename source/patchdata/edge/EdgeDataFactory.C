//
// File:	EdgeDataFactory.C
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Factory class for creating edge data objects
//

#ifndef included_pdat_EdgeDataFactory_C
#define included_pdat_EdgeDataFactory_C

#include "EdgeDataFactory.h"
#include "tbox/Arena.h"
#include "tbox/ArenaManager.h"
#include "tbox/Utilities.h"
#include "Box.h"
#include "EdgeData.h"
#include "EdgeGeometry.h"
#include "OuteredgeDataFactory.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

#ifdef DEBUG_NO_INLINE
#include "EdgeDataFactory.I"
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
EdgeDataFactory<DIM,TYPE>::EdgeDataFactory(
   int depth,
   const hier::IntVector<DIM>& ghosts,
   bool fine_boundary_represents_var)
:  hier::PatchDataFactory<DIM>(),
   d_depth(depth),
   d_ghosts(ghosts),
   d_fine_boundary_represents_var(fine_boundary_represents_var)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(depth > 0);
   assert(ghosts.min() >= 0);
#endif
}

template<int DIM, class TYPE>
EdgeDataFactory<DIM,TYPE>::~EdgeDataFactory()
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
EdgeDataFactory<DIM,TYPE>::cloneFactory()
{
   return(new EdgeDataFactory<DIM,TYPE>(d_depth, 
                                          d_ghosts,
                                          d_fine_boundary_represents_var));
}

/*
*************************************************************************
*									*
* Allocate the concrete edge data classes.  If no arena is specified,	*
* then the standard memory arena is used.				*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
tbox::Pointer< hier::PatchData<DIM> >
EdgeDataFactory<DIM,TYPE>::allocate(
   const hier::Box<DIM>& box, tbox::Pointer<tbox::Arena> pool) const
{
   if (pool.isNull()) {
      pool = tbox::ArenaManager::getManager()->getStandardAllocator();
   }

   hier::PatchData<DIM> *patchdata =
      new (pool) EdgeData<DIM,TYPE>(box, d_depth, d_ghosts, pool);
   return(tbox::Pointer< hier::PatchData<DIM> >(patchdata, pool));
}

/*
*************************************************************************
*									*
* Return the box geometry type for edge data objects.			*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
tbox::Pointer< hier::BoxGeometry<DIM> >
EdgeDataFactory<DIM,TYPE>::getBoxGeometry(const hier::Box<DIM>& box) const
{
   hier::BoxGeometry<DIM> *boxgeometry = new EdgeGeometry<DIM>(box, d_ghosts);
   return(tbox::Pointer< hier::BoxGeometry<DIM> >(boxgeometry));
}

/*
*************************************************************************
*									*
* Get and set the default ghost cell widths for the edge data objects	*
* created with this factory.						*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
const hier::IntVector<DIM>&
EdgeDataFactory<DIM,TYPE>::getDefaultGhostCellWidth() const
{
   return(d_ghosts);
}

template<int DIM, class TYPE>
void EdgeDataFactory<DIM,TYPE>::setDefaultGhostCellWidth(
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
size_t EdgeDataFactory<DIM,TYPE>::getSizeOfMemory(const hier::Box<DIM>& box) const
{
   const size_t obj =
      tbox::Arena::align(sizeof(EdgeData<DIM,TYPE>));
   const size_t data =
      EdgeData<DIM,TYPE>::getSizeOfData(box, d_depth, d_ghosts);
   return(obj+data);
}

/*
*************************************************************************
*									*
* Determine whether this is a valid copy operation to/from EdgeData     *
* between the supplied datatype.                                        *
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
bool EdgeDataFactory<DIM,TYPE>::validCopyTo(
   const tbox::Pointer<hier::PatchDataFactory<DIM> >& dst_pdf) const
{

   bool valid_copy = false;

   /*
    * Valid options are EdgeData and OuteredgeData.
    */
   if (!valid_copy) {
      tbox::Pointer< EdgeDataFactory<DIM,TYPE> > edf = dst_pdf;
      if (!edf.isNull()) {
         valid_copy = true;
      }
   }

   if (!valid_copy) {
      tbox::Pointer< OuteredgeDataFactory<DIM,TYPE> > oedf = dst_pdf;
      if (!oedf.isNull()) {
         valid_copy = true;
      }
   }

   return(valid_copy);
}


}
}
#endif
