//
// File:	CellDataFactory.C
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Factory class for creating cell data objects
//

#ifndef included_pdat_CellDataFactory_C
#define included_pdat_CellDataFactory_C

#include "CellDataFactory.h"
#include "tbox/ArenaManager.h"
#include "tbox/Utilities.h"
#include "CellData.h"
#include "CellGeometry.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

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
:  hier::PatchDataFactory<DIM>(),
   d_depth(depth),
   d_ghosts(ghosts)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(depth > 0);
   assert(ghosts.min() >= 0);
#endif
}

template<int DIM, class TYPE>
CellDataFactory<DIM,TYPE>::~CellDataFactory()
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
CellDataFactory<DIM,TYPE>::cloneFactory()
{
   return(new CellDataFactory<DIM,TYPE>(d_depth, d_ghosts));
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
      new (pool) CellData<DIM,TYPE>(box, d_depth, d_ghosts, pool);
   return(tbox::Pointer< hier::PatchData<DIM> >(patchdata, pool));
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
   hier::BoxGeometry<DIM> *boxgeometry = new CellGeometry<DIM>(box, d_ghosts);
   return(tbox::Pointer< hier::BoxGeometry<DIM> >(boxgeometry));
}

/*
*************************************************************************
*									*
* Get and set the default ghost cell widths for the cell data objects	*
* created with this factory.						*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
const hier::IntVector<DIM>&
CellDataFactory<DIM,TYPE>::getDefaultGhostCellWidth() const
{
   return(d_ghosts);
}

template<int DIM, class TYPE>
void CellDataFactory<DIM,TYPE>::setDefaultGhostCellWidth(
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
size_t CellDataFactory<DIM,TYPE>::getSizeOfMemory(const hier::Box<DIM>& box) const
{
   const size_t obj =
      tbox::Arena::align(sizeof(CellData<DIM,TYPE>));
   const size_t data =
      CellData<DIM,TYPE>::getSizeOfData(box, d_depth, d_ghosts);
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
