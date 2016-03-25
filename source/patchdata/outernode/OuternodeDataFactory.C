//
// File:	OuternodeDataFactory.C
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Release:	$Name$
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Factory class for creating outernode data objects
//

#ifndef included_pdat_OuternodeDataFactory_C
#define included_pdat_OuternodeDataFactory_C

#include "OuternodeDataFactory.h"
#include "tbox/Arena.h"
#include "tbox/ArenaManager.h"
#include "tbox/Utilities.h"
#include "Box.h"
#include "NodeDataFactory.h"
#include "OuternodeData.h"
#include "OuternodeGeometry.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

#ifdef DEBUG_NO_INLINE
#include "OuternodeDataFactory.I"
#endif

namespace SAMRAI {
    namespace pdat {

/*
*************************************************************************
*                                                                       *
* The constructor simply caches the default depth of the patch data.    *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
OuternodeDataFactory<DIM,TYPE>::OuternodeDataFactory(int depth)
:  hier::PatchDataFactory<DIM>(),
   d_depth(depth)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(depth > 0);
#endif

   d_no_ghosts = new hier::IntVector<DIM>(0);
}

template <int DIM, class TYPE>
OuternodeDataFactory<DIM,TYPE>::~OuternodeDataFactory()
{
   if(d_no_ghosts)
      delete d_no_ghosts;
   d_no_ghosts = NULL;
}

/*
*************************************************************************
*                                                                       *
* Clone the factory and copy the default parameters to the new factory. *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
tbox::Pointer<hier::PatchDataFactory<DIM> >
OuternodeDataFactory<DIM,TYPE>::cloneFactory()
{
   return(new OuternodeDataFactory<DIM,TYPE>(d_depth));
}

/*
*************************************************************************
*                                                                       *
* Allocate the concrete outernode data classes.  If no arena is         *
* specified, then the standard memory arena is used.                    *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
tbox::Pointer<hier::PatchData<DIM> > OuternodeDataFactory<DIM,TYPE>::allocate(
   const hier::Box<DIM>& box,
   tbox::Pointer<tbox::Arena> pool) const
{
   if (pool.isNull()) {
      pool = tbox::ArenaManager::getManager()->getStandardAllocator();
   }

   hier::PatchData<DIM> *patchdata =
      new (pool) OuternodeData<DIM,TYPE>(box, d_depth, pool);
   return(tbox::Pointer<hier::PatchData<DIM> >(patchdata, pool));
}

/*
*************************************************************************
*                                                                       *
* Return the box geometry type for outernode data objects.              *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
tbox::Pointer<hier::BoxGeometry<DIM> >
OuternodeDataFactory<DIM,TYPE>::getBoxGeometry(const hier::Box<DIM>& box) const
{
   hier::BoxGeometry<DIM> *boxgeometry = new OuternodeGeometry<DIM>(box, 0);
   return(tbox::Pointer<hier::BoxGeometry<DIM> >(boxgeometry));
}

/*
*************************************************************************
*                                                                       *
* Get and set the default ghost cell widths for the outernode data      *
* objects created with this factory.                                    *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
const hier::IntVector<DIM>&
OuternodeDataFactory<DIM,TYPE>::getDefaultGhostCellWidth() const
{
   return(*d_no_ghosts);
}

template <int DIM, class TYPE>
void OuternodeDataFactory<DIM,TYPE>::setDefaultGhostCellWidth(
   const hier::IntVector<DIM>& ghosts)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(ghosts == hier::IntVector<DIM>(0));
#endif
   (void) ghosts;
}

/*
*************************************************************************
*                                                                       *
* Calculate the amount of memory needed to allocate the data object.    *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
size_t OuternodeDataFactory<DIM,TYPE>::getSizeOfMemory(
   const hier::Box<DIM>& box) const
{
   const size_t obj = tbox::Arena::align(sizeof(OuternodeData<DIM,TYPE>));
   const size_t data = OuternodeData<DIM,TYPE>::getSizeOfData(box,
								d_depth);
   return(obj+data);
}

/*
*************************************************************************
*									*
* Determine whether this is a valid copy operation to/from OuternodeData*
* between the supplied datatype.                                        *
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
bool OuternodeDataFactory<DIM,TYPE>::validCopyTo(
   const tbox::Pointer<hier::PatchDataFactory<DIM> >& dst_pdf) const
{

   bool valid_copy = false;

   /*
    * Valid options are NodeData and OuternodeData.
    */
   if (!valid_copy) {
      tbox::Pointer< NodeDataFactory<DIM,TYPE> > ndf = dst_pdf;
      if (!ndf.isNull()) {
         valid_copy = true;
      }
   }

   if (!valid_copy) {
      tbox::Pointer< OuternodeDataFactory<DIM,TYPE> > ondf = dst_pdf;
      if (!ondf.isNull()) {
         valid_copy = true;
      }
   }

   return(valid_copy);
}


}
}
#endif

