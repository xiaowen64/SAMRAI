//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-1/source/patchdata/outernode/OuternodeDataFactory.C $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Release:	$Name$
// Revision:	$LastChangedRevision: 1846 $
// Modified:	$LastChangedDate: 2008-01-11 09:51:05 -0800 (Fri, 11 Jan 2008) $
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
#include "Patch.h"

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
   :  hier::PatchDataFactory<DIM>(hier::IntVector<DIM>(0)),
      d_depth(depth),
      d_no_ghosts(hier::IntVector<DIM>(0))
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(depth > 0);
#endif
}

template <int DIM, class TYPE>
OuternodeDataFactory<DIM,TYPE>::~OuternodeDataFactory()
{
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
OuternodeDataFactory<DIM,TYPE>::cloneFactory(const hier::IntVector<DIM>& ghosts)
{
   NULL_USE(ghosts);
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

template<int DIM, class TYPE>
tbox::Pointer< hier::PatchData<DIM> >
OuternodeDataFactory<DIM,TYPE>::allocate(const hier::Patch<DIM>& patch,
                                         tbox::Pointer<tbox::Arena> pool) const
{
   return (allocate(patch.getBox(), pool));
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

