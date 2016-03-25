//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/outeredge/OuteredgeDataFactory.C $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1973 $
// Modified:	$LastChangedDate: 2008-02-11 16:39:15 -0800 (Mon, 11 Feb 2008) $
// Description: Factory class for creating outeredge data objects
//

#ifndef included_pdat_OuteredgeDataFactory_C
#define included_pdat_OuteredgeDataFactory_C

#include "EdgeDataFactory.h"
#include "tbox/Arena.h"
#include "tbox/ArenaManager.h"
#include "tbox/Utilities.h"
#include "Box.h"
#include "OuteredgeData.h"
#include "OuteredgeDataFactory.h"
#include "OuteredgeGeometry.h"
#include "Patch.h"


#ifdef DEBUG_NO_INLINE
#include "OuteredgeDataFactory.I"
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
OuteredgeDataFactory<DIM,TYPE>::OuteredgeDataFactory(int depth)
:  hier::PatchDataFactory<DIM>(hier::IntVector<DIM>(0)),
   d_depth(depth),
   d_no_ghosts(hier::IntVector<DIM>(0))
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(depth > 0);
#endif
}

template <int DIM, class TYPE>
OuteredgeDataFactory<DIM,TYPE>::~OuteredgeDataFactory()
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
OuteredgeDataFactory<DIM,TYPE>::cloneFactory(const hier::IntVector<DIM>& ghosts)
{
   NULL_USE(ghosts);
   return(new OuteredgeDataFactory<DIM,TYPE>(d_depth));
}

/*
*************************************************************************
*                                                                       *
* Allocate the concrete outeredge data classes.  If no arena is         *
* specified, then the standard memory arena is used.                    *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
tbox::Pointer<hier::PatchData<DIM> > OuteredgeDataFactory<DIM,TYPE>::allocate(
   const hier::Box<DIM>& box,
   tbox::Pointer<tbox::Arena> pool) const
{
   if (pool.isNull()) {
      pool = tbox::ArenaManager::getManager()->getStandardAllocator();
   }

   hier::PatchData<DIM> *patchdata =
      new (pool) OuteredgeData<DIM,TYPE>(box, d_depth, pool);
   return(tbox::Pointer<hier::PatchData<DIM> >(patchdata, pool));
}

template<int DIM, class TYPE>
tbox::Pointer< hier::PatchData<DIM> >
OuteredgeDataFactory<DIM,TYPE>::allocate(const hier::Patch<DIM>& patch,
                                         tbox::Pointer<tbox::Arena> pool) const
{
   return (allocate(patch.getBox(), pool));
}

/*
*************************************************************************
*                                                                       *
* Return the box geometry type for outeredge data objects.              *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
tbox::Pointer<hier::BoxGeometry<DIM> >
OuteredgeDataFactory<DIM,TYPE>::getBoxGeometry(const hier::Box<DIM>& box) const
{
   hier::BoxGeometry<DIM>* boxgeometry = new OuteredgeGeometry<DIM>(box, 0);
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
size_t OuteredgeDataFactory<DIM,TYPE>::getSizeOfMemory(
   const hier::Box<DIM>& box) const
{
   const size_t obj = tbox::Arena::align(sizeof(OuteredgeData<DIM,TYPE>));
   const size_t data = OuteredgeData<DIM,TYPE>::getSizeOfData(box,
                                                              d_depth);
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
bool OuteredgeDataFactory<DIM,TYPE>::validCopyTo(
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

