//
// File:	PatchSideDataOpsInteger.C
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Operations for integer side-centered patch data.
//

#ifndef included_math_PatchSideDataOpsInteger_C
#define included_math_PatchSideDataOpsInteger_C

#include "PatchSideDataOpsInteger.h"
#include "SideGeometry.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

namespace SAMRAI {
    namespace math {

template<int DIM>  PatchSideDataOpsInteger<DIM>::PatchSideDataOpsInteger()
{
}

template<int DIM>  PatchSideDataOpsInteger<DIM>::~PatchSideDataOpsInteger()
{
}

/*
*************************************************************************
*                                                                       *
* Compute the number of data entries on a patch in the given box.       *
*                                                                       *
*************************************************************************
*/

template<int DIM> int PatchSideDataOpsInteger<DIM>::numberOfEntries(
   const tbox::Pointer< pdat::SideData<DIM,int> >& data,
   const hier::Box<DIM>& box) const
{
   int retval = 0;
   const hier::Box<DIM> ibox = box * data->getGhostBox();
   const hier::IntVector<DIM>& directions = data->getDirectionVector();
   for (int d = 0; d < DIM; d++) {
      if (directions(d)) {
         const hier::Box<DIM> dbox = pdat::SideGeometry<DIM>::toSideBox(ibox, d);
         retval += (dbox.size() * data->getDepth());
      }
   }
   return( retval );
}

/*
*************************************************************************
*                                                                       *
* General operations for integer side-centered patch data.              *
*                                                                       *
*************************************************************************
*/

template<int DIM> void PatchSideDataOpsInteger<DIM>::swapData(
   tbox::Pointer< hier::Patch<DIM> > patch,
   const int data1_id,
   const int data2_id) const
{
   tbox::Pointer< pdat::SideData<DIM,int> > d1 = patch->getPatchData(data1_id);
   tbox::Pointer< pdat::SideData<DIM,int> > d2 = patch->getPatchData(data2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d1.isNull() && !d2.isNull());
   assert(d1->getDepth() && d2->getDepth());
   assert(d1->getBox() == d2->getBox());
   assert(d1->getDirectionVector() == d2->getDirectionVector()); 
   assert(d1->getGhostBox() == d2->getGhostBox());
#endif
   patch->setPatchData( data1_id, d2 );
   patch->setPatchData( data2_id, d1 );
}

template<int DIM> void PatchSideDataOpsInteger<DIM>::printData(
   const tbox::Pointer< pdat::SideData<DIM,int> >& data,
   const hier::Box<DIM>& box,
   ostream& s) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!data.isNull());
#endif
   s << "Data box = " << box << endl;
   data->print(box, s);
   s << "\n";
}

template<int DIM> void PatchSideDataOpsInteger<DIM>::copyData(
   tbox::Pointer< pdat::SideData<DIM,int> >& dst,
   const tbox::Pointer< pdat::SideData<DIM,int> >& src,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!dst.isNull() && !src.isNull());
   assert(dst->getDirectionVector() == src->getDirectionVector());
#endif
   const hier::IntVector<DIM>& directions = dst->getDirectionVector();
   for (int d = 0; d < DIM; d++) {
      if (directions(d)) {
         const hier::Box<DIM> side_box = pdat::SideGeometry<DIM>::toSideBox(box, d);
         (dst->getArrayData(d)).copy(src->getArrayData(d), side_box);
      }
   }
}

template<int DIM> void PatchSideDataOpsInteger<DIM>::setToScalar(
   tbox::Pointer< pdat::SideData<DIM,int> >& dst,
   const int& alpha,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!dst.isNull());
#endif
   dst->fillAll(alpha, box);
}

template<int DIM> void PatchSideDataOpsInteger<DIM>::abs(
   tbox::Pointer< pdat::SideData<DIM,int> >& dst,
   const tbox::Pointer< pdat::SideData<DIM,int> >& src,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!dst.isNull() && !src.isNull());
   assert(dst->getDirectionVector() == src->getDirectionVector());
#endif
   const hier::IntVector<DIM>& directions = dst->getDirectionVector();
   for (int d = 0; d < DIM; d++) {
      if (directions(d)) {
         const hier::Box<DIM> side_box = pdat::SideGeometry<DIM>::toSideBox(box, d);
         d_array_ops.abs(dst->getArrayData(d),
                         src->getArrayData(d),
                         side_box);
      }
   }
}


}
}
#endif
