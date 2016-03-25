//
// File:	PatchCellDataOpsInteger.C
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Operations for integer cell-centered patch data.
//

#ifndef included_math_PatchCellDataOpsInteger_C
#define included_math_PatchCellDataOpsInteger_C

#include "PatchCellDataOpsInteger.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

namespace SAMRAI {
    namespace math {

template<int DIM>  PatchCellDataOpsInteger<DIM>::PatchCellDataOpsInteger()
{
}

template<int DIM>  PatchCellDataOpsInteger<DIM>::~PatchCellDataOpsInteger()
{
}

/*
*************************************************************************
*                                                                       *
* Compute the number of data entries on a patch in the given box.       *
*                                                                       *
*************************************************************************
*/

template<int DIM> int PatchCellDataOpsInteger<DIM>::numberOfEntries(
   const tbox::Pointer< pdat::CellData<DIM,int> >& data,
   const hier::Box<DIM>& box) const
{
   const hier::Box<DIM> ibox = box * data->getGhostBox();
   int retval = ibox.size() * data->getDepth();
   return( retval );
}

/*
*************************************************************************
*                                                                       *
* General operations for integer cell-centered patch data.              *
*                                                                       *
*************************************************************************
*/

template<int DIM> void PatchCellDataOpsInteger<DIM>::swapData(
   tbox::Pointer< hier::Patch<DIM> > patch,
   const int data1_id,
   const int data2_id) const
{
   tbox::Pointer< pdat::CellData<DIM,int> > d1 = patch->getPatchData(data1_id);
   tbox::Pointer< pdat::CellData<DIM,int> > d2 = patch->getPatchData(data2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d1.isNull() && !d2.isNull());
   assert(d1->getDepth() && d2->getDepth());
   assert(d1->getBox() == d2->getBox());
   assert(d1->getGhostBox() == d2->getGhostBox());
#endif
   patch->setPatchData( data1_id, d2 );
   patch->setPatchData( data2_id, d1 );
}

template<int DIM> void PatchCellDataOpsInteger<DIM>::printData(
   const tbox::Pointer< pdat::CellData<DIM,int> >& data,
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

template<int DIM> void PatchCellDataOpsInteger<DIM>::copyData(
   tbox::Pointer< pdat::CellData<DIM,int> >& dst,
   const tbox::Pointer< pdat::CellData<DIM,int> >& src,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!dst.isNull() && !src.isNull());
#endif
   (dst->getArrayData()).copy(src->getArrayData(), box);
}

template<int DIM> void PatchCellDataOpsInteger<DIM>::setToScalar(
   tbox::Pointer< pdat::CellData<DIM,int> >& dst,
   const int& alpha,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!dst.isNull());
#endif
   dst->fillAll(alpha, box);
}

template<int DIM> void PatchCellDataOpsInteger<DIM>::abs(
   tbox::Pointer< pdat::CellData<DIM,int> >& dst,
   const tbox::Pointer< pdat::CellData<DIM,int> >& src,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!dst.isNull() && !src.isNull());
#endif
   d_array_ops.abs(dst->getArrayData(),
                   src->getArrayData(),
                   box);
}

}
}
#endif
