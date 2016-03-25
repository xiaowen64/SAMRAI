//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/mathops/node/PatchNodeDataOpsInteger.C $
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Operations for integer node-centered patch data.
//

#ifndef included_math_PatchNodeDataOpsInteger_C
#define included_math_PatchNodeDataOpsInteger_C

#include "PatchNodeDataOpsInteger.h"
#include "NodeGeometry.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#include "tbox/Utilities.h"
#endif

namespace SAMRAI {
    namespace math {

template<int DIM>  PatchNodeDataOpsInteger<DIM>::PatchNodeDataOpsInteger()
{
}

template<int DIM>  PatchNodeDataOpsInteger<DIM>::~PatchNodeDataOpsInteger()
{
}

/*
*************************************************************************
*                                                                       *
* Compute the number of data entries on a patch in the given box.       *
*                                                                       *
*************************************************************************
*/

template<int DIM> int PatchNodeDataOpsInteger<DIM>::numberOfEntries(
   const tbox::Pointer< pdat::NodeData<DIM,int> >& data,
   const hier::Box<DIM>& box) const
{
   const hier::Box<DIM> ibox =
                   pdat::NodeGeometry<DIM>::toNodeBox(box * data->getGhostBox());
   int retval = ibox.size() * data->getDepth();
   return( retval );
}

/*
*************************************************************************
*                                                                       *
* General operations for integer node-centered patch data.              *
*                                                                       *
*************************************************************************
*/

template<int DIM> void PatchNodeDataOpsInteger<DIM>::swapData(
   tbox::Pointer< hier::Patch<DIM> > patch,
   const int data1_id,
   const int data2_id) const
{
   tbox::Pointer< pdat::NodeData<DIM,int> > d1 = patch->getPatchData(data1_id);
   tbox::Pointer< pdat::NodeData<DIM,int> > d2 = patch->getPatchData(data2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d1.isNull() && !d2.isNull());
   TBOX_ASSERT(d1->getDepth() && d2->getDepth());
   TBOX_ASSERT(d1->getBox() == d2->getBox());
   TBOX_ASSERT(d1->getGhostBox() == d2->getGhostBox());
#endif
   patch->setPatchData( data1_id, d2 );
   patch->setPatchData( data2_id, d1 );
}

template<int DIM> void PatchNodeDataOpsInteger<DIM>::printData(
   const tbox::Pointer< pdat::NodeData<DIM,int> >& data,
   const hier::Box<DIM>& box,
   std::ostream& s) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
#endif
   s << "Data box = " << box << std::endl;
   data->print(box, s);
   s << "\n";
}

template<int DIM> void PatchNodeDataOpsInteger<DIM>::copyData(
   tbox::Pointer< pdat::NodeData<DIM,int> >& dst,
   const tbox::Pointer< pdat::NodeData<DIM,int> >& src,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
#endif
   const hier::Box<DIM> node_box = pdat::NodeGeometry<DIM>::toNodeBox(box);
   (dst->getArrayData()).copy(src->getArrayData(), node_box);
}

template<int DIM> void PatchNodeDataOpsInteger<DIM>::setToScalar(
   tbox::Pointer< pdat::NodeData<DIM,int> >& dst,
   const int& alpha,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull());
#endif
   dst->fillAll(alpha, box);
}

template<int DIM> void PatchNodeDataOpsInteger<DIM>::abs(
   tbox::Pointer< pdat::NodeData<DIM,int> >& dst,
   const tbox::Pointer< pdat::NodeData<DIM,int> >& src,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
#endif
   const hier::Box<DIM> node_box = pdat::NodeGeometry<DIM>::toNodeBox(box);
   d_array_ops.abs(dst->getArrayData(),
                   src->getArrayData(),
                   node_box);
}

}
}
#endif
