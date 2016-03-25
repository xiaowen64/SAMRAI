//
// File:	PatchEdgeDataOpsComplex.C
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Operations for complex edge-centered patch data.
//

#ifndef included_math_PatchEdgeDataOpsComplex_C
#define included_math_PatchEdgeDataOpsComplex_C

#include "PatchEdgeDataOpsComplex.h"
#include "EdgeGeometry.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

namespace SAMRAI {
    namespace math {

template<int DIM>  PatchEdgeDataOpsComplex<DIM>::PatchEdgeDataOpsComplex()
{
}

template<int DIM>  PatchEdgeDataOpsComplex<DIM>::~PatchEdgeDataOpsComplex()
{
}

/*
*************************************************************************
*                                                                       *
* General operations for complex edge-centered patch data.              *
*                                                                       *
*************************************************************************
*/

template<int DIM> void PatchEdgeDataOpsComplex<DIM>::swapData(
   tbox::Pointer< hier::Patch<DIM> > patch,
   const int data1_id,
   const int data2_id) const
{
   tbox::Pointer< pdat::EdgeData<DIM,dcomplex> > d1 = patch->getPatchData(data1_id);
   tbox::Pointer< pdat::EdgeData<DIM,dcomplex> > d2 = patch->getPatchData(data2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d1.isNull() && !d2.isNull());
   assert(d1->getDepth() && d2->getDepth());
   assert(d1->getBox() == d2->getBox());
   assert(d1->getGhostBox() == d2->getGhostBox());
#endif
   patch->setPatchData( data1_id, d2 );
   patch->setPatchData( data2_id, d1 );
}

template<int DIM> void PatchEdgeDataOpsComplex<DIM>::printData(
   const tbox::Pointer< pdat::EdgeData<DIM,dcomplex> >& data,
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

template<int DIM> void PatchEdgeDataOpsComplex<DIM>::copyData(
   tbox::Pointer< pdat::EdgeData<DIM,dcomplex> >& dst,
   const tbox::Pointer< pdat::EdgeData<DIM,dcomplex> >& src,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!dst.isNull() && !src.isNull());
#endif
   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> edge_box = pdat::EdgeGeometry<DIM>::toEdgeBox(box, d);
      (dst->getArrayData(d)).copy(src->getArrayData(d), edge_box);
   }
}

template<int DIM> void PatchEdgeDataOpsComplex<DIM>::setToScalar(
   tbox::Pointer< pdat::EdgeData<DIM,dcomplex> >& dst,
   const dcomplex& alpha,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!dst.isNull());
#endif
   dst->fillAll(alpha, box);
}

}
}
#endif
