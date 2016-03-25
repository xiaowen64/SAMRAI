//
// File:	PatchEdgeDataMiscellaneousOpsReal.C
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Templated miscellaneous operations for real edge-centered data.
//

#ifndef included_math_PatchEdgeDataMiscellaneousOpsReal_C
#define included_math_PatchEdgeDataMiscellaneousOpsReal_C

#include "PatchEdgeDataMiscellaneousOpsReal.h"
#include "EdgeGeometry.h"
#include "tbox/Utilities.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

namespace SAMRAI {
    namespace math {

template<int DIM, class TYPE>
PatchEdgeDataMiscellaneousOpsReal<DIM,TYPE>::PatchEdgeDataMiscellaneousOpsReal()
{
}

template<int DIM, class TYPE>
PatchEdgeDataMiscellaneousOpsReal<DIM,TYPE>::~PatchEdgeDataMiscellaneousOpsReal()
{
}

/*
*************************************************************************
*                                                                       *
* The const constructor and assignment operator are not actually used   *
* but are defined here for compilers that require an implementation for *
* every declaration.                                                    *
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
PatchEdgeDataMiscellaneousOpsReal<DIM,TYPE>::PatchEdgeDataMiscellaneousOpsReal(
   const PatchEdgeDataMiscellaneousOpsReal<DIM,TYPE>& foo)
{
   (void) foo;  // not implemented (but needed by some compilers)
}

template<int DIM, class TYPE>
void PatchEdgeDataMiscellaneousOpsReal<DIM,TYPE>::operator=(
   const PatchEdgeDataMiscellaneousOpsReal<DIM,TYPE>& foo)
{
   (void) foo;  // not implemented (but needed by some compilers)
}

/*
*************************************************************************
*                                                                       *
* Templated miscellaneous opertions for real edge-centered data.        * 
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
int PatchEdgeDataMiscellaneousOpsReal<DIM,TYPE>::computeConstrProdPos(
   const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& data1,
   const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& data2,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::EdgeData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!data1.isNull() && !data2.isNull());
#endif
   int retval = 1;
   if (cvol.isNull()) {
      for (int d = 0; d < DIM; d++) {
         const hier::Box<DIM> edge_box = pdat::EdgeGeometry<DIM>::toEdgeBox(box, d);
         retval = tbox::Utilities::imin(retval,
                     d_array_ops.computeConstrProdPos(
                        data1->getArrayData(d),
                        data2->getArrayData(d),
                        edge_box) );
      }   
   } else {
      for (int d = 0; d < DIM; d++) {
         const hier::Box<DIM> edge_box = pdat::EdgeGeometry<DIM>::toEdgeBox(box, d);
         retval = tbox::Utilities::imin(retval,
                     d_array_ops.computeConstrProdPosWithControlVolume(
                        data1->getArrayData(d),
                        data2->getArrayData(d),
                        cvol->getArrayData(d),
                        edge_box) );
      }
   }
   return( retval );
}

template<int DIM, class TYPE>
void PatchEdgeDataMiscellaneousOpsReal<DIM,TYPE>::compareToScalar(
   tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& src,
   const TYPE& alpha,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::EdgeData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!dst.isNull() && !src.isNull());
#endif
   if (cvol.isNull()) {
      for (int d = 0; d < DIM; d++) {
         const hier::Box<DIM> edge_box = pdat::EdgeGeometry<DIM>::toEdgeBox(box, d);
         d_array_ops.compareToScalar(dst->getArrayData(d),
                                     src->getArrayData(d),
                                     alpha,
                                     edge_box);
      }
   } else {
      for (int d = 0; d < DIM; d++) {
         const hier::Box<DIM> edge_box = pdat::EdgeGeometry<DIM>::toEdgeBox(box, d);
         d_array_ops.compareToScalarWithControlVolume(dst->getArrayData(d),
                                                      src->getArrayData(d),
                                                      alpha,
                                                      cvol->getArrayData(d),
                                                      edge_box);
      }
   }
}

template<int DIM, class TYPE>
int PatchEdgeDataMiscellaneousOpsReal<DIM,TYPE>::testReciprocal(
   tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& src,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::EdgeData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!dst.isNull() && !src.isNull());
#endif
   int retval = 1;
   if (cvol.isNull()) {
      for (int d = 0; d < DIM; d++) {
         const hier::Box<DIM> edge_box = pdat::EdgeGeometry<DIM>::toEdgeBox(box, d);
         retval = tbox::Utilities::imin(retval,
                     d_array_ops.testReciprocal(
                        dst->getArrayData(d),
                        src->getArrayData(d),
                        edge_box) );
      }
   } else {
      for (int d = 0; d < DIM; d++) {
         const hier::Box<DIM> edge_box = pdat::EdgeGeometry<DIM>::toEdgeBox(box, d);
         retval = tbox::Utilities::imin(retval,
                     d_array_ops.testReciprocalWithControlVolume(
                        dst->getArrayData(d),
                        src->getArrayData(d),
                        cvol->getArrayData(d),
                        edge_box) );
      }
   }
   return( retval );
}

template<int DIM, class TYPE>
TYPE PatchEdgeDataMiscellaneousOpsReal<DIM,TYPE>::maxPointwiseDivide(
   const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& numer,
   const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& denom,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!numer.isNull() && !denom.isNull());
#endif
   TYPE retval = 0.0;
   for ( int d = 0; d < DIM; d++ ) {
      const hier::Box<DIM> edge_box = pdat::EdgeGeometry<DIM>::toEdgeBox(box, d);
      TYPE dirval = d_array_ops.maxPointwiseDivide(numer->getArrayData(d),
						   denom->getArrayData(d),
						   edge_box);
      retval = tbox::Utilities::dmax(retval, dirval);
   }
   return( retval );
}

template <int DIM, class TYPE>
TYPE PatchEdgeDataMiscellaneousOpsReal<DIM,TYPE>::minPointwiseDivide(
   const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& numer,
   const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& denom,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!numer.isNull() && !denom.isNull());
#endif
   TYPE retval = 0.0;
   for ( int d = 0; d < DIM; d++ ) {
      const hier::Box<DIM> edge_box = pdat::EdgeGeometry<DIM>::toEdgeBox(box, d);
      TYPE dirval = d_array_ops.minPointwiseDivide(numer->getArrayData(d),
						   denom->getArrayData(d),
						   edge_box);
      retval = tbox::Utilities::dmin(retval, dirval);
   }
   return( retval );
}

}
}
#endif
