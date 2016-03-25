//
// File:	PatchCellDataMiscellaneousOpsReal.C
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Templated miscellaneous operations for real cell-centered data.
//

#ifndef included_math_PatchCellDataMiscellaneousOpsReal_C
#define included_math_PatchCellDataMiscellaneousOpsReal_C

#include "PatchCellDataMiscellaneousOpsReal.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

namespace SAMRAI {
    namespace math {

template<int DIM, class TYPE>
PatchCellDataMiscellaneousOpsReal<DIM,TYPE>::PatchCellDataMiscellaneousOpsReal()
{
}

template<int DIM, class TYPE>
PatchCellDataMiscellaneousOpsReal<DIM,TYPE>::~PatchCellDataMiscellaneousOpsReal()
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
PatchCellDataMiscellaneousOpsReal<DIM,TYPE>::PatchCellDataMiscellaneousOpsReal(
   const PatchCellDataMiscellaneousOpsReal<DIM,TYPE>& foo)
{
   (void) foo;  // not implemented (but needed by some compilers)
}

template<int DIM, class TYPE>
void PatchCellDataMiscellaneousOpsReal<DIM,TYPE>::operator=(
   const PatchCellDataMiscellaneousOpsReal<DIM,TYPE>& foo)
{
   (void) foo;  // not implemented (but needed by some compilers)
}

/*
*************************************************************************
*                                                                       *
* Templated miscellaneous operations for real cell-centered data.       *
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
int PatchCellDataMiscellaneousOpsReal<DIM,TYPE>::computeConstrProdPos(
   const tbox::Pointer< pdat::CellData<DIM,TYPE> >& data1,
   const tbox::Pointer< pdat::CellData<DIM,TYPE> >& data2,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::CellData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!data1.isNull() && !data2.isNull());
#endif
   int retval;
   if (cvol.isNull()) {
      retval = d_array_ops.computeConstrProdPos(data1->getArrayData(),
                                                data2->getArrayData(),
                                                box);
   } else {
      retval = d_array_ops.computeConstrProdPosWithControlVolume(
                           data1->getArrayData(),
                           data2->getArrayData(),
                           cvol->getArrayData(),
                           box);
   }
   return( retval );
}

template<int DIM, class TYPE>
void PatchCellDataMiscellaneousOpsReal<DIM,TYPE>::compareToScalar(
   tbox::Pointer< pdat::CellData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::CellData<DIM,TYPE> >& src,
   const TYPE& alpha,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::CellData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!dst.isNull() && !src.isNull());
#endif
   if (cvol.isNull()) {
      d_array_ops.compareToScalar(dst->getArrayData(),
                                  src->getArrayData(),
                                  alpha,
                                  box);
   } else {
      d_array_ops.compareToScalarWithControlVolume(dst->getArrayData(),
                                                   src->getArrayData(),
                                                   alpha,
                                                   cvol->getArrayData(),
                                                   box);
   }
}

template<int DIM, class TYPE>
int PatchCellDataMiscellaneousOpsReal<DIM,TYPE>::testReciprocal(
   tbox::Pointer< pdat::CellData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::CellData<DIM,TYPE> >& src,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::CellData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!dst.isNull() && !src.isNull());
#endif
   int retval;
   if (cvol.isNull()) {
      retval = d_array_ops.testReciprocal(dst->getArrayData(),
                                          src->getArrayData(),
                                          box);
   } else {
      retval = d_array_ops.testReciprocalWithControlVolume(
                           dst->getArrayData(),
                           src->getArrayData(),
                           cvol->getArrayData(),
                           box);
   }
   return( retval );
}

template<int DIM, class TYPE>
TYPE PatchCellDataMiscellaneousOpsReal<DIM,TYPE>::maxPointwiseDivide(
   const tbox::Pointer< pdat::CellData<DIM,TYPE> >& numer,
   const tbox::Pointer< pdat::CellData<DIM,TYPE> >& denom,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!numer.isNull() && !denom.isNull());
#endif
   TYPE retval;
   retval = d_array_ops.maxPointwiseDivide(numer->getArrayData(),
					   denom->getArrayData(),
					   box);
   return( retval );
}


template <int DIM, class TYPE>
TYPE PatchCellDataMiscellaneousOpsReal<DIM,TYPE>::minPointwiseDivide(
   const tbox::Pointer< pdat::CellData<DIM,TYPE> >& numer,
   const tbox::Pointer< pdat::CellData<DIM,TYPE> >& denom,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!numer.isNull() && !denom.isNull());
#endif
   TYPE retval;
   retval = d_array_ops.minPointwiseDivide(numer->getArrayData(),
					   denom->getArrayData(),
					   box);
   return( retval );
}

}
}
#endif
