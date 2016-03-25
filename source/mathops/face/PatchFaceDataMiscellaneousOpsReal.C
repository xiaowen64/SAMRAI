//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-0/source/mathops/face/PatchFaceDataMiscellaneousOpsReal.C $
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2249 $
// Modified:	$LastChangedDate: 2008-07-03 08:17:20 -0700 (Thu, 03 Jul 2008) $
// Description:	Templated miscellaneous operations for real face-centered data.
//

#ifndef included_math_PatchFaceDataMiscellaneousOpsReal_C
#define included_math_PatchFaceDataMiscellaneousOpsReal_C

#include "PatchFaceDataMiscellaneousOpsReal.h"
#include "tbox/MathUtilities.h"
#include "tbox/Utilities.h"
#include "FaceGeometry.h"

namespace SAMRAI {
    namespace math {

template<int DIM, class TYPE>
PatchFaceDataMiscellaneousOpsReal<DIM,TYPE>::PatchFaceDataMiscellaneousOpsReal()
{
}

template<int DIM, class TYPE>
PatchFaceDataMiscellaneousOpsReal<DIM,TYPE>::~PatchFaceDataMiscellaneousOpsReal()
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
PatchFaceDataMiscellaneousOpsReal<DIM,TYPE>::PatchFaceDataMiscellaneousOpsReal(
   const PatchFaceDataMiscellaneousOpsReal<DIM,TYPE>& foo)
{
   NULL_USE(foo);
}

template<int DIM, class TYPE>
void PatchFaceDataMiscellaneousOpsReal<DIM,TYPE>::operator=(
   const PatchFaceDataMiscellaneousOpsReal<DIM,TYPE>& foo)
{
   NULL_USE(foo);
}

/*
*************************************************************************
*                                                                       *
* Templated miscellaneous opertions for real face-centered data.        * 
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
int PatchFaceDataMiscellaneousOpsReal<DIM,TYPE>::computeConstrProdPos(
   const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& data1,
   const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& data2,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::FaceData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data1.isNull() && !data2.isNull());
#endif
   int retval = 1;
   if (cvol.isNull()) {
      for (int d = 0; d < DIM; d++) {
         const hier::Box<DIM> face_box = 
            pdat::FaceGeometry<DIM>::toFaceBox(box, d);
         retval = tbox::MathUtilities<int>::Min( retval,
                     d_array_ops.computeConstrProdPos(
                        data1->getArrayData(d),
                        data2->getArrayData(d),
                        face_box) );
      }   
   } else {
      for (int d = 0; d < DIM; d++) {
         const hier::Box<DIM> face_box = 
            pdat::FaceGeometry<DIM>::toFaceBox(box, d);
         retval = tbox::MathUtilities<int>::Min( retval,
                     d_array_ops.computeConstrProdPosWithControlVolume(
                        data1->getArrayData(d),
                        data2->getArrayData(d),
                        cvol->getArrayData(d),
                        face_box) );
      }
   }
   return( retval );
}

template<int DIM, class TYPE>
void PatchFaceDataMiscellaneousOpsReal<DIM,TYPE>::compareToScalar(
   tbox::Pointer< pdat::FaceData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& src,
   const TYPE& alpha,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::FaceData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
#endif
   if (cvol.isNull()) {
      for (int d = 0; d < DIM; d++) {
         const hier::Box<DIM> face_box = pdat::FaceGeometry<DIM>::toFaceBox(box, d);
         d_array_ops.compareToScalar(dst->getArrayData(d),
                                     src->getArrayData(d),
                                     alpha,
                                     face_box);
      }
   } else {
      for (int d = 0; d < DIM; d++) {
         const hier::Box<DIM> face_box = pdat::FaceGeometry<DIM>::toFaceBox(box, d);
         d_array_ops.compareToScalarWithControlVolume(dst->getArrayData(d),
                                                      src->getArrayData(d),
                                                      alpha,
                                                      cvol->getArrayData(d),
                                                      face_box);
      }
   }
}

template<int DIM, class TYPE>
int PatchFaceDataMiscellaneousOpsReal<DIM,TYPE>::testReciprocal(
   tbox::Pointer< pdat::FaceData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& src,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::FaceData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
#endif
   int retval = 1;
   if (cvol.isNull()) {
      for (int d = 0; d < DIM; d++) {
         const hier::Box<DIM> face_box = 
            pdat::FaceGeometry<DIM>::toFaceBox(box, d);
         retval = tbox::MathUtilities<int>::Min( retval,
                     d_array_ops.testReciprocal(
                        dst->getArrayData(d),
                        src->getArrayData(d),
                        face_box) );
      }
   } else {
      for (int d = 0; d < DIM; d++) {
         const hier::Box<DIM> face_box = 
            pdat::FaceGeometry<DIM>::toFaceBox(box, d);
         retval = tbox::MathUtilities<int>::Min( retval,
                     d_array_ops.testReciprocalWithControlVolume(
                        dst->getArrayData(d),
                        src->getArrayData(d),
                        cvol->getArrayData(d),
                        face_box) );
      }
   }
   return( retval );
}

template<int DIM, class TYPE>
TYPE PatchFaceDataMiscellaneousOpsReal<DIM,TYPE>::maxPointwiseDivide(
   const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& numer,
   const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& denom,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!numer.isNull() && !denom.isNull());
#endif
   TYPE retval = 0.0;
   for ( int d = 0; d < DIM; d++ ) {
      const hier::Box<DIM> face_box = 
         pdat::FaceGeometry<DIM>::toFaceBox(box, d);
      TYPE dirval = d_array_ops.maxPointwiseDivide(numer->getArrayData(d),
						   denom->getArrayData(d),
						   face_box);
      retval = tbox::MathUtilities<TYPE>::Max(retval, dirval);
   }
   return( retval );
}

template <int DIM, class TYPE>
TYPE PatchFaceDataMiscellaneousOpsReal<DIM,TYPE>::minPointwiseDivide(
   const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& numer,
   const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& denom,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!numer.isNull() && !denom.isNull());
#endif
   TYPE retval = 0.0;
   for ( int d = 0; d < DIM; d++ ) {
      const hier::Box<DIM> face_box = pdat::FaceGeometry<DIM>::toFaceBox(box, d);
      TYPE dirval = d_array_ops.minPointwiseDivide(numer->getArrayData(d),
						   denom->getArrayData(d),
						   face_box);
      retval = tbox::MathUtilities<TYPE>::Min(retval, dirval);
   }
   return( retval );
}

}
}
#endif
