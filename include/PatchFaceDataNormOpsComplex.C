//
// File:	PatchFaceDataNormOpsComplex.C
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Norm operations for complex face-centered patch data.
//

#ifndef included_math_PatchFaceDataNormOpsComplex_C
#define included_math_PatchFaceDataNormOpsComplex_C

#include "PatchFaceDataNormOpsComplex.h"
#include "FaceGeometry.h"
#include "tbox/Utilities.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

namespace SAMRAI {
    namespace math {

template<int DIM>  PatchFaceDataNormOpsComplex<DIM>::PatchFaceDataNormOpsComplex()
{
}

template<int DIM>  PatchFaceDataNormOpsComplex<DIM>::~PatchFaceDataNormOpsComplex()
{
}

/*
*************************************************************************
*                                                                       *
* Compute the number of data entries on a patch in the given box.       *
*                                                                       *
*************************************************************************
*/

template<int DIM> int PatchFaceDataNormOpsComplex<DIM>::numberOfEntries(
   const tbox::Pointer< pdat::FaceData<DIM,dcomplex> >& data,
   const hier::Box<DIM>& box) const
{
   int retval = 0;
   const hier::Box<DIM> ibox = box * data->getGhostBox();
   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> dbox = pdat::FaceGeometry<DIM>::toFaceBox(ibox, d);
      retval += (dbox.size() * data->getDepth()); 
   }
   return( retval );
}

/*
*************************************************************************
*                                                                       *
* Norm operations for complex face-centered data.                       *
*                                                                       *
*************************************************************************
*/

template<int DIM> double PatchFaceDataNormOpsComplex<DIM>::sumControlVolumes(
   const tbox::Pointer< pdat::FaceData<DIM,dcomplex> >& data,
   const tbox::Pointer< pdat::FaceData<DIM,double> >& cvol,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!data.isNull() && !cvol.isNull());
#endif
   double retval = 0.0;
   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> face_box = pdat::FaceGeometry<DIM>::toFaceBox(box, d);
      retval += d_array_ops.sumControlVolumes(data->getArrayData(d),
                                              cvol->getArrayData(d),
                                              face_box);
   }
   return( retval );
}

template<int DIM> void PatchFaceDataNormOpsComplex<DIM>::abs(
   tbox::Pointer< pdat::FaceData<DIM,double> >& dst,
   const tbox::Pointer< pdat::FaceData<DIM,dcomplex> >& src,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!dst.isNull() && !src.isNull());
#endif
   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> face_box = pdat::FaceGeometry<DIM>::toFaceBox(box, d);
      d_array_ops.abs(dst->getArrayData(d),
                      src->getArrayData(d),
                      face_box);
   }
}

template<int DIM> double PatchFaceDataNormOpsComplex<DIM>::L1Norm(
   const tbox::Pointer< pdat::FaceData<DIM,dcomplex> >& data,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::FaceData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!data.isNull());
#endif
   double retval = 0.0;
   if (cvol.isNull()) {
      for (int d = 0; d < DIM; d++) {
         const hier::Box<DIM> face_box = pdat::FaceGeometry<DIM>::toFaceBox(box, d);
         retval += d_array_ops.L1Norm(data->getArrayData(d), face_box);
      }
   } else {
      for (int d = 0; d < DIM; d++) {
         const hier::Box<DIM> face_box = pdat::FaceGeometry<DIM>::toFaceBox(box, d);
         retval += d_array_ops.L1NormWithControlVolume(data->getArrayData(d),
                                                       cvol->getArrayData(d),
                                                       face_box);
      }
   }
   return( retval );
}

template<int DIM> double PatchFaceDataNormOpsComplex<DIM>::L2Norm(
   const tbox::Pointer< pdat::FaceData<DIM,dcomplex> >& data,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::FaceData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!data.isNull());
#endif
   double retval = 0.0;
   if (cvol.isNull()) {
      for (int d = 0; d < DIM; d++) {
         const hier::Box<DIM> face_box = pdat::FaceGeometry<DIM>::toFaceBox(box, d);
         double aval = d_array_ops.L2Norm(data->getArrayData(d), face_box);
         retval += aval * aval;
      }
   } else {
      for (int d = 0; d < DIM; d++) {
         const hier::Box<DIM> face_box = pdat::FaceGeometry<DIM>::toFaceBox(box, d);
         double aval = d_array_ops.L2NormWithControlVolume(
                                   data->getArrayData(d), 
                                   cvol->getArrayData(d),
                                   face_box);
         retval += aval * aval;
      }
   }
   return( sqrt(retval) );
}

template<int DIM> double PatchFaceDataNormOpsComplex<DIM>::weightedL2Norm(
   const tbox::Pointer< pdat::FaceData<DIM,dcomplex> >& data,
   const tbox::Pointer< pdat::FaceData<DIM,dcomplex> >& weight,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::FaceData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!data.isNull() && !weight.isNull());
#endif
   double retval = 0.0;
   if (cvol.isNull()) {
      for (int d = 0; d < DIM; d++) {
         const hier::Box<DIM> face_box = pdat::FaceGeometry<DIM>::toFaceBox(box, d);
         double aval = d_array_ops.weightedL2Norm(data->getArrayData(d),
                                                  weight->getArrayData(d),
                                                  face_box);
         retval += aval * aval;
      } 
   } else {
      for (int d = 0; d < DIM; d++) {
         const hier::Box<DIM> face_box = pdat::FaceGeometry<DIM>::toFaceBox(box, d);
         double aval = d_array_ops.weightedL2NormWithControlVolume(
                                   data->getArrayData(d),
                                   weight->getArrayData(d),
                                   cvol->getArrayData(d),
                                   face_box);
         retval += aval * aval;
      }
   }
   return( sqrt(retval) );
}

template<int DIM> double PatchFaceDataNormOpsComplex<DIM>::RMSNorm(
   const tbox::Pointer< pdat::FaceData<DIM,dcomplex> >& data,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::FaceData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!data.isNull());
#endif
   double retval = L2Norm(data, box, cvol);
   if (cvol.isNull()) {
      retval /= sqrt( (double)numberOfEntries(data, box) );
   } else {
      retval /= sqrt( sumControlVolumes(data, cvol, box) );
   }
   return( retval );
}

template<int DIM> double PatchFaceDataNormOpsComplex<DIM>::weightedRMSNorm(
   const tbox::Pointer< pdat::FaceData<DIM,dcomplex> >& data,
   const tbox::Pointer< pdat::FaceData<DIM,dcomplex> >& weight,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::FaceData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!data.isNull() && !weight.isNull());
#endif
   double retval = weightedL2Norm(data, weight, box, cvol);
   if (cvol.isNull()) {
      retval /= sqrt( (double)numberOfEntries(data, box) );
   } else {
      retval /= sqrt( sumControlVolumes(data, cvol, box) );
   }
   return( retval );
}

template<int DIM> double PatchFaceDataNormOpsComplex<DIM>::maxNorm(
   const tbox::Pointer< pdat::FaceData<DIM,dcomplex> >& data,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::FaceData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!data.isNull());
#endif
   double retval = 0.0;
   if (cvol.isNull()) {
      for (int d = 0; d < DIM; d++) {
         const hier::Box<DIM> face_box = pdat::FaceGeometry<DIM>::toFaceBox(box, d);
         retval = tbox::Utilities::dmax(retval, 
                     d_array_ops.maxNorm(data->getArrayData(d), face_box) );
      }
   } else {
      for (int d = 0; d < DIM; d++) {
         const hier::Box<DIM> face_box = pdat::FaceGeometry<DIM>::toFaceBox(box, d);
         retval = tbox::Utilities::dmax(retval, 
                     d_array_ops.maxNormWithControlVolume(
                        data->getArrayData(d), cvol->getArrayData(d), face_box) );
      }
   }
   return( retval );
}

template<int DIM> dcomplex PatchFaceDataNormOpsComplex<DIM>::dot(
   const tbox::Pointer< pdat::FaceData<DIM,dcomplex> >& data1,
   const tbox::Pointer< pdat::FaceData<DIM,dcomplex> >& data2,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::FaceData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!data1.isNull() && !data2.isNull());
#endif
   dcomplex retval = dcomplex(0.0,0.0);
   if (cvol.isNull()) {
      for (int d = 0; d < DIM; d++) {
         const hier::Box<DIM> face_box = pdat::FaceGeometry<DIM>::toFaceBox(box, d);
         retval += d_array_ops.dot(data1->getArrayData(d),
                                   data2->getArrayData(d),
                                   face_box);
      }
   } else {
      for (int d = 0; d < DIM; d++) {
         const hier::Box<DIM> face_box = pdat::FaceGeometry<DIM>::toFaceBox(box, d);
         retval += d_array_ops.dotWithControlVolume(
                               data1->getArrayData(d),
                               data2->getArrayData(d),
                               cvol->getArrayData(d),
                               face_box);
      }
   }
   return( retval );
}

template<int DIM> dcomplex PatchFaceDataNormOpsComplex<DIM>::integral(
   const tbox::Pointer< pdat::FaceData<DIM,dcomplex> >& data,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::FaceData<DIM,double> > vol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!data.isNull());
#endif
   dcomplex retval = dcomplex(0.0,0.0);

   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> face_box = pdat::FaceGeometry<DIM>::toFaceBox(box, d);
      retval += d_array_ops.integral(data->getArrayData(d),
                                     vol->getArrayData(d),
                                     face_box);
   }

   return( retval );
}

}
}
#endif
