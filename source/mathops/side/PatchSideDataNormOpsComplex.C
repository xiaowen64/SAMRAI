//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/mathops/side/PatchSideDataNormOpsComplex.C $
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	Norm operations for complex side-centered patch data.
//

#ifndef included_math_PatchSideDataNormOpsComplex_C
#define included_math_PatchSideDataNormOpsComplex_C

#include "PatchSideDataNormOpsComplex.h"
#include "SideGeometry.h"
#include "tbox/MathUtilities.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#include "tbox/Utilities.h"
#endif

namespace SAMRAI {
    namespace math {

template<int DIM>  PatchSideDataNormOpsComplex<DIM>::PatchSideDataNormOpsComplex()
{
}

template<int DIM>  PatchSideDataNormOpsComplex<DIM>::~PatchSideDataNormOpsComplex()
{
}

/*
*************************************************************************
*                                                                       *
* Compute the number of data entries on a patch in the given box.       *
*                                                                       *
*************************************************************************
*/

template<int DIM> int PatchSideDataNormOpsComplex<DIM>::numberOfEntries(
   const tbox::Pointer< pdat::SideData<DIM,dcomplex> >& data,
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
* Norm operations for complex side-centered data.                       *
*                                                                       *
*************************************************************************
*/

template<int DIM> double PatchSideDataNormOpsComplex<DIM>::sumControlVolumes(
   const tbox::Pointer< pdat::SideData<DIM,dcomplex> >& data,
   const tbox::Pointer< pdat::SideData<DIM,double> >& cvol,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull() && !cvol.isNull());
#endif
   double retval = 0.0;
   const hier::IntVector<DIM>& directions = data->getDirectionVector();
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(directions ==
          hier::IntVector<DIM>::min(directions, cvol->getDirectionVector()));
#endif
   for (int d = 0; d < DIM; d++) {
      if (directions(d)) {
         const hier::Box<DIM> side_box = pdat::SideGeometry<DIM>::toSideBox(box, d);
         retval += d_array_ops.sumControlVolumes(data->getArrayData(d),
                                                 cvol->getArrayData(d),
                                                 side_box);
      }
   }
   return( retval );
}

template<int DIM> void PatchSideDataNormOpsComplex<DIM>::abs(
   tbox::Pointer< pdat::SideData<DIM,double> >& dst,
   const tbox::Pointer< pdat::SideData<DIM,dcomplex> >& src,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
   TBOX_ASSERT(dst->getDirectionVector() == src->getDirectionVector());
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

template<int DIM> double PatchSideDataNormOpsComplex<DIM>::L1Norm(
   const tbox::Pointer< pdat::SideData<DIM,dcomplex> >& data,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::SideData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
#endif
   double retval = 0.0;
   const hier::IntVector<DIM>& directions = data->getDirectionVector();
   if (cvol.isNull()) {
      for (int d = 0; d < DIM; d++) {
         if (directions(d)) { 
            const hier::Box<DIM> side_box = pdat::SideGeometry<DIM>::toSideBox(box, d);
            retval += d_array_ops.L1Norm(data->getArrayData(d), side_box);
         }
      }
   } else {
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(directions ==
             hier::IntVector<DIM>::min(directions, cvol->getDirectionVector()));
#endif
      for (int d = 0; d < DIM; d++) {
         if (directions(d)) {
            const hier::Box<DIM> side_box = pdat::SideGeometry<DIM>::toSideBox(box, d);
            retval += d_array_ops.L1NormWithControlVolume(data->getArrayData(d),
                                                          cvol->getArrayData(d),
                                                          side_box);
         }
      }
   }
   return( retval );
}

template<int DIM> double PatchSideDataNormOpsComplex<DIM>::L2Norm(
   const tbox::Pointer< pdat::SideData<DIM,dcomplex> >& data,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::SideData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
#endif
   double retval = 0.0;
   const hier::IntVector<DIM>& directions = data->getDirectionVector();
   if (cvol.isNull()) {
      for (int d = 0; d < DIM; d++) {
         if (directions(d)) {
            const hier::Box<DIM> side_box = pdat::SideGeometry<DIM>::toSideBox(box, d);
            double aval = d_array_ops.L2Norm(data->getArrayData(d), side_box);
            retval += aval * aval;
         }
      }
   } else {
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(directions ==
             hier::IntVector<DIM>::min(directions, cvol->getDirectionVector()));
#endif
      for (int d = 0; d < DIM; d++) {
         if (directions(d)) {
            const hier::Box<DIM> side_box = pdat::SideGeometry<DIM>::toSideBox(box, d);
            double aval = d_array_ops.L2NormWithControlVolume(
                                      data->getArrayData(d), 
                                      cvol->getArrayData(d),
                                      side_box);
            retval += aval * aval;
         }
      }
   }
   return( sqrt(retval) );
}

template<int DIM> double PatchSideDataNormOpsComplex<DIM>::weightedL2Norm(
   const tbox::Pointer< pdat::SideData<DIM,dcomplex> >& data,
   const tbox::Pointer< pdat::SideData<DIM,dcomplex> >& weight,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::SideData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull() && !weight.isNull());
#endif
   double retval = 0.0;
   const hier::IntVector<DIM>& directions = data->getDirectionVector();
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(directions ==
          hier::IntVector<DIM>::min(directions, weight->getDirectionVector()));
#endif
   if (cvol.isNull()) {
      for (int d = 0; d < DIM; d++) {
         if (directions(d)) {
            const hier::Box<DIM> side_box = pdat::SideGeometry<DIM>::toSideBox(box, d);
            double aval = d_array_ops.weightedL2Norm(data->getArrayData(d),
                                                     weight->getArrayData(d),
                                                     side_box);
            retval += aval * aval;
         }
      } 
   } else {
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(directions ==
             hier::IntVector<DIM>::min(directions, cvol->getDirectionVector()));
#endif
      for (int d = 0; d < DIM; d++) {
         if (directions(d)) {
            const hier::Box<DIM> side_box = pdat::SideGeometry<DIM>::toSideBox(box, d);
            double aval = d_array_ops.weightedL2NormWithControlVolume(
                                      data->getArrayData(d),
                                      weight->getArrayData(d),
                                      cvol->getArrayData(d),
                                      side_box);
            retval += aval * aval;
         }
      }
   }
   return( sqrt(retval) );
}

template<int DIM> double PatchSideDataNormOpsComplex<DIM>::RMSNorm(
   const tbox::Pointer< pdat::SideData<DIM,dcomplex> >& data,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::SideData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
#endif
   double retval = L2Norm(data, box, cvol);
   if (cvol.isNull()) {
      retval /= sqrt( (double)numberOfEntries(data, box) );
   } else {
      retval /= sqrt( sumControlVolumes(data, cvol, box) );
   }
   return( retval );
}

template<int DIM> double PatchSideDataNormOpsComplex<DIM>::weightedRMSNorm(
   const tbox::Pointer< pdat::SideData<DIM,dcomplex> >& data,
   const tbox::Pointer< pdat::SideData<DIM,dcomplex> >& weight,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::SideData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull() && !weight.isNull());
#endif
   double retval = weightedL2Norm(data, weight, box, cvol);
   if (cvol.isNull()) {
      retval /= sqrt( (double)numberOfEntries(data, box) );
   } else {
      retval /= sqrt( sumControlVolumes(data, cvol, box) );
   }
   return( retval );
}

template<int DIM> double PatchSideDataNormOpsComplex<DIM>::maxNorm(
   const tbox::Pointer< pdat::SideData<DIM,dcomplex> >& data,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::SideData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
#endif
   double retval = 0.0;
   const hier::IntVector<DIM>& directions = data->getDirectionVector();
   if (cvol.isNull()) {
      for (int d = 0; d < DIM; d++) {
         if (directions(d)) {
            const hier::Box<DIM> side_box = 
               pdat::SideGeometry<DIM>::toSideBox(box, d);
            retval = tbox::MathUtilities<double>::Max(retval, 
                        d_array_ops.maxNorm(data->getArrayData(d), side_box) );
         }
      }
   } else {
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(directions ==
             hier::IntVector<DIM>::min(directions, cvol->getDirectionVector()));
#endif
      for (int d = 0; d < DIM; d++) {
         if (directions(d)) {
            const hier::Box<DIM> side_box = 
               pdat::SideGeometry<DIM>::toSideBox(box, d);
            retval = tbox::MathUtilities<double>::Max(retval, 
                     d_array_ops.maxNormWithControlVolume(
                                 data->getArrayData(d), 
                                 cvol->getArrayData(d), side_box) );
         }
      }
   }
   return( retval );
}

template<int DIM> dcomplex PatchSideDataNormOpsComplex<DIM>::dot(
   const tbox::Pointer< pdat::SideData<DIM,dcomplex> >& data1,
   const tbox::Pointer< pdat::SideData<DIM,dcomplex> >& data2,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::SideData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data1.isNull() && !data2.isNull());
   TBOX_ASSERT(data1->getDirectionVector() == data2->getDirectionVector()); 
#endif
   dcomplex retval = dcomplex(0.0,0.0);
   const hier::IntVector<DIM>& directions = data1->getDirectionVector();
   if (cvol.isNull()) {
      for (int d = 0; d < DIM; d++) {
         if (directions(d)) {
            const hier::Box<DIM> side_box = pdat::SideGeometry<DIM>::toSideBox(box, d);
            retval += d_array_ops.dot(data1->getArrayData(d),
                                      data2->getArrayData(d),
                                      side_box);
         }
      }
   } else {
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(directions ==
             hier::IntVector<DIM>::min(directions, cvol->getDirectionVector()));
#endif
      for (int d = 0; d < DIM; d++) {
         if (directions(d)) {
            const hier::Box<DIM> side_box = pdat::SideGeometry<DIM>::toSideBox(box, d);
            retval += d_array_ops.dotWithControlVolume(
                                  data1->getArrayData(d),
                                  data2->getArrayData(d),
                                  cvol->getArrayData(d),
                                  side_box);
         }
      }
   }
   return( retval );
}

template<int DIM> dcomplex PatchSideDataNormOpsComplex<DIM>::integral(
   const tbox::Pointer< pdat::SideData<DIM,dcomplex> >& data,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::SideData<DIM,double> > vol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
#endif
   dcomplex retval = dcomplex(0.0,0.0);
   const hier::IntVector<DIM>& directions = data->getDirectionVector();

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(directions ==
          hier::IntVector<DIM>::min(directions, vol->getDirectionVector()));
#endif
   for (int d = 0; d < DIM; d++) {
      if (directions(d)) {
         const hier::Box<DIM> side_box = pdat::SideGeometry<DIM>::toSideBox(box, d);
         retval += d_array_ops.integral(
                               data->getArrayData(d),
                               vol->getArrayData(d),
                               side_box);
      }
   }

   return( retval );
}

}
}
#endif
