//
// File:	PatchSideDataNormOpsReal.C
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Templated norm operations for real side-centered patch data.
//

#ifndef included_math_PatchSideDataNormOpsReal_C
#define included_math_PatchSideDataNormOpsReal_C

#include "PatchSideDataNormOpsReal.h"
#include "SideGeometry.h"
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
PatchSideDataNormOpsReal<DIM,TYPE>::PatchSideDataNormOpsReal()
{
}

template<int DIM, class TYPE>
PatchSideDataNormOpsReal<DIM,TYPE>::~PatchSideDataNormOpsReal()
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
PatchSideDataNormOpsReal<DIM,TYPE>::PatchSideDataNormOpsReal(
   const PatchSideDataNormOpsReal<DIM,TYPE>& foo)
{
   (void) foo;  // not implemented (but needed by some compilers)
}

template<int DIM, class TYPE>
void PatchSideDataNormOpsReal<DIM,TYPE>::operator=(
   const PatchSideDataNormOpsReal<DIM,TYPE>& foo)
{
   (void) foo;  // not implemented (but needed by some compilers)
}

/*
*************************************************************************
*                                                                       *
* Compute the number of data entries on a patch in the given box.       *
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
int PatchSideDataNormOpsReal<DIM,TYPE>::numberOfEntries(
   const tbox::Pointer< pdat::SideData<DIM,TYPE> >& data,
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
* Templated norm operations for real side-centered data.                *
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
double PatchSideDataNormOpsReal<DIM,TYPE>::sumControlVolumes(
   const tbox::Pointer< pdat::SideData<DIM,TYPE> >& data,
   const tbox::Pointer< pdat::SideData<DIM,double> >& cvol,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!data.isNull() && !cvol.isNull());
#endif
   double retval = 0.0;
   const hier::IntVector<DIM>& directions = data->getDirectionVector();
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(directions ==
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

template<int DIM, class TYPE>
void PatchSideDataNormOpsReal<DIM,TYPE>::abs(
   tbox::Pointer< pdat::SideData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::SideData<DIM,TYPE> >& src,
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

template<int DIM, class TYPE>
double PatchSideDataNormOpsReal<DIM,TYPE>::L1Norm(
   const tbox::Pointer< pdat::SideData<DIM,TYPE> >& data,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::SideData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!data.isNull());
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
      assert(directions ==
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

template<int DIM, class TYPE>
double PatchSideDataNormOpsReal<DIM,TYPE>::L2Norm(
   const tbox::Pointer< pdat::SideData<DIM,TYPE> >& data,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::SideData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!data.isNull());
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
      assert(directions ==
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

template<int DIM, class TYPE>
double PatchSideDataNormOpsReal<DIM,TYPE>::weightedL2Norm(
   const tbox::Pointer< pdat::SideData<DIM,TYPE> >& data,
   const tbox::Pointer< pdat::SideData<DIM,TYPE> >& weight,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::SideData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!data.isNull() && !weight.isNull());
#endif
   double retval = 0.0;
   const hier::IntVector<DIM>& directions = data->getDirectionVector();
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(directions ==
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
      assert(directions ==
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

template<int DIM, class TYPE>
double PatchSideDataNormOpsReal<DIM,TYPE>::RMSNorm(
   const tbox::Pointer< pdat::SideData<DIM,TYPE> >& data,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::SideData<DIM,double> > cvol) const
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

template<int DIM, class TYPE>
double PatchSideDataNormOpsReal<DIM,TYPE>::weightedRMSNorm(
   const tbox::Pointer< pdat::SideData<DIM,TYPE> >& data,
   const tbox::Pointer< pdat::SideData<DIM,TYPE> >& weight,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::SideData<DIM,double> > cvol) const
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

template<int DIM, class TYPE>
double PatchSideDataNormOpsReal<DIM,TYPE>::maxNorm(
   const tbox::Pointer< pdat::SideData<DIM,TYPE> >& data,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::SideData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!data.isNull());
#endif
   double retval = 0.0;
   const hier::IntVector<DIM>& directions = data->getDirectionVector();
   if (cvol.isNull()) {
      for (int d = 0; d < DIM; d++) {
         if (directions(d)) {
            const hier::Box<DIM> side_box = pdat::SideGeometry<DIM>::toSideBox(box, d);
            retval = tbox::Utilities::dmax(retval,
                        d_array_ops.maxNorm(data->getArrayData(d), side_box) );
         }
      }
   } else {
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(directions ==
             hier::IntVector<DIM>::min(directions, cvol->getDirectionVector()));
#endif
      for (int d = 0; d < DIM; d++) {
         if (directions(d)) {
            const hier::Box<DIM> side_box = pdat::SideGeometry<DIM>::toSideBox(box, d);
            retval = tbox::Utilities::dmax(retval,
                        d_array_ops.maxNormWithControlVolume(
                           data->getArrayData(d), 
                           cvol->getArrayData(d),
                           side_box));
         }
      }
   }
   return( retval );
}

template<int DIM, class TYPE>
TYPE PatchSideDataNormOpsReal<DIM,TYPE>::dot(
   const tbox::Pointer< pdat::SideData<DIM,TYPE> >& data1,
   const tbox::Pointer< pdat::SideData<DIM,TYPE> >& data2,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::SideData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!data1.isNull() && !data2.isNull());
   assert(data1->getDirectionVector() == data2->getDirectionVector());
#endif
   TYPE retval = 0.0;
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
      assert(directions ==
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

template<int DIM, class TYPE>
TYPE PatchSideDataNormOpsReal<DIM,TYPE>::integral(
   const tbox::Pointer< pdat::SideData<DIM,TYPE> >& data,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::SideData<DIM,double> > vol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!data.isNull());
#endif
   TYPE retval = 0.0;
   const hier::IntVector<DIM>& directions = data->getDirectionVector();

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(directions ==
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
