//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-0/source/mathops/array/ArrayDataMiscellaneousOpsReal.C $
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2249 $
// Modified:	$LastChangedDate: 2008-07-03 08:17:20 -0700 (Thu, 03 Jul 2008) $
// Description:	Miscellaneous templated operations for real array data
//

#ifndef included_math_ArrayDataMiscellaneousOpsReal_C
#define included_math_ArrayDataMiscellaneousOpsReal_C

#include "ArrayDataMiscellaneousOpsReal.h"

#include "tbox/MathUtilities.h"
#include "tbox/Utilities.h"

namespace SAMRAI {
    namespace math {

template<int DIM, class TYPE>
ArrayDataMiscellaneousOpsReal<DIM,TYPE>::ArrayDataMiscellaneousOpsReal()
{
}

template<int DIM, class TYPE>
ArrayDataMiscellaneousOpsReal<DIM,TYPE>::~ArrayDataMiscellaneousOpsReal()
{
}

/*
*************************************************************************
*									*
* The const constructor and assignment operator are not actually used	*
* but are defined here for compilers that require an implementation for	*
* every declaration.							*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
ArrayDataMiscellaneousOpsReal<DIM,TYPE>::ArrayDataMiscellaneousOpsReal(
   const ArrayDataMiscellaneousOpsReal<DIM,TYPE>& foo)
{
   NULL_USE(foo);
}

template<int DIM, class TYPE>
void ArrayDataMiscellaneousOpsReal<DIM,TYPE>::operator=(
   const ArrayDataMiscellaneousOpsReal<DIM,TYPE>& foo)
{
   NULL_USE(foo);
}

/*
*************************************************************************
*									*
* General templated miscellaneous operations for array data.            *
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
int 
ArrayDataMiscellaneousOpsReal<DIM,TYPE>::computeConstrProdPosWithControlVolume(
   const pdat::ArrayData<DIM,TYPE>& data1,
   const pdat::ArrayData<DIM,TYPE>& data2,
   const pdat::ArrayData<DIM,double>& cvol,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(data1.getDepth() == data2.getDepth());
#endif

   int test = 1;

   const hier::Box<DIM> d1_box = data1.getBox();
   const hier::Box<DIM> d2_box = data2.getBox();
   const hier::Box<DIM> cv_box = cvol.getBox();
   const hier::Box<DIM> ibox = box * d1_box * d2_box * cv_box;

   if (!ibox.empty()) {
      const int ddepth  = data1.getDepth();
      const int cvdepth = cvol.getDepth();
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT((ddepth == cvdepth) || (cvdepth == 1));
#endif

      int box_w[DIM];
      int d1_w[DIM];
      int d2_w[DIM];
      int cv_w[DIM];
      int dim_counter[DIM];
      for (int i = 0; i < DIM; i++) {
         box_w[i] = ibox.numberCells(i);
         d1_w[i] = d1_box.numberCells(i);
         d2_w[i] = d2_box.numberCells(i);
         cv_w[i] = cv_box.numberCells(i);
         dim_counter[i] = 0;
      } 

      const int d1_offset = data1.getOffset();
      const int d2_offset = data2.getOffset();
      const int cv_offset = ( (cvdepth == 1) ? 0 : cvol.getOffset() );

      const int num_d0_blocks = ibox.size() / box_w[0];

      int d1_begin = d1_box.offset(ibox.lower());
      int d2_begin = d2_box.offset(ibox.lower());
      int cv_begin = cv_box.offset(ibox.lower());

      const TYPE* dd1   = data1.getPointer();
      const TYPE* dd2   = data2.getPointer();
      const double* cvd = cvol.getPointer();

      for (int d = 0; d < ddepth; d++) {

         int d1_counter = d1_begin;
         int d2_counter = d2_begin;
         int cv_counter = cv_begin;

         int d1_b[DIM];
         int d2_b[DIM];
         int cv_b[DIM];
         for (int nd = 0; nd < DIM; nd++) {
            d1_b[nd] = d1_counter;
            d2_b[nd] = d2_counter;
            cv_b[nd] = cv_counter;
         }

         for (int nb = 0; nb < num_d0_blocks; nb++) {

            for (int i0 = 0; i0 < box_w[0]; i0++) {
               if (cvd[cv_counter+i0] > 0.0) {
                  if ( tbox::MathUtilities<TYPE>::Abs(dd2[d2_counter+i0]) > 0.0
                       && (dd1[d1_counter+i0] * dd2[d2_counter+i0] <= 0.0) 
                     ) {
                     test = 0;
                  }
               }
            }

            int dim_jump = 0;

            for (int j = 1; j < DIM; j++) {
               if (dim_counter[j] < box_w[j]-1) {
                  ++dim_counter[j];
                  dim_jump = j;
                  break;
               } else {
                  dim_counter[j] = 0;
               }
            }

            if (dim_jump > 0) {
               int d1_step = 1;
               int d2_step = 1;
               int cv_step = 1;
               for (int k = 0; k < dim_jump; k++) {
                  d1_step *= d1_w[k];
                  d2_step *= d2_w[k];
                  cv_step *= cv_w[k];
               }
               d1_counter = d1_b[dim_jump-1] + d1_step;
               d2_counter = d2_b[dim_jump-1] + d1_step;
               cv_counter = cv_b[dim_jump-1] + cv_step;

               for (int m = 0; m < dim_jump; m++) {
                  d1_b[m] = d1_counter;
                  d2_b[m] = d2_counter;
                  cv_b[m] = cv_counter;
               }
            }

         }

         d1_begin += d1_offset;
         d2_begin += d2_offset;
         cv_begin += cv_offset;
      }

   }

   return( test );
}

template<int DIM, class TYPE>
int 
ArrayDataMiscellaneousOpsReal<DIM,TYPE>::computeConstrProdPos(
   const pdat::ArrayData<DIM,TYPE>& data1,
   const pdat::ArrayData<DIM,TYPE>& data2,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(data1.getDepth() == data2.getDepth());
#endif

   int test = 1;

   const hier::Box<DIM> d1_box = data1.getBox();
   const hier::Box<DIM> d2_box = data2.getBox();
   const hier::Box<DIM> ibox = box * d1_box * d2_box;

   if (!ibox.empty()) {
      const int ddepth  = data1.getDepth();

      int box_w[DIM];
      int d1_w[DIM];
      int d2_w[DIM];
      int dim_counter[DIM];
      for (int i = 0; i < DIM; i++) {
         box_w[i] = ibox.numberCells(i);
         d1_w[i] = d1_box.numberCells(i);
         d2_w[i] = d2_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const int d1_offset = data1.getOffset();
      const int d2_offset = data2.getOffset();

      const int num_d0_blocks = ibox.size() / box_w[0];

      int d1_begin = d1_box.offset(ibox.lower());
      int d2_begin = d2_box.offset(ibox.lower());

      const TYPE* dd1   = data1.getPointer();
      const TYPE* dd2   = data2.getPointer();

      for (int d = 0; d < ddepth; d++) {

         int d1_counter = d1_begin;
         int d2_counter = d2_begin;

         int d1_b[DIM];
         int d2_b[DIM];
         for (int nd = 0; nd < DIM; nd++) {
            d1_b[nd] = d1_counter;
            d2_b[nd] = d2_counter;
         }

         for (int nb = 0; nb < num_d0_blocks; nb++) {

            for (int i0 = 0; i0 < box_w[0]; i0++) {
               if ( tbox::MathUtilities<TYPE>::Abs(dd2[d2_counter+i0]) > 0.0
                    && (dd1[d1_counter+i0] * dd2[d2_counter+i0] <= 0.0) ) {
                  test = 0;
               }
            }

            int dim_jump = 0;

            for (int j = 1; j < DIM; j++) {
               if (dim_counter[j] < box_w[j]-1) {
                  ++dim_counter[j];
                  dim_jump = j;
                  break;
               } else {
                  dim_counter[j] = 0;
               }
            }

            if (dim_jump > 0) {
               int d1_step = 1;
               int d2_step = 1;
               for (int k = 0; k < dim_jump; k++) {
                  d1_step *= d1_w[k];
                  d2_step *= d2_w[k];
               }
               d1_counter = d1_b[dim_jump-1] + d1_step;
               d2_counter = d2_b[dim_jump-1] + d1_step;

               for (int m = 0; m < dim_jump; m++) {
                  d1_b[m] = d1_counter;
                  d2_b[m] = d2_counter;
               }
            }

         }

         d1_begin += d1_offset;
         d2_begin += d2_offset;
      }

   }

   return( test );
}

template<int DIM, class TYPE>
void 
ArrayDataMiscellaneousOpsReal<DIM,TYPE>::compareToScalarWithControlVolume(
   pdat::ArrayData<DIM,TYPE>& dst,
   const pdat::ArrayData<DIM,TYPE>& src,
   const TYPE& alpha,
   const pdat::ArrayData<DIM,double>& cvol,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(dst.getDepth() == src.getDepth());
#endif

   const hier::Box<DIM> d_box = dst.getBox();
   const hier::Box<DIM> s_box = src.getBox();
   const hier::Box<DIM> cv_box = cvol.getBox();
   const hier::Box<DIM> ibox = box * d_box * s_box * cv_box;

   if (!ibox.empty()) {
      const int ddepth  = dst.getDepth();
      const int cvdepth = cvol.getDepth();
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT((ddepth == cvdepth) || (cvdepth == 1));
#endif

      int box_w[DIM];
      int d_w[DIM];
      int s_w[DIM];
      int cv_w[DIM];
      int dim_counter[DIM];
      for (int i = 0; i < DIM; i++) {
         box_w[i] = ibox.numberCells(i);
         d_w[i] = d_box.numberCells(i);
         s_w[i] = s_box.numberCells(i);
         cv_w[i] = cv_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const int d_offset = dst.getOffset();
      const int s_offset = src.getOffset();
      const int cv_offset = ( (cvdepth == 1) ? 0 : cvol.getOffset() );

      const int num_d0_blocks = ibox.size() / box_w[0];

      int d_begin = d_box.offset(ibox.lower());
      int s_begin = s_box.offset(ibox.lower());
      int cv_begin = cv_box.offset(ibox.lower());

      TYPE* dd    = dst.getPointer();
      const TYPE* sd    = src.getPointer();
      const double* cvd = cvol.getPointer();

      for (int d = 0; d < ddepth; d++) {

         int d_counter = d_begin;
         int s_counter = s_begin;
         int cv_counter = cv_begin;

         int d_b[DIM];
         int s_b[DIM];
         int cv_b[DIM];
         for (int nd = 0; nd < DIM; nd++) {
            d_b[nd] = d_counter;
            s_b[nd] = s_counter;
            cv_b[nd] = cv_counter;
         }

         for (int nb = 0; nb < num_d0_blocks; nb++) {

            for (int i0 = 0; i0 < box_w[0]; i0++) {

               if (cvd[cv_counter+i0] > 0.0) {
                  dd[d_counter+i0] = ( 
                  (tbox::MathUtilities<TYPE>::Abs(sd[s_counter+i0]) >= alpha) 
                      ?  1.0F : 0.0F );
               }
            }

            int dim_jump = 0;

            for (int j = 1; j < DIM; j++) {
               if (dim_counter[j] < box_w[j]-1) {
                  ++dim_counter[j];
                  dim_jump = j;
                  break;
               } else {
                  dim_counter[j] = 0;
               }
            }

            if (dim_jump > 0) {
               int d_step = 1;
               int s_step = 1;
               int cv_step = 1;
               for (int k = 0; k < dim_jump; k++) {
                  d_step *= d_w[k];
                  s_step *= s_w[k];
                  cv_step *= cv_w[k];
               }
               d_counter = d_b[dim_jump-1] + d_step;
               s_counter = s_b[dim_jump-1] + s_step;
               cv_counter = cv_b[dim_jump-1] + cv_step;

               for (int m = 0; m < dim_jump; m++) {
                  d_b[m] = d_counter;
                  s_b[m] = s_counter;
                  cv_b[m] = cv_counter;
               }
            }
         }

         d_begin += d_offset;
         s_begin += s_offset;
         cv_begin += cv_offset;
      }

   }
}

template<int DIM, class TYPE>
void 
ArrayDataMiscellaneousOpsReal<DIM,TYPE>::compareToScalar(
   pdat::ArrayData<DIM,TYPE>& dst,
   const pdat::ArrayData<DIM,TYPE>& src,
   const TYPE& alpha,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(dst.getDepth() == src.getDepth());
#endif

   const hier::Box<DIM> d_box = dst.getBox();
   const hier::Box<DIM> s_box = src.getBox();
   const hier::Box<DIM> ibox = box * d_box * s_box;

   if (!ibox.empty()) {

      const int ddepth  = dst.getDepth();

      int box_w[DIM];
      int d_w[DIM];
      int s_w[DIM];
      int dim_counter[DIM];
      for (int i = 0; i < DIM; i++) {
         box_w[i] = ibox.numberCells(i);
         d_w[i] = d_box.numberCells(i);
         s_w[i] = s_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const int d_offset = dst.getOffset();
      const int s_offset = src.getOffset();

      const int num_d0_blocks = ibox.size() / box_w[0];

      int d_begin = d_box.offset(ibox.lower());
      int s_begin = s_box.offset(ibox.lower());

      TYPE* dd    = dst.getPointer();
      const TYPE* sd    = src.getPointer();

      for (int d = 0; d < ddepth; d++) {

         int d_counter = d_begin;
         int s_counter = s_begin;

         int d_b[DIM];
         int s_b[DIM];
         for (int nd = 0; nd < DIM; nd++) {
            d_b[nd] = d_counter;
            s_b[nd] = s_counter;
         }

         for (int nb = 0; nb < num_d0_blocks; nb++) {

            for (int i0 = 0; i0 < box_w[0]; i0++) {
               dd[d_counter+i0] = ( 
                  (tbox::MathUtilities<TYPE>::Abs(sd[s_counter+i0]) >= alpha) 
                      ?  1.0F : 0.0F );
            }

            int dim_jump = 0;

            for (int j = 1; j < DIM; j++) {
               if (dim_counter[j] < box_w[j]-1) {
                  ++dim_counter[j];
                  dim_jump = j;
                  break;
               } else {
                  dim_counter[j] = 0;
               }
            }

            if (dim_jump > 0) {
               int d_step = 1;
               int s_step = 1;
               for (int k = 0; k < dim_jump; k++) {
                  d_step *= d_w[k];
                  s_step *= s_w[k];
               }
               d_counter = d_b[dim_jump-1] + d_step;
               s_counter = s_b[dim_jump-1] + s_step;

               for (int m = 0; m < dim_jump; m++) {
                  d_b[m] = d_counter;
                  s_b[m] = s_counter;
               }
            }
         }

         d_begin += d_offset;
         s_begin += s_offset;
      }

   }
}

template<int DIM, class TYPE>
int 
ArrayDataMiscellaneousOpsReal<DIM,TYPE>::testReciprocalWithControlVolume(
   pdat::ArrayData<DIM,TYPE>& dst,
   const pdat::ArrayData<DIM,TYPE>& src,
   const pdat::ArrayData<DIM,double>& cvol,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(dst.getDepth() == src.getDepth());
#endif

// Ignore Intel warning about floating point comparisons
#ifdef __INTEL_COMPILER
#pragma warning (disable:1572)
#endif

   int test = 1;

   const hier::Box<DIM> d_box = dst.getBox();
   const hier::Box<DIM> s_box = src.getBox();
   const hier::Box<DIM> cv_box = cvol.getBox();
   const hier::Box<DIM> ibox = box * d_box * s_box * cv_box;

   if (!ibox.empty()) {
      const int ddepth  = dst.getDepth();
      const int cvdepth = cvol.getDepth();
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT((ddepth == cvdepth) || (cvdepth == 1));
#endif

      int box_w[DIM];
      int d_w[DIM];
      int s_w[DIM];
      int cv_w[DIM];
      int dim_counter[DIM];
      for (int i = 0; i < DIM; i++) {
         box_w[i] = ibox.numberCells(i);
         d_w[i] = d_box.numberCells(i);
         s_w[i] = s_box.numberCells(i);
         cv_w[i] = cv_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const int d_offset = dst.getOffset();
      const int s_offset = src.getOffset();
      const int cv_offset = ( (cvdepth == 1) ? 0 : cvol.getOffset() );

      const int num_d0_blocks = ibox.size() / box_w[0];

      int d_begin = d_box.offset(ibox.lower());
      int s_begin = s_box.offset(ibox.lower());
      int cv_begin = cv_box.offset(ibox.lower());

      TYPE* dd    = dst.getPointer();
      const TYPE* sd    = src.getPointer();
      const double* cvd = cvol.getPointer();

      for (int d = 0; d < ddepth; d++) {

         int d_counter = d_begin;
         int s_counter = s_begin;
         int cv_counter = cv_begin;

         int d_b[DIM];
         int s_b[DIM];
         int cv_b[DIM];
         for (int nd = 0; nd < DIM; nd++) {
            d_b[nd] = d_counter;
            s_b[nd] = s_counter;
            cv_b[nd] = cv_counter;
         }

         for (int nb = 0; nb < num_d0_blocks; nb++) {

            for (int i0 = 0; i0 < box_w[0]; i0++) {
               if (cvd[cv_counter+i0] > 0.0) {
                  if ( sd[s_counter+i0] == 0.0 ) {
                     test = 0;
                     dd[d_counter+i0] = 0.0;
                  } else {
                     dd[d_counter+i0] = 1.0F / sd[s_counter+i0];
                  }
               }
            }

            int dim_jump = 0;

            for (int j = 1; j < DIM; j++) {
               if (dim_counter[j] < box_w[j]-1) {
                  ++dim_counter[j];
                  dim_jump = j;
                  break;
               } else {
                  dim_counter[j] = 0;
               }
            }
            if (dim_jump > 0) {
               int d_step = 1;
               int s_step = 1;
               int cv_step = 1;
               for (int k = 0; k < dim_jump; k++) {
                  d_step *= d_w[k];
                  s_step *= s_w[k];
                  cv_step *= cv_w[k];
               }
               d_counter = d_b[dim_jump-1] + d_step;
               s_counter = s_b[dim_jump-1] + s_step;
               cv_counter = cv_b[dim_jump-1] + cv_step;

               for (int m = 0; m < dim_jump; m++) {
                  d_b[m] = d_counter;
                  s_b[m] = s_counter;
                  cv_b[m] = cv_counter;
               }
            }
         }

         d_begin += d_offset;
         s_begin += s_offset;
         cv_begin += cv_offset;
      }

   }

   return( test );
}

template<int DIM, class TYPE>
int 
ArrayDataMiscellaneousOpsReal<DIM,TYPE>::testReciprocal(
   pdat::ArrayData<DIM,TYPE>& dst,
   const pdat::ArrayData<DIM,TYPE>& src,
   const hier::Box<DIM>& box) const
{
// Ignore Intel warning about floating point comparisons
#ifdef __INTEL_COMPILER
#pragma warning (disable:1572)
#endif

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(dst.getDepth() == src.getDepth());
#endif

   int test = 1;

   const hier::Box<DIM> d_box = dst.getBox();
   const hier::Box<DIM> s_box = src.getBox();
   const hier::Box<DIM> ibox = box * d_box * s_box;

   if (!ibox.empty()) {
      const int ddepth  = dst.getDepth();

      int box_w[DIM];
      int d_w[DIM];
      int s_w[DIM];
      int dim_counter[DIM];
      for (int i = 0; i < DIM; i++) {
         box_w[i] = ibox.numberCells(i);
         d_w[i] = d_box.numberCells(i);
         s_w[i] = s_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const int d_offset = dst.getOffset();
      const int s_offset = src.getOffset();

      const int num_d0_blocks = ibox.size() / box_w[0];

      int d_begin = d_box.offset(ibox.lower());
      int s_begin = s_box.offset(ibox.lower());

      TYPE* dd    = dst.getPointer();
      const TYPE* sd    = src.getPointer();

      for (int d = 0; d < ddepth; d++) {

         int d_counter = d_begin;
         int s_counter = s_begin;

         int d_b[DIM];
         int s_b[DIM];
         for (int nd = 0; nd < DIM; nd++) {
            d_b[nd] = d_counter;
            s_b[nd] = s_counter;
         }

         for (int nb = 0; nb < num_d0_blocks; nb++) {

            for (int i0 = 0; i0 < box_w[0]; i0++) {
               if ( sd[s_counter+i0] == 0.0 ) {
                  test = 0;
                  dd[d_counter+i0] = 0.0F;
               } else {
                  dd[d_counter+i0] = 1.0F / sd[s_counter+i0];
               }
            }

            int dim_jump = 0;

            for (int j = 1; j < DIM; j++) {
               if (dim_counter[j] < box_w[j]-1) {
                  ++dim_counter[j];
                  dim_jump = j;
                  break;
               } else {
                  dim_counter[j] = 0;
               }
            }
            if (dim_jump > 0) {
               int d_step = 1;
               int s_step = 1;
               for (int k = 0; k < dim_jump; k++) {
                  d_step *= d_w[k];
                  s_step *= s_w[k];
               }
               d_counter = d_b[dim_jump-1] + d_step;
               s_counter = s_b[dim_jump-1] + s_step;

               for (int m = 0; m < dim_jump; m++) {
                  d_b[m] = d_counter;
                  s_b[m] = s_counter;
               }
            }
         }

         d_begin += d_offset;
         s_begin += s_offset;
      }

   }

   return( test );
}

template<int DIM, class TYPE>
TYPE 
ArrayDataMiscellaneousOpsReal<DIM,TYPE>::maxPointwiseDivide(
   const pdat::ArrayData<DIM,TYPE>& numer,
   const pdat::ArrayData<DIM,TYPE>& denom,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(denom.getDepth() == numer.getDepth());
#endif

   TYPE max = 0.0, quot;

   const hier::Box<DIM> n_box = numer.getBox();
   const hier::Box<DIM> d_box = denom.getBox();
   const hier::Box<DIM> ibox = box * d_box * n_box;

   if (!ibox.empty()) {
      const int ddepth  = denom.getDepth();

      int box_w[DIM];
      int n_w[DIM];
      int d_w[DIM];
      int dim_counter[DIM];
      for (int i = 0; i < DIM; i++) {
         box_w[i] = ibox.numberCells(i);
         n_w[i] = n_box.numberCells(i);
         d_w[i] = d_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const int n_offset = numer.getOffset();
      const int d_offset = denom.getOffset();

      const int num_d0_blocks = ibox.size() / box_w[0];

      int n_begin = n_box.offset(ibox.lower());
      int d_begin = d_box.offset(ibox.lower());

      const TYPE* nd    = numer.getPointer();
      const TYPE* dd    = denom.getPointer();

      for (int d = 0; d < ddepth; d++) {

         int n_counter = n_begin;
         int d_counter = d_begin;

         int n_b[DIM];
         int d_b[DIM];
         for (int nm = 0; nm < DIM; nm++) {
            n_b[nm] = n_counter;
            d_b[nm] = d_counter;
         }

         for (int nb = 0; nb < num_d0_blocks; nb++) {

            for (int i0 = 0; i0 < box_w[0]; i0++) {
               if ( dd[d_counter+i0] == 0.0 ) {
                  quot = tbox::MathUtilities<TYPE>::Abs(nd[n_counter+i0]);
               } else {
                  quot = tbox::MathUtilities<TYPE>::Abs(nd[n_counter+i0] /
                         dd[d_counter+i0]);
               }
               if ( max < quot ) max = quot;
            }
            int dim_jump = 0;

            for (int j = 1; j < DIM; j++) {
               if (dim_counter[j] < box_w[j]-1) {
                  ++dim_counter[j];
                  dim_jump = j;
                  break;
               } else {
                  dim_counter[j] = 0;
               }
            }

            if (dim_jump > 0) {
               int n_step = 1;
               int d_step = 1;
               for (int k = 0; k < dim_jump; k++) {
                  n_step *= n_w[k];
                  d_step *= d_w[k];
               }
               n_counter = n_b[dim_jump-1] + n_step;
               d_counter = d_b[dim_jump-1] + d_step;

               for (int m = 0; m < dim_jump; m++) {
                  n_b[m] = n_counter;
                  d_b[m] = d_counter;
               }
            }
         }

         n_begin += n_offset;
         d_begin += d_offset;
      }

   }

   return( max );
}


template <int DIM, class TYPE>
TYPE 
ArrayDataMiscellaneousOpsReal<DIM,TYPE>::minPointwiseDivide(
   const pdat::ArrayData<DIM,TYPE>& numer,
   const pdat::ArrayData<DIM,TYPE>& denom,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(denom.getDepth() == numer.getDepth());
#endif

   TYPE min  = tbox::MathUtilities<TYPE>::getMax();
   TYPE quot = tbox::MathUtilities<TYPE>::getMax();

   const hier::Box<DIM> n_box = numer.getBox();
   const hier::Box<DIM> d_box = denom.getBox();
   const hier::Box<DIM> ibox = box * d_box * n_box;

   if (!ibox.empty()) {
      const int ddepth  = denom.getDepth();

      int box_w[DIM];
      int n_w[DIM];
      int d_w[DIM];
      int dim_counter[DIM];
      for (int i = 0; i < DIM; i++) {
         box_w[i] = ibox.numberCells(i);
         n_w[i] = n_box.numberCells(i);
         d_w[i] = d_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const int n_offset = numer.getOffset();
      const int d_offset = denom.getOffset();

      const int num_d0_blocks = ibox.size() / box_w[0];

      int n_begin = n_box.offset(ibox.lower());
      int d_begin = d_box.offset(ibox.lower());

      const TYPE* nd    = numer.getPointer();
      const TYPE* dd    = denom.getPointer();

      for (int d = 0; d < ddepth; d++) {

         int n_counter = n_begin;
         int d_counter = d_begin;

         int n_b[DIM];
         int d_b[DIM];
         for (int nm = 0; nm < DIM; nm++) {
            n_b[nm] = n_counter;
            d_b[nm] = d_counter;
         }

         for (int nb = 0; nb < num_d0_blocks; nb++) {

            for (int i0 = 0; i0 < box_w[0]; i0++) {
               if ( dd[d_counter+i0] != 0.0 ) {
                  quot = nd[n_counter+i0]/dd[d_counter+i0];
               }
               if ( quot < min ) min = quot;
            }
            int dim_jump = 0;

            for (int j = 1; j < DIM; j++) {
               if (dim_counter[j] < box_w[j]-1) {
                  ++dim_counter[j];
                  dim_jump = j;
                  break;
               } else {
                  dim_counter[j] = 0;
               }
            }

            if (dim_jump > 0) {
               int n_step = 1;
               int d_step = 1;
               for (int k = 0; k < dim_jump; k++) {
                  n_step *= n_w[k];
                  d_step *= d_w[k];
               }
               n_counter = n_b[dim_jump-1] + n_step;
               d_counter = d_b[dim_jump-1] + d_step;

               for (int m = 0; m < dim_jump; m++) {
                  n_b[m] = n_counter;
                  d_b[m] = d_counter;
               }
            }
         }

         n_begin += n_offset;
         d_begin += d_offset;
      }

   }

   return( min );
}

}
}
#endif
