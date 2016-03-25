//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/mathops/array/ArrayDataBasicOps.C $
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Basic templated opertions for array data.
//

#ifndef included_math_ArrayDataBasicOps_C
#define included_math_ArrayDataBasicOps_C

#include "ArrayDataBasicOps.h"

#include "tbox/MathUtilities.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#include "tbox/Utilities.h"
#endif



namespace SAMRAI {
    namespace math {

template<int DIM, class TYPE>
ArrayDataBasicOps<DIM,TYPE>::ArrayDataBasicOps()
{
}

template<int DIM, class TYPE>
ArrayDataBasicOps<DIM,TYPE>::~ArrayDataBasicOps()
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
ArrayDataBasicOps<DIM,TYPE>::ArrayDataBasicOps(
   const ArrayDataBasicOps<DIM,TYPE>& foo)
{
   (void) foo;	// not implemented (but needed by some compilers)
}

template<int DIM, class TYPE>
void ArrayDataBasicOps<DIM,TYPE>::operator=(
   const ArrayDataBasicOps<DIM,TYPE>& foo)
{
   (void) foo;	// not implemented (but needed by some compilers)
}

/*
*************************************************************************
*									*
* General templated operations for array data.                          *
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
void ArrayDataBasicOps<DIM,TYPE>::scale(
   pdat::ArrayData<DIM,TYPE>& dst,
   const TYPE& alpha,
   const pdat::ArrayData<DIM,TYPE>& src,
   const hier::Box<DIM>& box) const
{
   if (alpha == tbox::MathUtilities<TYPE>::getZero()) {
      dst.fillAll(alpha, box);
   } else if (alpha == tbox::MathUtilities<TYPE>::getOne()) {
      dst.copy(src, box);
   } else {
      const int ddepth = dst.getDepth();
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(ddepth == src.getDepth());
#endif
      const hier::Box<DIM> dst_box = dst.getBox();
      const hier::Box<DIM> src_box = src.getBox();
      const hier::Box<DIM> ibox = box * dst_box * src_box;

      if (!ibox.empty()) {

         int box_w[DIM];
         int dst_w[DIM];
         int src_w[DIM];
         int dim_counter[DIM];
         for (int i = 0; i < DIM; i++) {
            box_w[i] = ibox.numberCells(i);
            dst_w[i] = dst_box.numberCells(i);
            src_w[i] = src_box.numberCells(i);
            dim_counter[i] = 0;
         }

         const int dst_offset = dst.getOffset();
         const int src_offset = src.getOffset();

         const int num_d0_blocks = ibox.size() / box_w[0];

         int dst_begin = dst_box.offset(ibox.lower());
         int src_begin = src_box.offset(ibox.lower());
                            
         TYPE* dd       = dst.getPointer();
         const TYPE* sd = src.getPointer();

         for (int d = 0; d < ddepth; d++) {

            int dst_counter = dst_begin;
            int src_counter = src_begin;

            int dst_b[DIM];
            int src_b[DIM];
            for (int nd = 0; nd < DIM; nd++) {
               dst_b[nd] = dst_counter;
               src_b[nd] = src_counter;
            }

            /*
             * Loop over each contiguous block of data.
             */
            for (int nb = 0; nb < num_d0_blocks; nb++) {

               for (int i0 = 0; i0 < box_w[0]; i0++) {
                  dd[dst_counter+i0] = alpha * sd[src_counter+i0];
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
                  int dst_step = 1;
                  int src_step = 1;
                  for (int k = 0; k < dim_jump; k++) {
                     dst_step *= dst_w[k];
                     src_step *= src_w[k];
                  }
                  dst_counter = dst_b[dim_jump-1] + dst_step;
                  src_counter = src_b[dim_jump-1] + src_step;

                  for (int m = 0; m < dim_jump; m++) {
                     dst_b[m] = dst_counter;
                     src_b[m] = src_counter;
                  }

               }
            }

            dst_begin += dst_offset;
            src_begin += src_offset;

         }
      }
   }
}

template<int DIM, class TYPE>
void ArrayDataBasicOps<DIM,TYPE>::addScalar(
   pdat::ArrayData<DIM,TYPE>& dst, 
   const pdat::ArrayData<DIM,TYPE>& src,
   const TYPE& alpha, 
   const hier::Box<DIM>& box) const
{
   if (alpha == tbox::MathUtilities<TYPE>::getZero()) {
      dst.copy(src, box);
   } else {
      const int ddepth = dst.getDepth();
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(ddepth == src.getDepth());
#endif
      const hier::Box<DIM> dst_box = dst.getBox();
      const hier::Box<DIM> src_box = src.getBox();
      const hier::Box<DIM> ibox = box * dst_box * src_box;

      if (!ibox.empty()) {

         int box_w[DIM];
         int dst_w[DIM];
         int src_w[DIM];
         int dim_counter[DIM];
         for (int i = 0; i < DIM; i++) {
            box_w[i] = ibox.numberCells(i);
            dst_w[i] = dst_box.numberCells(i);
            src_w[i] = src_box.numberCells(i);
            dim_counter[i] = 0;
         }

         const int dst_offset = dst.getOffset();
         const int src_offset = src.getOffset();

         const int num_d0_blocks = ibox.size() / box_w[0];

         int dst_begin = dst_box.offset(ibox.lower());
         int src_begin = src_box.offset(ibox.lower());

         TYPE* dd       = dst.getPointer();
         const TYPE* sd = src.getPointer();

         for (int d = 0; d < ddepth; d++) {

            int dst_counter = dst_begin;
            int src_counter = src_begin;

            int dst_b[DIM];
            int src_b[DIM];
            for (int nd = 0; nd < DIM; nd++) {
               dst_b[nd] = dst_counter;
               src_b[nd] = src_counter;
            }

            /*
             * Loop over each contiguous block of data.
             */
            for (int nb = 0; nb < num_d0_blocks; nb++) {

               for (int i0 = 0; i0 < box_w[0]; i0++) {
                  dd[dst_counter+i0] = alpha + sd[src_counter+i0];
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
                  int dst_step = 1;
                  int src_step = 1;
                  for (int k = 0; k < dim_jump; k++) {
                     dst_step *= dst_w[k];
                     src_step *= src_w[k];
                  }
                  dst_counter = dst_b[dim_jump-1] + dst_step;
                  src_counter = src_b[dim_jump-1] + src_step;

                  for (int m = 0; m < dim_jump; m++) {
                     dst_b[m] = dst_counter;
                     src_b[m] = src_counter;
                  }

               }
            }

            dst_begin += dst_offset;
            src_begin += src_offset;

         }
      }
   }
}

template<int DIM, class TYPE>
void ArrayDataBasicOps<DIM,TYPE>::add(
   pdat::ArrayData<DIM,TYPE>& dst,
   const pdat::ArrayData<DIM,TYPE>& src1,
   const pdat::ArrayData<DIM,TYPE>& src2,
   const hier::Box<DIM>& box) const
{
   const int ddepth = dst.getDepth();
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(ddepth == src1.getDepth() && ddepth == src2.getDepth());
#endif

   const hier::Box<DIM> dst_box  = dst.getBox();
   const hier::Box<DIM> src1_box = src1.getBox();
   const hier::Box<DIM> src2_box = src2.getBox();
   const hier::Box<DIM> ibox = box * dst_box * src1_box * src2_box;

   if (!ibox.empty()) {

      int box_w[DIM];
      int dst_w[DIM];
      int src1_w[DIM];
      int src2_w[DIM];
      int dim_counter[DIM];
      for (int i = 0; i < DIM; i++) {
         box_w[i] = ibox.numberCells(i);
         dst_w[i] = dst_box.numberCells(i);
         src1_w[i] = src1_box.numberCells(i);
         src2_w[i] = src2_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const int dst_offset = dst.getOffset();
      const int src1_offset = src1.getOffset();
      const int src2_offset = src2.getOffset();

      const int num_d0_blocks = ibox.size() / box_w[0];

      int dst_begin = dst_box.offset(ibox.lower());
      int src1_begin = src1_box.offset(ibox.lower());
      int src2_begin = src2_box.offset(ibox.lower());

      TYPE* dd       = dst.getPointer();
      const TYPE* s1d = src1.getPointer();
      const TYPE* s2d = src2.getPointer();

      for (int d = 0; d < ddepth; d++) {

         int dst_counter = dst_begin;
         int src1_counter = src1_begin;
         int src2_counter = src2_begin;

         int dst_b[DIM];
         int src1_b[DIM];
         int src2_b[DIM];
         for (int nd = 0; nd < DIM; nd++) {
            dst_b[nd] = dst_counter;
            src1_b[nd] = src1_counter;
            src2_b[nd] = src2_counter;
         }

         /*
          * Loop over each contiguous block of data.
          */
         for (int nb = 0; nb < num_d0_blocks; nb++) {

            for (int i0 = 0; i0 < box_w[0]; i0++) {
               dd[dst_counter+i0] = s1d[src1_counter+i0] +
                                    s2d[src2_counter+i0];
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
               int dst_step = 1;
               int src1_step = 1;
               int src2_step = 1;
               for (int k = 0; k < dim_jump; k++) {
                  dst_step *= dst_w[k];
                  src1_step *= src1_w[k];
                  src2_step *= src2_w[k];
               }
               dst_counter = dst_b[dim_jump-1] + dst_step;
               src1_counter = src1_b[dim_jump-1] + src1_step;
               src2_counter = src2_b[dim_jump-1] + src2_step;

               for (int m = 0; m < dim_jump; m++) {
                  dst_b[m] = dst_counter;
                  src1_b[m] = src1_counter;
                  src2_b[m] = src2_counter;
               }

            }
         }

         dst_begin += dst_offset;
         src1_begin += src1_offset;
         src2_begin += src2_offset;

      }
   }
}

template<int DIM, class TYPE>
void ArrayDataBasicOps<DIM,TYPE>::subtract(
   pdat::ArrayData<DIM,TYPE>& dst,
   const pdat::ArrayData<DIM,TYPE>& src1,
   const pdat::ArrayData<DIM,TYPE>& src2,
   const hier::Box<DIM>& box) const
{
   const int ddepth = dst.getDepth();
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(ddepth == src1.getDepth() && ddepth == src2.getDepth());
#endif

   const hier::Box<DIM> dst_box  = dst.getBox();
   const hier::Box<DIM> src1_box = src1.getBox();
   const hier::Box<DIM> src2_box = src2.getBox();
   const hier::Box<DIM> ibox = box * dst_box * src1_box * src2_box;

   if (!ibox.empty()) {

      int box_w[DIM];
      int dst_w[DIM];
      int src1_w[DIM];
      int src2_w[DIM];
      int dim_counter[DIM];
      for (int i = 0; i < DIM; i++) {
         box_w[i] = ibox.numberCells(i);
         dst_w[i] = dst_box.numberCells(i);
         src1_w[i] = src1_box.numberCells(i);
         src2_w[i] = src2_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const int dst_offset = dst.getOffset();
      const int src1_offset = src1.getOffset();
      const int src2_offset = src2.getOffset();

      const int num_d0_blocks = ibox.size() / box_w[0];

      int dst_begin = dst_box.offset(ibox.lower());
      int src1_begin = src1_box.offset(ibox.lower());
      int src2_begin = src2_box.offset(ibox.lower());

      TYPE* dd       = dst.getPointer();
      const TYPE* s1d = src1.getPointer();
      const TYPE* s2d = src2.getPointer();

      for (int d = 0; d < ddepth; d++) {

         int dst_counter = dst_begin;
         int src1_counter = src1_begin;
         int src2_counter = src2_begin;

         int dst_b[DIM];
         int src1_b[DIM];
         int src2_b[DIM];
         for (int nd = 0; nd < DIM; nd++) {
            dst_b[nd] = dst_counter;
            src1_b[nd] = src1_counter;
            src2_b[nd] = src2_counter;
         }

         /*
          * Loop over each contiguous block of data.
          */
         for (int nb = 0; nb < num_d0_blocks; nb++) {

            for (int i0 = 0; i0 < box_w[0]; i0++) {
               dd[dst_counter+i0] = s1d[src1_counter+i0] -
                                    s2d[src2_counter+i0];
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
               int dst_step = 1;
               int src1_step = 1;
               int src2_step = 1;
               for (int k = 0; k < dim_jump; k++) {
                  dst_step *= dst_w[k];
                  src1_step *= src1_w[k];
                  src2_step *= src2_w[k];
               }
               dst_counter = dst_b[dim_jump-1] + dst_step;
               src1_counter = src1_b[dim_jump-1] + src1_step;
               src2_counter = src2_b[dim_jump-1] + src2_step;

               for (int m = 0; m < dim_jump; m++) {
                  dst_b[m] = dst_counter;
                  src1_b[m] = src1_counter;
                  src2_b[m] = src2_counter;
               }

            }
         }

         dst_begin += dst_offset;
         src1_begin += src1_offset;
         src2_begin += src2_offset;

      }
   }
}

template<int DIM, class TYPE>
void ArrayDataBasicOps<DIM,TYPE>::multiply(
   pdat::ArrayData<DIM,TYPE>& dst,
   const pdat::ArrayData<DIM,TYPE>& src1,
   const pdat::ArrayData<DIM,TYPE>& src2,
   const hier::Box<DIM>& box) const
{
   const int ddepth = dst.getDepth();
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(ddepth == src1.getDepth() && ddepth == src2.getDepth());
#endif

   const hier::Box<DIM> dst_box  = dst.getBox();
   const hier::Box<DIM> src1_box = src1.getBox();
   const hier::Box<DIM> src2_box = src2.getBox();
   const hier::Box<DIM> ibox = box * dst_box * src1_box * src2_box;

   if (!ibox.empty()) {

      int box_w[DIM];
      int dst_w[DIM];
      int src1_w[DIM];
      int src2_w[DIM];
      int dim_counter[DIM];
      for (int i = 0; i < DIM; i++) {
         box_w[i] = ibox.numberCells(i);
         dst_w[i] = dst_box.numberCells(i);
         src1_w[i] = src1_box.numberCells(i);
         src2_w[i] = src2_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const int dst_offset = dst.getOffset();
      const int src1_offset = src1.getOffset();
      const int src2_offset = src2.getOffset();

      const int num_d0_blocks = ibox.size() / box_w[0];

      int dst_begin = dst_box.offset(ibox.lower());
      int src1_begin = src1_box.offset(ibox.lower());
      int src2_begin = src2_box.offset(ibox.lower());

      TYPE* dd       = dst.getPointer();
      const TYPE* s1d = src1.getPointer();
      const TYPE* s2d = src2.getPointer();

      for (int d = 0; d < ddepth; d++) {

         int dst_counter = dst_begin;
         int src1_counter = src1_begin;
         int src2_counter = src2_begin;

         int dst_b[DIM];
         int src1_b[DIM];
         int src2_b[DIM];
         for (int nd = 0; nd < DIM; nd++) {
            dst_b[nd] = dst_counter;
            src1_b[nd] = src1_counter;
            src2_b[nd] = src2_counter;
         }

         /*
          * Loop over each contiguous block of data.
          */
         for (int nb = 0; nb < num_d0_blocks; nb++) {

            for (int i0 = 0; i0 < box_w[0]; i0++) {
               dd[dst_counter+i0] = s1d[src1_counter+i0] *
                                    s2d[src2_counter+i0];
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
               int dst_step = 1;
               int src1_step = 1;
               int src2_step = 1;
               for (int k = 0; k < dim_jump; k++) {
                  dst_step *= dst_w[k];
                  src1_step *= src1_w[k];
                  src2_step *= src2_w[k];
               }
               dst_counter = dst_b[dim_jump-1] + dst_step;
               src1_counter = src1_b[dim_jump-1] + src1_step;
               src2_counter = src2_b[dim_jump-1] + src2_step;

               for (int m = 0; m < dim_jump; m++) {
                  dst_b[m] = dst_counter;
                  src1_b[m] = src1_counter;
                  src2_b[m] = src2_counter;
               }

            }
         }

         dst_begin += dst_offset;
         src1_begin += src1_offset;
         src2_begin += src2_offset;

      }
   }
}

template<int DIM, class TYPE>
void ArrayDataBasicOps<DIM,TYPE>::divide(
   pdat::ArrayData<DIM,TYPE>& dst,
   const pdat::ArrayData<DIM,TYPE>& src1,
   const pdat::ArrayData<DIM,TYPE>& src2,
   const hier::Box<DIM>& box) const
{
   const int ddepth = dst.getDepth();
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(ddepth == src1.getDepth() && ddepth == src2.getDepth());
#endif

   const hier::Box<DIM> dst_box  = dst.getBox();
   const hier::Box<DIM> src1_box = src1.getBox();
   const hier::Box<DIM> src2_box = src2.getBox();
   const hier::Box<DIM> ibox = box * dst_box * src1_box * src2_box;

   if (!ibox.empty()) {

      int box_w[DIM];
      int dst_w[DIM];
      int src1_w[DIM];
      int src2_w[DIM];
      int dim_counter[DIM];
      for (int i = 0; i < DIM; i++) {
         box_w[i] = ibox.numberCells(i);
         dst_w[i] = dst_box.numberCells(i);
         src1_w[i] = src1_box.numberCells(i);
         src2_w[i] = src2_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const int dst_offset = dst.getOffset();
      const int src1_offset = src1.getOffset();
      const int src2_offset = src2.getOffset();

      const int num_d0_blocks = ibox.size() / box_w[0];

      int dst_begin = dst_box.offset(ibox.lower());
      int src1_begin = src1_box.offset(ibox.lower());
      int src2_begin = src2_box.offset(ibox.lower());

      TYPE* dd       = dst.getPointer();
      const TYPE* s1d = src1.getPointer();
      const TYPE* s2d = src2.getPointer();

      for (int d = 0; d < ddepth; d++) {

         int dst_counter = dst_begin;
         int src1_counter = src1_begin;
         int src2_counter = src2_begin;

         int dst_b[DIM];
         int src1_b[DIM];
         int src2_b[DIM];
         for (int nd = 0; nd < DIM; nd++) {
            dst_b[nd] = dst_counter;
            src1_b[nd] = src1_counter;
            src2_b[nd] = src2_counter;
         }

         for (int nb = 0; nb < num_d0_blocks; nb++) {
                                   
            for (int i0 = 0; i0 < box_w[0]; i0++) {
               dd[dst_counter+i0] = s1d[src1_counter+i0] /
                                    s2d[src2_counter+i0];
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
               int dst_step = 1;
               int src1_step = 1;
               int src2_step = 1;
               for (int k = 0; k < dim_jump; k++) {
                  dst_step *= dst_w[k];
                  src1_step *= src1_w[k];
                  src2_step *= src2_w[k];
               }
               dst_counter = dst_b[dim_jump-1] + dst_step;
               src1_counter = src1_b[dim_jump-1] + src1_step;
               src2_counter = src2_b[dim_jump-1] + src2_step;

               for (int m = 0; m < dim_jump; m++) {
                  dst_b[m] = dst_counter;
                  src1_b[m] = src1_counter;
                  src2_b[m] = src2_counter;
               }

            }
         }

         dst_begin += dst_offset;
         src1_begin += src1_offset;
         src2_begin += src2_offset;

      }
   }
}

template<int DIM, class TYPE>
void ArrayDataBasicOps<DIM,TYPE>::reciprocal(
   pdat::ArrayData<DIM,TYPE>& dst,
   const pdat::ArrayData<DIM,TYPE>& src,
   const hier::Box<DIM>& box) const
{
   const int ddepth = dst.getDepth();
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(ddepth == src.getDepth());
#endif

   const hier::Box<DIM> dst_box  = dst.getBox();
   const hier::Box<DIM> src_box = src.getBox();
   const hier::Box<DIM> ibox = box * dst_box * src_box;

   if (!ibox.empty()) {

      int box_w[DIM];
      int dst_w[DIM];
      int src_w[DIM];
      int dim_counter[DIM];
      for (int i = 0; i < DIM; i++) {
         box_w[i] = ibox.numberCells(i);
         dst_w[i] = dst_box.numberCells(i);
         src_w[i] = src_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const int dst_offset = dst.getOffset();
      const int src_offset = src.getOffset();

      const int num_d0_blocks = ibox.size() / box_w[0];

      int dst_begin = dst_box.offset(ibox.lower());
      int src_begin = src_box.offset(ibox.lower());

      TYPE* dd       = dst.getPointer();
      const TYPE* sd = src.getPointer();

      for (int d = 0; d < ddepth; d++) {

         int dst_counter = dst_begin;
         int src_counter = src_begin;

         int dst_b[DIM];
         int src_b[DIM];
         for (int nd = 0; nd < DIM; nd++) {
            dst_b[nd] = dst_counter;
            src_b[nd] = src_counter;
         }

         /*
          * Loop over each contiguous block of data.
          */
         for (int nb = 0; nb < num_d0_blocks; nb++) {

            for (int i0 = 0; i0 < box_w[0]; i0++) {
               dd[dst_counter+i0] = 
                  tbox::MathUtilities<TYPE>::getOne() / sd[src_counter+i0];
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
               int dst_step = 1;
               int src_step = 1;
               for (int k = 0; k < dim_jump; k++) {
                  dst_step *= dst_w[k];
                  src_step *= src_w[k];
               }
               dst_counter = dst_b[dim_jump-1] + dst_step;
               src_counter = src_b[dim_jump-1] + src_step;

               for (int m = 0; m < dim_jump; m++) {
                  dst_b[m] = dst_counter;
                  src_b[m] = src_counter;
               }

            }
         }

         dst_begin += dst_offset;
         src_begin += src_offset;

      }
   }
}

template<int DIM, class TYPE>
void ArrayDataBasicOps<DIM,TYPE>::linearSum(
   pdat::ArrayData<DIM,TYPE>& dst,
   const TYPE& alpha,
   const pdat::ArrayData<DIM,TYPE>& src1,
   const TYPE& beta,
   const pdat::ArrayData<DIM,TYPE>& src2,
   const hier::Box<DIM>& box) const
{
   const int ddepth = dst.getDepth();
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(ddepth == src1.getDepth() && ddepth == src2.getDepth());
#endif
   if (alpha == tbox::MathUtilities<TYPE>::getZero()) {
      if (beta == tbox::MathUtilities<TYPE>::getZero()) {
         dst.fillAll(tbox::MathUtilities<TYPE>::getZero(), box);
      } else { 
         scale(dst, beta, src2, box);
      }
   } else if (beta == tbox::MathUtilities<TYPE>::getZero()) {
      scale(dst, alpha, src1, box);
   } else if (alpha == tbox::MathUtilities<TYPE>::getOne()) {
      axpy(dst, beta, src2, src1, box);
   } else if (beta == tbox::MathUtilities<TYPE>::getOne()) {
      axpy(dst, alpha, src1, src2, box);
   } else if (alpha == -tbox::MathUtilities<TYPE>::getOne()) {
      axmy(dst, beta, src2, src1, box);
   } else if (beta == -tbox::MathUtilities<TYPE>::getOne()) {
      axmy(dst, alpha, src1, src2, box);
   } else {

      const hier::Box<DIM> dst_box  = dst.getBox();
      const hier::Box<DIM> src1_box = src1.getBox();
      const hier::Box<DIM> src2_box = src2.getBox();
      const hier::Box<DIM> ibox = box * dst_box * src1_box * src2_box;

      if (!ibox.empty()) {

         int box_w[DIM];
         int dst_w[DIM];
         int src1_w[DIM];
         int src2_w[DIM];
         int dim_counter[DIM];
         for (int i = 0; i < DIM; i++) {
            box_w[i] = ibox.numberCells(i);
            dst_w[i] = dst_box.numberCells(i);
            src1_w[i] = src1_box.numberCells(i);
            src2_w[i] = src2_box.numberCells(i);
            dim_counter[i] = 0;
         }

         const int dst_offset = dst.getOffset();
         const int src1_offset = src1.getOffset();
         const int src2_offset = src2.getOffset();

         const int num_d0_blocks = ibox.size() / box_w[0];

         int dst_begin = dst_box.offset(ibox.lower());
         int src1_begin = src1_box.offset(ibox.lower());
         int src2_begin = src2_box.offset(ibox.lower());

         TYPE* dd       = dst.getPointer();
         const TYPE* s1d = src1.getPointer();
         const TYPE* s2d = src2.getPointer();

         for (int d = 0; d < ddepth; d++) {

            int dst_counter = dst_begin;
            int src1_counter = src1_begin;
            int src2_counter = src2_begin;

            int dst_b[DIM];
            int src1_b[DIM];
            int src2_b[DIM];
            for (int nd = 0; nd < DIM; nd++) {
               dst_b[nd] = dst_counter;
               src1_b[nd] = src1_counter;
               src2_b[nd] = src2_counter;
            }

            /*
             * Loop over each contiguous block of data.
             */
            for (int nb = 0; nb < num_d0_blocks; nb++) {

               for (int i0 = 0; i0 < box_w[0]; i0++) {
                  dd[dst_counter+i0] = alpha * s1d[src1_counter+i0] +
                                       beta * s2d[src2_counter+i0];
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
                  int dst_step = 1;
                  int src1_step = 1;
                  int src2_step = 1;
                  for (int k = 0; k < dim_jump; k++) {
                     dst_step *= dst_w[k];
                     src1_step *= src1_w[k];
                     src2_step *= src2_w[k];
                  }
                  dst_counter = dst_b[dim_jump-1] + dst_step;
                  src1_counter = src1_b[dim_jump-1] + src1_step;
                  src2_counter = src2_b[dim_jump-1] + src2_step;

                  for (int m = 0; m < dim_jump; m++) {
                     dst_b[m] = dst_counter;
                     src1_b[m] = src1_counter;
                     src2_b[m] = src2_counter;
                  }

               }
            }

            dst_begin += dst_offset;
            src1_begin += src1_offset;
            src2_begin += src2_offset;

         }
      }
   }
}

template<int DIM, class TYPE>
void ArrayDataBasicOps<DIM,TYPE>::axpy(
   pdat::ArrayData<DIM,TYPE>& dst,
   const TYPE& alpha,
   const pdat::ArrayData<DIM,TYPE>& src1,
   const pdat::ArrayData<DIM,TYPE>& src2,
   const hier::Box<DIM>& box) const
{
   if (alpha == tbox::MathUtilities<TYPE>::getZero()) {
      dst.copy(src2, box);
   } else if (alpha == tbox::MathUtilities<TYPE>::getOne()) {
      add(dst, src1, src2, box);
   } else if (alpha == -tbox::MathUtilities<TYPE>::getOne()) {
      subtract(dst, src2, src1, box);
   } else {
      const int ddepth = dst.getDepth();
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(ddepth == src1.getDepth() && ddepth == src2.getDepth());
#endif
      const hier::Box<DIM> dst_box  = dst.getBox();
      const hier::Box<DIM> src1_box = src1.getBox();
      const hier::Box<DIM> src2_box = src2.getBox();
      const hier::Box<DIM> ibox = box * dst_box * src1_box * src2_box;

      if (!ibox.empty()) {

         int box_w[DIM];
         int dst_w[DIM];
         int src1_w[DIM];
         int src2_w[DIM];
         int dim_counter[DIM];
         for (int i = 0; i < DIM; i++) {
            box_w[i] = ibox.numberCells(i);
            dst_w[i] = dst_box.numberCells(i);
            src1_w[i] = src1_box.numberCells(i);
            src2_w[i] = src2_box.numberCells(i);
            dim_counter[i] = 0;
         }

         const int dst_offset = dst.getOffset();
         const int src1_offset = src1.getOffset();
         const int src2_offset = src2.getOffset();

         const int num_d0_blocks = ibox.size() / box_w[0];

         int dst_begin = dst_box.offset(ibox.lower());
         int src1_begin = src1_box.offset(ibox.lower());
         int src2_begin = src2_box.offset(ibox.lower());

         TYPE* dd       = dst.getPointer();
         const TYPE* s1d = src1.getPointer();
         const TYPE* s2d = src2.getPointer();

         for (int d = 0; d < ddepth; d++) {

            int dst_counter = dst_begin;
            int src1_counter = src1_begin;
            int src2_counter = src2_begin;

            int dst_b[DIM];
            int src1_b[DIM];
            int src2_b[DIM];
            for (int nd = 0; nd < DIM; nd++) {
               dst_b[nd] = dst_counter;
               src1_b[nd] = src1_counter;
               src2_b[nd] = src2_counter;
            }

            /*
             * Loop over each contiguous block of data.
             */
            for (int nb = 0; nb < num_d0_blocks; nb++) {

               for (int i0 = 0; i0 < box_w[0]; i0++) {
                  dd[dst_counter+i0] = alpha * s1d[src1_counter+i0] +
                                       s2d[src2_counter+i0];
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
                  int dst_step = 1;
                  int src1_step = 1;
                  int src2_step = 1;
                  for (int k = 0; k < dim_jump; k++) {
                     dst_step *= dst_w[k];
                     src1_step *= src1_w[k];
                     src2_step *= src2_w[k];
                  }
                  dst_counter = dst_b[dim_jump-1] + dst_step;
                  src1_counter = src1_b[dim_jump-1] + src1_step;
                  src2_counter = src2_b[dim_jump-1] + src2_step;

                  for (int m = 0; m < dim_jump; m++) {
                     dst_b[m] = dst_counter;
                     src1_b[m] = src1_counter;
                     src2_b[m] = src2_counter;
                  }

               }
            }

            dst_begin += dst_offset;
            src1_begin += src1_offset;
            src2_begin += src2_offset;

         }
      }
   }
}

template<int DIM, class TYPE>
void ArrayDataBasicOps<DIM,TYPE>::axmy(
   pdat::ArrayData<DIM,TYPE>& dst,
   const TYPE& alpha,
   const pdat::ArrayData<DIM,TYPE>& src1,
   const pdat::ArrayData<DIM,TYPE>& src2,
   const hier::Box<DIM>& box) const
{
   if (alpha == tbox::MathUtilities<TYPE>::getZero()) {
      scale(dst, -tbox::MathUtilities<TYPE>::getOne(), src2, box);
   } else if (alpha == tbox::MathUtilities<TYPE>::getOne()) {
      subtract(dst, src1, src2, box);      
   } else {
      const int ddepth = dst.getDepth();
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(ddepth == src1.getDepth() && ddepth == src2.getDepth());
#endif
      const hier::Box<DIM> dst_box  = dst.getBox();
      const hier::Box<DIM> src1_box = src1.getBox();
      const hier::Box<DIM> src2_box = src2.getBox();
      const hier::Box<DIM> ibox = box * dst_box * src1_box * src2_box;

      if (!ibox.empty()) {

         int box_w[DIM];
         int dst_w[DIM];
         int src1_w[DIM];
         int src2_w[DIM];
         int dim_counter[DIM];
         for (int i = 0; i < DIM; i++) {
            box_w[i] = ibox.numberCells(i);
            dst_w[i] = dst_box.numberCells(i);
            src1_w[i] = src1_box.numberCells(i);
            src2_w[i] = src2_box.numberCells(i);
            dim_counter[i] = 0;
         }

         const int dst_offset = dst.getOffset();
         const int src1_offset = src1.getOffset();
         const int src2_offset = src2.getOffset();

         const int num_d0_blocks = ibox.size() / box_w[0];

         int dst_begin = dst_box.offset(ibox.lower());
         int src1_begin = src1_box.offset(ibox.lower());
         int src2_begin = src2_box.offset(ibox.lower());

         TYPE* dd       = dst.getPointer();
         const TYPE* s1d = src1.getPointer();
         const TYPE* s2d = src2.getPointer();

         for (int d = 0; d < ddepth; d++) {

            int dst_counter = dst_begin;
            int src1_counter = src1_begin;
            int src2_counter = src2_begin;

            int dst_b[DIM];
            int src1_b[DIM];
            int src2_b[DIM];
            for (int nd = 0; nd < DIM; nd++) {
               dst_b[nd] = dst_counter;
               src1_b[nd] = src1_counter;
               src2_b[nd] = src2_counter;
            }

            /*
             * Loop over each contiguous block of data.
             */
            for (int nb = 0; nb < num_d0_blocks; nb++) {

               for (int i0 = 0; i0 < box_w[0]; i0++) {
                  dd[dst_counter+i0] = alpha * s1d[src1_counter+i0] -
                                       s2d[src2_counter+i0];
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
                  int dst_step = 1;
                  int src1_step = 1;
                  int src2_step = 1;
                  for (int k = 0; k < dim_jump; k++) {
                     dst_step *= dst_w[k];
                     src1_step *= src1_w[k];
                     src2_step *= src2_w[k];
                  }
                  dst_counter = dst_b[dim_jump-1] + dst_step;
                  src1_counter = src1_b[dim_jump-1] + src1_step;
                  src2_counter = src2_b[dim_jump-1] + src2_step;

                  for (int m = 0; m < dim_jump; m++) {
                     dst_b[m] = dst_counter;
                     src1_b[m] = src1_counter;
                     src2_b[m] = src2_counter;
                  }

               }
            }

            dst_begin += dst_offset;
            src1_begin += src1_offset;
            src2_begin += src2_offset;

         }
      }
   }
}

template<int DIM, class TYPE> 
TYPE ArrayDataBasicOps<DIM,TYPE>::min(
   const pdat::ArrayData<DIM,TYPE>& data,
   const hier::Box<DIM>& box) const
{
   TYPE minval = tbox::MathUtilities<TYPE>::getMax();

   const hier::Box<DIM> d_box = data.getBox();
   const hier::Box<DIM> ibox = box * d_box;

   if (!ibox.empty()) {
      int box_w[DIM];
      int d_w[DIM];
      int dim_counter[DIM];
      for (int i = 0; i < DIM; i++) {
         box_w[i] = ibox.numberCells(i);
         d_w[i] = d_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const int d_offset = data.getOffset();

      int d_begin = d_box.offset(ibox.lower());

      const int num_d0_blocks = ibox.size() / box_w[0];

      const int ddepth = data.getDepth();

      const TYPE* dd    = data.getPointer();

      for (int d = 0; d < ddepth; d++) {

         int d_counter = d_begin;

         int d_b[DIM];
         for (int nd = 0; nd < DIM; nd++) {
            d_b[nd] = d_counter;
         }

         for (int nb = 0; nb < num_d0_blocks; nb++) {

            for (int i0 = 0; i0 < box_w[0]; i0++) {
               minval = 
                  tbox::MathUtilities<TYPE>::Min( minval, dd[d_counter+i0] );
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
               for (int k = 0; k < dim_jump; k++) {
                  d_step *= d_w[k];
               }
               d_counter = d_b[dim_jump-1] + d_step;

               for (int m = 0; m < dim_jump; m++) {
                  d_b[m] = d_counter;
               }
            }
         }

         d_begin += d_offset;
      }
   }

   return(minval);
}

template<int DIM, class TYPE>
TYPE ArrayDataBasicOps<DIM,TYPE>::max(
   const pdat::ArrayData<DIM,TYPE>& data,
   const hier::Box<DIM>& box) const
{
   TYPE maxval = -(tbox::MathUtilities<TYPE>::getMax());

   const hier::Box<DIM> d_box = data.getBox();
   const hier::Box<DIM> ibox = box * d_box;

   if (!ibox.empty()) {

      int box_w[DIM];
      int d_w[DIM];
      int dim_counter[DIM];
      for (int i = 0; i < DIM; i++) {
         box_w[i] = ibox.numberCells(i);
         d_w[i] = d_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const int d_offset = data.getOffset();

      int d_begin = d_box.offset(ibox.lower());

      const int num_d0_blocks = ibox.size() / box_w[0];

      const int ddepth = data.getDepth();

      const TYPE* dd    = data.getPointer();

      for (int d = 0; d < ddepth; d++) {

         int d_counter = d_begin;

         int d_b[DIM];
         for (int nd = 0; nd < DIM; nd++) {
            d_b[nd] = d_counter;
         }

         for (int nb = 0; nb < num_d0_blocks; nb++) {

            for (int i0 = 0; i0 < box_w[0]; i0++) {
               maxval = 
                  tbox::MathUtilities<TYPE>::Max( maxval, dd[d_counter+i0] );
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
               for (int k = 0; k < dim_jump; k++) {
                  d_step *= d_w[k];
               }
               d_counter = d_b[dim_jump-1] + d_step;

               for (int m = 0; m < dim_jump; m++) {
                  d_b[m] = d_counter;
               }
            }
         }

         d_begin += d_offset;
      }
   }
   
   return maxval;
}


template<int DIM, class TYPE>
void ArrayDataBasicOps<DIM,TYPE>::setRandomValues(
   pdat::ArrayData<DIM,TYPE>& dst,
   const TYPE& width,
   const TYPE& low,
   const hier::Box<DIM>& box) const
{
   const hier::Box<DIM> d_box = dst.getBox();
   const hier::Box<DIM> ibox = box * d_box;

   if (!ibox.empty()) {

      int box_w[DIM];
      int d_w[DIM];
      int dim_counter[DIM];
      for (int i = 0; i < DIM; i++) {
         box_w[i] = ibox.numberCells(i);
         d_w[i] = d_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const int d_offset = dst.getOffset();

      int d_begin = d_box.offset(ibox.lower());

      const int num_d0_blocks = ibox.size() / box_w[0];

      const int ddepth = dst.getDepth();

      TYPE* dd    = dst.getPointer();

      for (int d = 0; d < ddepth; d++) {

         int d_counter = d_begin;

         int d_b[DIM];
         for (int nd = 0; nd < DIM; nd++) {
            d_b[nd] = d_counter;
         }

         for (int nb = 0; nb < num_d0_blocks; nb++) {

            for (int i0 = 0; i0 < box_w[0]; i0++) {
               dd[d_counter+i0] =  tbox::MathUtilities<TYPE>::Rand(low,width);
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
               for (int k = 0; k < dim_jump; k++) {
                  d_step *= d_w[k];
               }
               d_counter = d_b[dim_jump-1] + d_step;

               for (int m = 0; m < dim_jump; m++) {
                  d_b[m] = d_counter;
               }
            }
         }

         d_begin += d_offset;
      }
   }
}

}
}
#endif
