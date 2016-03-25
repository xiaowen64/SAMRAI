//
// File:        ArrayDataNormOpsInteger.C
// Package:     SAMRAI mathops
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Norm operations for integer data arrays.
//

#ifndef included_tbox_ArrayDataNormOpsInteger_C
#define included_tbox_ArrayDataNormOpsInteger_C

#include "ArrayDataNormOpsInteger.h"
#include "tbox/Utilities.h"
#include <stdlib.h>
#include <math.h>
#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

namespace SAMRAI {
    namespace math {

template<int DIM>  ArrayDataNormOpsInteger<DIM>::ArrayDataNormOpsInteger()
{
}

template<int DIM>  ArrayDataNormOpsInteger<DIM>::~ArrayDataNormOpsInteger()
{
}

/*
*************************************************************************
*                                                                       *
* Norm operations for integer array data.				*
*                                                                       *
*************************************************************************
*/

template<int DIM> void ArrayDataNormOpsInteger<DIM>::abs(
   pdat::ArrayData<DIM,int>& dst,
   const pdat::ArrayData<DIM,int>& src,
   const hier::Box<DIM>& box) const 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(dst.getDepth() == src.getDepth());
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

      int* dd       = dst.getPointer();
      const int* sd = src.getPointer();

      const int ddepth = dst.getDepth();
      for (int d = 0; d < ddepth; d++) {

         int dst_counter = dst_begin;
         int src_counter = src_begin;

         int dst_b[DIM];
         int src_b[DIM];
         for (int nd = 0; nd < DIM; nd++) {
            dst_b[nd] = dst_counter;
            src_b[nd] = src_counter;
         }

         for (int nb = 0; nb < num_d0_blocks; nb++) {

            for (int i0 = 0; i0 < box_w[0]; i0++) {
               dd[dst_counter+i0] = tbox::Utilities::iabs(sd[src_counter+i0]);
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
}
#endif
