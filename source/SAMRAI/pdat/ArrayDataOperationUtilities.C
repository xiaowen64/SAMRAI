/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2016 Lawrence Livermore National Security, LLC
 * Description:   Templated array data looping operations supporting patch data types
 *
 ************************************************************************/

#ifndef included_pdat_ArrayDataOperationUtilities_C
#define included_pdat_ArrayDataOperationUtilities_C

#include "SAMRAI/pdat/ArrayDataOperationUtilities.h"
#include "SAMRAI/pdat/ArrayData.h"
#include "SAMRAI/tbox/Utilities.h"
#if defined(HAVE_RAJA)
#include "SAMRAI/tbox/RAJA_API.h"
#endif

namespace SAMRAI {
namespace pdat {

/*
 *************************************************************************
 *
 * Function that performs specified operation involving source and
 * destination array data objects and puts result in destination array
 * data object using explicit dimension-generic looping constructs.
 *
 *************************************************************************
 */

template<class TYPE, class OP>
void ArrayDataOperationUtilities<TYPE, OP>::doArrayDataOperationOnBox(
   ArrayData<TYPE>& dst,
   const ArrayData<TYPE>& src,
   const hier::Box& opbox,
   const hier::IntVector& src_shift,
   unsigned int dst_start_depth,
   unsigned int src_start_depth,
   unsigned int num_depth,
   const OP& op)
{
   TBOX_ASSERT_OBJDIM_EQUALITY4(dst, src, opbox, src_shift);
   TBOX_ASSERT((dst_start_depth + num_depth <= dst.getDepth()));
   TBOX_ASSERT((src_start_depth + num_depth <= src.getDepth()));

   const tbox::Dimension& dim(dst.getDim());

   TYPE * const dst_ptr = dst.getPointer();
   const TYPE * const src_ptr = src.getPointer();

   const hier::Box& dst_box(dst.getBox());
   const hier::Box& src_box(src.getBox());

   int box_w[SAMRAI::MAX_DIM_VAL];
   int dst_w[SAMRAI::MAX_DIM_VAL];
   int src_w[SAMRAI::MAX_DIM_VAL];
   int dim_counter[SAMRAI::MAX_DIM_VAL];
   for (tbox::Dimension::dir_t i = 0; i < dim.getValue(); ++i) {
      box_w[i] = opbox.numberCells(i);
      dst_w[i] = dst_box.numberCells(i);
      src_w[i] = src_box.numberCells(i);
      dim_counter[i] = 0;
   }

   const size_t dst_offset = dst.getOffset();
   const size_t src_offset = src.getOffset();

   /*
    * Data on the opbox can be decomposed into a set of
    * contiguous array sections representing data in a straight line
    * in the 0 coordinate direction.
    *
    * num_d0_blocks is the number of such array sections.
    * dst_begin, src_begin are the array indices for the first
    * data items in each array section to be copied.
    */

   const int num_d0_blocks = static_cast<int>(opbox.size() / box_w[0]);

   size_t dst_begin = dst_box.offset(opbox.lower())
      + dst_start_depth * dst_offset;
   size_t src_begin = src_box.offset(opbox.lower() - src_shift)
      + src_start_depth * src_offset;

   /*
    * Loop over the depth sections of the data arrays.
    */

   for (unsigned int d = 0; d < num_depth; ++d) {

#if defined(HAVE_RAJA)
     switch (dim.getValue()) {
       case 1: {
         //tbox::ArrayView<1, const TYPE> source(src, d);
         tbox::ArrayView<1, const TYPE> source(src, d);
         tbox::ArrayView<1, TYPE> dest(dst, d);

         tbox::for_all1< tbox::policy::parallel >(opbox, [=] SAMRAI_DEVICE (int i) {
             op(dest(i), source(i));
         });
         } break;
               
       case 2: {
         // tbox::ArrayView<2, const TYPE> source(src, d);
         tbox::ArrayView<2, const TYPE> source(src, d);
         tbox::ArrayView<2, TYPE> dest(dst, d);

         tbox::for_all2<tbox::policy::parallel>(opbox, [=] SAMRAI_DEVICE (int j, int i) {
             op(dest(i, j), source(i, j));
         });
         } break;

       case 3: {

         // tbox::ArrayView<3, const TYPE> source(src, d);
         tbox::ArrayView<3, const TYPE> source(src, d);
         tbox::ArrayView<3, TYPE> dest(dst, d);

         tbox::for_all3<tbox::policy::parallel>(opbox, [=] SAMRAI_DEVICE (int k, int j, int i) {
             op(dest(i,j,k), source(i,j,k));
         });
         } break;

       default: tbox::perr << "Dimension > 3 not supported with RAJA parallel kernels" << std::endl;
                TBOX_ERROR("Aborting in " __FILE__);
                break;
     }

#else

      size_t dst_counter = dst_begin;
      size_t src_counter = src_begin;

      size_t dst_b[SAMRAI::MAX_DIM_VAL];
      size_t src_b[SAMRAI::MAX_DIM_VAL];
      for (tbox::Dimension::dir_t nd = 0; nd < dim.getValue(); ++nd) {
         dst_b[nd] = dst_counter;
         src_b[nd] = src_counter;
      }

      /*
       * Loop over each contiguous block of data.
       */

      for (int nb = 0; nb < num_d0_blocks; ++nb) {

         for (int i0 = 0; i0 < box_w[0]; ++i0) {
            op(dst_ptr[dst_counter + i0], src_ptr[src_counter + i0]);
         }
         int dim_jump = 0;

         /*
          * After each contiguous block is copied, calculate the
          * beginning array index for the next block.
          */

         for (tbox::Dimension::dir_t j = 1; j < dim.getValue(); ++j) {
            if (dim_counter[j] < box_w[j] - 1) {
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
            for (int k = 0; k < dim_jump; ++k) {
               dst_step *= dst_w[k];
               src_step *= src_w[k];
            }
            dst_counter = dst_b[dim_jump - 1] + dst_step;
            src_counter = src_b[dim_jump - 1] + src_step;

            for (int m = 0; m < dim_jump; ++m) {
               dst_b[m] = dst_counter;
               src_b[m] = src_counter;
            }

         }  // if dim_jump > 0

      }  // nb loop over contiguous data blocks

      /*
       * After copy is complete on a full box for one depth index,
       * advance by the offset values.
       */

      dst_begin += dst_offset;
      src_begin += src_offset;
#endif

   }  // d loop over depth indices

}

/*
 *************************************************************************
 *
 * Function that performs specified operation involving source and
 * destination data pointers and puts result in destination array
 * data object using explicit dimension-generic looping constructs.
 *
 *************************************************************************
 */

template<class TYPE, class OP>
void ArrayDataOperationUtilities<TYPE, OP>::doArrayDataBufferOperationOnBox(
   const ArrayData<TYPE>& arraydata,
   const TYPE* buffer,
   const hier::Box& opbox,
   bool src_is_buffer,
   const OP& op)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(arraydata, opbox);
   TBOX_ASSERT(buffer != 0);
   TBOX_ASSERT(opbox.isSpatiallyEqual((opbox * arraydata.getBox())));

   const tbox::Dimension& dim(arraydata.getDim());

   TYPE * const dst_ptr =
      (src_is_buffer ? const_cast<TYPE *>(arraydata.getPointer())
       : const_cast<TYPE *>(buffer));
   const TYPE * const src_ptr =
      (src_is_buffer ? buffer : arraydata.getPointer());

   const hier::Box& array_d_box(arraydata.getBox());
   const unsigned int array_d_depth = arraydata.getDepth();

   int box_w[SAMRAI::MAX_DIM_VAL];
   int dat_w[SAMRAI::MAX_DIM_VAL];
   int dim_counter[SAMRAI::MAX_DIM_VAL];
   for (tbox::Dimension::dir_t i = 0; i < dim.getValue(); ++i) {
      box_w[i] = opbox.numberCells(i);
      dat_w[i] = array_d_box.numberCells(i);
      dim_counter[i] = 0;
   }

   const size_t dat_offset = arraydata.getOffset();
   const size_t buf_offset = box_w[0];

   /*
    * Data on the opbox can be decomposed into a set of
    * contiguous array sections representing data in a straight line
    * in the 0 coordinate direction.
    *
    * num_d0_blocks is the number of such array sections.
    * dat_begin, buf_begin are the array indices for the first
    * data items in each array section to be copied.
    */

   const int num_d0_blocks = static_cast<int>(opbox.size() / box_w[0]);

   size_t dat_begin = array_d_box.offset(opbox.lower());
   size_t buf_begin = 0;

   /*
    * Loop over the depth sections of the data arrays.
    */

   for (unsigned int d = 0; d < array_d_depth; ++d) {

#if defined(HAVE_RAJA)
     switch (dim.getValue()) {
       case 1: {
         tbox::ArrayView<1, TYPE> dest(src_is_buffer ? 
             tbox::ArrayView<1, TYPE>(const_cast<ArrayData<TYPE>&>(arraydata), d) : 
             tbox::ArrayView<1, TYPE>(const_cast<TYPE *>(buffer), arraydata.getBox()));

         // tbox::ArrayView<1, const TYPE> source(src_is_buffer ?
         tbox::ArrayView<1, const TYPE> source(src_is_buffer ?
             tbox::ArrayView<1, const TYPE>(buffer, arraydata.getBox(), d) :
             tbox::ArrayView<1, const TYPE>(arraydata, d));

         tbox::for_all1< tbox::policy::parallel >(opbox, [=] SAMRAI_DEVICE (int i) {
             op(dest(i), source(i));
         });
         } break;

       case 2: {
         tbox::ArrayView<2, TYPE> dest(src_is_buffer ?
             tbox::ArrayView<2, TYPE>(const_cast<ArrayData<TYPE>&>(arraydata), d) :
             tbox::ArrayView<2, TYPE>(const_cast<TYPE *>(buffer), arraydata.getBox(), d));

         //tbox::ArrayView<2, const TYPE> source(src_is_buffer ?
         tbox::ArrayView<2, const TYPE> source(src_is_buffer ?
             tbox::ArrayView<2, const TYPE>(buffer, arraydata.getBox(), d) :
             tbox::ArrayView<2, const TYPE>(arraydata, d));

         tbox::for_all2<tbox::policy::parallel>(opbox, [=] SAMRAI_DEVICE (int j, int i) {
             op(dest(i, j), source(i, j));
         });
         } break;

       case 3: {
         tbox::ArrayView<3, TYPE> dest(src_is_buffer ? 
             tbox::ArrayView<3, TYPE>(const_cast<ArrayData<TYPE>&>(arraydata), d) :
             tbox::ArrayView<3, TYPE>(const_cast<TYPE *>(buffer), arraydata.getBox(), d));

         // tbox::ArrayView<3, const TYPE> source(src_is_buffer ? 
         tbox::ArrayView<3, const TYPE> source(src_is_buffer ? 
             tbox::ArrayView<3, const TYPE>(buffer, arraydata.getBox(), d) :
             tbox::ArrayView<3, const TYPE>(arraydata, d));

         tbox::for_all3<tbox::policy::parallel>(opbox, [=] SAMRAI_DEVICE (int k, int j, int i) {
             op(dest(i, j, k), source(i, j, k));
         });
         } break;

       default: tbox::perr << "Dimension > 3 not supported with RAJA parallel kernels" << std::endl;
                TBOX_ERROR("Aborting in " __FILE__);
                break;
     }
#else
      size_t dat_counter = dat_begin;
      size_t buf_counter = buf_begin;

      size_t& dst_counter = (src_is_buffer ? dat_counter : buf_counter);
      size_t& src_counter = (src_is_buffer ? buf_counter : dat_counter);

      int dat_b[SAMRAI::MAX_DIM_VAL];
      for (tbox::Dimension::dir_t nd = 0; nd < dim.getValue(); ++nd) {
         dat_b[nd] = static_cast<int>(dat_counter);
      }

      /*
       * Loop over each contiguous block of data.
       */

      for (int nb = 0; nb < num_d0_blocks; ++nb) {

         for (int i0 = 0; i0 < box_w[0]; ++i0) {
            op(dst_ptr[dst_counter + i0], src_ptr[src_counter + i0]);
         }
         int dim_jump = 0;

         /*
          * After each contiguous block is packed, calculate the
          * beginning array index for the next block.
          */

         for (int j = 1; j < dim.getValue(); ++j) {
            if (dim_counter[j] < box_w[j] - 1) {
               ++dim_counter[j];
               dim_jump = j;
               break;
            } else {
               dim_counter[j] = 0;
            }
         }

         if (dim_jump > 0) {

            int dat_step = 1;
            for (int k = 0; k < dim_jump; ++k) {
               dat_step *= dat_w[k];
            }
            dat_counter = dat_b[dim_jump - 1] + dat_step;

            for (int m = 0; m < dim_jump; ++m) {
               dat_b[m] = static_cast<int>(dat_counter);
            }

         }  // if dim_jump > 0

         buf_counter += buf_offset;

      }  // nb loop over contiguous data blocks

      /*
       * After packing is complete on a full box for one depth index,
       * advance by the offset value.
       */

      dat_begin += dat_offset;
      buf_begin = buf_counter;
#endif

   }  // d loop over depth indices

}

}
}
#endif
