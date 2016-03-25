//
// File:	ArrayData.C
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Templated array data structure supporting patch data types
//

#ifndef included_tbox_ArrayData_C
#define included_tbox_ArrayData_C

#include "tbox/AbstractStream.h"
#include "BoxList.h"
#include "ArrayData.h"
#include "tbox/ArenaManager.h"
#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"

#define PDAT_ARRAYDATA_VERSION 1

#ifdef DEBUG_NO_INLINE
#include "ArrayData.I"
#endif

namespace SAMRAI {
    namespace pdat {

/*
*************************************************************************
*									*
* The default constructor creates an object that absolutely should      *
* not be used until it is initialized using initializeArray().          *
*									*
*************************************************************************
*/
template<int DIM, class TYPE>
ArrayData<DIM,TYPE>::ArrayData() :
   d_array(0, isStandardType()),
   d_box( hier::Index<DIM>(-1), hier::Index<DIM>(-2) ),
   d_depth(0),
   d_offset(0)
{
   return;
}

/*
*************************************************************************
*									*
* The main constructor allocates data for the given box and depth from	*
* the specified memory pool.  It does not initialize the memory.  The	*
* destructor automatically deallocates memory via the array destructor.	*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
ArrayData<DIM,TYPE>::ArrayData(
   const hier::Box<DIM>& box,
   const int depth,
   tbox::Pointer<tbox::Arena> pool)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(depth > 0);
#endif
   if (pool.isNull()) {
      pool = tbox::ArenaManager::getManager()->getStandardAllocator();
   }
   d_depth  = depth;
   d_offset = box.size();
   d_box    = box;
   d_array  = tbox::Array<TYPE>(d_depth * d_offset, pool, isStandardType());
#ifdef DEBUG_INITIALIZE_UNDEFINED
   undefineData();
#endif
}

template<int DIM, class TYPE>
ArrayData<DIM,TYPE>::~ArrayData()
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
ArrayData<DIM,TYPE>::ArrayData(const ArrayData<DIM,TYPE>& foo)
{
   NULL_USE(foo);
}

template<int DIM, class TYPE>
void ArrayData<DIM,TYPE>::operator=(const ArrayData<DIM,TYPE>& foo)
{
   NULL_USE(foo);
}

/*
*************************************************************************
*									*
* Initialize the array using the specified box, depth, and memory pool.	*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
void ArrayData<DIM,TYPE>::initializeArray(
   const hier::Box<DIM>& box,
   const int depth,
   const tbox::Pointer<tbox::Arena>& pool)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(depth > 0);
#endif
   d_depth  = depth;
   d_offset = box.size();
   d_box    = box;
   d_array  = tbox::Array<TYPE>(depth * d_offset, pool, isStandardType());
#ifdef DEBUG_INITIALIZE_UNDEFINED
   undefineData();
#endif
}

/*
*************************************************************************
*									*
* Copy data between two array data objects on a specified box domain.	*
* Don't use C++ indexing member functions, since compilers are probably	*
* too stupid to do strength reduction on the loops to get performance.	*
*									*
* If the source box, destination box, and copy box are the same and the	*
* source and destination have the same depth, then perform a fast copy	*
* of all data.								*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
void ArrayData<DIM,TYPE>::copy(
   const ArrayData<DIM,TYPE>& src,
   const hier::Box<DIM>& box)
{

   /*
    * Try to do a fast copy of data if all the data aligns perfectly
    */

   if ((d_depth == src.d_depth) && (d_box == src.d_box) && (box == d_box)) {
      TYPE *const dst_ptr = d_array.getPointer();
      const TYPE *const src_ptr = src.d_array.getPointer();
      const int n = d_offset * d_depth;
      for (int i = 0; i < n; i++) {
         dst_ptr[i] = src_ptr[i];
      }

   /*
    * Otherwise, do a painful copy using explicit looping constructs
    */

   } else {
      const hier::Box<DIM> copybox = box * d_box * src.d_box;
      if (!copybox.empty()) {
         TYPE *const dst_ptr = d_array.getPointer();
         const TYPE *const src_ptr = src.d_array.getPointer();
         const int depth = (d_depth < src.d_depth ? d_depth : src.d_depth);

         int box_w[DIM];
         int dst_w[DIM];
         int src_w[DIM];
         int dim_counter[DIM];
         for (int i = 0; i < DIM; i++) {
            box_w[i] = copybox.numberCells(i);
            dst_w[i] = d_box.numberCells(i);
            src_w[i] = src.d_box.numberCells(i);
            dim_counter[i] = 0;
         }

         const int dst_offset = d_offset;
         const int src_offset = src.d_offset;

         /*
          * Data on the copybox can be decomposed into a set of
          * contiguous arrays that represent data in a straight line
          * in the 0 direction.  num_d0_blocks is the number of such arrays.
          */
         const int num_d0_blocks = copybox.size() / box_w[0];

         /*
          * Find the array indices for the first item of data to be copied
          */
         int dst_begin = d_box.offset(copybox.lower());
         int src_begin = src.d_box.offset(copybox.lower());

         for (int d = 0; d < depth; d++) {

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
                  dst_ptr[dst_counter+i0] = src_ptr[src_counter+i0];
               }
               int dim_jump = 0;

               /*
                * After each contiguous block is copied, calculate the
                * beginning array index for the next block.
                */
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

            /*
             * After copy is complete on a full box for one depth index,
             * advance by the offset values.
             */
            dst_begin += dst_offset;
            src_begin += src_offset;

         }
      }
   }
}

/*
*************************************************************************
*									*
* Copy data between two ArrayData objects on different source and	*
* destination domains.  Don't use C++ indexing member functions, since	*
* compilers are probably too stupid to do strength reduction on the	*
* loops to get performance.						*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
void ArrayData<DIM,TYPE>::copy(
   const ArrayData<DIM,TYPE>& src,
   const hier::Box<DIM>& box,
   const hier::IntVector<DIM>& offset)
{
   const hier::Box<DIM> copybox = box * d_box * hier::Box<DIM>::shift(src.d_box, offset);
   if (!copybox.empty()) {
      TYPE *const dst_ptr = d_array.getPointer();
      const TYPE *const src_ptr = src.d_array.getPointer();

      const int depth = (d_depth < src.d_depth ? d_depth : src.d_depth);

      int box_w[DIM];
      int dst_w[DIM];
      int src_w[DIM];
      int dim_counter[DIM];
      for (int i = 0; i < DIM; i++) {
         box_w[i] = copybox.numberCells(i);
         dst_w[i] = d_box.numberCells(i);
         src_w[i] = src.d_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const int dst_offset = d_offset;
      const int src_offset = src.d_offset;

      const int num_d0_blocks = copybox.size() / box_w[0];

      int dst_begin = d_box.offset(copybox.lower());
      int src_begin = src.d_box.offset(copybox.lower()-offset);

      for (int d = 0; d < depth; d++) {

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
               dst_ptr[dst_counter+i0] = src_ptr[src_counter+i0];
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

/*
*************************************************************************
*									*
* Copy over the boxlist by calling the single-box copy for each box in	*
* the boxlist.								*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
void ArrayData<DIM,TYPE>::copy(
   const ArrayData<DIM,TYPE>& src,
   const hier::BoxList<DIM>& boxes,
   const hier::IntVector<DIM>& offset)
{
   for (typename hier::BoxList<DIM>::Iterator b(boxes); b; b++) {
      this->copy(src, b(), offset);
   }
}

/*
*************************************************************************
*									*
* Copy data between two array data objects on a specified box domain.	*
* Don't use C++ indexing member functions, since compilers are probably	*
* too stupid to do strength reduction on the loops to get performance.	*
*									*
* If the source box, destination box, and copy box are the same and the	*
* source and destination have the same depth, then perform a fast copy	*
* of all data.								*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
void ArrayData<DIM,TYPE>::copyDepth(int dst_depth,
				      const ArrayData<DIM,TYPE>& src,
				      int src_depth,
				      const hier::Box<DIM>& box)
{
   TYPE *const dst_ptr = d_array.getPointer();
   const TYPE *const src_ptr = src.d_array.getPointer();

   /*
    * Try to do a fast copy of data if all the data aligns perfectly
    */

   if ((d_box == src.d_box) && (box == d_box)) {
     TYPE *const dst_ptr_d = dst_ptr + dst_depth*d_offset;
     const TYPE *const src_ptr_d = src_ptr + src_depth*d_offset;
      for (int i = 0; i < d_offset; i++) {
         dst_ptr_d[i] = src_ptr_d[i];
      }

   /*
    * Otherwise, do a painful copy using explicit looping constructs
    */

   } else {
      const hier::Box<DIM> copybox = box * d_box * src.d_box;
      if (!copybox.empty()) {

         int box_w[DIM];
         int dst_w[DIM];
         int src_w[DIM];
         int dim_counter[DIM];
         for (int i = 0; i < DIM; i++) {
            box_w[i] = copybox.numberCells(i);
            dst_w[i] = d_box.numberCells(i);
            src_w[i] = src.d_box.numberCells(i);
            dim_counter[i] = 0;
         }

         const int dst_offset = d_offset;
         const int src_offset = src.d_offset;

         const int num_d0_blocks = copybox.size() / box_w[0];

         int dst_counter = d_box.offset(copybox.lower()) +
                         dst_depth*dst_offset;
         int src_counter = src.d_box.offset(copybox.lower()) +
                         src_depth*src_offset;

         int dst_b[DIM];
         int src_b[DIM];
         for (int nd = 0; nd < DIM; nd++) {
            dst_b[nd] = dst_counter;
            src_b[nd] = src_counter;
         }

         for (int nb = 0; nb < num_d0_blocks; nb++) {

            for (int i0 = 0; i0 < box_w[0]; i0++) {
               dst_ptr[dst_counter+i0] = src_ptr[src_counter+i0];
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
      }
   }
}

/*
*************************************************************************
*									*
* Pack data into the message stream.  Both packing routines add one	*
* level of copy into a temporary buffer to reduce the number of calls	*
* to the abstract stream packing routines.  These definitions will only	*
* work for the standard built-in types of bool, char, double, float,	*
* and int.								*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
void ArrayData<DIM,TYPE>::packStream(
   tbox::AbstractStream& stream,
   const hier::Box<DIM>& dest_box,
   const hier::IntVector<DIM>& source_offset) const
{
   const int size = d_depth * dest_box.size();
   tbox::Pointer<tbox::Arena> scratch =
      tbox::ArenaManager::getManager()->getScratchAllocator();
   tbox::Array<TYPE> buffer(size, scratch);

   packBuffer(buffer.getPointer(), hier::Box<DIM>::shift(dest_box, -source_offset));
   stream.pack(buffer.getPointer(), size);
}

template<int DIM, class TYPE>
void ArrayData<DIM,TYPE>::packStream(
   tbox::AbstractStream& stream,
   const hier::BoxList<DIM>& dest_boxes,
   const hier::IntVector<DIM>& source_offset) const
{
   const int size = d_depth * dest_boxes.getTotalSizeOfBoxes();
   tbox::Pointer<tbox::Arena> scratch =
      tbox::ArenaManager::getManager()->getScratchAllocator();
   tbox::Array<TYPE> buffer(size, scratch);

   int ptr = 0;
   for (typename hier::BoxList<DIM>::Iterator b(dest_boxes); b; b++) {
      packBuffer(buffer.getPointer(ptr), hier::Box<DIM>::shift(b(), -source_offset));
      ptr += d_depth * b().size();
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(ptr == size);
#endif
   stream.pack(buffer.getPointer(), size);
}

/*
*************************************************************************
*									*
* Unpack data from the message stream.  Both unpacking routines add one	*
* level of copy into a temporary buffer to reduce the number of calls	*
* to the abstract stream packing routines.  These definitions will only	*
* work for the standard built-in types of bool, char, double, float,	*
* and int.								*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
void ArrayData<DIM,TYPE>::unpackStream(
   tbox::AbstractStream& stream, 
   const hier::Box<DIM>& dest_box,
   const hier::IntVector<DIM>& source_offset)
{
   (void) source_offset;

   const int size = d_depth * dest_box.size();
   tbox::Pointer<tbox::Arena> scratch =
      tbox::ArenaManager::getManager()->getScratchAllocator();
   tbox::Array<TYPE> buffer(size, scratch);

   stream.unpack(buffer.getPointer(), size);
   unpackBuffer(buffer.getPointer(), dest_box);
}

template<int DIM, class TYPE>
void ArrayData<DIM,TYPE>::unpackStream(
   tbox::AbstractStream& stream,
   const hier::BoxList<DIM>& dest_boxes,
   const hier::IntVector<DIM>& source_offset)
{
   (void) source_offset;

   const int size = d_depth * dest_boxes.getTotalSizeOfBoxes();
   tbox::Pointer<tbox::Arena> scratch =
      tbox::ArenaManager::getManager()->getScratchAllocator();
   tbox::Array<TYPE> buffer(size, scratch);

   stream.unpack(buffer.getPointer(), size);

   int ptr = 0;
   for (typename hier::BoxList<DIM>::Iterator b(dest_boxes); b; b++) {
      unpackBuffer(buffer.getPointer(ptr), b());
      ptr += d_depth * b().size();
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(ptr == size);
#endif
}

/*
*************************************************************************
*									*
* Pack data on the specified box (for all components) into the buffer.	*
* Do not use C++ indexing because most compilers are too stupid to do a	*
* good job optimizing the loops.					*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
void ArrayData<DIM,TYPE>::packBuffer(
   TYPE *buffer, const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((box * d_box) == box);
#endif
   TYPE *const buf_ptr = buffer;
   const TYPE *const dat_ptr = d_array.getPointer();

   int box_w[DIM];
   int dat_w[DIM];
   int dim_counter[DIM];
   for (int i = 0; i < DIM; i++) {
      box_w[i] = box.numberCells(i);
      dat_w[i] = d_box.numberCells(i);
      dim_counter[i] = 0;
   }

   const int dat_offset = d_offset;

   /*
    * Data on the box can be decomposed into a set of
    * contiguous arrays that represent data in a straight line
    * in the 0 direction.  num_d0_blocks is the number of such arrays.
    */
   const int num_d0_blocks = box.size() / box_w[0];

   /*
    * Find the array index for the first item of data to be packed
    */
   int dat_begin = d_box.offset(box.lower());
   int buf_begin = 0;

   for (int d = 0; d < d_depth; d++) {

      int dat_counter = dat_begin;
      int buf_counter = buf_begin;

      int dat_b[DIM];
      for (int nd = 0; nd < DIM; nd++) {
         dat_b[nd] = dat_counter;
      }

      /*
       * Loop over each contiguous block of data.
       */
      for (int nb = 0; nb < num_d0_blocks; nb++) {

         for (int i0 = 0; i0 < box_w[0]; i0++) {
            buf_ptr[buf_counter+i0] = dat_ptr[dat_counter+i0];
         }

         /*
          * After each contiguous block is packed, calculate the
          * beginning array index for the next block.
          */
         buf_counter += box_w[0];

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
            int dat_step = 1;
            for (int k = 0; k < dim_jump; k++) {
               dat_step *= dat_w[k];
            }
            dat_counter = dat_b[dim_jump-1] + dat_step;

            for (int m = 0; m < dim_jump; m++) {
               dat_b[m] = dat_counter;
            }
         }
      }

      /*
       * After packing is complete on a full box for one depth index,
       * advance by the offset value.
       */
      dat_begin += dat_offset;
      buf_begin = buf_counter;
   }
}

/*
*************************************************************************
*									*
* Unpack data from the specified box (for all components) into the	*
* buffer.  Do not use C++ indexing because most compilers are too	*
* stupid to do a good job optimizing the loops.				*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
void ArrayData<DIM,TYPE>::unpackBuffer(
   const TYPE *buffer, const hier::Box<DIM>& box)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((box * d_box) == box);
#endif
   const TYPE *const buf_ptr = buffer;
   TYPE *const dat_ptr = d_array.getPointer();

   int box_w[DIM];
   int dat_w[DIM];
   int dim_counter[DIM];
   for (int i = 0; i < DIM; i++) {
      box_w[i] = box.numberCells(i);
      dat_w[i] = d_box.numberCells(i);
      dim_counter[i] = 0;
   }

   const int dat_offset = d_offset;

   const int num_d0_blocks = box.size() / box_w[0];

   int dat_begin = d_box.offset(box.lower());
   int buf_begin = 0;

   for (int d = 0; d < d_depth; d++) {

      int dat_counter = dat_begin;
      int buf_counter = buf_begin;

      int dat_b[DIM];
      for (int nd = 0; nd < DIM; nd++) {
         dat_b[nd] = dat_counter;
      }

      for (int nb = 0; nb < num_d0_blocks; nb++) {

         for (int i0 = 0; i0 < box_w[0]; i0++) {
            dat_ptr[dat_counter+i0] = buf_ptr[buf_counter+i0];
         }

         buf_counter += box_w[0];

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
            int dat_step = 1;
            for (int k = 0; k < dim_jump; k++) {
               dat_step *= dat_w[k];
            }
            dat_counter = dat_b[dim_jump-1] + dat_step;

            for (int m = 0; m < dim_jump; m++) {
               dat_b[m] = dat_counter;
            }

         }
      }

      dat_begin += dat_offset;
      buf_begin = buf_counter;
   }
}

/*
*************************************************************************
*									*
* Fill all or portions of the array with the specified data value.	*
* The templated TYPE must define the assignment operator.		*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
void ArrayData<DIM,TYPE>::fillAll(const TYPE& t)
{
   if ( ! d_box.empty() ) {
      TYPE *ptr = d_array.getPointer();
      const int n = d_depth * d_offset;
      for (int i = 0; i < n; i++) {
         ptr[i] = t;
      }
   }
}

template<int DIM, class TYPE>
void ArrayData<DIM,TYPE>::fillAll(const TYPE& t, const hier::Box<DIM>& box)
{
   for (int d = 0; d < d_depth; d++) {
      fill(t, box, d);
   }
}

template<int DIM, class TYPE>
void ArrayData<DIM,TYPE>::fill(const TYPE& t, const int d)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((d >= 0) && (d < d_depth));
#endif
   if ( ! d_box.empty() ) {
      TYPE *ptr = d_array.getPointer(d * d_offset);
      const int n = d_offset;
      for (int i = 0; i < n; i++) {
         ptr[i] = t;
      }
   }
}

template<int DIM, class TYPE>
void ArrayData<DIM,TYPE>::fill(
   const TYPE& t, const hier::Box<DIM>& box, const int d)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((d >= 0) && (d < d_depth));
#endif
   const hier::Box<DIM> ispace = d_box * box;

   if ( ! ispace.empty() ) {

      int box_w[DIM];
      int dst_w[DIM];
      int dim_counter[DIM];
      for (int i = 0; i < DIM; i++) {
         box_w[i] = ispace.numberCells(i);
         dst_w[i] = d_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const int num_d0_blocks = ispace.size() / box_w[0];

      int dst_counter = d_box.offset(ispace.lower()) + d*d_offset;

      int dst_b[DIM];
      for (int nd = 0; nd < DIM; nd++) {
         dst_b[nd] = dst_counter;
      }

      TYPE *const dst_ptr = d_array.getPointer();

      for (int nb = 0; nb < num_d0_blocks; nb++) {

         for (int i0 = 0; i0 < box_w[0]; i0++) {
            dst_ptr[dst_counter+i0] = t;
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
            for (int k = 0; k < dim_jump; k++) {
               dst_step *= dst_w[k];
            }
            dst_counter = dst_b[dim_jump-1] + dst_step;
            
            for (int m = 0; m < dim_jump; m++) {
               dst_b[m] = dst_counter;
            }
         }
      }
   }
}

/*
*************************************************************************
*									*
* Checks to make sure that class and restart file version numbers are   *
* equal.  If so, reads in d_depth, d_offset, and d_box from the 	*
* database.  Then calls getSpecializedFromDatabase() to read in the	*
* actual data.								*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
void ArrayData<DIM,TYPE>::getFromDatabase(
     tbox::Pointer<tbox::Database> database)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!database.isNull());
#endif

   int ver =  database->getInteger("PDAT_ARRAYDATA_VERSION");
   if (ver != PDAT_ARRAYDATA_VERSION) {
      TBOX_ERROR("ArrayData<DIM>::getFromDatabase error...\n"
          << " : Restart file version different than class version" << endl);
   }

   d_depth = database->getInteger("d_depth");
   d_offset = database->getInteger("d_offset");
   d_box = database->getDatabaseBox("d_box");

   getSpecializedFromDatabase(database);
}

/*
*************************************************************************
*									*
* Writes out the class version number, d_depth, d_offset, and d_box	*
* to the database.  Then calls putSpecializedToDatabase() to write	*
* in the actual data.  							*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
void ArrayData<DIM,TYPE>::putToDatabase(
   tbox::Pointer<tbox::Database> database,
   bool data_only)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!database.isNull());
#endif
  
   if (!data_only) { 
      database->putInteger("PDAT_ARRAYDATA_VERSION",PDAT_ARRAYDATA_VERSION);

      database->putInteger("d_depth",d_depth);
      database->putInteger("d_offset",d_offset);
      database->putDatabaseBox("d_box",d_box);
   }

   putSpecializedToDatabase(database);
}


template<int DIM, class TYPE>
void ArrayData<DIM,TYPE>::putSpecializedToDatabase(
                 tbox::Pointer<tbox::Database> database)
{
   database->putArray("d_array",d_array);
}

template<int DIM, class TYPE> 
void ArrayData<DIM,TYPE>::getSpecializedFromDatabase(
                 tbox::Pointer<tbox::Database> database)
{
   database->getArray("d_array", d_array);
}

template<int DIM, class TYPE>
void ArrayData<DIM,TYPE>::undefineData()
{
   fillAll(  tbox::MathUtilities<TYPE>::getUndefined() );
}

}
}

#endif
