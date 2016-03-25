//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/array/ArrayData.C $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1777 $
// Modified:	$LastChangedDate: 2007-12-13 16:51:06 -0800 (Thu, 13 Dec 2007) $
// Description:	Templated array data structure supporting patch data types
//

#ifndef included_pdat_ArrayData_C
#define included_pdat_ArrayData_C

#include "tbox/AbstractStream.h"
#include "tbox/ArenaManager.h"
#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"
#include "BoxList.h"
#include "ArrayData.h"
#include "ArrayDataOperationUtilities.h"
#include "CopyOperation.h"
#include "SumOperation.h"


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
   int depth,
   tbox::Pointer<tbox::Arena> pool)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(depth > 0);
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
   int depth,
   const tbox::Pointer<tbox::Arena> pool)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(depth > 0);
#endif
   tbox::Pointer<tbox::Arena> mem_pool = pool;
   if (mem_pool.isNull()) {
      mem_pool = tbox::ArenaManager::getManager()->getStandardAllocator();
   }
   d_depth  = depth;
   d_offset = box.size();
   d_box    = box;
   d_array  = tbox::Array<TYPE>(depth * d_offset, mem_pool, isStandardType());
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
void ArrayData<DIM,TYPE>::copy(const ArrayData<DIM,TYPE>& src,
                               const hier::Box<DIM>& box)
{
   
   CopyOperation<TYPE> copyop;

   /*
    * Do a fast copy of data if all data aligns with copy region
    */

   if ( (d_depth == src.d_depth) && 
        (d_box == src.d_box) && 
        (box == d_box)) {

      TYPE* const dst_ptr = d_array.getPointer();
      const TYPE* const src_ptr = src.d_array.getPointer();
      const int n = d_offset * d_depth;
      for (int i = 0; i < n; i++) {
         copyop(dst_ptr[i], src_ptr[i]);
      }

   } else {

      const hier::Box<DIM> copybox = box * d_box * src.d_box;

      if (!copybox.empty()) {

         const int dst_start_depth = 0;
         const int src_start_depth = 0;
         const int num_depth = (d_depth < src.d_depth ? d_depth : src.d_depth);
         const hier::IntVector<DIM> src_shift(0);

         ArrayDataOperationUtilities< DIM, TYPE, CopyOperation<TYPE> >::
            doArrayDataOperationOnBox(*this,
                                      src,
                                      copybox,
                                      src_shift,
                                      dst_start_depth,
                                      src_start_depth,
                                      num_depth,
                                      copyop);

      } 

   }

}

/*
*************************************************************************
*									*
* Copy data from source ArrayData object to this (destination)          *
* ArrayData object on given box domain.                                 *
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
void ArrayData<DIM,TYPE>::copy(const ArrayData<DIM,TYPE>& src,
                               const hier::Box<DIM>& box,
                               const hier::IntVector<DIM>& src_shift)
{

   if ( src_shift == hier::IntVector<DIM>(0) ) {

      copy(src, box);

   } else {

      const hier::Box<DIM> copybox = 
         box * d_box * hier::Box<DIM>::shift(src.d_box, src_shift);

      if (!copybox.empty()) {

         const int dst_start_depth = 0;
         const int src_start_depth = 0;
         const int num_depth = (d_depth < src.d_depth ? d_depth : src.d_depth);

         CopyOperation<TYPE> copyop;

         ArrayDataOperationUtilities< DIM, TYPE, CopyOperation<TYPE> >::
            doArrayDataOperationOnBox(*this,
                                      src,
                                      copybox,
                                      src_shift,
                                      dst_start_depth,
                                      src_start_depth,
                                      num_depth,
                                      copyop);

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
   const hier::IntVector<DIM>& src_shift)
{
   for (typename hier::BoxList<DIM>::Iterator b(boxes); b; b++) {
      this->copy(src, b(), src_shift);
   }
}

/*
*************************************************************************
*									*
* Copy data between two array data objects on a specified box domain.	*
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
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( (0 <= dst_depth) && (dst_depth <= d_depth) );
   TBOX_ASSERT( (0 <= src_depth) && (src_depth <= src.d_depth) );
#endif
   
   CopyOperation<TYPE> copyop;

   /*
    * Do a fast copy of data if all data aligns with copy region
    */

   if ( (d_box == src.d_box) && (box == d_box) ) {

     TYPE* const dst_ptr = d_array.getPointer();
     const TYPE* const src_ptr = src.d_array.getPointer();

     TYPE* const dst_ptr_d = dst_ptr + dst_depth*d_offset;
     const TYPE* const src_ptr_d = src_ptr + src_depth*d_offset;
      for (int i = 0; i < d_offset; i++) {
         copyop(dst_ptr_d[i], src_ptr_d[i]);
      }

   } else {

      const hier::Box<DIM> copybox = box * d_box * src.d_box;

      if (!copybox.empty()) {

         const int dst_start_depth = dst_depth;
         const int src_start_depth = src_depth;
         const int num_depth = 1;
         const hier::IntVector<DIM> src_shift(0);

         ArrayDataOperationUtilities< DIM, TYPE, CopyOperation<TYPE> >::
            doArrayDataOperationOnBox(*this,
                                      src,
                                      copybox,
                                      src_shift,
                                      dst_start_depth,
                                      src_start_depth,
                                      num_depth,
                                      copyop); 

      }

   }

}

/*
*************************************************************************
*									*
* Add data from source ArrayData object to this (destination)           *
* ArrayData object on given box region.                                 *
*									*
* If the source box, destination box, and copy box are the same and the	*
* source and destination have the same depth, then perform a fast sum	*
* on all data rather than performing explicit looping operations.       *
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
void ArrayData<DIM,TYPE>::sum(const ArrayData<DIM,TYPE>& src,
                              const hier::Box<DIM>& box)
{

   SumOperation<TYPE> sumop;

   /*
    * Do a fast copy and add if all data aligns with copy region
    */

   if ( (d_depth == src.d_depth) && 
        (d_box == src.d_box) && 
        (box == d_box)) {

      TYPE* const dst_ptr = d_array.getPointer();
      const TYPE* const src_ptr = src.d_array.getPointer();
      const int n = d_offset * d_depth;
      for (int i = 0; i < n; i++) {
         sumop(dst_ptr[i], src_ptr[i]);
      }

   } else {

      const hier::Box<DIM> copybox = box * d_box * src.d_box;

      if (!copybox.empty()) {

         const int dst_start_depth = 0;
         const int src_start_depth = 0;
         const int num_depth = (d_depth < src.d_depth ? d_depth : src.d_depth);
         const hier::IntVector<DIM> src_shift(0);

         ArrayDataOperationUtilities< DIM, TYPE, SumOperation<TYPE> >::
            doArrayDataOperationOnBox(*this,
                                      src,
                                      copybox,
                                      src_shift,
                                      dst_start_depth,
                                      src_start_depth,
                                      num_depth,
                                      sumop);

      } 

   }

}

/*
*************************************************************************
*									*
* Add data from source ArrayData object to this (destination)           *
* ArrayData object on region described by given box and offset.         *
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
void ArrayData<DIM,TYPE>::sum(const ArrayData<DIM,TYPE>& src,
                              const hier::Box<DIM>& box,
                              const hier::IntVector<DIM>& src_shift)
{

   if ( src_shift == hier::IntVector<DIM>(0) ) {

      sum(src, box);

   } else {

      const hier::Box<DIM> copybox = 
         box * d_box * hier::Box<DIM>::shift(src.d_box, src_shift);

      if (!copybox.empty()) {

         const int dst_start_depth = 0;
         const int src_start_depth = 0;
         const int num_depth = (d_depth < src.d_depth ? d_depth : src.d_depth);

         SumOperation<TYPE> sumop;

         ArrayDataOperationUtilities< DIM, TYPE, SumOperation<TYPE> >::
            doArrayDataOperationOnBox(*this,
                                      src,
                                      copybox,
                                      src_shift,
                                      dst_start_depth,
                                      src_start_depth,
                                      num_depth,
                                      sumop);

      }

   }

}

/*
*************************************************************************
*									*
* Add data from source ArrayData object to this (destination)           *
* ArrayData object on regions described by given boxes and offset.      *
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
void ArrayData<DIM,TYPE>::sum(
   const ArrayData<DIM,TYPE>& src,
   const hier::BoxList<DIM>& boxes,
   const hier::IntVector<DIM>& src_shift)
{
   for (typename hier::BoxList<DIM>::Iterator b(boxes); b; b++) {
      this->sum(src, b(), src_shift);
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
   const hier::IntVector<DIM>& src_shift) const
{

   const int size = d_depth * dest_box.size();
   tbox::Pointer<tbox::Arena> scratch =
      tbox::ArenaManager::getManager()->getScratchAllocator();
   tbox::Array<TYPE> buffer(size, scratch);

   packBuffer(buffer.getPointer(), 
              hier::Box<DIM>::shift(dest_box, -src_shift));

   stream.pack(buffer.getPointer(), size);

}

template<int DIM, class TYPE>
void ArrayData<DIM,TYPE>::packStream(
   tbox::AbstractStream& stream,
   const hier::BoxList<DIM>& dest_boxes,
   const hier::IntVector<DIM>& src_shift) const
{

   const int size = d_depth * dest_boxes.getTotalSizeOfBoxes();
   tbox::Pointer<tbox::Arena> scratch =
      tbox::ArenaManager::getManager()->getScratchAllocator();
   tbox::Array<TYPE> buffer(size, scratch);

   int ptr = 0;
   for (typename hier::BoxList<DIM>::Iterator b(dest_boxes); b; b++) {
      packBuffer(buffer.getPointer(ptr), 
                 hier::Box<DIM>::shift(b(), -src_shift));
      ptr += d_depth * b().size();
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(ptr == size);
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
   const hier::IntVector<DIM>& src_shift)
{

   (void) src_shift;

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
   const hier::IntVector<DIM>& src_shift)
{

   (void) src_shift;

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
   TBOX_ASSERT(ptr == size);
#endif

}

/*
*************************************************************************
*									*
* Unpack data from the message stream and add to this array data object.*
* Both unpacking routines add one level of copy into a temporary buffer *
* to reduce the number of calls	to the abstract stream packing routines.*
* These definitions will only work for the standard built-in types of   *
* bool, char, double, float, and int.					*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
void ArrayData<DIM,TYPE>::unpackStreamAndSum(
   tbox::AbstractStream& stream, 
   const hier::Box<DIM>& dest_box,
   const hier::IntVector<DIM>& src_shift)
{

   (void) src_shift;

   const int size = d_depth * dest_box.size();
   tbox::Pointer<tbox::Arena> scratch =
      tbox::ArenaManager::getManager()->getScratchAllocator();
   tbox::Array<TYPE> buffer(size, scratch);

   stream.unpack(buffer.getPointer(), size);
   unpackBufferAndSum(buffer.getPointer(), dest_box);

}

template<int DIM, class TYPE>
void ArrayData<DIM,TYPE>::unpackStreamAndSum(
   tbox::AbstractStream& stream,
   const hier::BoxList<DIM>& dest_boxes,
   const hier::IntVector<DIM>& src_shift)
{

   (void) src_shift;

   const int size = d_depth * dest_boxes.getTotalSizeOfBoxes();
   tbox::Pointer<tbox::Arena> scratch =
      tbox::ArenaManager::getManager()->getScratchAllocator();
   tbox::Array<TYPE> buffer(size, scratch);

   stream.unpack(buffer.getPointer(), size);

   int ptr = 0;
   for (typename hier::BoxList<DIM>::Iterator b(dest_boxes); b; b++) {
      unpackBufferAndSum(buffer.getPointer(ptr), b());
      ptr += d_depth * b().size();
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(ptr == size);
#endif

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
      TYPE* ptr = d_array.getPointer();
      const int n = d_depth * d_offset;
      for (int i = 0; i < n; i++) {
         ptr[i] = t;
      }
   }
}

template<int DIM, class TYPE>
void ArrayData<DIM,TYPE>::fillAll(const TYPE& t, 
                                  const hier::Box<DIM>& box)
{
   for (int d = 0; d < d_depth; d++) {
      fill(t, box, d);
   }
}

template<int DIM, class TYPE>
void ArrayData<DIM,TYPE>::fill(const TYPE& t, 
                               const int d)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((d >= 0) && (d < d_depth));
#endif
   if ( ! d_box.empty() ) {
      TYPE* ptr = d_array.getPointer(d * d_offset);
      const int n = d_offset;
      for (int i = 0; i < n; i++) {
         ptr[i] = t;
      }
   }
}

template<int DIM, class TYPE>
void ArrayData<DIM,TYPE>::fill(const TYPE& t, 
                               const hier::Box<DIM>& box, 
                               const int d)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((d >= 0) && (d < d_depth));
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

      TYPE* const dst_ptr = d_array.getPointer();

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
   TBOX_ASSERT(!database.isNull());
#endif

   int ver =  database->getInteger("PDAT_ARRAYDATA_VERSION");
   if (ver != PDAT_ARRAYDATA_VERSION) {
      TBOX_ERROR("ArrayData<DIM>::getFromDatabase error...\n"
          << " : Restart file version different than class version" << std::endl);
   }

   d_depth = database->getInteger("d_depth");
   d_offset = database->getInteger("d_offset");
   d_box = database->getDatabaseBox("d_box");

   getSpecializedFromDatabase(database);
}

/*
*************************************************************************
*									*
* Write out the class version number, d_depth, d_offset, and d_box	*
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
   TBOX_ASSERT(!database.isNull());
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
   database->putArray("d_array", d_array);
}

template<int DIM, class TYPE> 
void ArrayData<DIM,TYPE>::getSpecializedFromDatabase(
                 tbox::Pointer<tbox::Database> database)
{
   database->getArray("d_array", d_array);
}

/*
*************************************************************************
*									*
* Set all array data to undefined values.                               * 
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
void ArrayData<DIM,TYPE>::undefineData()
{
   fillAll(  tbox::MathUtilities<TYPE>::getSignalingNaN() );
}

/*
*************************************************************************
*									*
* Private member functions to pack and unpack data on the specified box *
* (for all components) into/from the buffer.	                        *
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
void ArrayData<DIM,TYPE>::packBuffer(TYPE* buffer, 
                                     const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((box * d_box) == box);
#endif

   bool src_is_buffer = false;

   CopyOperation<TYPE> copyop;

   ArrayDataOperationUtilities< DIM, TYPE, CopyOperation<TYPE> >::
      doArrayDataBufferOperationOnBox(*this,
                                      buffer,
                                      box,
                                      src_is_buffer,
                                      copyop);

}

template<int DIM, class TYPE>
void ArrayData<DIM,TYPE>::unpackBuffer(const TYPE* buffer,
                                       const hier::Box<DIM>& box)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((box * d_box) == box);
#endif

   bool src_is_buffer = true;

   CopyOperation<TYPE> copyop;

   ArrayDataOperationUtilities< DIM, TYPE, CopyOperation<TYPE> >::
      doArrayDataBufferOperationOnBox(*this,
                                      buffer,
                                      box,
                                      src_is_buffer,
                                      copyop);

}

/*
*************************************************************************
*									*
* Private member function to unpack data on the specified box           *
* (all components) from the buffer and add to this array data object.	*
*									*
*************************************************************************
*/

template<int DIM, class TYPE>
void ArrayData<DIM,TYPE>::unpackBufferAndSum(const TYPE* buffer,
                                             const hier::Box<DIM>& box)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((box * d_box) == box);
#endif

   bool src_is_buffer = true;

   SumOperation<TYPE> sumop;

   ArrayDataOperationUtilities< DIM, TYPE, SumOperation<TYPE> >::
      doArrayDataBufferOperationOnBox(*this,
                                      buffer,
                                      box,
                                      src_is_buffer,
                                      sumop);

}


}
}

#endif
