//
// File:	$URL$
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Templated array data looping operations supporting patch data types
//

#ifndef included_pdat_ArrayDataOperationUtilities
#define included_pdat_ArrayDataOperationUtilities

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_tbox_Box
#include "Box.h"
#endif
#ifndef included_hier_IntVector
#include "IntVector.h"
#endif

namespace SAMRAI {
    namespace pdat {

template<int DIM, class TYPE> class ArrayData;

/*!
 * @brief Struct ArrayDataOperationUtilities<DIM, TYPE, OP> provides
 * generic looping operations for all array-based patch data types.
 * The operations are templated on spatial dimension, data type,
 * and the operation that will be performed on individual array
 * elements in the innermost loop. 
 *
 * @see ArrayData
 */

template<int DIM, class TYPE, class OP> 
class ArrayDataOperationUtilities
{
public:
   /*!
    * Perform operation on a subset of data components of source and destination 
    * array data objects and put results in destination array data object.  
    * 
    * @param dst    Reference to destination array data object.
    * @param src    Const reference to source array data object.
    * @param opbox  Const reference to Box indicating index space region of operation.
    * @param src_shift  Const reference to IntVector indicating shift required to put
    *                   source index space region into destination index space region.
    * @param dst_start_depth  Integer specifying starting depth component of 
    *                         operation in destination array. 
    * @param src_start_depth  Integer specifying starting depth component of 
    *                         operation in source array. 
    * @param num_depth  Integer number of depth components on which to perform operation.
    * @param op  Const reference to object that performs operations on individual
    *            data array elements. 
    *
    * When assertion checking is active, assertion will result when any of the depth-related
    * arguments are out-of-bounds for the given array data objects.
    */
   static void doArrayDataOperationOnBox(ArrayData<DIM,TYPE>& dst,
                                         const ArrayData<DIM,TYPE>& src,
                                         const hier::Box<DIM>& opbox,
                                         const hier::IntVector<DIM>& src_shift,
                                         int dst_start_depth,
                                         int src_start_depth,
                                         int num_depth,
                                         const OP& op);

   /*!
    * Perform operation on all data components of array data object and corresponding
    * buffer data, putting results in either the array data object or buffer.
    *
    * @param arraydata   Const reference to array data object.
    * @param buffer      Const pointer to first element in buffer.
    * @param opbox       Const reference to Box indicating operation region in
    *                    index space of array data object.
    * @param src_is_buffer  Boolean value indicating whether buffer is source data
    *                       for operation; if true results will be placed in array
    *                       data object, otherwise results will go in buffer.
    * @param op  Const reference to object that performs operations on individual
    *            data array elements.
    *
    * When assertion checking is active, assertion will result when buffer pointer
    * is null or operation box is not the same as the array data box.
    */
   static void doArrayDataBufferOperationOnBox(const ArrayData<DIM,TYPE>& arraydata,
                                               const TYPE* buffer,
                                               const hier::Box<DIM>& opbox,
                                               bool src_is_buffer,
                                               const OP& op);
private:
   // the following are not implemented:
   ArrayDataOperationUtilities();
   ~ArrayDataOperationUtilities();
   ArrayDataOperationUtilities(const ArrayDataOperationUtilities<DIM,TYPE,OP>&);
   void operator=(const ArrayDataOperationUtilities<DIM,TYPE,OP>&);

};

}
}

#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "ArrayDataOperationUtilities.C"
#endif
