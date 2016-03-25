//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/array/ArrayData.h $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2172 $
// Modified:	$LastChangedDate: 2008-05-02 11:02:08 -0700 (Fri, 02 May 2008) $
// Description:	Templated array data structure supporting patch data types
//

#ifndef included_pdat_ArrayData
#define included_pdat_ArrayData

#include "SAMRAI_config.h"

#include <typeinfo>

#include "Box.h"
#include "BoxList.h"
#include "Index.h"
#include "IntVector.h"
#include "ArrayDataIterator.h"
#include "tbox/AbstractStream.h"
#include "tbox/Arena.h"
#include "tbox/Array.h"
#include "tbox/Complex.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

#ifndef NULL
#define NULL (0)
#endif

namespace SAMRAI {
    namespace pdat {

/*!
 * @brief Class ArrayData<DIM, TYPE> is a basic templated array structure defined
 * over the index space of a box (with a specified depth) that provides
 * the support for the various standard array-based patch data subclasses.
 *
 * The data storage is in (i,...,k,d) order, where i,...,k indicates
 * spatial indices and the d indicates the component at that location.
 * Memory allocation is in column-major ordering (e.g., Fortran style)
 * so that the leftmost index runs fastest in memory.
 *
 * The data type TYPE must define a default constructor (that takes no
 * arguments) and also the assignment operator.  Note that a number of
 * functions only work for standard built-in types (bool, char, double,
 * float, and int).  To use this class with other user-defined types,
 * many of these functions will need to be specialized, especially those
 * that deal with message packing and unpacking.
 */

template<int DIM, class TYPE>
class ArrayData
{
public:

   /*!
    * Static member function that returns tru when the amount of buffer space in a 
    * message stream can be estimated from box only.  For built-in types (bool, char, 
    * double, float, int, and dcomplex), this routine returns true.  For other 
    * data types (template paramters) that may require special handling, 
    * a different implementation must be provided.
    */
   static bool canEstimateStreamSizeFromBox();

   /*!
    * Static member function that returns the amount of memory space needed to 
    * store data of given depth on a box.  
    *
    * Note that this function is only defined for the standard data types: 
    * bool, char, double, float, int, and dcomplex.  It must be provided for other
    * template parameter types.
    *
    * @return size_t value indicating the amount of memory space needed for the data.
    *
    * @param box   Const reference to box object describing the spatial extents
    *              of the array data index region of interest.
    * @param depth Integer number of data values at each spatial location in
    *              the array.
    */
   static size_t getSizeOfData(const hier::Box<DIM>& box, 
                               int depth);

   /*!
    * The default constructor creates an empty array data object.
    * The initializeArray() member function must be called before the
    * array can be used.
    */
   ArrayData();

   /*!
    * Construct an array data object.  
    * 
    * @param box   Const reference to box object describing the spatial extents 
    *              of the index space associated with the array data object.
    * @param depth Integer number of data values at each spatial location in 
    *              the array.  
    * @param pool  Optional pointer to memory pool for allocating array data storage.
    *              Default is null indicating that standard allocator will be used.
    */
   ArrayData(
      const hier::Box<DIM>& box,
      int depth,
      tbox::Pointer<tbox::Arena> pool = tbox::Pointer<tbox::Arena>(NULL));

   /*!
    * The destructor for an array data object releases all memory allocated
    * for the array elements.
    */
   ~ArrayData();

   /*!
    * Initialize the array data.  This routine is normally called to allocate
    * storage when the array data object is created using the default constructor.
    * 
    * @param box   Const reference to box object describing the spatial extents
    *              of the index space associated with the array data object.
    * @param depth Integer number of data values at each spatial location in
    *              the array.
    * @param pool  Optional pointer to memory pool for allocating array data storage.
    *              Default is null indicating that standard allocator will be used.
    */
   void initializeArray(
      const hier::Box<DIM>& box,
      int depth,
      const tbox::Pointer<tbox::Arena> pool = tbox::Pointer<tbox::Arena>(NULL));

   /*!
    * @brief Returns true when the array has been properly initialized
    * and storage has been allocated; otherwise, return false.
    *
    * Note: Only arrays that have been initialized can do anything useful. 
    * Initialize an uninitialized array by calling the initializeArray() method.
    */
   bool isInitialized() const;

   /*!
    * Set the array data to an ``undefined'' state appropriate for the data type.  
    * For example, for float and double, this means setting data to signaling NaNs 
    * that cause a floating point exception when used in a numerical expression
    * without being set to valid values.
    */
   void undefineData();

   /*!
    * Return the box over which the array is defined.
    */
   const hier::Box<DIM>& getBox() const;

   /*!
    * Return the depth (e.g., the number of data values at each spatial
    * location) of this array.
    */
   int getDepth() const;

   /*!
    * Return the offset (e.g., the number of data values for each 
    * depth component) of this array.
    */
   int getOffset() const;

   /*!
    * Get a non-const pointer to the beginning of the given depth 
    * component of this data array. 
    */
   TYPE* getPointer(const int d = 0);

   /*!
    * Get a const pointer to the beginning of the given depth 
    * component of this data array.
    */
   const TYPE* getPointer(const int d = 0) const;

   /*!
    * Return reference to value in this array associated with the given
    * box index and depth component.
    */
   TYPE& operator()(const hier::Index<DIM>& i, 
                    const int d);

   /*!
    * Return const reference to value in this array associated with the given
    * box index and depth component.
    */
   const TYPE& operator()(const hier::Index<DIM>& i, 
                          const int d) const;

   /*!
    * Copy data from the source array data object to this array data object 
    * on the specified index space region.
    *
    * Note that this routine assumes that the source and destination
    * box regions require no shifting to make them consistent.  This routine 
    * will intersect the specified box with the source and destination boxes 
    * to find the region of intersection.
    * 
    * @param src   Const reference to source array data object.
    * @param box   Const reference to box object describing the spatial extents
    *              of the index space region over which to perform the copy operation.
    *              Note: the box is in either the source or destination index space
    *                    (which are assumed to be the same).
    */
   void copy(const ArrayData<DIM,TYPE>& src, 
             const hier::Box<DIM>& box);

   /*!
    * Copy data from the source array data object to this array data object
    * on the specified index space region.  
    *
    * Note that this routine assumes that the source array box region must
    * be shifted to be consistent with the destination (this) array box region.
    * This routine will intersect the specified box with the destination box and
    * shifted source box to find the region of intersection.
    *
    * @param src   Const reference to source array data object.
    * @param box   Const reference to box object describing the spatial extents
    *              of the index space region over which to perform the copy operation.
    *              Note: the box is in the destination index space.
    * @param src_shift Const reference to shift vector used to put the source
    *              array data box into the index space region of this array data object.
    */
   void copy(const ArrayData<DIM,TYPE>& src,
             const hier::Box<DIM>& box,
             const hier::IntVector<DIM>& src_shift);

   /*!
    * Copy data from the source array data object to this array data object
    * on the specified index space regions.
    *
    * Note that this routine assumes that the source array box region must
    * be shifted to be consistent with the destination (this) array box region.
    * This routine will intersect the specified boxes with the destination box and
    * shifted source box to find the regions of intersection.
    *
    * @param src   Const reference to source array data object.
    * @param boxes Const reference to box list describing the spatial extents
    *              of the index space regions over which to perform the copy operation.
    *              Note: the boxes are in the destination index space.
    * @param src_shift Const reference to shift vector used to put the source
    *              array data box into the index space region of this array data object.
    */
   void copy(const ArrayData<DIM,TYPE>& src,
             const hier::BoxList<DIM>& boxes,
             const hier::IntVector<DIM>& src_shift);

   /*!
    * Copy given source depth of source array data object to given destination 
    * depth of this array data object on the specified index space region.
    *
    * Note that this routine assumes that the source and destination
    * box regions require no shifting to make them consistent.  This routine
    * will intersect the specified box with the source and destination boxes
    * to find the region of intersection.
    *
    * @param dst_depth Integer depth of destination array.
    * @param src       Const reference to source array data object.
    * @param src_depth Integer depth of source array.
    * @param box       Const reference to box object describing the spatial extents
    *                  of the index space region over which to perform the copy operation.
    *                  Note: the box is in either the source or destination index space
    *                        (which are assumed to be the same).
    */
   void copyDepth(int dst_depth,
		  const ArrayData<DIM,TYPE>& src,
		  int src_depth,
		  const hier::Box<DIM>& box);

   /*!
    * Add data from the source array data object to this array data object
    * on the specified index space region.
    *
    * Note that this routine assumes that the source and destination
    * box regions require no shifting to make them consistent.  This routine
    * will intersect the specified box with the source and destination boxes
    * to find the region of intersection.
    *
    * @param src   Const reference to source array data object.
    * @param box   Const reference to box object describing the spatial extents
    *              of the index space region over which to perform the sum operation.
    *              Note: the box is in either the source or destination index space
    *                    (which are assumed to be the same).
    */
   void sum(const ArrayData<DIM,TYPE>& src,
             const hier::Box<DIM>& box);

   /*!
    * Add data from the source array data object to this array data object
    * on the specified index space region.
    *
    * Note that this routine assumes that the source array box region must
    * be shifted to be consistent with the destination (this) array box region.
    * This routine will intersect the specified box with the destination box and
    * shifted source box to find the region of intersection.
    *
    * @param src   Const reference to source array data object.
    * @param box   Const reference to box object describing the spatial extents
    *              of the index space region over which to perform the sum operation.
    *              Note: the box is in the destination index space.
    * @param src_shift Const reference to shift vector used to put the source
    *              array data box into the index space region of this array data object.
    */
   void sum(const ArrayData<DIM,TYPE>& src,
             const hier::Box<DIM>& box,
             const hier::IntVector<DIM>& src_shift);

   /*!
    * Add data from the source array data object to this array data object
    * on the specified index space regions.
    *
    * Note that this routine assumes that the source array box region must
    * be shifted to be consistent with the destination (this) array box region.
    * This routine will intersect the specified boxes with the destination box and
    * shifted source box to find the regions of intersection.
    *
    * @param src   Const reference to source array data object.
    * @param boxes Const reference to box list describing the spatial extents
    *              of the index space regions over which to perform the sum operation.
    *              Note: the boxes are in the destination index space.
    * @param src_shift Const reference to shift vector used to put the source
    *              array data box into the index space region of this array data object.
    */
   void sum(const ArrayData<DIM,TYPE>& src,
             const hier::BoxList<DIM>& boxes,
             const hier::IntVector<DIM>& src_shift);

   /*!
    * Calculate the number of bytes needed to stream the data living
    * in the specified box domains.  This routine is only defined for
    * the built-in types of bool, char, double, float, int, and dcomplex.  For
    * all other types, a specialized implementation must be provided.
    * 
    * @param boxes Const reference to box list describing the spatial extents
    *              of the index space regions of interest.
    *              Note: the boxes are assumed to be in the index space of this
    *              array data object.
    * @param src_shift Const reference to vector used to shift the given
    *              boxes into the index space region of this array data object.
    *              Note: this argument is currently ignored.
    */
   int getDataStreamSize(const hier::BoxList<DIM>& boxes,
                         const hier::IntVector<DIM>& src_shift) const;

   /*!
    * Pack data living on the specified index region into the stream.
    * 
    * Note that this routine assumes that the given box region must
    * be shifted to be consistent with the source (this) array box region.
    * 
    * @param stream Reference to stream into which to pack data.
    * @param dest_box Const reference to box describing the spatial extent
    *              of the destination index space region of interest.
    * @param src_shift Const reference to vector used to shift the given
    *              box into the index space region of this (source) array data 
    *              object.
    * 
    * Note: The shifted box must lie completely within the index space of this
    * array data object.  When assertion checking is active, the routine will
    * abort if the box is not contained in the index space of this array.
    */
   void packStream(tbox::AbstractStream& stream,
                   const hier::Box<DIM>& dest_box,
                   const hier::IntVector<DIM>& src_shift) const;

   /*!
    * Pack data living on the specified index regions into the stream.
    *  
    * Note that this routine assumes that the given box regions must
    * be shifted to be consistent with the source (this) array box region.
    *
    * @param stream Reference to stream into which to pack data.
    * @param dest_boxes Const reference to boxes describing the spatial extents
    *              of the destination index space regions of interest.
    * @param src_shift Const reference to vector used to shift the given
    *              boxes into the index space region of this (source) array data
    *              object.
    *
    * Note: The shifted boxes must lie completely within the index space of this
    * array.  If compiled with assertions enabled, the routine will abort if
    * the shifted boxes are not contained in the index space of this array.
    */
   void packStream(tbox::AbstractStream& stream,
                   const hier::BoxList<DIM>& dest_boxes,
                   const hier::IntVector<DIM>& src_shift) const;

   /*!
    * Unpack data from the stream into the index region specified.
    * 
    * @param stream Reference to stream from which to unpack data.
    * @param dest_box Const reference to box describing the spatial extent
    *              of the destination index space region of interest.
    * @param src_offset Const reference to vector used to offset
    *              box into the index space region of some (source) array data
    *              object. Currently, this argument is ignored.
    *
    * Note: The given box must lie completely within the index space of this
    * array data object.  When assertion checking is active, the routine will
    * abort if the box is not contained in the index space of this array.
    */
   void unpackStream(tbox::AbstractStream& stream, 
                     const hier::Box<DIM>& dest_box,
                     const hier::IntVector<DIM>& src_offset);

   /*!
    * Unpack data from the stream into the index regions specified.
    *
    * @param stream Reference to stream from which to unpack data.
    * @param dest_boxes Const reference to box list describing the spatial extents
    *              of the destination index space regions of interest.
    * @param src_offset Const reference to vector used to offset the given
    *              boxes into the index space region of some (source) array data
    *              object. Currently, this argument is ignored.
    *
    * Note: The given boxes must lie completely within the index space of this
    * array data object.  When assertion checking is active, the routine will 
    * abort if some box is not contained in the index space of this array.
    */
   void unpackStream(tbox::AbstractStream& stream, 
                     const hier::BoxList<DIM>& dest_boxes,
                     const hier::IntVector<DIM>& src_offset);

   /*!
    * Unpack data from the stream and add to the array in the index region specified.
    *
    * @param stream Reference to stream from which to unpack data.
    * @param dest_box Const reference to box describing the spatial extent
    *              of the destination index space region of interest.
    * @param src_offset Const reference to vector used to offset the given
    *              box into the index space region of some (source) array data
    *              object. Currently, this argument is ignored.
    *
    * Note: The given box must lie completely within the index space of this
    * array data object.  When assertion checking is active, the routine will
    * abort if the box is not contained in the index space of this array.
    */
   void unpackStreamAndSum(tbox::AbstractStream& stream,
                           const hier::Box<DIM>& dest_box,
                           const hier::IntVector<DIM>& src_offset);

   /*!
    * Unpack data from the stream and ad to the array in the index region specified.
    *
    * @param stream Reference to stream from which to unpack data.
    * @param dest_boxes Const reference to box list describing the spatial extents
    *              of the destination index space regions of interest.
    * @param src_offset Const reference to vector used to offset the given
    *              boxes into the index space region of some (source) array data
    *              object. Currently, this argument is ignored.
    *
    * Note: The given boxes must lie completely within the index space of this
    * array.  If compiled with assertions enabled, the routine will abort if
    * some box is not contained in the index space of this array.
    */
   void unpackStreamAndSum(tbox::AbstractStream& stream,
                           const hier::BoxList<DIM>& dest_boxes,
                           const hier::IntVector<DIM>& src_offset);

   /*!
    * Fill all array values with value t.
    */
   void fillAll(const TYPE& t);

   /*!
    * Fill all array values within the box with value t.
    */
   void fillAll(const TYPE& t, 
                const hier::Box<DIM>& box);

   /*!
    * Fill all array values associated with depth component d with the value t.
    */
   void fill(const TYPE& t, const int d = 0);

   /*!
    * Fill all array values associated with depth component d 
    * within the box with the value t.
    */
   void fill(const TYPE& t, 
             const hier::Box<DIM>& box, 
             const int d = 0);

   /*!
    * Check to make sure that the class version and restart file
    * version are equal.  If so, read in data from database.  This
    * routine calls getSpecializedFromDatabase() to read in the 
    * proper data type.
    *
    * Assertions:  database must be a non-null pointer.
    */
   void getFromDatabase(tbox::Pointer<tbox::Database> database);

   /*!
    * Write out array data object data to database.  This
    * routine calls putSpecializedToDatabase() to read in the 
    * proper data type.  The default behavior (boolean argument is
    * false) is to put all data members in database.  Otherwise, only
    * the array contents are written out.
    *
    * Assertions:  database must be a non-null pointer.
    */
   void putToDatabase(tbox::Pointer<tbox::Database> database,
                      bool data_only = false);

   /*!
    * Use specialized template method to get the correct behavior 
    * when reading in the array of data.
    */
   void getSpecializedFromDatabase(tbox::Pointer<tbox::Database> database);

   /*!
    * Use specialized template method to get the correct behavior 
    * when writing out the array of data.
    */
   void putSpecializedToDatabase(tbox::Pointer<tbox::Database> database);

   /*!
    * The array data iterator iterates over the elements of a box
    * associated with an ArrayData object.  This typedef is
    * convenient link to the ArrayDataIterator<DIM> class.
    */
   typedef ArrayDataIterator<DIM> Iterator;

private:
   ArrayData(const ArrayData<DIM,TYPE>&);	// not implemented
   void operator=(const ArrayData<DIM,TYPE>&);	// not implemented

   /*
    * Private member functions to pack/unpack data to/from buffer.
    *
    * Note: box of this array data object must completely contain given box.
    */
   void packBuffer(TYPE* buffer, 
                   const hier::Box<DIM>& box) const;
   void unpackBuffer(const TYPE* buffer, 
                     const hier::Box<DIM>& box);

   /*
    * Private member functions to unpack data from buffer and add to
    * this array data object.
    *
    * Note: box of this array data object must completely contain given box.
    */
   void unpackBufferAndSum(const TYPE* buffer,
                           const hier::Box<DIM>& box);

   tbox::Array<TYPE> d_array;
   hier::Box<DIM> d_box;
   int d_depth;
   int d_offset;
};

}
}

#ifndef DEBUG_NO_INLINE
#include "ArrayData.I"
#endif

#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "ArrayData.C"
#endif
