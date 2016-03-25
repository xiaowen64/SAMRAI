//
// File:	ArrayData.h
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Templated array data structure supporting patch data types
//

#ifndef included_tbox_ArrayData
#define included_tbox_ArrayData

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_hier_Box
#include "Box.h"
#endif
#ifndef included_hier_BoxList
#include "BoxList.h"
#endif
#ifndef included_hier_Index
#include "Index.h"
#endif
#ifndef included_hier_IntVector
#include "IntVector.h"
#endif
#ifndef included_tbox_ArrayDataIterator
#include "ArrayDataIterator.h"
#endif
#ifndef included_tbox_AbstractStream
#include "tbox/AbstractStream.h"
#endif
#ifndef included_tbox_Arena
#include "tbox/Arena.h"
#endif
#ifndef included_tbox_Array
#include "tbox/Array.h"
#endif
#ifndef included_tbox_Complex
#include "tbox/Complex.h"
#endif
#ifndef included_tbox_Database
#include "tbox/Database.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif

#ifndef NULL
#define NULL (0)
#endif

namespace SAMRAI {
    namespace pdat {

/**
 * Class ArrayData<DIM> is a basic templated array structure defined
 * over the index space of a box (with a specified depth) that provides
 * the support for the various flavors of patch data subclasses.
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
   /**
    * The no-arguments constructor creates an empty array data object.
    * The initializeArray() member function must be called before the
    * array can be used.
    */
   ArrayData();

   /**
    * The constructor for an array data object.  The box describes the
    * spatial extents of the index space and the depth gives the number
    * of data components for each spatial location in the array.  tbox::Array
    * memory is allocated from the specified memory pool.
    */
   ArrayData(
      const hier::Box<DIM>& box,
      const int depth,
      tbox::Pointer<tbox::Arena> pool = tbox::Pointer<tbox::Arena>(NULL));

   /**
    * The destructor for an array data object releases all memory allocated
    * for the array elements.
    */
   ~ArrayData();

   /**
    * Return the box over which the array is defined.
    */
   const hier::Box<DIM>& getBox() const;

   /**
    * Initialize the array data for the specified box and depth.
    * This routine is normally called to initialize an array that
    * was created by the no-arguments constructor.
    */
   void initializeArray(
      const hier::Box<DIM>& box,
      const int depth,
      const tbox::Pointer<tbox::Arena>& pool);

   /*!
    * @brief Returns whether array has been properly initialized
    * proper depth and box.
    *
    * Only initialized arrays can be used.  Uninitialized arrays
    * should not be used until initializeArray() is called.
    */
   bool isInitialized() const;

   /**
    * Set the array data to some ``undefined'' state.  For floats and
    * doubles, this means setting data to signaling NaNs that cause
    * a floating point exception when used in a numerical expression.
    */
   void undefineData();

   /**
    * Return the amount of memory space needed to allocate a box of the
    * specified size and depth.  This function will need to be redefined
    * for data types other than the standard bool, char, double, float,
    * and int types.
    */
   static size_t getSizeOfData(const hier::Box<DIM>& box, const int depth);

   /**
    * Return true if the TYPE is a standard fundamental data type
    */
   bool isStandardType();

   /**
    * Copy data from the source array data on the specified index domain.
    * This routine will intersect the specified box against the source and
    * destination boxes to find the region of intersection.
    */
   void copy(const ArrayData<DIM,TYPE>& src, const hier::Box<DIM>& box);

   /**
    * Copy data shifted in the index space from the source into the
    * destination.  The boxes are in the destination index space, and
    * the source data is taken from the boxes shifted into the source
    * index space.  This routine will intersect the specified box against
    * the source and destination boxes to find the region of intersection.
    */
   void copy(
      const ArrayData<DIM,TYPE>& src,
      const hier::Box<DIM>& box,
      const hier::IntVector<DIM>& offset);

   /**
    * Copy data shifted in the index space from the source into the
    * destination.  The boxes are in the destination index space, and
    * the source data is taken from the boxes shifted into the source
    * index space.
    */
   void copy(
      const ArrayData<DIM,TYPE>& src,
      const hier::BoxList<DIM>& boxes,
      const hier::IntVector<DIM>& offset);

   /*!
    * @brief Copy one depth of a source array to a destination depth
    * of another.
    *
    * This routine will intersect the specified box against the source and
    * destination boxes to find the region of intersection.
    */
   void copyDepth(int dst_depth,
		  const ArrayData<DIM,TYPE>& src,
		  int src_depth,
		  const hier::Box<DIM>& box);

   /**
    * Return whether the amount of buffer space in the message stream can be
    * estimated using the box alone (i.e., without specific type information).
    * For built-in types (bool, char, double, float, int, and dcomplex), this 
    * routine returns true.  For other user-defined data types that may require
    * special handling, the user MUST define a different implementation.
    */
   static bool canEstimateStreamSizeFromBox();

   /**
    * Calculate the number of bytes needed to stream the data lying
    * in the specified box domains.  This routine is only defined for
    * the built-in types of bool, char, double, float, and int.  For
    * all other types, the user must define a specialized implementation.
    */
   int getDataStreamSize(const hier::BoxList<DIM>& dest_boxes,
                         const hier::IntVector<DIM>& source_offset) const;

   /**
    * Pack data lying on the specified index set into the output stream.
    * The shifted box must lie completely within the index space of the
    * array.  If compiled with assertions enabled, the packing routine
    * will abort if the box is not contained within the index space of
    * the array.
    */
   void packStream(
      tbox::AbstractStream& stream,
      const hier::Box<DIM>& dest_box,
      const hier::IntVector<DIM>& source_offset) const;

   /**
    * Pack data lying on the specified index set into the output stream.
    * The shifted boxes in the box list must lie completely within the
    * index space of the array.  If compiled with assertions enabled, the
    * packing routine will abort if the boxes are not contained within the
    * index space of the array.
    */
   void packStream(
      tbox::AbstractStream& stream,
      const hier::BoxList<DIM>& dest_boxes,
      const hier::IntVector<DIM>& source_offset) const;

   /**
    * Unpack data from the message stream into the index set lying under
    * the specified index set.  The box must lie completely within the
    * index space of the array.  If compiled with assertions enabled, the
    * unpacking routine will abort if the box is not contained within the
    * index space of the array.
    */
   void unpackStream(tbox::AbstractStream& stream, 
                     const hier::Box<DIM>& dest_box,
                     const hier::IntVector<DIM>& source_offset);

   /**
    * Unpack data from the message stream into the index set lying under
    * the specified index set.  The boxes in the box list must lie completely
    * within the index space of the array.  If compiled with assertions
    * enabled, the unpacking routine will abort if the boxes are not
    * contained within the index space of the array.
    */
   void unpackStream(tbox::AbstractStream& stream, 
                     const hier::BoxList<DIM>& boxes,
                     const hier::IntVector<DIM>& source_offset);

   /**
    * Return the depth (e.g., the number of components in each spatial
    * location) of the array.
    */
   int getDepth() const;

   /**
    * Return the offset (e.g., the number of data values for each 
    * depth component) of the array.
    */
   int getOffset() const;

   /**
    * Get a pointer to the beginning of a particular component of the
    * patch data array.
    */
   TYPE *getPointer(const int d = 0);

   /**
    * Get a const pointer to the beginning of a particular component
    * of the patch data array.
    */
   const TYPE *getPointer(const int d = 0) const;

   /**
    * hier::Index into the array using an index and the component.
    */
   TYPE& operator()(const hier::Index<DIM>& i, const int d);

   /**
    * hier::Index into the array (via a const reference) using an index
    * and the component.
    */
   const TYPE& operator()(const hier::Index<DIM>& i, const int d) const;

   /**
    * Fill all components with value t.
    */
   void fillAll(const TYPE& t);

   /**
    * Fill all components within the box with value t.
    */
   void fillAll(const TYPE& t, const hier::Box<DIM>& box);

   /**
    * Fill all values of component d with the value t.
    */
   void fill(const TYPE& t, const int d = 0);

   /**
    * Fill all values of component d within the box with the value t.
    */
   void fill(const TYPE& t, const hier::Box<DIM>& box, const int d = 0);

   /**
    * Check to make sure that the class version and restart file
    * version are equal.  If so, read in data from database.  This
    * routine calls getSpecializedFromDatabase() to read in the 
    * proper data type.
    *
    * Assertions:  database must be a non-null pointer.
    */
   void getFromDatabase(tbox::Pointer<tbox::Database> database);

   /**
    * Writes out data for class representation to database.  This
    * routine calls putSpecializedToDatabase() to read in the 
    * proper data type.  The default behavior (boolean argument is
    * false) is to put all data members in database.  Otherwise, only
    * the array contents are written out.
    *
    * Assertions:  database must be a non-null pointer.
    */
   void putToDatabase(tbox::Pointer<tbox::Database> database,
                      bool data_only = false);

   /**
    * Use specialized template method to get the correct behavior 
    * when reading in the array of data.
    */
   void getSpecializedFromDatabase(tbox::Pointer<tbox::Database> database);

   /**
    * Use specialized template method to get the correct behavior 
    * when writing out the array of data.
    */
   void putSpecializedToDatabase(tbox::Pointer<tbox::Database> database);

   /**
    * The array data iterator iterates over the elements of a box
    * associated with an ArrayData object.  This typedef is
    * convenient link to the ArrayDataIterator<DIM> class.
    */
   typedef ArrayDataIterator<DIM> Iterator;

private:
   void packBuffer(TYPE *buffer, const hier::Box<DIM>& box) const;
   void unpackBuffer(const TYPE *buffer, const hier::Box<DIM>& box);

   ArrayData(const ArrayData<DIM,TYPE>&);	// not implemented
   void operator=(const ArrayData<DIM,TYPE>&);	// not implemented

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
