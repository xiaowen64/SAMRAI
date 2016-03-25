//
// File:	Array.h
// Package:	SAMRAI toolbox for memory management
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	A simple array template class
//

#ifndef included_tbox_Array
#define included_tbox_Array

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_tbox_ReferenceCounter
#include "tbox/ReferenceCounter.h"
#endif


namespace SAMRAI {
   namespace tbox {


class Arena;
template <class TYPE> class Pointer;

/**
 * Class Array<TYPE> defines a smart pointer to an array of TYPE
 * and manages all reference counting and deallocation of the array (even
 * if the data was originally allocated from an arena).  When the reference
 * count on a Array<TYPE> object goes to zero, the array objects are
 * automatically deallocated.  The array class frees the user from deleting
 * and tracking aliases for object arrays.
 *
 * A block with references count and arena pointer is allocated for all
 * non-empty arrays.  These reference counted blocks are freed at the end
 * of the lifetime of the array.
 *
 * Class TYPE must define a default constructor and an assignment operator
 * for the resize member function.
 *
 * @see tbox::ReferenceCounter
 * @see tbox::Pointer
 */

template <class TYPE>
class Array
{
public:
   /**
    * Create an array of zero elements.
    */
   Array();

   /**
    * Create an array of ``n'' elements.  The storage for the objects
    * is allocated via the standard ``new'' operator.   The
    * ``standard_type'' argument should only be true when TYPE is
    * fundamental data type such as int, double, etc.  In those cases
    * it is not necessary for the argument to be true but if so it will
    * improve efficiency.
    */
   Array(const int n, const bool standard_type = false);

   /**
    * Allocate an array of ``n'' elements using the memory arena.  The
    * storage for the objects is allocated from the specified arena; on
    * deallocation, the storage will be returned to the arena.  The
    * ``standard_type'' argument should only be true when TYPE is 
    * fundamental data type such as int, double, etc.  In those cases
    * it is not necessary for the argument to be true but if so it will
    * improve efficiency. 
    */
   Array(const int n, const Pointer<Arena>& pool,
              const bool standard_type = false);

   /**
    * Const constructor for the array.  This creates an alias to the
    * right hand side and increments the reference count.  The
    * ``standard_type'' argument should only be true when TYPE is
    * fundamental data type such as int, double, etc.  In those cases
    * it is not necessary for the argument to be true but if so it will
    * improve efficiency.
    *
    * CAUTION: invoking resizeArray() forces a deep copy.
    * Upon return, two objects that formerly were aliases to the
    * same underlying data will point to separate data.  For this
    * reason, it is best to pass a Array by reference, instead
    * of by value.
    */
   Array(const Array<TYPE>& rhs);

   /**
    * Destructor for the array.  If the reference count for the array data
    * has gone to zero, then the array data is deallocated from the memory
    * arena from which it was allocated.
    */
   ~Array();

   /**
    * Array assignment.  The assignment operator copies a pointer to the
    * array data and increments the reference count.  Both array objects refer 
    * to the same data, and changes to individual array entry values in one will 
    * be reflected in the other array.  However, this assignment operation DOES NOT 
    * involve a "deep copy" (see the resizeArray() routines below). Thus, changes 
    * to one Array object container will not necessarily be reflected in the 
    * other container.
    */
   Array<TYPE>& operator=(const Array<TYPE>& rhs);

   /**
    * Non-const array subscripting.  Return a reference the object at array
    * index ``i'' (between 0 and N-1, where N is the number of elements in
    * the array.
    */
   TYPE& operator[](const int i);

   /**
    * Const array subscripting.  Return a const reference to the object
    * at array index ``i'' (between 0 and N-1, where N is the number of
    * elements in the array.
    */
   const TYPE& operator[](const int i) const;

   /**
    * Test whether the array is NULL (has any elements).
    */
   bool isNull() const;

   /**
    * Set the length of the array to zero.  If the reference count for
    * the objects has dropped to zero, then the array data is deallocated.
    */
   void setNull();

   /**
    * Return a non-const pointer to the i-th object.  The index must be
    * between 0 and N-1, where N is the number of elements in the array.
    */
   TYPE *getPointer(const int i = 0);

   /**
    * Return a const pointer to the i-th object.  The index must be
    * between 0 and N-1, where N is the number of elements in the array.
    */
   const TYPE *getPointer(const int i = 0) const;

   /**
    * Return the number of elements in the array.
    */
   int getSize() const;

   /**
    * Return the number of elements in the array.  Identical to getSize(),
    * but this method is common to several container classes.
    */
   int size() const;

   /**
    * Resize the array by allocating new array storage and copying from the
    * old array into the new; i.e., a "deep" copy.  Space for the new array 
    * is allocated via the standard ``new'' operator. 
    */
   void resizeArray(const int n);

   /**
    * Resize the array by allocating new array storage from the provided 
    * memory pool and copying from the old array into the new; i.e., a "deep" 
    * copy.  
    */
   void resizeArray(const int n, const Pointer<Arena>& pool);

private:
   TYPE *allocateObjects(const int n, Arena *arena);
   void deleteObjects();

   TYPE *d_objects;
   ReferenceCounter *d_counter;
   int d_elements;

   bool d_standard_type;
};


}
}

#ifndef DEBUG_NO_INLINE
#include "tbox/Array.I"
#endif
#endif
