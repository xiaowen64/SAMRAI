//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/array/ArrayDataIterator.h $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1818 $
// Modified:	$LastChangedDate: 2007-12-20 15:50:44 -0800 (Thu, 20 Dec 2007) $
// Description:	Iterator for array patch data types
//

#ifndef included_pdat_ArrayDataIterator
#define included_pdat_ArrayDataIterator

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_hier_Box
#include "Box.h"
#endif
#ifndef included_hier_Index
#include "Index.h"
#endif

namespace SAMRAI {
    namespace pdat {

/**
 * Class ArrayDataIterator is an iterator that provides methods for
 * stepping through the index space associated with a pdat_ArrayData object.
 * The indices are enumerated in column-major (e.g., Fortran) order.
 * The iterator should be used as follows:

   \verbatim
   hier::Box<DIM> box;
   ...
   for (ArrayDataIterator<DIM> c(box); c; c++) {
      // use index c of the box
   }
   \endverbatim

 * Note that the array data iterator may not compile to efficient code,
 * depending on your compiler.  Many compilers are not smart enough to
 * optimize the looping constructs and indexing operations.
 *
 * @see pdat::ArrayData
 * @see hier::Index
 */

template<int DIM> class ArrayDataIterator
{
public:
   /**
    * Default constructor for the array data iterator.  The iterator must
    * be initialized before it can be used to iterate over a box.
    */
   ArrayDataIterator();

   /**
    * Constructor for the array data iterator.  The iterator will enumerate
    * the indices in the argument box.
    */
   ArrayDataIterator(const hier::Box<DIM>& box);

   /**
    * Copy constructor for the array data iterator
    */
   ArrayDataIterator(const ArrayDataIterator<DIM>& iterator);

   /**
    * Assignment operator for the array data iterator.
    */
   ArrayDataIterator<DIM>& operator=(const ArrayDataIterator<DIM>& iterator);

   /**
    * Destructor for the array data iterator.
    */
   ~ArrayDataIterator<DIM>();

   /**
    * Extract the index corresponding to the iterator position in the box.
    */
   const hier::Index<DIM>& operator*() const;

   /**
    * Extract the index corresponding to the iterator position in the box.
    */
   const hier::Index<DIM>& operator()() const;

   /**
    * Return true if the iterator points to a valid index within the box.
    */
   operator bool() const;

#ifndef LACKS_BOOL_VOID_RESOLUTION
   /**
    * Return a non-NULL if the iterator points to a valid index within the box.
    */
   operator const void*() const;
#endif

   /**
    * Return whether the iterator points to a valid index within the box.
    * This operator mimics the !p operation applied to a pointer p.
    */
   bool operator!() const;

   /**
    * Increment the iterator to point to the next index in the box.
    */
   void operator++(int);

   /**
    * Test two iterators for equality (same index value).
    */
   bool operator==(const ArrayDataIterator<DIM>& iterator) const;

   /**
    * Test two iterators for inequality (different index values).
    */
   bool operator!=(const ArrayDataIterator<DIM>& iterator) const;

private:
   hier::Index<DIM> d_index;
   hier::Box<DIM> d_box;
};

}
}
#ifndef DEBUG_NO_INLINE
#include "ArrayDataIterator.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "ArrayDataIterator.C"
#endif
