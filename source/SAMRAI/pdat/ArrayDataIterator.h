/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Iterator for array patch data types 
 *
 ************************************************************************/

#ifndef included_pdat_ArrayDataIterator
#define included_pdat_ArrayDataIterator

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/Index.h"

namespace SAMRAI {
namespace pdat {

/**
 * Class ArrayDataIterator is an iterator that provides methods for
 * stepping through the index space associated with a pdat_ArrayData object.
 * The indices are enumerated in column-major (e.g., Fortran) order.
 * The iterator should be used as follows:
 *
 * \verbatim
 * hier::Box box;
 * ...
 * for (ArrayDataIterator c(box); c; c++) {
 *    // use index c of the box
 * }
 * \endverbatim
 *
 * Note that the array data iterator may not compile to efficient code,
 * depending on your compiler.  Many compilers are not smart enough to
 * optimize the looping constructs and indexing operations.
 *
 * @see pdat::ArrayData
 * @see hier::Index
 */

class ArrayDataIterator
{
public:
   /**
    * Constructor for the array data iterator.  The iterator will enumerate
    * the indices in the argument box.
    */
   explicit ArrayDataIterator(
      const hier::Box& box);

   /**
    * Copy constructor for the array data iterator
    */
   ArrayDataIterator(
      const ArrayDataIterator& iterator);

   /**
    * Assignment operator for the array data iterator.
    */
   ArrayDataIterator&
   operator = (
      const ArrayDataIterator& iterator);

   /**
    * Destructor for the array data iterator.
    */
   ~ArrayDataIterator();

   /**
    * Extract the index corresponding to the iterator position in the box.
    */
   const hier::Index&
   operator * () const;

   /**
    * Extract the index corresponding to the iterator position in the box.
    */
   const hier::Index&
   operator () () const;

   /**
    * Return true if the iterator points to a valid index within the box.
    */
   operator bool () const;

#ifndef LACKS_BOOL_VOID_RESOLUTION
   /**
    * Return a non-NULL if the iterator points to a valid index within the box.
    */
   operator const void
   * () const;
#endif

   /**
    * Return whether the iterator points to a valid index within the box.
    * This operator mimics the !p operation applied to a pointer p.
    */
   bool
   operator ! () const;

   /**
    * Increment the iterator to point to the next index in the box.
    */
   void
   operator ++ (
      int);

   /**
    * Test two iterators for equality (same index value).
    */
   bool
   operator == (
      const ArrayDataIterator& iterator) const;

   /**
    * Test two iterators for inequality (different index values).
    */
   bool
   operator != (
      const ArrayDataIterator& iterator) const;

private:
   hier::Index d_index;
   hier::Box d_box;
};

}
}
#ifdef SAMRAI_INLINE
#include "SAMRAI/pdat/ArrayDataIterator.I"
#endif
#endif
