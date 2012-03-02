/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Iterator for side centered patch data types
 *
 ************************************************************************/

#ifndef included_pdat_SideIterator
#define included_pdat_SideIterator

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/pdat/SideGeometry.h"
#include "SAMRAI/pdat/SideIndex.h"

namespace SAMRAI {
namespace pdat {

/**
 * Class SideIterator is an iterator that provides methods for
 * stepping through the index space associated with a side centered box.
 * The indices are enumerated in column-major (e.g., Fortran) order.
 * The iterator should be used as follows:
 * \verbatim
 * hier::Box box;
 * ...
 * for (SideIterator c(box, axis); c; c++) {
 *    // use index c of the box
 * }
 * \endverbatim
 * Note that the side iterator may not compile to efficient code, depending
 * on your compiler.  Many compilers are not smart enough to optimize the
 * looping constructs and indexing operations.
 *
 * @see pdat::SideData
 * @see pdat::SideGeometry
 * @see pdat::SideIndex
 */

class SideIterator
{
public:
   /**
    * Constructor for the side iterator.  The iterator will enumerate
    * the indices in the argument box.
    */
   SideIterator(
      const hier::Box& box,
      const int axis);

   /**
    * Copy constructor for the side iterator
    */
   SideIterator(
      const SideIterator& iterator);

   /**
    * Assignment operator for the side iterator.
    */
   SideIterator&
   operator = (
      const SideIterator& iterator);

   /**
    * Destructor for the side iterator.
    */
   ~SideIterator();

   /**
    * Extract the side index corresponding to the iterator position in the box.
    */
   const SideIndex&
   operator * () const;

   /**
    * Extract the side index corresponding to the iterator position in the box.
    */
   const SideIndex&
   operator () () const;

   /**
    * Return true if the iterator points to a valid index within the box.
    */
   operator bool () const;

#ifndef LACKS_BOOL_VOID_RESOLUTION
   /**
    * Return a non-NULL if the iterator points to a valid index within the box.
    */
   operator const void * () const;
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
      const SideIterator& iterator) const;

   /**
    * Test two iterators for inequality (different index values).
    */
   bool
   operator != (
      const SideIterator& iterator) const;

private:
   SideIndex d_index;
   hier::Box d_box;
};

}
}
#ifdef SAMRAI_INLINE
#include "SAMRAI/pdat/SideIterator.I"
#endif
#endif
