//
// File:	EdgeIterator.h
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Iterator for edge centered patch data types
//

#ifndef included_pdat_EdgeIterator
#define included_pdat_EdgeIterator

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_hier_Box
#include "Box.h"
#endif
#ifndef included_pdat_EdgeGeometry
#include "EdgeGeometry.h"
#endif
#ifndef included_pdat_EdgeIndex
#include "EdgeIndex.h"
#endif

namespace SAMRAI {
    namespace pdat {

/**
 * Class EdgeIterator<DIM> is an iterator that provides methods for
 * stepping through the index space associated with a edge centered box.
 * The indices are enumerated in column-major (e.g., Fortran) order.
 * The iterator should be used as follows:
   \verbatim
   hier::Box<DIM> box;
   ...
   for (EdgeIterator<DIM> c(box, axis); c; c++) {
      // use index c of the box
   }
   \endverbatim
 * Note that the edge iterator may not compile to efficient code, depending
 * on your compiler.  Many compilers are not smart enough to optimize the
 * looping constructs and indexing operations.
 *
 * @see pdat::EdgeData
 * @see pdat::EdgeGeometry
 * @see pdat::EdgeIndex
 */

template<int DIM> class EdgeIterator
{
public:
   /**
    * Default constructor for the edge iterator.  The iterator must
    * be initialized before it can be used to iterate over a box.
    */
   EdgeIterator();

   /**
    * Constructor for the edge iterator.  The iterator will enumerate
    * the indices in the argument box.
    */
   EdgeIterator(const hier::Box<DIM>& box, const int axis);

   /**
    * Copy constructor for the edge iterator
    */
   EdgeIterator(const EdgeIterator<DIM>& iterator);

   /**
    * Assignment operator for the edge iterator.
    */
   EdgeIterator<DIM>& operator=(const EdgeIterator<DIM>& iterator);

   /**
    * Destructor for the edge iterator.
    */
   ~EdgeIterator<DIM>();

   /**
    * Extract the edge index corresponding to the iterator position in the box.
    */
   const EdgeIndex<DIM>& operator*() const;

   /**
    * Extract the edge index corresponding to the iterator position in the box.
    */
   const EdgeIndex<DIM>& operator()() const;

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
   bool operator==(const EdgeIterator<DIM>& iterator) const;

   /**
    * Test two iterators for inequality (different index values).
    */
   bool operator!=(const EdgeIterator<DIM>& iterator) const;

private:
   hier::Box<DIM> d_box;
   EdgeIndex<DIM> d_index;
};

}
}
#ifndef DEBUG_NO_INLINE
#include "EdgeIterator.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "EdgeIterator.C"
#endif
