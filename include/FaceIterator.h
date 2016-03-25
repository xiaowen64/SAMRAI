//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/face/FaceIterator.h $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Iterator for face centered patch data types
//

#ifndef included_pdat_FaceIterator
#define included_pdat_FaceIterator

#include "SAMRAI_config.h"
#include "Box.h"
#include "FaceGeometry.h"
#include "FaceIndex.h"

namespace SAMRAI {
    namespace pdat {

/**
 * Class FaceIterator<DIM> is an iterator that provides methods for
 * stepping through the index space associated with a face centered box.
 * The indices are enumerated in column-major (e.g., Fortran) order.
 * The iterator should be used as follows:
   \verbatim
   hier::Box<DIM> box;
   ...
   for (FaceIterator<DIM> c(box, axis); c; c++) {
      // use index c of the box
   }
   \endverbatim
 * Note that the face iterator may not compile to efficient code, depending
 * on your compiler.  Many compilers are not smart enough to optimize the
 * looping constructs and indexing operations.
 *
 * @see pdat::FaceData
 * @see pdat::FaceGeometry
 * @see pdat::FaceIndex
 */

template<int DIM> class FaceIterator
{
public:
   /**
    * Default constructor for the face iterator.  The iterator must
    * be initialized before it can be used to iterate over a box.
    */
   FaceIterator();

   /**
    * Constructor for the face iterator.  The iterator will enumerate
    * the indices in the argument box.
    */
   FaceIterator(const hier::Box<DIM>& box, const int axis);

   /**
    * Copy constructor for the face iterator
    */
   FaceIterator(const FaceIterator<DIM>& iterator);

   /**
    * Assignment operator for the face iterator.
    */
   FaceIterator<DIM>& operator=(const FaceIterator<DIM>& iterator);

   /**
    * Destructor for the face iterator.
    */
   ~FaceIterator<DIM>();

   /**
    * Extract the face index corresponding to the iterator position in the box.
    */
   const FaceIndex<DIM>& operator*() const;

   /**
    * Extract the face index corresponding to the iterator position in the box.
    */
   const FaceIndex<DIM>& operator()() const;

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
   bool operator==(const FaceIterator<DIM>& iterator) const;

   /**
    * Test two iterators for inequality (different index values).
    */
   bool operator!=(const FaceIterator<DIM>& iterator) const;

private:
   FaceIndex<DIM> d_index;
   hier::Box<DIM> d_box;
};

}
}
#ifndef DEBUG_NO_INLINE
#include "FaceIterator.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "FaceIterator.C"
#endif
