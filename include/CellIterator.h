//
// File:	CellIterator.h
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Iterator for cell centered patch data types
//

#ifndef included_pdat_CellIterator
#define included_pdat_CellIterator

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_hier_Box
#include "Box.h"
#endif
#ifndef included_pdat_CellGeometry
#include "CellGeometry.h"
#endif
#ifndef included_pdat_CellIndex
#include "CellIndex.h"
#endif

namespace SAMRAI {
    namespace pdat {

/**
 * Class CellIterator<DIM> is an iterator that provides methods for
 * stepping through the index space associated with a cell centered box.
 * The indices are enumerated in column-major (e.g., Fortran) order.
 * The iterator should be used as follows:
   \verbatim
   hier::Box<DIM> box;
   ...
   for (CellIterator<DIM> c(box); c; c++) {
      // use index c of the box
   }
   \endverbatim
 * Note that the cell iterator may not compile to efficient code, depending
 * on your compiler.  Many compilers are not smart enough to optimize the
 * looping constructs and indexing operations.
 *
 * @see pdat::CellData
 * @see pdat::CellGeometry
 * @see pdat::CellIndex
 */

template<int DIM> class CellIterator
{
public:
   /**
    * Default constructor for the cell iterator.  The iterator must
    * be initialized before it can be used to iterate over a box.
    */
   CellIterator();

   /**
    * Constructor for the cell iterator.  The iterator will enumerate
    * the indices in the argument box.
    */
   CellIterator(const hier::Box<DIM>& box);

   /**
    * Copy constructor for the cell iterator
    */
   CellIterator(const CellIterator<DIM>& iterator);

   /**
    * Assignment operator for the cell iterator.
    */
   CellIterator<DIM>& operator=(const CellIterator<DIM>& iterator);

   /**
    * Destructor for the cell iterator.
    */
   ~CellIterator<DIM>();

   /**
    * Extract the cell index corresponding to the iterator position in the box.
    */
   const CellIndex<DIM>& operator*() const;

   /**
    * Extract the cell index corresponding to the iterator position in the box.
    */
   const CellIndex<DIM>& operator()() const;

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
   bool operator==(const CellIterator<DIM>& iterator) const;

   /**
    * Test two iterators for inequality (different index values).
    */
   bool operator!=(const CellIterator<DIM>& iterator) const;

private:
   CellIndex<DIM> d_index;
   hier::Box<DIM> d_box;
};

}
}
#ifndef DEBUG_NO_INLINE
#include "CellIterator.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "CellIterator.C"
#endif
