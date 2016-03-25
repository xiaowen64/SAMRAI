//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/node/NodeIterator.h $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Iterator for node centered patch data types
//

#ifndef included_pdat_NodeIterator
#define included_pdat_NodeIterator

#include "SAMRAI_config.h"
#include "Box.h"
#include "NodeGeometry.h"
#include "NodeIndex.h"

namespace SAMRAI {
    namespace pdat {

/**
 * Class NodeIterator<DIM> is an iterator that provides methods for
 * stepping through the index space associated with a node centered box.
 * The indices are enumerated in column-major (e.g., Fortran) order.
 * The iterator should be used as follows:
   \verbatim
   hier::Box<DIM> box;
   ...
   for (NodeIterator<DIM> c(box); c; c++) {
      // use index c of the box
   }
   \endverbatim
 * Note that the node iterator may not compile to efficient code, depending
 * on your compiler.  Many compilers are not smart enough to optimize the
 * looping constructs and indexing operations.
 *
 * @see pdat::NodeData
 * @see pdat::NodeGeometry
 * @see pdat::NodeIndex
 */

template<int DIM> class NodeIterator
{
public:
   /**
    * Default constructor for the node iterator.  The iterator must
    * be initialized before it can be used to iterate over a box.
    */
   NodeIterator();

   /**
    * Constructor for the node iterator.  The iterator will enumerate
    * the indices in the argument box.
    */
   NodeIterator(const hier::Box<DIM>& box);

   /**
    * Copy constructor for the node iterator
    */
   NodeIterator(const NodeIterator<DIM>& iterator);

   /**
    * Assignment operator for the node iterator.
    */
   NodeIterator<DIM>& operator=(const NodeIterator<DIM>& iterator);

   /**
    * Destructor for the node iterator.
    */
   ~NodeIterator<DIM>();

   /**
    * Extract the node index corresponding to the iterator position in the box.
    */
   const NodeIndex<DIM>& operator*() const;

   /**
    * Extract the node index corresponding to the iterator position in the box.
    */
   const NodeIndex<DIM>& operator()() const;

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
   bool operator==(const NodeIterator<DIM>& iterator) const;

   /**
    * Test two iterators for inequality (different index values).
    */
   bool operator!=(const NodeIterator<DIM>& iterator) const;

private:
   NodeIndex<DIM> d_index;
   hier::Box<DIM> d_box;
};

}
}
#ifndef DEBUG_NO_INLINE
#include "NodeIterator.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "NodeIterator.C"
#endif
