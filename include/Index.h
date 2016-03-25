//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/boxes/Index.h $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2147 $
// Modified:	$LastChangedDate: 2008-04-23 16:48:12 -0700 (Wed, 23 Apr 2008) $
// Description:	Interface for the AMR Index object
//

#ifndef included_hier_Index
#define included_hier_Index

#include "SAMRAI_config.h"

#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif

#include "tbox/Array.h"
#include "IntVector.h"


namespace SAMRAI {
   namespace hier {

/**
 * Class Index<DIM> implements a simple n-dimensional integer vector
 * in the AMR index space.  Index is used as lower and upper bounds when
 * creating a box and also when iterating over the cells in a box.  An index
 * is essentially an integer vector but it carries along the notion of indexing
 * into AMR's abstract index space.
 *
 * @see hier::Box
 * @see hier::BoxIterator
 * @see hier::IntVector 
 */

template<int DIM> class Index : public IntVector<DIM>
{
public:
   /**
    * The default constructor for Index creates an uninitialized index.
    */
   Index();

   /**
    * Construct an index with all components equal to the argument.
    */
   Index(const int i);

   /**
    * Construct a two-dimensional index with the value (i,j).
    */
   Index(const int i, const int j);

   /**
    * Construct a three-dimensional index with the value (i,j,k).
    */
   Index(const int i, const int j, const int k);

   /**
    * Construct an n-dimensional index with the values copied 
    * from the integer tbox::Array i of size n.
    */
   Index(const tbox::Array<int> i);

   /**
    * The copy constructor creates an index equal to the argument.
    */
   Index(const Index<DIM>& rhs);

   /**
    * Construct an index equal to the argument IntVector.
    */
   Index(const IntVector<DIM>& rhs);

   /**
    * The assignment operator sets the index equal to the argument.
    */
   Index<DIM>& operator=(const Index<DIM>& rhs);

   /**
    * The assignment operator sets the index equal to the argument IntVector.
    */
   Index<DIM>& operator=(const IntVector<DIM>& rhs);

   /**
    * The index destructor does nothing interesting.
    */
   ~Index();

   /**
    * Plus-equals operator for an index and an integer vector.
    */
   Index<DIM>& operator+=(const IntVector<DIM>& rhs);

   /**
    * Plus operator for an index and an integer vector.
    */
   Index<DIM> operator+(const IntVector<DIM>& rhs) const;

   /**
    * Plus-equals operator for an index and an integer.
    */
   Index<DIM>& operator+=(const int rhs);

   /**
    * Plus operator for an index and an integer.
    */
   Index<DIM> operator+(const int rhs) const;

   /**
    * Minus-equals operator for an index and an integer vector.
    */
   Index<DIM>& operator-=(const IntVector<DIM>& rhs);

   /**
    * Minus operator for an index and an integer vector.
    */
   Index<DIM> operator-(const IntVector<DIM>& rhs) const;

   /**
    * Minus-equals operator for an index and an integer.
    */
   Index<DIM>& operator-=(const int rhs);

   /**
    * Minus operator for an index and an integer.
    */
   Index<DIM> operator-(const int rhs) const;

   /**
    * Times-equals operator for an index and an integer vector.
    */
   Index<DIM>& operator*=(const IntVector<DIM>& rhs);

   /**
    * Times operator for an index and an integer vector.
    */
   Index<DIM> operator*(const IntVector<DIM>& rhs) const;

   /**
    * Times-equals operator for an index and an integer.
    */
   Index<DIM>& operator*=(const int rhs);

   /**
    * Times operator for an index and an integer.
    */
   Index<DIM> operator*(const int rhs) const;

   /**
    * Assign-quotient operator for an index and an integer vector.
    */
   Index<DIM>& operator/=(const IntVector<DIM>& rhs);

   /**
    * Quotient operator for an index and an integer vector.
    */
   Index<DIM> operator/(const IntVector<DIM>& rhs) const;

   /**
    * Assign-quotient operator for an index and an integer.
    */
   Index<DIM>& operator/=(const int rhs);

   /**
    * Quotient operator for an index and an integer.
    */
   Index<DIM> operator/(const int rhs) const;
};


}
}

#ifndef DEBUG_NO_INLINE
#include "Index.I"
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "Index.C"
#endif

#endif

