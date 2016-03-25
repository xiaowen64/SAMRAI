//
// File:	Index.h
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 601 $
// Modified:	$Date: 2005-09-06 11:23:15 -0700 (Tue, 06 Sep 2005) $
// Description:	Interface for the AMR Index object
//

#ifndef included_hier_Index
#define included_hier_Index

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_iostream
#define included_iostream
#include <iostream>
using namespace std;
#endif

#ifndef included_tbox_Array
#include "tbox/Array.h"
#endif
#ifndef included_hier_IntVector
#include "IntVector.h"
#endif


namespace SAMRAI {
   namespace hier {

/**
 * Class Index<DIM> implements a simple n-dimensional integer vector
 * in the AMR index space.  Index is used as lower and upper bounds when
 * creating a box and also when iterating over the cells in a box.  An index
 * is essentially an integer vector but it carries along the notion of indexing
 * into AMR's abstract index space.
 *
 * Class Index<DIM> is translated into classes Index1, Index2,
 * and Index3 after being passed through a preprocessor.
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

#if INCLUDE_DEPRECATED < 2
   /**
    * Construct a two-dimensional index with the value (i,j).
    */
   Index(const int i, const int j);

   /**
    * Construct a three-dimensional index with the value (i,j,k).
    */
   Index(const int i, const int j, const int k);
#endif

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
    * The assignment operator sets the index equal to the argument.
    */
   Index<DIM>& operator=(const Index<DIM>& rhs);

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

