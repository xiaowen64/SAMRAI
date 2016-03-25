//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/boxes/IntVector.h $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2143 $
// Modified:	$LastChangedDate: 2008-04-23 09:37:15 -0700 (Wed, 23 Apr 2008) $
// Description:	A N-dimensional integer vector
//

#ifndef included_hier_IntVector
#define included_hier_IntVector

#include "SAMRAI_config.h"

#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif

#include "tbox/Array.h"

namespace SAMRAI {
   namespace hier {

/**
 * Class IntVector implements a simple N-dimensional integer
 * vector.  This class is the base class for most of the simple indexing
 * classes.
 *
 */

template<int DIM> class IntVector
{
public:
   /**
    * The default constructor creates an uninitialized vector.
    */
   IntVector();

   /**
    * Construct an integer vector with all components equal to the argument.
    */
   IntVector(const int i);


   /**
    * Construct a two-dimensional integer vector with the value (i,j).
    * Provided for 2D and 3D only.
    */
   IntVector(const int i, const int j);

   /**
    * Construct a three-dimensional integer vector with the value (i,j,k).
    * Provided for 2D and 3D only.
    */
   IntVector(const int i, const int j, const int k);

   /**
    * Construct a n-dimensional integer vector with the value with
    * values provided by the array.
    */
   IntVector(const tbox::Array<int>& a);

   /**
    * Construct an integer vector equal to the argument.
    */
   IntVector(const IntVector<DIM>& rhs);

   /**
    * The assignment operator sets the integer vector equal to the argument.
    */
   IntVector& operator=(const IntVector<DIM>& rhs);

   /**
    * The integer vector destructor does nothing interesting.
    */
   virtual ~IntVector();

   /**
    * Return the specified component of the vector.  No bounds checking.
    */
   int& operator()(const int i);

   /**
    * Return the specified component of the vector as a const integer.
    * No bounds checking.
    */
   const int& operator()(const int i) const;

   /**
    * Return a pointer to the beginning of the vector data values.
    */
   operator int*();

   /**
    * Return a const pointer to the beginning of the vector data values.
    */
   operator const int*() const;

   /**
    * Plus-equals operator for two integer vectors.
    */
   IntVector<DIM>& operator+=(const IntVector<DIM>& rhs);

   /**
    * Plus operator for two integer vectors.
    */
   IntVector<DIM> operator+(const IntVector<DIM>& rhs) const;

   /**
    * Plus-equals operator for an integer vector and an integer.
    */
   IntVector<DIM>& operator+=(const int rhs);

   /**
    * Plus operator for an integer vector and an integer.
    */
   IntVector<DIM> operator+(const int rhs) const;

   /**
    * Minus-equals operator for two integer vectors.
    */
   IntVector<DIM>& operator-=(const IntVector<DIM>& rhs);

   /**
    * Minus operator for two integer vectors.
    */
   IntVector<DIM> operator-(const IntVector<DIM>& rhs) const;

   /**
    * Minus-equals operator for an integer vector and an integer.
    */
   IntVector<DIM>& operator-=(const int rhs);

   /**
    * Minus operator for an integer vector and an integer.
    */
   IntVector<DIM> operator-(const int rhs) const;

   /**
    * Times-equals operator for two integer vectors.
    */
   IntVector<DIM>& operator*=(const IntVector<DIM>& rhs);

   /**
    * Times operator for two integer vectors.
    */
   IntVector<DIM> operator*(const IntVector<DIM>& rhs) const;

   /**
    * Times-equals operator for an integer vector and an integer.
    */
   IntVector<DIM>& operator*=(const int rhs);

   /**
    * Times operator for an integer vector and an integer.
    */
   IntVector<DIM> operator*(const int rhs) const;

   /**
    * Assign-quotient operator for two integer vectors.
    */
   IntVector<DIM>& operator/=(const IntVector<DIM>& rhs);

   /**
    * Quotient operator for two integer vectors.
    */
   IntVector<DIM> operator/(const IntVector<DIM>& rhs) const;

   /**
    * Assign-quotient operator for an integer vector and an integer.
    */
   IntVector<DIM>& operator/=(const int rhs);

   /**
    * Quotient operator for an integer vector and an integer.
    */
   IntVector<DIM> operator/(const int rhs) const;

   /**
    * Unary minus to negate an integer vector.
    */
   IntVector<DIM> operator-() const;

   /**
    * Returns true if two vector objects are equal.  All components
    * must be the same for equality.
    */
   bool operator==(const IntVector<DIM>& rhs) const;

   /**
    * Returns true if two vector objects are not equal.  Any of
    * the components may be different for inequality.
    */
   bool operator!=(const IntVector<DIM>& rhs) const;

   /**
    * Returns true if each integer in vector is less than
    * corresponding integer in comparison vector.
    */
   bool operator<(const IntVector<DIM>& rhs) const;

   /**
    * Returns true if each integer in vector is less or equal to
    * corresponding integer in comparison vector.
    */
   bool operator<=(const IntVector<DIM>& rhs) const;

   /**
    * Returns true if each integer in vector is greater than
    * corresponding integer in comparison vector.
    */
   bool operator>(const IntVector<DIM>& rhs) const;

   /**
    * Returns true if each integer in vector is greater or equal to
    * corresponding integer in comparison vector.
    */
   bool operator>=(const IntVector<DIM>& rhs) const;

   /**
    * Return the component-wise minimum of two integer vector objects.
    */
   void min(const IntVector<DIM>& rhs);

   /**
    * Return the minimum entry in an integer vector.
    */
   int min() const;

   /**
    * Return the component-wise maximum of two integer vector objects.
    */
   void max(const IntVector<DIM>& rhs);

   /**
    * Return the maximum entry in an integer vector.
    */
   int max() const;

   /**
    * Utility function to take the minimum of two integer vector objects.
    */
   static IntVector<DIM> min(const IntVector<DIM>& a,
                              const IntVector<DIM>& b);

   /**
    * Utility function to take the maximum of two integer vector objects.
    */
   static IntVector<DIM> max(const IntVector<DIM>& a,
                              const IntVector<DIM>& b);

   /**
    * Return the product of the entries in the integer vector.
    */
   int getProduct() const;

   /**
    * Read an integer vector from an input stream.  The format for
    * the input is (i0,...,in) for an n-dimensional vector.
    */
    template<int DIMENSION> 
       friend std::istream& operator>> (std::istream& s, IntVector<DIMENSION>& rhs);

   /**
    * Write an integer vector into an output stream.  The format for
    * the output is (i0,...,in) for an n-dimensional vector.
    */
    template<int DIMENSION> 
       friend std::ostream& operator<< (std::ostream& s, 
				   const IntVector<DIMENSION>& rhs);

private:
   int d_vector[DIM];
};

}
}

#ifndef DEBUG_NO_INLINE
#include "IntVector.I"
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "IntVector.C"
#endif

#endif

