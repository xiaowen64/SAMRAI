/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   A N-dimensional integer vector 
 *
 ************************************************************************/

#ifndef included_hier_IntVector
#define included_hier_IntVector

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/Dimension.h"

#include <vector>
#include <iostream>

namespace SAMRAI {

namespace hier {

/**
 * Class IntVector implements a simple N-dimensional integer
 * vector.  This class is the base class for most of the simple indexing
 * classes.
 *
 */

class IntVector
{
public:
   /**
    * Creates an uninitialized vector.
    */
   explicit IntVector(
      const tbox::Dimension& dim);

   /**
    * Construct an integer vector with all components equal to the argument.
    */
   explicit IntVector(
      const tbox::Dimension& dim,
      const int i);

   /**
    * Construct a n-dimensional integer vector with the value with
    * values provided by the array.
    *
    * Dimension inferred from array size.
    */
   explicit IntVector(
      const tbox::Array<int>& a);

   /**
    * Construct a n-dimensional integer vector with the value with
    * values provided by the array.
    *
    */
   explicit IntVector(
      const tbox::Dimension& dim,
      const int array[]);

   /**
    * Construct an integer vector equal to the argument.
    */
   IntVector(
      const IntVector& rhs);

   /**
    * The assignment operator sets the integer vector equal to the argument.
    *
    * An assignment to an uninitialized Index is allowed but assigning
    * from an uninitialized Index will result in an assert.
    */
   IntVector&
   operator = (
      const IntVector& rhs);

   /**
    * The integer vector destructor does nothing interesting.
    */
   virtual ~IntVector();

   /**
    * Return the specified component of the vector.  No bounds checking.
    */
   int&
   operator [] (
      const int i);

   /**
    * Return the specified component of the vector as a const integer.
    * No bounds checking.
    */
   const int&
   operator [] (
      const int i) const;

   /**
    * Return the specified component of the vector.  No bounds checking.
    */
   int&
   operator () (
      const int i);

   /**
    * Return the specified component of the vector as a const integer.
    * No bounds checking.
    */
   const int&
   operator () (
      const int i) const;

   /**
    * Plus-equals operator for two integer vectors.
    */
   IntVector&
   operator += (
      const IntVector& rhs);

   /**
    * Plus operator for two integer vectors.
    */
   IntVector
   operator + (
      const IntVector& rhs) const;

   /**
    * Plus-equals operator for an integer vector and an integer.
    */
   IntVector&
   operator += (
      const int rhs);

   /**
    * Plus operator for an integer vector and an integer.
    */
   IntVector
   operator + (
      const int rhs) const;

   /**
    * Minus-equals operator for two integer vectors.
    */
   IntVector&
   operator -= (
      const IntVector& rhs);

   /**
    * Minus operator for two integer vectors.
    */
   IntVector
   operator - (
      const IntVector& rhs) const;

   /**
    * Minus-equals operator for an integer vector and an integer.
    */
   IntVector&
   operator -= (
      const int rhs);

   /**
    * Minus operator for an integer vector and an integer.
    */
   IntVector
   operator - (
      const int rhs) const;

   /**
    * Times-equals operator for two integer vectors.
    */
   IntVector&
   operator *= (
      const IntVector& rhs);

   /**
    * Times operator for two integer vectors.
    */
   IntVector
   operator * (
      const IntVector& rhs) const;

   /**
    * Times-equals operator for an integer vector and an integer.
    */
   IntVector&
   operator *= (
      const int rhs);

   /**
    * Times operator for an integer vector and an integer.
    */
   IntVector
   operator * (
      const int rhs) const;

   /**
    * Assign-quotient operator for two integer vectors.
    */
   IntVector&
   operator /= (
      const IntVector& rhs);

   /**
    * Component-wise ceiling quotient (integer divide with rounding up).
    */
   void
   ceiling(
      const IntVector& denominator);

   /**
    * Component-wise ceiling quotient (integer divide with rounding up).
    */
   static IntVector
   ceiling(
      const IntVector& numerator,
      const IntVector& denominator);

   /**
    * Quotient operator for two integer vectors.
    */
   IntVector
   operator / (
      const IntVector& rhs) const;

   /**
    * Assign-quotient operator for an integer vector and an integer.
    */
   IntVector&
   operator /= (
      const int rhs);

   /**
    * Quotient operator for an integer vector and an integer.
    */
   IntVector
   operator / (
      const int rhs) const;

   /**
    * Unary minus to negate an integer vector.
    */
   IntVector
   operator - () const;

   /**
    * Returns true if all components are equal to a given integer.
    */
   bool
   operator == (
      int rhs) const;

   /**
    * Returns true if some components are not equal to a given integer.
    */
   bool
   operator != (
      int rhs) const;

   /**
    * Returns true if two vector objects are equal.  All components
    * must be the same for equality.
    */
   bool
   operator == (
      const IntVector& rhs) const;

   /**
    * Returns true if two vector objects are not equal.  Any of
    * the components may be different for inequality.
    */
   bool
   operator != (
      const IntVector& rhs) const;

   /**
    * Returns true if each integer in vector is less than
    * corresponding integer in comparison vector.
    */
   bool
   operator < (
      const IntVector& rhs) const;

   /**
    * Returns true if each integer in vector is less or equal to
    * corresponding integer in comparison vector.
    */
   bool
   operator <= (
      const IntVector& rhs) const;

   /**
    * Returns true if each integer in vector is greater than
    * corresponding integer in comparison vector.
    */
   bool
   operator > (
      const IntVector& rhs) const;

   /**
    * Returns true if each integer in vector is greater or equal to
    * corresponding integer in comparison vector.
    */
   bool
   operator >= (
      const IntVector& rhs) const;

   /**
    * Return the component-wise minimum of two integer vector objects.
    */
   void
   min(
      const IntVector& rhs);

   /**
    * Return the minimum entry in an integer vector.
    */
   int
   min() const;

   /**
    * Return the component-wise maximum of two integer vector objects.
    */
   void
   max(
      const IntVector& rhs);

   /**
    * Return the maximum entry in an integer vector.
    */
   int
   max() const;

   /**
    * Utility function to take the minimum of two integer vector objects.
    */
   static IntVector
   min(
      const IntVector& a,
      const IntVector& b);

   /**
    * Utility function to take the maximum of two integer vector objects.
    */
   static IntVector
   max(
      const IntVector& a,
      const IntVector& b);

   /**
    * Return the product of the entries in the integer vector.
    */
   int
   getProduct() const;

   /**
    * Store the object state to the specified database
    * with the provided name.
    *
    */
   virtual void
   putToDatabase(
      tbox::Database& database,
      const std::string& name) const;

   /**
    * Restores the object state from the specified database
    * with the provided name.
    *
    */
   virtual void
   getFromDatabase(
      tbox::Database& database,
      const std::string& name);

   /**
    * Return the dimension of this object.
    */
   const tbox::Dimension&
   getDim() const;

   /*!
    * @brief Return an IntVector of zeros of the specified dimension.
    *
    * Can be used to avoid object creation overheads.
    */
   static const IntVector&
   getZero(
      const tbox::Dimension& dim);

   /*!
    * @brief Return an IntVector of ones of the specified dimension.
    *
    * Can be used to avoid object creation overheads.
    */
   static const IntVector&
   getOne(
      const tbox::Dimension& dim);

   /**
    * Read an integer vector from an input stream.  The format for
    * the input is (i0,...,in) for an n-dimensional vector.
    */
   friend std::istream&
   operator >> (
      std::istream& s,
      IntVector& rhs);

   /**
    * Write an integer vector into an output stream.  The format for
    * the output is (i0,...,in) for an n-dimensional vector.
    */
   friend std::ostream&
   operator << (
      std::ostream& s,
      const IntVector& rhs);

   friend class std::vector<IntVector>;

protected:
   /**
    * Default ctor for IntVector is protected to disallow normal use.
    * This is needed by the poorly designed STL container library.
    *
    *
    */
   IntVector();

private:
   /*!
    * @brief Initialize static objects and register shutdown routine.
    *
    * Only called by StartupShutdownManager.
    *
    */
   static void
   initializeCallback();

   /*!
    * @brief Method registered with ShutdownRegister to cleanup statics.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   finalizeCallback();

   tbox::Dimension d_dim;

   int d_vector[SAMRAI::tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

   static IntVector* s_zeros[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
   static IntVector* s_ones[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

   static tbox::StartupShutdownManager::Handler
   s_initialize_finalize_handler;

};

}
}

#ifdef SAMRAI_INLINE
#include "SAMRAI/hier/IntVector.I"
#endif

#endif
