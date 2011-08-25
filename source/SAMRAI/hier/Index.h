/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Interface for the AMR Index object
 *
 ************************************************************************/

#ifndef included_hier_Index
#define included_hier_Index

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/tbox/Array.h"

namespace SAMRAI {
namespace hier {

/**
 * Class Index implements a simple n-dimensional integer vector in the
 * AMR index space.  Index is used as lower and upper bounds when
 * creating a box and also when iterating over the cells in a box.  An
 * index is essentially an integer vector but it carries along the
 * notion of indexing into AMR's abstract index space.
 *
 * @see hier::Box
 * @see hier::BoxIterator
 * @see hier::IntVector
 */

class Index:public IntVector
{
public:
   /**
    * Creates an uninitialized vector.
    */
   explicit Index(
      const tbox::Dimension& dim);

   /**
    * Construct an index with all components equal to the argument.
    */
   explicit Index(
      const tbox::Dimension& dim,
      const int i);

   /**
    * Construct a two-dimensional index with the value (i,j).
    */
   explicit Index(
      const int i,
      const int j);

   /**
    * Construct a three-dimensional index with the value (i,j,k).
    */
   explicit Index(
      const int i,
      const int j,
      const int k);

   /**
    * Construct an n-dimensional index with the values copied
    * from the integer tbox::Array i of size n.
    */
   explicit Index(
      const tbox::Array<int>& i);

   /**
    * The copy constructor creates an index equal to the argument.
    */
   Index(
      const Index& rhs);

   /**
    * Construct an index equal to the argument IntVector.
    */
   Index(
      const IntVector& rhs);

   /**
    * Construct an index equal to the argument array.
    */
   explicit Index(
      const tbox::Dimension& dim,
      const int array[]);

   /**
    * The assignment operator sets the index equal to the argument.
    *
    * An assignment to an uninitialized Index is allowed but assigning
    * from an uninitialized Index will result in an assert.
    */
   Index&
   operator = (
      const Index& rhs);

   /**
    * The assignment operator sets the index equal to the argument IntVector.
    */
   Index&
   operator = (
      const IntVector& rhs);

   /**
    * The index destructor does nothing interesting.
    */
   ~Index();

   /**
    * Plus-equals operator for an index and an integer vector.
    */
   Index&
   operator += (
      const IntVector& rhs);

   /**
    * Plus operator for an index and an integer vector.
    */
   Index
   operator + (
      const IntVector& rhs) const;

   /**
    * Plus-equals operator for an index and an integer.
    */
   Index&
   operator += (
      const int rhs);

   /**
    * Plus operator for an index and an integer.
    */
   Index
   operator + (
      const int rhs) const;

   /**
    * Minus-equals operator for an index and an integer vector.
    */
   Index&
   operator -= (
      const IntVector& rhs);

   /**
    * Minus operator for an index and an integer vector.
    */
   Index
   operator - (
      const IntVector& rhs) const;

   /**
    * Minus-equals operator for an index and an integer.
    */
   Index&
   operator -= (
      const int rhs);

   /**
    * Minus operator for an index and an integer.
    */
   Index
   operator - (
      const int rhs) const;

   /**
    * Times-equals operator for an index and an integer vector.
    */
   Index&
   operator *= (
      const IntVector& rhs);

   /**
    * Times operator for an index and an integer vector.
    */
   Index
   operator * (
      const IntVector& rhs) const;

   /**
    * Times-equals operator for an index and an integer.
    */
   Index&
   operator *= (
      const int rhs);

   /**
    * Times operator for an index and an integer.
    */
   Index
   operator * (
      const int rhs) const;

   /**
    * Assign-quotient operator for an index and an integer vector.
    */
   Index&
   operator /= (
      const IntVector& rhs);

   /**
    * Quotient operator for an index and an integer vector.
    */
   Index
   operator / (
      const IntVector& rhs) const;

   /**
    * Assign-quotient operator for an index and an integer.
    */
   Index&
   operator /= (
      const int rhs);

   /**
    * Quotient operator for an index and an integer.
    */
   Index
   operator / (
      const int rhs) const;

   /*!
    * @brief Coarsen the Index by a given ratio.
    *
    * For positive indices, this is the same as dividing by the ratio.
    */
   Index&
   coarsen(
      const IntVector& ratio);

   /*!
    * @brief Return an Index of zeros of the specified dimension.
    *
    * Can be used to avoid object creation overheads.
    */
   static const Index&
   getZeroIndex(
      const tbox::Dimension& dim);

   /*!
    * @brief Return an Index of ones of the specified dimension.
    *
    * Can be used to avoid object creation overheads.
    */
   static const Index&
   getOneIndex(
      const tbox::Dimension& dim);

   /*!
    * @brief Return an Index with minimum index values for the
    * specified dimension.
    *
    * Can be used to avoid object creation overheads.
    */
   static const Index&
   getMinIndex(
      const tbox::Dimension& dim);

   /*!
    * @brief Return an Index with maximum index values for the
    * specified dimension.
    *
    * Can be used to avoid object creation overheads.
    */
   static const Index&
   getMaxIndex(
      const tbox::Dimension& dim);

   /*!
    * @brief Coarsen an Index by a given ratio.
    *
    * For positive indices, this is the same as dividing by the ratio.
    */
   static Index
   coarsen(
      const Index& index,
      const IntVector& ratio);

private:
   friend class std::vector<Index>;

   /**
    * The default constructor for Index creates an uninitialized index.
    */
   Index();

   static int
   coarsen(
      const int index,
      const int ratio);

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
    *
    */
   static void
   finalizeCallback();

   static Index* s_zeros[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
   static Index* s_ones[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

   static Index* s_maxs[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
   static Index* s_mins[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

   static tbox::StartupShutdownManager::Handler
      s_initialize_finalize_handler;

};

}
}

#ifdef SAMRAI_INLINE
#include "SAMRAI/hier/Index.I"
#endif

#endif
