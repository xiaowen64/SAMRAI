/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Interface for the AMR Index object
 *
 ************************************************************************/

#ifndef included_hier_Index
#define included_hier_Index

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/Utilities.h"

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
    * @brief Creates an uninitialized vector.
    */
   explicit Index(
      const tbox::Dimension& dim);

   /**
    * @brief Construct an index with all components equal to the argument.
    */
   Index(
      const tbox::Dimension& dim,
      const int i);

   /**
    * @brief Construct a two-dimensional index with the value (i,j).
    */
   Index(
      const int i,
      const int j);

   /**
    * @brief Construct a three-dimensional index with the value (i,j,k).
    */
   Index(
      const int i,
      const int j,
      const int k);

   /**
    * @brief Construct an n-dimensional index with the values copied
    *        from the integer tbox::Array i of size n.
    */
   explicit Index(
      const tbox::Array<int>& i);

   /**
    * @brief The copy constructor creates an index equal to the argument.
    */
   Index(
      const Index& rhs);

   /**
    * @brief Construct an index equal to the argument IntVector.
    */
   explicit Index(
      const IntVector& rhs);

   /**
    * @brief Construct an index equal to the argument array.
    */
   Index(
      const tbox::Dimension& dim,
      const int array[]);

   /**
    * @brief The assignment operator sets the index equal to the argument.
    *
    * An assignment to an uninitialized Index is allowed but assigning
    * from an uninitialized Index will result in an assert.
    *
    * @pre getDim() == rhs.getDim()
    */
   Index&
   operator = (
      const Index& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      IntVector::operator = (rhs);
      return *this;
   }

   /**
    * @brief The assignment operator sets the index equal to the argument
    *        IntVector.
    *
    * @pre getDim() == rhs.getDim()
    */
   Index&
   operator = (
      const IntVector& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      IntVector::operator = (rhs);
      return *this;
   }

   /**
    * @brief The index destructor does nothing interesting.
    */
   virtual ~Index();

   /**
    * @brief Plus-equals operator for an index and an integer vector.
    *
    * @pre getDim() == rhs.getDim()
    */
   Index&
   operator += (
      const IntVector& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      IntVector::operator += (rhs);
      return *this;
   }

   /**
    * @brief Plus operator for an index and an integer vector.
    *
    * @pre getDim() == rhs.getDim()
    */
   Index
   operator + (
      const IntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      Index tmp = *this;
      tmp += rhs;
      return tmp;
   }

   /**
    * @brief Plus-equals operator for an index and an integer.
    */
   Index&
   operator += (
      const int rhs)
   {
      IntVector::operator += (rhs);
      return *this;
   }

   /**
    * @brief Plus operator for an index and an integer.
    */
   Index
   operator + (
      const int rhs) const
   {
      Index tmp = *this;
      tmp += rhs;
      return tmp;
   }

   /**
    * @brief Minus-equals operator for an index and an integer vector.
    *
    * @pre getDim() == rhs.getDim()
    */
   Index&
   operator -= (
      const IntVector& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      IntVector::operator -= (rhs);
      return *this;
   }

   /**
    * @brief Minus operator for an index and an integer vector.
    *
    * @pre getDim() == rhs.getDim()
    */
   Index
   operator - (
      const IntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      Index tmp = *this;
      tmp -= rhs;
      return tmp;
   }

   /**
    * @brief Minus-equals operator for an index and an integer.
    */
   Index&
   operator -= (
      const int rhs)
   {
      IntVector::operator -= (rhs);
      return *this;
   }

   /**
    * @brief Minus operator for an index and an integer.
    */
   Index
   operator - (
      const int rhs) const
   {
      Index tmp = *this;
      tmp -= rhs;
      return tmp;
   }

   /**
    * @brief Times-equals operator for an index and an integer vector.
    *
    * @pre getDim() == rhs.getDim()
    */
   Index&
   operator *= (
      const IntVector& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      IntVector::operator *= (rhs);
      return *this;
   }

   /**
    * @brief Times operator for an index and an integer vector.
    *
    * @pre getDim() == rhs.getDim()
    */
   Index
   operator * (
      const IntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      Index tmp = *this;
      tmp *= rhs;
      return tmp;
   }

   /**
    * @brief Times-equals operator for an index and an integer.
    */
   Index&
   operator *= (
      const int rhs)
   {
      IntVector::operator *= (rhs);
      return *this;
   }

   /**
    * @brief Times operator for an index and an integer.
    */
   Index
   operator * (
      const int rhs) const
   {
      Index tmp = *this;
      tmp *= rhs;
      return tmp;
   }

   /**
    * @brief Assign-quotient operator for an index and an integer vector.
    *
    * @pre getDim() == rhs.getDim()
    */
   Index&
   operator /= (
      const IntVector& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      IntVector::operator /= (rhs);
      return *this;
   }

   /**
    * @brief Quotient operator for an index and an integer vector.
    *
    * @pre getDim() == rhs.getDim()
    */
   Index
   operator / (
      const IntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      Index tmp = *this;
      tmp /= rhs;
      return tmp;
   }

   /**
    * @brief Assign-quotient operator for an index and an integer.
    */
   Index&
   operator /= (
      const int rhs)
   {
      IntVector::operator /= (rhs);
      return *this;
   }

   /**
    * @brief Quotient operator for an index and an integer.
    */
   Index
   operator / (
      const int rhs) const
   {
      Index tmp = *this;
      tmp /= rhs;
      return tmp;
   }

   /*!
    * @brief Coarsen the Index by a given ratio.
    *
    * For positive indices, this is the same as dividing by the ratio.
    *
    * @pre getDim() == ratio.getDim()
    */
   Index&
   coarsen(
      const IntVector& ratio)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, ratio);
      for (int d = 0; d < getDim().getValue(); ++d) {
         (*this)(d) = coarsen((*this)(d), ratio(d));
      }
      return *this;
   }

   /*!
    * @brief Return an Index of zeros of the specified dimension.
    *
    * Can be used to avoid object creation overheads.
    */
   static const Index&
   getZeroIndex(
      const tbox::Dimension& dim)
   {
      return *(s_zeros[dim.getValue() - 1]);
   }

   /*!
    * @brief Return an Index of ones of the specified dimension.
    *
    * Can be used to avoid object creation overheads.
    */
   static const Index&
   getOneIndex(
      const tbox::Dimension& dim)
   {
      return *(s_ones[dim.getValue() - 1]);
   }

   /*!
    * @brief Return an Index with minimum index values for the
    * specified dimension.
    *
    * Can be used to avoid object creation overheads.
    */
   static const Index&
   getMinIndex(
      const tbox::Dimension& dim)
   {
      return *(s_mins[dim.getValue() - 1]);
   }

   /*!
    * @brief Return an Index with maximum index values for the
    * specified dimension.
    *
    * Can be used to avoid object creation overheads.
    */
   static const Index&
   getMaxIndex(
      const tbox::Dimension& dim)
   {
      return *(s_maxs[dim.getValue() - 1]);
   }

   /*!
    * @brief Coarsen an Index by a given ratio.
    *
    * For positive indices, this is the same as dividing by the ratio.
    *
    * @pre index.getDim() == ratio.getDim()
    */
   static Index
   coarsen(
      const Index& index,
      const IntVector& ratio)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(index, ratio);
      tbox::Dimension dim(index.getDim());
      Index tmp(dim);
      for (int d = 0; d < dim.getValue(); ++d) {
         tmp(d) = coarsen(index(d), ratio(d));
      }
      return tmp;
   }

private:
   /*
    * Unimplemented default constructor
    */
   Index();

   static int
   coarsen(
      const int index,
      const int ratio)
   {
      return index < 0 ? (index + 1) / ratio - 1 : index / ratio;
   }

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

   static Index* s_zeros[SAMRAI::MAX_DIM_VAL];
   static Index* s_ones[SAMRAI::MAX_DIM_VAL];

   static Index* s_maxs[SAMRAI::MAX_DIM_VAL];
   static Index* s_mins[SAMRAI::MAX_DIM_VAL];

   static tbox::StartupShutdownManager::Handler
      s_initialize_finalize_handler;

};

}
}

#endif
