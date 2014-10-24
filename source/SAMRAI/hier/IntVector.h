/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2014 Lawrence Livermore National Security, LLC
 * Description:   A N-dimensional integer vector
 *
 ************************************************************************/

#ifndef included_hier_IntVector
#define included_hier_IntVector

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/BlockId.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Utilities.h"

#include <vector>
#include <iostream>

namespace SAMRAI {

namespace hier {

class Index;

/*!
 * @brief Simple integer vector class with size based on a dimension value
 *
 * Class IntVector implements a vector of integers that, depending on the usage
 * context, has a length equal to the number of spatial dimensions being used
 * or its length is the number of dimension multiplied by a certain number of
 * blocks.
 *
 * The number of blocks associated with an IntVector is called the "block
 * size".  When used in the context of a single-block problem, the block size
 * is always 1.  In a multiblock context, the block size of an IntVector
 * may be either 1 or the number of blocks being used in the problem.
 */

class IntVector
{
public:
   typedef tbox::Dimension::dir_t dir_t;

   /*!
    * @brief Creates an uninitialized IntVector of block size 1.
    *
    * @param dim
    */
   explicit IntVector(
      const tbox::Dimension& dim);

   /*!
    * @brief Creates an uninitialized IntVector of a given block size.
    *
    * @pre block_size >=1
    *
    * @param block_size
    * @param dim
    */
   IntVector(
      int block_size,
      const tbox::Dimension& dim);

   /*!
    * @brief Construct an IntVector with all components equal to the
    * value argument.
    *
    * @pre block_size >=1
    *
    * @param dim
    * @param value
    * @param block_size
    */
   IntVector(
      const tbox::Dimension& dim,
      int value,
      int block_size = 1);

   /*!
    * @brief Construct an IntVector with the values provided by
    * an STL vector of ints.
    *
    * The dimension of the constructed IntVector will be the size of the
    * vec argument.  If block_size has a value greater than 1, then the
    * IntVector will be constructed with the values held by vec duplicated for
    * every block.
    *
    * @pre vec.size() >= 1
    *
    * @param vec Vector of integers with a size equal to the desired dimension
    * @param block_size
    */
   IntVector(
      const std::vector<int>& vec,
      int block_size = 1);

   /*!
    * @brief Construct an IntVector with values provided by a raw array.
    *
    * This constructor assumes that the given array contains a number of values
    * equal to the dimension value.  As this constructor can do no error-
    * checking that the array is properly allocated and initialized, it is
    * up to the calling code to ensure that the array argument is valid.
    *
    * If block_size has a value greater than 1, then the IntVector
    * will be constructed with the values held by array duplicated for every
    * block.
    *
    * @param dim
    * @param array  Array of ints that should be allocated and initialized
    *               at a length equal to dim.getValue()
    * @param block_size 
    */
   IntVector(
      const tbox::Dimension& dim,
      const int array[],
      int block_size = 1);

   /*!
    * @brief Copy constructor.
    *
    * @pre rhs.getBlockSize() >= 1
    */
   IntVector(
      const IntVector& rhs);

   /*!
    * @brief Construct an IntVector from another IntVector.
    *
    * The main use case for this constructor is to use an IntVector of block
    * size 1 to construct an IntVector of a larger block size.  When used
    * in this way, the constructed IntVector will duplicate the contents of 
    * the argument for every block.
    * 
    * If block_size is equal to the rhs argument's block size, then this
    * constructor is equivalent to the copy constructor.
    *
    * @pre block_size >=1
    * @pre (rhs.getBlockSize() == block_size || rhs.getBlockSize() == 1)
    *
    * @param rhs
    * @param block_size
    */
   IntVector(
      const IntVector& rhs,
      int block_size);

   /*!
    * @brief Construct an IntVector from an Index.
    *
    * The constructed IntVector will have the same dimension value as the
    * Index.  If block_size is greater than 1, the values held by the Index
    * argument will be duplicated for every block.
    * 
    * @param rhs
    * @param block_size
    *
    * @pre block_size >=1
    *
    * @param rhs
    * @param block_size
    */
   IntVector(
      const Index& rhs,
      int block_size = 1);

   /*!
    * @brief The assignment operator sets the IntVector equal to the
    *        argument.
    *
    * @pre getDim() == rhs.getDim()
    */
   IntVector&
   operator = (
      const IntVector& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      d_block_size = rhs.d_block_size;
      d_vector = rhs.d_vector;

      return *this;
   }

   /*!
    * @brief Assignment operator assigning the values of an Index to an
    * IntVector.
    *
    * The assigned IntVector will have a block size of 1.
    *
    * @pre getDim() == rhs.getDim()
    */
   IntVector&
   operator = (
      const Index& rhs);

   /*!
    * @brief The IntVector destructor does nothing interesting.
    */
   virtual ~IntVector();

   /*!
    * @brief Return the block size of this IntVector
    */
   int getBlockSize() const
   {
      return d_block_size;
   }

   /*!
    * @brief Return an IntVector of block size 1 extracted from a possibly
    * larger IntVector
    *
    * This constructs an IntVector of block size 1 using the int values
    * associated with a given block.  The constructed IntVector is returned
    * by value.
    *
    * @pre block_id.getBlockValue() < getBlockSize()
    *
    * @param block_id  BlockId indicates which block is associated with
    *                  the desired integer values.
    *
    * @return A constructed IntVector of block size 1.
    */
   IntVector getBlockVector(const BlockId& block_id) const
   {
      TBOX_ASSERT(block_id.getBlockValue() < d_block_size);
      IntVector block_vec(d_dim,
                          &(d_vector[block_id.getBlockValue()*d_dim.getValue()]));

      return block_vec; 
   }

   /*!
    * @brief Return the specified component of the vector.
    *
    * @pre (i >= 0) && (i < getDim().getValue())
    * @pre getBlockSize() == 1
    */
   int&
   operator [] (
      const int i)
   {
      TBOX_ASSERT(i >= 0 && i < d_dim.getValue());
      TBOX_ASSERT(d_block_size == 1);
      return d_vector[i];
   }

   /*!
    * @brief Return the specified component of the vector as a const integer
    * reference.
    *
    * @pre (i >= 0) && (i < getDim().getValue())
    * @pre getBlockSize() == 1
    */
   const int&
   operator [] (
      const int i) const
   {
      TBOX_ASSERT(i >= 0 && i < d_dim.getValue());
      TBOX_ASSERT(d_block_size == 1);
      return d_vector[i];
   }

   /*!
    * @brief Return the specified component of the vector.
    *
    * @pre (i >= 0) && (i < getDim().getValue())
    * @pre getBlockSize() == 1
    */
   int&
   operator () (
      const int i)
   {
      TBOX_ASSERT(i >= 0 && i < d_dim.getValue());
      TBOX_ASSERT(d_block_size == 1);
      return d_vector[i];
   }

   /*!
    * @brief Return the specified component of the vector as a const integer
    * reference.
    *
    * @pre (i >= 0) && (i < getDim().getValue())
    * @pre getBlockSize() == 1
    */
   const int&
   operator () (
      const int i) const
   {
      TBOX_ASSERT(i >= 0 && i < d_dim.getValue());
      TBOX_ASSERT(d_block_size == 1);
      return d_vector[i];
   }

   /*!
    * @brief Return the specified component of the vector.
    *
    * @pre (b >= 0) && (b < getBlockSize())
    * @pre (i >= 0) && (i < getDim().getValue())
    *
    * The desired component is specified by the pair of the block number b
    * and the dimensional index i.
    */
   int&
   operator () (
      const int b, 
      const int i)
   {
      TBOX_ASSERT(b >= 0 && b < d_block_size);
      TBOX_ASSERT(i >= 0 && i < d_dim.getValue());
      return d_vector[b*d_dim.getValue() + i];
   }

   /*!
    * @brief Return the specified component of the vector as a const integer
    * reference 
    *
    * @pre (b >= 0) && (b < getBlockSize())
    * @pre (i >= 0) && (i < getDim().getValue())
    *
    * The desired component is specified by the pair of the block number b
    * and the dimensional index i.
    */
   const int&
   operator () (
      const int b, 
      const int i) const
   {
      TBOX_ASSERT(b >= 0 && b < d_block_size);
      TBOX_ASSERT(i >= 0 && i < d_dim.getValue());
      return d_vector[b*d_dim.getValue() + i];
   }

   //@{
   /*! @brief Arithmetic operators and related methods
    *
    * The arithmetic operators and related methods such as ceilingDivide
    * allow both full component-wise operations on IntVectors that have equal
    * block sizes, as well as operations between multiblock IntVectors and
    * IntVectors of block size 1.
    *
    * When an operator is invoked on a multiblock IntVector with a right-hand-
    * side IntVector of block size 1, the values within the right-hand-side are
    * applied to every block within the "this" IntVector.
    *
    * Assertion failures occur if operations are called using two multiblock
    * IntVectors that have different block sizes, or if the right-hand-side is
    * a multiblock IntVector while "this" has a block size of 1.
    *
    * There also are some operators that take a single integer as a right-
    * hand-side, in which case that integer value is applied in the operation
    * to every component of "this".
    */

   /*!
    * @brief Plus-equals operator for two integer vectors.
    *
    * @pre getDim() == rhs.getDim()
    * @pre getBlockSize() == rhs.getBlockSize() || rhs.getBlockSize() == 1
    */
   IntVector&
   operator += (
      const IntVector& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(d_block_size == rhs.d_block_size || rhs.d_block_size == 1);
      if (rhs.d_block_size == 1 && d_block_size != 1) {
         for (int b = 0; b < d_block_size; ++b) {
            int offset = b*d_dim.getValue();
            for (int i = 0; i < d_dim.getValue(); ++i) {
               d_vector[offset + i] += rhs.d_vector[i];
            }
         }
      } else {
         int length = d_block_size * d_dim.getValue();
         for (int i = 0; i < length; ++i) {
            d_vector[i] += rhs.d_vector[i];
         }
      }
      return *this;
   }

   /*!
    * @brief Plus operator for two integer vectors.
    *
    * @pre getDim() == rhs.getDim()
    */
   IntVector
   operator + (
      const IntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      IntVector tmp(*this);
      tmp += rhs;
      return tmp;
   }

   /*!
    * @brief Plus-equals operator for an integer vector and an integer.
    */
   IntVector&
   operator += (
      const int rhs)
   {
      int length = d_block_size * d_dim.getValue();
      for (int i = 0; i < length; ++i) {
         d_vector[i] += rhs;
      }
      return *this;
   }

   /*!
    * @brief Plus operator for an integer vector and an integer.
    */
   IntVector
   operator + (
      const int rhs) const
   {
      IntVector tmp(*this);
      tmp += rhs;
      return tmp;
   }

   /*!
    * @brief Minus-equals operator for two integer vectors.
    *
    * @pre getDim() == rhs.getDim()
    * @pre getBlockSize() == rhs.getBlockSize() || rhs.getBlockSize() == 1
    *
    */
   IntVector&
   operator -= (
      const IntVector& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(d_block_size == rhs.d_block_size || rhs.d_block_size == 1);
      if (rhs.d_block_size == 1 && d_block_size != 1) {
         for (int b = 0; b < d_block_size; ++b) {
            int offset = b*d_dim.getValue();
            for (int i = 0; i < d_dim.getValue(); ++i) {
               d_vector[offset + i] -= rhs.d_vector[i];
            }
         }
      } else {
         int length = d_block_size * d_dim.getValue();
         for (int i = 0; i < length; ++i) {
            d_vector[i] -= rhs.d_vector[i];
         }
      }
      return *this;
   }

   /*!
    * @brief Minus operator for two integer vectors.
    *
    * @pre getDim() == rhs.getDim()
    */
   IntVector
   operator - (
      const IntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      IntVector tmp(*this);
      tmp -= rhs;
      return tmp;
   }

   /*!
    * @brief Minus-equals operator for an integer vector and an integer.
    */
   IntVector&
   operator -= (
      const int rhs)
   {
      int length = d_block_size * d_dim.getValue();
      for (int i = 0; i < length; ++i) {
         d_vector[i] -= rhs;
      }
      return *this;
   }

   /*!
    * @brief Minus operator for an integer vector and an integer.
    */
   IntVector
   operator - (
      const int rhs) const
   {
      IntVector tmp(*this);
      tmp -= rhs;
      return tmp;
   }

   /*!
    * @brief Times-equals operator for two integer vectors.
    *
    * @pre getDim() == rhs.getDim()
    * @pre getBlockSize() == rhs.getBlockSize() || rhs.getBlockSize() == 1
    *
    */
   IntVector&
   operator *= (
      const IntVector& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(d_block_size == rhs.d_block_size || rhs.d_block_size == 1);
      if (rhs.d_block_size == 1 && d_block_size != 1) {
         for (int b = 0; b < d_block_size; ++b) {
            int offset = b*d_dim.getValue();
            for (int i = 0; i < d_dim.getValue(); ++i) {
               d_vector[offset + i] *= rhs.d_vector[i];
            }
         }
      } else {
         int length = d_block_size * d_dim.getValue();
         for (int i = 0; i < length; ++i) {
            d_vector[i] *= rhs.d_vector[i];
         }
      }
      return *this;
   }

   /*!
    * @brief Times operator for two integer vectors.
    *
    * @pre getDim() == rhs.getDim()
    */
   IntVector
   operator * (
      const IntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      IntVector tmp(*this);
      tmp *= rhs;
      return tmp;
   }

   /*!
    * @brief Times-equals operator for an integer vector and an integer.
    */
   IntVector&
   operator *= (
      const int rhs)
   {
      int length = d_block_size * d_dim.getValue();
      for (int i = 0; i < length; ++i) {
         d_vector[i] *= rhs;
      }
      return *this;
   }

   /*!
    * @brief Times operator for an integer vector and an integer.
    */
   IntVector
   operator * (
      const int rhs) const
   {
      IntVector tmp(*this);
      tmp *= rhs;
      return tmp;
   }

   /*!
    * @brief Quotient-equals operator for two integer vectors.
    *
    * @pre getDim() == rhs.getDim()
    * @pre getBlockSize() == rhs.getBlockSize() || rhs.getBlockSize() == 1
    */
   IntVector&
   operator /= (
      const IntVector& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(d_block_size == rhs.d_block_size || rhs.d_block_size == 1);
      if (rhs.d_block_size == 1 && d_block_size != 1) {
         for (int b = 0; b < d_block_size; ++b) {
            int offset = b*d_dim.getValue();
            for (int i = 0; i < d_dim.getValue(); ++i) {
               d_vector[offset + i] /= rhs.d_vector[i];
            }
         }
      } else {
         int length = d_block_size * d_dim.getValue();
         for (int i = 0; i < length; ++i) {
            d_vector[i] /= rhs.d_vector[i];
         }
      }
      return *this;
   }

   /*!
    * @brief Component-wise ceilingDivide quotient (integer divide with
    *        rounding up).
    *
    * @pre getDim() == denominator.getDim()
    * @pre getBlockSize() == denominator.getBlockSize() || denominator.getBlockSize() == 1
    */
   void
   ceilingDivide(
      const IntVector& denominator)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, denominator);
      TBOX_ASSERT(d_block_size == denominator.d_block_size ||
                  denominator.d_block_size == 1);

      /*
       * This is the formula for integer divide, rounding away from
       * zero.  It is meant as an extension of the ceilingDivide quotient of
       * 2 positive integers.
       *
       * The ceilingDivide is the integer divide plus 0, -1 or 1 representing
       * the results of rounding.
       * - Add zero if there's no remainder to round.
       * - Round remainder to 1 if numerator and denominator has same sign.
       * - Round remainder to -1 if numerator and denominator has opposite sign.
       */
      if (denominator.d_block_size == 1 && d_block_size != 1) {
         for (int b = 0; b < d_block_size; ++b) {
            int offset = b*d_dim.getValue();
            for (int i = 0; i < d_dim.getValue(); ++i) {
               d_vector[offset + i] = (d_vector[offset + i] / denominator[i]) +
               ((d_vector[offset + i] % denominator[i]) ?
                  ((d_vector[offset + i] > 0) == (denominator[i] > 0) ? 1 : -1) : 0);
            }
         }
      } else {
         int length = d_block_size * d_dim.getValue();
         for (int i = 0; i < length; ++i) {
            d_vector[i] = (d_vector[i] / denominator.d_vector[i]) +
               ((d_vector[i] % denominator.d_vector[i]) ?
               ((d_vector[i] > 0) == (denominator.d_vector[i] > 0) ? 1 : -1) : 0);
         }
      }
   }

   /*!
    * @brief Component-wise ceilingDivide quotient (integer divide with
    *        rounding up).
    *
    * @pre numerator.getDim() == denominator.getDim()
    * @pre numerator.getBlockSize() == denominator.getBlockSize() || denominator.getBlockSize() == 1
    */
   static IntVector
   ceilingDivide(
      const IntVector& numerator,
      const IntVector& denominator)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(numerator, denominator);
      TBOX_ASSERT(numerator.d_block_size == denominator.d_block_size ||
                  denominator.d_block_size == 1);
      IntVector rval(numerator);
      rval.ceilingDivide(denominator);
      return rval;
   }

   /*!
    * @brief Quotient operator for two integer vectors.
    *
    * @pre getDim() == rhs.getDim()
    */
   IntVector
   operator / (
      const IntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      IntVector tmp(*this);
      tmp /= rhs;
      return tmp;
   }

   /*!
    * @brief Quotient-equals operator for an integer vector and an integer.
    */
   IntVector&
   operator /= (
      const int rhs)
   {
      int length = d_block_size * d_dim.getValue();
      for (int i = 0; i < length; ++i) {
         d_vector[i] /= rhs;
      }
      return *this;
   }

   /*!
    * @brief Quotient operator for an integer vector and an integer.
    */
   IntVector
   operator / (
      const int rhs) const
   {
      IntVector tmp(*this);
      tmp /= rhs;
      return tmp;
   }

   /*!
    * @brief Unary minus to negate an integer vector.
    */
   IntVector
   operator - () const
   {
      IntVector tmp(*this);
      tmp *= -1; 
      return tmp;
   }

   /*!
    * @brief Returns true if all components are equal to a given integer.
    */
   bool
   operator == (
      int rhs) const
   {
      bool result = true;
      int length = d_block_size * d_dim.getValue();
      for (int i = 0; result && (i < length); ++i) {
         result = d_vector[i] == rhs;
      }
      return result;
   }

   /*!
    * @brief Returns true if some components are not equal to a given integer.
    */
   bool
   operator != (
      int rhs) const
   {
      return !(*this == rhs);
   }

   /*!
    * @brief Returns true if two vector objects are equal.  All components
    *        must be the same for equality.
    *
    * @pre getDim() == rhs.getDim()
    * @pre getBlockSize() == rhs.getBlockSize() || rhs.getBlockSize() == 1
    */
   bool
   operator == (
      const IntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(d_block_size == rhs.d_block_size || rhs.d_block_size == 1);
      bool result = true;
      if (rhs.d_block_size == 1 && d_block_size != 1) {
         for (int b = 0; result && b < d_block_size; ++b) {
            int offset = b*d_dim.getValue();
            for (int i = 0; result && (i < d_dim.getValue()); ++i) {
               result = result && (d_vector[offset + i] == rhs.d_vector[i]);
            }
         }
      } else {
         int length = d_block_size * d_dim.getValue();
         for (int i = 0; result && (i < length); ++i) {
            result = result && (d_vector[i] == rhs.d_vector[i]);
         }
      }
      return result;
   }

   /*!
    * @brief Returns true if two vector objects are not equal.  Any of
    *        the components may be different for inequality.
    *
    * @pre getDim() == rhs.getDim()
    */
   bool
   operator != (
      const IntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      return !(*this == rhs);
   }

   /*!
    * @brief Returns true if each integer in vector is less than
    *        corresponding integer in comparison vector.
    *
    * @pre getDim() == rhs.getDim()
    * @pre getBlockSize() == rhs.getBlockSize() || rhs.getBlockSize() == 1
    */
   bool
   operator < (
      const IntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(d_block_size == rhs.d_block_size || rhs.d_block_size == 1);
      bool result = true;
      if (rhs.d_block_size == 1 && d_block_size != 1) {
         for (int b = 0; result && b < d_block_size; ++b) {
            int offset = b*d_dim.getValue();
            for (int i = 0; result && (i < d_dim.getValue()); ++i) {
               result = result && (d_vector[offset + i] < rhs.d_vector[i]);
            }
         }
      } else {
         int length = d_block_size * d_dim.getValue();
         for (int i = 0; result && (i < length); ++i) {
            result = result && (d_vector[i] < rhs.d_vector[i]);
         }
      }
      return result;
   }

   /*!
    * @brief Returns true if each integer in vector is less or equal to
    *        corresponding integer in comparison vector.
    *
    * @pre getDim() == rhs.getDim()
    * @pre getBlockSize() == rhs.getBlockSize() || rhs.getBlockSize() == 1
    */
   bool
   operator <= (
      const IntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(d_block_size == rhs.d_block_size || rhs.d_block_size == 1);
      bool result = true;
      if (rhs.d_block_size == 1 && d_block_size != 1) {
         for (int b = 0; result && b < d_block_size; ++b) {
            int offset = b*d_dim.getValue();
            for (int i = 0; result && (i < d_dim.getValue()); ++i) {
               result = result && (d_vector[offset + i] <= rhs.d_vector[i]);
            }
         }
      } else {
         int length = d_block_size * d_dim.getValue();
         for (int i = 0; result && (i < length); ++i) {
            result = result && (d_vector[i] <= rhs.d_vector[i]);
         }
      }
      return result;
   }

   /*!
    * @brief Returns true if each integer in vector is greater than
    *        corresponding integer in comparison vector.
    *
    * @pre getDim() == rhs.getDim()
    * @pre getBlockSize() == rhs.getBlockSize() || rhs.getBlockSize() == 1
    */
   bool
   operator > (
      const IntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(d_block_size == rhs.d_block_size || rhs.d_block_size == 1);
      bool result = true;
      if (rhs.d_block_size == 1 && d_block_size != 1) {
         for (int b = 0; result && b < d_block_size; ++b) {
            int offset = b*d_dim.getValue();
            for (int i = 0; result && (i < d_dim.getValue()); ++i) {
               result = result && (d_vector[offset + i] > rhs.d_vector[i]);
            }
         }
      } else {
         int length = d_block_size * d_dim.getValue();
         for (int i = 0; result && (i < length); ++i) {
            result = result && (d_vector[i] > rhs.d_vector[i]);
         }
      }
      return result;
   }

   /*!
    * @brief Returns true if each integer in vector is greater or equal to
    *        corresponding integer in comparison vector.
    *
    * @pre getDim() == rhs.getDim()
    * @pre getBlockSize() == rhs.getBlockSize() || rhs.getBlockSize() == 1
    */
   bool
   operator >= (
      const IntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(d_block_size == rhs.d_block_size || rhs.d_block_size == 1);
      bool result = true;
      if (rhs.d_block_size == 1 && d_block_size != 1) {
         for (int b = 0; result && b < d_block_size; ++b) {
            int offset = b*d_dim.getValue();
            for (int i = 0; result && (i < d_dim.getValue()); ++i) {
               result = result && (d_vector[offset + i] >= rhs.d_vector[i]);
            }
         }
      } else {
         int length = d_block_size * d_dim.getValue();
         for (int i = 0; result && (i < length); ++i) {
            result = result && (d_vector[i] >= rhs.d_vector[i]);
         }
      }
      return result;
   }

   /*!
    * @brief Return the component-wise minimum of two integer vector objects.
    *
    * @pre getDim() == rhs.getDim()
    * @pre getBlockSize() == rhs.getBlockSize() || rhs.getBlockSize() == 1
    */
   void
   min(
      const IntVector& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(d_block_size == rhs.d_block_size || rhs.d_block_size == 1);
      if (rhs.d_block_size == 1 && d_block_size != 1) {
         for (int b = 0; b < d_block_size; ++b) {
            int offset = b*d_dim.getValue();
            for (int i = 0; i < d_dim.getValue(); ++i) {
               if (rhs.d_vector[i] < d_vector[offset + i]) {
                  d_vector[offset + i] = rhs.d_vector[i];
               }
            }
         }
      } else {
         int length = d_block_size * d_dim.getValue();
         for (int i = 0; i < length; ++i) {
            if (rhs.d_vector[i] < d_vector[i]) {
               d_vector[i] = rhs.d_vector[i];
            }
         }
      }
   }

   /*!
    * @brief Return the minimum entry in an integer vector.
    */
   int
   min() const
   {
      int min = d_vector[0];
      int length = d_block_size * d_dim.getValue();
      for (int i = 0; i < length; ++i) {
         if (d_vector[i] < min) {
            min = d_vector[i];
         }
      }
      return min;
   }

   /*!
    * @brief Return the component-wise maximum of two integer vector objects.
    */
   void
   max(
      const IntVector& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(d_block_size == rhs.d_block_size || rhs.d_block_size == 1);
      if (rhs.d_block_size == 1 && d_block_size != 1) {
         for (int b = 0; b < d_block_size; ++b) {
            int offset = b*d_dim.getValue();
            for (int i = 0; i < d_dim.getValue(); ++i) {
               if (rhs.d_vector[i] > d_vector[offset + i]) {
                  d_vector[offset + i] = rhs.d_vector[i];
               }
            }
         }
      } else {
         int length = d_block_size * d_dim.getValue();
         for (int i = 0; i < length; ++i) {
            if (rhs.d_vector[i] > d_vector[i]) {
               d_vector[i] = rhs.d_vector[i];
            }
         }
      }
   }

   //@}

   /*!
    * @brief Return the maximum entry in an integer vector.
    */
   int
   max() const
   {
      int max = d_vector[0];
      int length = d_block_size * d_dim.getValue();
      for (int i = 0; i < length; ++i) {
         if (d_vector[i] > max) {
            max = d_vector[i];
         }
      }
      return max;
   }

   /*!
    * @brief Utility function to take the minimum of two integer vector
    *        objects.
    *
    * @pre a.getDim() == b.getDim()
    */
   static IntVector
   min(
      const IntVector& a,
      const IntVector& b)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(a, b);
      IntVector tmp(a);
      tmp.min(b);
      return tmp;
   }

   /*!
    * @brief Utility function to take the maximum of two integer vector
    *        objects.
    *
    * @pre a.getDim() == b.getDim()
    */
   static IntVector
   max(
      const IntVector& a,
      const IntVector& b)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(a, b);
      IntVector tmp(a);
      tmp.max(b);
      return tmp;
   }

   /*!
    * @brief Set all block-wise components of an IntVector.
    *
    * If this IntVector and the argument IntVector are of the same block size,
    * this is the equivalent of a copy operation.  If this IntVector has a
    * block size greater than 1 while the argument is of block size 1, then
    * the values of the argument are copied to each block-wise component of
    * this IntVector.
    *
    * An error will occur if the block sizes are unequal and the argument
    * IntVector does not have a block size of 1.
    *
    * @param vector  Input IntVector
    */
   void setAll(const IntVector& vector)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, vector);
      if (vector.d_block_size == d_block_size) {
         *this = vector;
      } else if (d_block_size > 1 && vector.d_block_size == 1) {
         for (int b = 0; b < d_block_size; ++b) {
            int offset = b*d_dim.getValue();
            for (int d = 0; d < d_dim.getValue(); ++d) {
               d_vector[offset + d] = vector[d];
            }
         }
      } else {
         TBOX_ERROR("IntVector::setAll() attempted with argument of non-compatible block_size." << std::endl);
      }
   }

   /*!
    * @brief Return the product of the entries in the integer vector.
    *
    */
   long int
   getProduct() const
   {
      TBOX_ASSERT(d_block_size == 1);
      long int prod = 1;
      for (int i = 0; i < getDim().getValue(); ++i) {
         prod *= d_vector[i];
      }
      return prod;
   }

   /*!
    * @brief Store the object state to the specified restart database
    *        with the provided name.
    *
    */
   virtual void
   putToRestart(
      tbox::Database& restart_db,
      const std::string& name) const;

   /*!
    * @brief Restores the object giving it the provided name and getting its
    *        state from the specified restart database.
    *
    */
   virtual void
   getFromRestart(
      tbox::Database& restart_db,
      const std::string& name);

   /*!
    * @brief Return the dimension of this object.
    */
   const tbox::Dimension&
   getDim() const
   {
      return d_dim;
   }

   /*!
    * @brief Return an IntVector of zeros of the specified dimension.
    *
    * Can be used to avoid object creation overheads.  The block size of the
    * returned IntVector is 1.
    */
   static const IntVector&
   getZero(
      const tbox::Dimension& dim)
   {
      return *(s_zeros[dim.getValue() - 1]);
   }

   /*!
    * @brief Return an IntVector of zeros with a potentially larger block size.
    *
    * Can be used to avoid object creation overheads.  The block size of the
    * returned IntVector is set by the static method setNumberBlocks().
    */
   static const IntVector&
   getMultiZero(
      const tbox::Dimension& dim)
   {
      return *(s_mb_zeros[dim.getValue() - 1]);
   }

   /*!
    * @brief Return an IntVector of ones of the specified dimension.
    *
    * Can be used to avoid object creation overheads.  The block size of the
    * returned IntVector is 1.
    */
   static const IntVector&
   getOne(
      const tbox::Dimension& dim)
   {
      return *(s_ones[dim.getValue() - 1]);
   }

   /*!
    * @brief Return an IntVector of ones of the specified dimension.
    *
    * Can be used to avoid object creation overheads.  The block size of the
    * returned IntVector is set by the static method setNumberBlocks().
    */
   static const IntVector&
   getMultiOne(
      const tbox::Dimension& dim)
   {
      return *(s_mb_ones[dim.getValue() - 1]);
   }

   /*!
    * Set a static number of blocks for use in static IntVector methods
    */
   static void
   setNumberBlocks(int block_size)
   {
      for (unsigned short d = 0; d < SAMRAI::MAX_DIM_VAL; ++d) {
         tbox::Dimension dim(d+1);
         *(s_mb_zeros[d]) = IntVector(dim, 0, block_size);
         *(s_mb_ones[d]) = IntVector(dim, 1, block_size);
      }
   }

   /*!
    * @brief Sort the given IntVector by value.
    *
    * For an IntVector with block size 1, set the ith entry of this to the
    * position of the ith smallest value in the given IntVector.
    *
    * If block size is greater than 1, each section of the IntVector
    * associated with a block is sorted independently.
    */
   void
   sortIntVector(
      const IntVector& values);

   /*!
    * @brief Read an integer vector from an input stream.  The format for
    *        the input is (i0,...,in) for an n-dimensional vector.
    */
   friend std::istream&
   operator >> (
      std::istream& s,
      IntVector& rhs);

   /*!
    * @brief Write an integer vector into an output stream.  The format for
    *        the output is (i0,...,in) for an n-dimensional vector.
    */
   friend std::ostream&
   operator << (
      std::ostream& s,
      const IntVector& rhs);

private:
   /*
    * Unimplemented default constructor
    */
   IntVector();

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

   int d_block_size;

   std::vector<int> d_vector;

   static IntVector* s_zeros[SAMRAI::MAX_DIM_VAL];
   static IntVector* s_ones[SAMRAI::MAX_DIM_VAL];
   static IntVector* s_mb_zeros[SAMRAI::MAX_DIM_VAL];
   static IntVector* s_mb_ones[SAMRAI::MAX_DIM_VAL];

   static tbox::StartupShutdownManager::Handler
      s_initialize_finalize_handler;

};

}
}

#endif
