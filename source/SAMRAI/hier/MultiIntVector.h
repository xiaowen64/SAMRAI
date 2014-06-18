/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2014 Lawrence Livermore National Security, LLC
 * Description:  Wrapper for IntVector to be used in multiblock context
 *
 ************************************************************************/

#ifndef included_hier_MultiIntVector 
#define included_hier_MultiIntVector 

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/BlockId.h"
#include "SAMRAI/hier/IntVector.h"

/*!
 * @brief A container of multiple IntVectors, intended to hold one for
 * each block.
 *
 * MultiIntVector is a container of IntVectors, with the intended usage
 * being that it will hold one IntVector for every block in the mesh.
 * It is useful for containing N-dimensional information that may vary on
 * different blocks.
 */

namespace SAMRAI {
namespace hier {
#if 0
class MultiIntVector
{
public:

   /*!
    * @brief Set a default maximum size for the object.
    *
    * This sets a static value which will be used as the default size for
    * future MultiIntVector objects.  It is recommended that this be called
    * as soon as the number of blocks in the mesh is known.
    *
    * @param number_blocks
    *
    * @pre number_blocks > 0
    */
   static void setNumberBlocks(int number_blocks) {
      TBOX_ASSERT(number_blocks > 0);
      s_num_blocks = number_blocks;
   }

   /*!
    * @brief Construct a MultiIntVector from vector of IntVector.
    *
    * The MultiIntVector will be constructed to contain the IntVectors in
    * the argument vector.  Each IntVector will be associated with a BlockId
    * according to the ordering of the argument vector.
    *
    * If the default size has not already been set via setNumberBlocks(),
    * it will be set to the size of the argument vector;
    *
    * @param vector
    *
    * @pre !vector.empty();
    * @pre s_num_blocks == 0 || s_num_blocks == d_vector.size()
    *
    */
//   explicit MultiIntVector(
//      const std::vector<IntVector>& vector);

   /*!
    * @brief Construct a MultiIntVector from an IntVector.
    *
    * The MultiIntVector will be constructed to hold the same IntVector for
    * all blocks.
    *
    * @param vector
    *
    * @pre s_num_blocks > 0
    */
//   explicit MultiIntVector(
//      const IntVector& vector); 

   /*!
    * @brief Construct a constant MultiIntVector from a dimension and value
    *
    * All entries in the MultiIntVector will have the given value.
    *
    * @param dim
    * @param value
    *
    * @pre s_num_blocks > 0
    */
//   MultiIntVector(
//      const tbox::Dimension& dim,
//      int value);

   /*!
    * @brief Copy constructor 
    *
    * @param copy_vector
    */
//   MultiIntVector(
//      const MultiIntVector& copy_vector);

   /*!
    * @brief Assignment operator
    *
    * @param rhs
    *
    * @pre getDim() == rhs.getDim()
    */
   MultiIntVector&
   operator = (
      const MultiIntVector& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      d_vector = rhs.d_vector;
      return *this;
   }

   /*!
    * @brief Destructor
    */
   virtual ~MultiIntVector();

   /*!
    * @brief Clear the object so that it holds no IntVector data.
    */
   void clear()
   {
      d_vector.clear();
   }

   /*!
    * @brief Check if the obect is empty.
    */
   bool empty() const
   {
      return d_vector.empty();
   }

   /*!
    * @brief Set the object to hold the contents of the argument vector
    *
    * The size of the argument vector must equal the default size.
    * Any data previously held by this object will be discarded.
    *
    * @param vector
    *
    * @pre s_num_blocks > 0
    * @pre s_num_blocks == d_vector.size()
    * @pre getDim() == vector[0].getDim()
    */
   void set(const std::vector<IntVector>& vector)
   {
      TBOX_ASSERT(s_num_blocks > 0);
      TBOX_ASSERT(s_num_blocks == d_vector.size());
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, vector[0]);
      d_vector = vector;
   }

   /*!
    * @brief Set the object to hold the same IntVector for every block.
    *
    * Any data previously held by this object will be discarded.
    *
    * @param vector
    *
    * @pre s_num_blocks > 0
    * @pre getDim() == vector.getDim()
    */
   void setAll(const IntVector& vector)
   {
      TBOX_ASSERT(s_num_blocks > 0);
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, vector);
      d_vector.resize(s_num_blocks, vector);
      for (int b = 0; b < d_vector.size(); ++b) {
         d_vector[b] = vector;
      }
   }

   /*!
    * @brief Get the dimension of the object.
    */
   const tbox::Dimension& getDim() const
   {
      return (d_dim);
   }

   /*!
    * @brief Check if the MultiIntVector has the value 1 in all entries.
    */
   bool isOne() const {
      for (int b = 0; b < d_vector.size(); ++b) {
         if (d_vector[b] != IntVector::getOne(d_dim)) {
            return false;
         }
      }
      return true;
   }

   /*!
    * @brief Check if the MultiIntVector has the value 0 in all entries.
    */
   bool isZero() const {
      for (int b = 0; b < d_vector.size(); ++b) {
         if (d_vector[b] != IntVector::getZero(d_dim)) {
            return false;
         }
      }
      return true;
   }

   /*!
    * @brief Access a const reference to an IntVector for a certain block.
    *
    * Provides access to the IntVector associated with the block identified
    * by the given BlockId.
    *
    * @param block_id
    */ 
   const IntVector& getBlockVector(const BlockId& block_id) const
   {
      TBOX_ASSERT(block_id.getBlockValue() < d_vector.size());
      return (d_vector[block_id.getBlockValue()]);
   }  


   /*!
    * @brief Set this object to the component-wise minimum of itself
    * and another MultiIntVector.
    *
    * Compare each entry in two MultiIntVectors and set this object's entry
    * to the minimum of the two values.
    *
    * @param rhs
    */
   void
   min(const MultiIntVector& rhs) 
   {
      TBOX_ASSERT(d_vector.size() == rhs.d_vector.size());
      for (int b = 0; b < d_vector.size(); ++b) {
         d_vector[b].min(rhs.d_vector[b]);
      }
   }

   /*!
    * @brief Return the minimum of all the values held in this MultiIntVector.
    */
   int
   min() const
   {
      int min = d_vector[0][0];
      for (int b = 0; b < d_vector.size(); ++b) {
         for (int i = 1; i < getDim().getValue(); ++i) {
            if (d_vector[b][i] > min) {
               min = d_vector[b][i];
            }
         }
      }
      return min;
   }


   /*!
    * @brief Utility function to compute the minimum of two MultiIntVectors.
    *
    * @param a
    * @param b
    * @return Component-wise minimum of a and b
    * 
    * @pre a.getDim() == b.getDim()
    */
   static MultiIntVector
   min(const MultiIntVector& a,
       const MultiIntVector& b)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(a, b);
      MultiIntVector tmp = a;
      tmp.min(b);
      return tmp;
   }

   /*!
    * @brief Set this object to the component-wise maximum of itself
    * and another MultiIntVector.
    *
    * Compare each entry in two MultiIntVectors and set this object's entry
    * to the maximum of the two values.
    *
    * @param rhs
    */
   void
   max(const MultiIntVector& rhs)
   {
      TBOX_ASSERT(d_vector.size() == rhs.d_vector.size());
      for (int b = 0; b < d_vector.size(); ++b) {
         d_vector[b].max(rhs.d_vector[b]);
      }
   }

   /*!
    * @brief Return the maximum of all the values held in this MultiIntVector.
    */
   int
   max() const
   {
      int max = d_vector[0][0];
      for (int b = 0; b < d_vector.size(); ++b) {
         for (int i = 1; i < getDim().getValue(); ++i) {
            if (d_vector[b][i] > max) {
               max = d_vector[b][i];
            }
         }
      }
      return max;
   }

   /*!
    * @brief Utility function to compute the maximum of two MultiIntVectors.
    *
    * @param a 
    * @param b
    * @return Component-wise maximum of a and b.
    *
    * @pre a.getDim() == b.getDim()
    */
   static MultiIntVector
   max(
      const MultiIntVector& a,
      const MultiIntVector& b)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(a, b);
      MultiIntVector tmp = a;
      tmp.max(b);
      return tmp;
   }

   /*!
    * @brief Plus operator for two MultiIntVectors.
    *
    * @param rhs
    * @return Component-wise sum of this and rhs
    *
    * @pre getDim() == rhs.getDim()
    * @pre d_vector.size() == rhs.d_vector.size()
    */
   MultiIntVector
   operator + (
      const MultiIntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(d_vector.size() == rhs.d_vector.size());
      MultiIntVector tmp = *this;
      tmp += rhs;
      return tmp;
   }

   /*!
    * @brief Plus-equals operator for two MultiIntVectors.
    *
    * @param rhs
    * @return Reference to this object after plus-equal operation.
    * 
    * @pre getDim() == rhs.getDim()
    */
   MultiIntVector&
   operator += (
      const MultiIntVector& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(d_vector.size() == rhs.d_vector.size());
      for (int b = 0; b < d_vector.size(); ++b) {
         d_vector[b] += rhs.d_vector[b];
      }
      return *this;
   }

   /*!
    * @brief Minus operator for two MultiIntVectors.
    *
    * @param rhs
    * @return Component-wise difference of this and rhs
    *
    * @pre getDim() == rhs.getDim()
    * @pre d_vector.size() == rhs.d_vector.size()
    */
   MultiIntVector
   operator - (
      const MultiIntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(d_vector.size() == rhs.d_vector.size());
      MultiIntVector tmp = *this;
      tmp -= rhs;
      return tmp;
   }

   /*!
    * @brief Minus-equals operator for two MultiIntVectors.
    *
    * @param rhs
    * @return Reference to this object after minus-equal operation.
    *
    * @pre getDim() == rhs.getDim()
    */
   MultiIntVector&
   operator -= (
      const MultiIntVector& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(d_vector.size() == rhs.d_vector.size());
      for (int b = 0; b < d_vector.size(); ++b) {
         d_vector[b] -= rhs.d_vector[b];
      }
      return *this;
   }

   /*!
    * @brief Unary minus to negate a MultiIntVector
    *
    * @return Component-wise negation of this object.
    */
   MultiIntVector
   operator - () const
   {
      MultiIntVector tmp(*this);
      for (int b = 0; b < d_vector.size(); ++b) {
         tmp.d_vector[b] = -d_vector[b];
      }
      return tmp;
   }

   /*!
    * @brief Times operator for two MultiIntVectors.
    *
    * @param rhs
    * @return Component-wise product of this and rhs
    *
    * @pre getDim() == rhs.getDim()
    * @pre d_vector.size() == rhs.d_vector.size()
    */
   MultiIntVector
   operator * (
      const MultiIntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(d_vector.size() == rhs.d_vector.size());
      MultiIntVector tmp = *this;
      tmp *= rhs;
      return tmp;
   }

   /*!
    * @brief Times-equals operator for two MultiIntVectors.
    *
    * @param rhs
    * @return Reference to this object after times-equal operation.
    *
    * @pre getDim() == rhs.getDim()
    */
   MultiIntVector&
   operator *= (
      const MultiIntVector& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(d_vector.size() == rhs.d_vector.size());
      for (int b = 0; b < d_vector.size(); ++b) {
         d_vector[b] *= rhs.d_vector[b];
      }
      return *this;
   }

   /*!
    * @brief Times operator for MultiIntVector and IntVector
    *
    * Each IntVector held by this MultiIntVector is multiplied by the same
    * right-hand-side IntVector.
    *
    * @param rhs
    * @return Copy of this MultiIntVector multiplied by rhs in every entry
    *
    * @pre getDim() == rhs.getDim()
    */
   MultiIntVector
   operator * (
      const IntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      MultiIntVector tmp = *this;
      for (int b = 0; b < d_vector.size(); ++b) {
         tmp.d_vector[b] *= rhs;
      }
      return tmp;
   }

   /*!
    * @brief Times operator for MultiIntVector and integer scalar
    *
    * Each entry of this MultiIntVector is multiplied by the same
    * right-hand-side scalar integer.
    *
    * @param rhs
    * @return Copy of this MultiIntVector multiplied by rhs in every entry
    */
   MultiIntVector
   operator * (
      const int& rhs) const
   {
      MultiIntVector tmp = *this;
      for (int b = 0; b < d_vector.size(); ++b) {
         tmp.d_vector[b] *= rhs;
      }
      return tmp;
   }

   /*!
    * @brief Division operator for two MultiIntVectors.
    *
    * @param rhs
    * @return Component-wise integer quotient of this and rhs
    *
    * @pre getDim() == rhs.getDim()
    * @pre d_vector.size() == rhs.d_vector.size()
    */
   MultiIntVector
   operator / (
      const MultiIntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(d_vector.size() == rhs.d_vector.size());
      MultiIntVector tmp = *this;
      tmp /= rhs;
      return tmp;
   }

   /*!
    * @brief Divide-equals operator for two MultiIntVectors.
    *
    * @param rhs
    * @return Reference to this object after divide-equal operation.
    *
    * @pre getDim() == rhs.getDim()
    */
   MultiIntVector&
   operator /= (
      const MultiIntVector rhs)
   {
      TBOX_ASSERT(d_vector.size() == rhs.d_vector.size());
      for (int b = 0; b < d_vector.size(); ++b) {
         d_vector[b] /= rhs.d_vector[b];
      }
      return *this;
   }

   /*!
    * @brief Division operator for MultiIntVector divided by IntVector
    *
    * Each IntVector held by this MultiIntVector is divided by the same
    * right-hand-side IntVector.
    *
    * @param rhs
    * @return Copy of this MultiIntVector divided by rhs in every entry
    *
    * @pre getDim() == rhs.getDim()
    */
   MultiIntVector
   operator / (
      const IntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      MultiIntVector tmp = *this;
      for (int b = 0; b < d_vector.size(); ++b) {
         tmp.d_vector[b] /= rhs;
      }
      return tmp;
   }

   /*!
    * @brief Division operator for MultiIntVector divided by integer scalar
    *
    * Each entry of this MultiIntVector is divided by the same
    * right-hand-side scalar integer.
    *
    * @param rhs
    * @return Copy of this MultiIntVector divided by rhs in every entry
    *
    * @pre rhs != 0
    */
   MultiIntVector
   operator / (
      const int& rhs) const
   {
      TBOX_ASSERT(rhs != 0);
      MultiIntVector tmp = *this;
      for (int b = 0; b < d_vector.size(); ++b) {
         tmp.d_vector[b] /= rhs;
      }
      return tmp;
   }

   /*!
    * @brief Internal component-wise ceilingDivide quotient (integer divide
    *        with rounding up).
    *
    * Each entry in this object will be set to the ceiling-divide quotient
    * of itself and the corresponding entry in rhs
    *
    * @param rhs
    *
    * @pre d_vector.size() == rhs.d_vector.size()
    */
   void
   ceilingDivide(const MultiIntVector& rhs)
   {
      TBOX_ASSERT(d_vector.size() == rhs.d_vector.size());
      for (int b = 0; b < d_vector.size(); ++b) {
         d_vector[b].ceilingDivide(rhs.d_vector[b]);
      }
   }

   /*!
    * @brief Utility function for ceilingDivide quotient (integer divide
    *        with rounding up).
    *
    * Each entry in this object will be set to the ceiling-divide quotient
    * of itself and the corresponding entry in rhs
    *
    * @param numerator
    * @param denominator
    * @return Ceiling-divide quotient of (numerator / denominator)
    */
   static MultiIntVector
   ceilingDivide(
      const MultiIntVector& numerator,
      const MultiIntVector& denominator)
   {
      MultiIntVector tmp_num(numerator);
      tmp_num.ceilingDivide(denominator);
      return tmp_num;
   }

   /*!
    * @brief Returns true if each integer in this object is less than or equal
    * to the corresponding integer in comparison vector.
    *
    * @param rhs
    *
    * @pre getDim() == rhs.getDim()
    * @pre d_vector.size() == rhs.d_vector.size()
    */
   bool
   operator <= (
      const MultiIntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(d_vector.size() == rhs.d_vector.size());
      bool result = true;
      for (int b = 0; result && (b < d_vector.size()); ++b) {
         result = result && (d_vector[b] <= rhs.d_vector[b]);
      }
      return result;
   }

   /*!
    * @brief Returns true if each integer in this object is less than the
    * corresponding integer in comparison vector.
    *
    * @pre getDim() == rhs.getDim()
    * @pre d_vector.size() == rhs.d_vector.size()
    */
   bool
   operator < (
      const MultiIntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(d_vector.size() == rhs.d_vector.size());
      bool result = true;
      for (int b = 0; result && (b < d_vector.size()); ++b) {
         result = result && (d_vector[b] < rhs.d_vector[b]);
      }
      return result;
   }

   /*!
    * @brief Returns true if each integer in this object is less than or
    * equal to the corresponding integer in comparison IntVector.
    *
    * Each IntVector held by this object will be compared to the same rhs
    * IntVector.
    *
    * @param rhs
    *
    * @pre getDim() == rhs.getDim()
    */
   bool
   operator <= (
      const IntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      bool result = true;
      for (int b = 0; result && (b < d_vector.size()); ++b) {
         result = result && (d_vector[b] <= rhs);
      }
      return result;
   }

   /*!
    * @brief Returns true if each integer in this object is less than the
    * corresponding integer in comparison IntVector.
    *
    * Each IntVector held by this object will be compared to the same rhs
    * IntVector.
    *
    * @pre getDim() == rhs.getDim()
    */
   bool
   operator < (
      const IntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      bool result = true;
      for (int b = 0; result && (b < d_vector.size()); ++b) {
         result = result && (d_vector[b] < rhs);
      }
      return result;
   }

   /*!
    * @brief Returns true if each integer in this object is greater than or
    * equal to the corresponding integer in comparison vector.
    *
    * @param rhs
    *
    * @pre getDim() == rhs.getDim()
    * @pre d_vector.size() == rhs.d_vector.size()
    */
   bool
   operator >= (
      const MultiIntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(d_vector.size() == rhs.d_vector.size());
      bool result = true;
      for (int b = 0; result && (b < d_vector.size()); ++b) {
         result = result && (d_vector[b] >= rhs.d_vector[b]);
      }
      return result;
   }

   /*!
    * @brief Returns true if each integer in this object is greater than the
    * corresponding integer in comparison vector.
    *
    * @pre getDim() == rhs.getDim()
    * @pre d_vector.size() == rhs.d_vector.size()
    */
   bool
   operator > (
      const MultiIntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(d_vector.size() == rhs.d_vector.size());
      bool result = true;
      for (int b = 0; result && (b < d_vector.size()); ++b) {
         result = result && (d_vector[b] > rhs.d_vector[b]);
      }
      return result;
   }

   /*!
    * @brief Returns true if each integer in this object is greater than or
    * equal to the corresponding integer in comparison IntVector.
    *
    * Each IntVector held by this object will be compared to the same rhs
    * IntVector.
    *
    * @param rhs
    *
    * @pre getDim() == rhs.getDim()
    */
   bool
   operator >= (
      const IntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      bool result = true;
      for (int b = 0; result && (b < d_vector.size()); ++b) {
         result = result && (d_vector[b] >= rhs);
      }
      return result;
   }

   /*!
    * @brief Returns true if each integer in this object is greater than
    * the corresponding integer in comparison IntVector.
    *
    * Each IntVector held by this object will be compared to the same rhs
    * IntVector.
    *
    * @param rhs
    *
    * @pre getDim() == rhs.getDim()
    */
   bool
   operator > (
      const IntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      bool result = true;
      for (int b = 0; result && (b < d_vector.size()); ++b) {
         result = result && (d_vector[b] > rhs);
      }
      return result;
   }

   /*!
    * @brief Returns true if two objects are equal.  All components
    *        must be the same for equality.
    *
    * @param rhs 
    *
    * @pre getDim() == rhs.getDim()
    */
   bool
   operator == (
      const MultiIntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      bool result = (d_vector.size() == rhs.d_vector.size());
      for (int b = 0; result && (b < d_vector.size()); ++b) {
         result = result && (d_vector[b] == rhs.d_vector[b]);
      }
      return result;
   }

   /*!
    * @brief Returns true if two vector objects are not equal.  Any of
    * the components may be different for inequality.
    *
    * @pre getDim() == rhs.getDim()
    */
   bool
   operator != (
      const MultiIntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      return !(*this == rhs);
   }

   /*!
    * @brief Read an integer vector from an input stream.  The format for
    *        the input is (i0,...,in) for an n-dimensional vector.
    */
   friend std::istream&
   operator >> (
      std::istream& s,
      MultiIntVector& rhs);

   /*!
    * @brief Write an integer vector into an output stream.  The format for
    *        the output is (i0,...,in) for an n-dimensional vector.
    */
   friend std::ostream&
   operator << (
      std::ostream& s,
      const MultiIntVector& rhs);



private:
   /**
    * Unimplemented default constructor.
    */
   MultiIntVector();

   static int s_num_blocks;

   std::vector<IntVector> d_vector;
   tbox::Dimension d_dim;

};

#endif
}
}


#endif
