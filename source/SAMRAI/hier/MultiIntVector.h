/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2014 Lawrence Livermore National Security, LLC
 * Description:   Ratio representing a portion of the AMR index space
 *
 ************************************************************************/

#ifndef included_hier_Ratio 
#define included_hier_Ratio 

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/BlockId.h"
#include "SAMRAI/hier/IntVector.h"

namespace SAMRAI {
namespace hier {

class MultiIntVector
{
public:

//   MultiIntVector(
//      const IntVector& ratio,
//      const BlockId& block_id);
   static void setNumberBlocks(int number_blocks) {
      s_max_blocks = number_blocks;
   }

   explicit MultiIntVector(
      const IntVector& ratio); 

   explicit MultiIntVector(
      const std::vector<IntVector>& ratio);

   MultiIntVector(
      const tbox::Dimension& dim,
      int value);

   MultiIntVector(
      const tbox::Dimension& dim,
      int value,
      int nblocks);

   /*!
    * @brief Copy constructor 
    */
   MultiIntVector(
      const MultiIntVector& box);

   virtual ~MultiIntVector();

   void clear()
   {
      d_ratio.clear();
   }

   bool empty() const
   {
      return d_ratio.empty();
   }

   void set(const std::vector<IntVector>& ratio)
   {
      clear();
      d_ratio = ratio;
      if (d_ratio.size() > s_max_blocks) {
         s_max_blocks = d_ratio.size();
      } 
   }

   void setAll(const IntVector& ratio)
   {
      TBOX_ASSERT(s_max_blocks >= 1);
      d_ratio.resize(s_max_blocks, ratio);
      for (int b = 0; b < d_ratio.size(); ++b) {
         d_ratio[b] = ratio;
      }
   }

   const tbox::Dimension& getDim() const
   {
      TBOX_ASSERT(d_ratio.size() > 0);

      return (d_ratio[0].getDim());
   }

   bool isOne(const tbox::Dimension& dim) const {
      for (int b = 0; b < d_ratio.size(); ++b) {
         if (d_ratio[b] != IntVector::getOne(dim)) {
            return false;
         }
      }
      return true;
   }

   bool isZero(const tbox::Dimension& dim) const {
      for (int b = 0; b < d_ratio.size(); ++b) {
         if (d_ratio[b] != IntVector::getZero(dim)) {
            return false;
         }
      }
      return true;
   }

   const IntVector& getBlockVector(const BlockId& block_id) const
   {
      TBOX_ASSERT(block_id.getBlockValue() < d_ratio.size());
      return (d_ratio[block_id.getBlockValue()]);
   }  


   /**
    * @brief Return the component-wise minimum of two integer vector objects.
    */
   void
   min(
      const MultiIntVector& rhs)
   {
      for (int b = 0; b < d_ratio.size(); ++b) {
         d_ratio[b].min(rhs.d_ratio[b]);
      }
   }


   /**
    * @brief Return the minimum entry in an integer vector.
    */
   int
   min() const
   {
      int min = d_ratio[0][0];
      for (int b = 0; b < d_ratio.size(); ++b) {
         for (int i = 1; i < getDim().getValue(); ++i) {
            if (d_ratio[b][i] > min) {
               min = d_ratio[b][i];
            }
         }
      }
      return min;
   }


   /**
    * @brief Utility function to take the minimum of two integer vector
    *        objects.
    *
    * @pre a.getDim() == b.getDim()
    */
   static MultiIntVector
   min(
      const MultiIntVector& a,
      const MultiIntVector& b)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(a, b);
      MultiIntVector tmp = a;
      tmp.min(b);
      return tmp;
   }

   /**
    * @brief Return the component-wise maximum of two integer vector objects.
    */
   void
   max(
      const MultiIntVector& rhs)
   {
      for (int b = 0; b < d_ratio.size(); ++b) {
         d_ratio[b].max(rhs.d_ratio[b]);
      }
   }


   /**
    * @brief Return the maximum entry in an integer vector.
    */
   int
   max() const
   {
      int max = d_ratio[0][0];
      for (int b = 0; b < d_ratio.size(); ++b) {
         for (int i = 1; i < getDim().getValue(); ++i) {
            if (d_ratio[b][i] > max) {
               max = d_ratio[b][i];
            }
         }
      }
      return max;
   }


   /**
    * @brief Utility function to take the maximum of two integer vector
    *        objects.
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



   /**
    * @brief Plus-equals operator for two integer vectors.
    *
    * @pre getDim() == rhs.getDim()
    */
   MultiIntVector&
   operator += (
      const MultiIntVector& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      for (int b = 0; b < d_ratio.size(); ++b) {
         d_ratio[b] += rhs.d_ratio[b];
      }
      return *this;
   }

   /**
    * @brief Plus operator for two integer vectors.
    *
    * @pre getDim() == rhs.getDim()
    */
   MultiIntVector
   operator + (
      const MultiIntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      MultiIntVector tmp = *this;
      tmp += rhs;
      return tmp;
   }


   /**
    * @brief Minus-equals operator for two integer vectors.
    *
    * @pre getDim() == rhs.getDim()
    */
   MultiIntVector&
   operator -= (
      const MultiIntVector& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      for (int b = 0; b < d_ratio.size(); ++b) {
         d_ratio[b] -= rhs.d_ratio[b];
      }
      return *this;
   }

   /**
    * @brief Minus operator for two integer vectors.
    *
    * @pre getDim() == rhs.getDim()
    */
   MultiIntVector
   operator - (
      const MultiIntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      MultiIntVector tmp = *this;
      tmp -= rhs;
      return tmp;
   }

   /**
    * @brief Unary minus to negate an integer vector.
    */
   MultiIntVector
   operator - () const
   {
      MultiIntVector tmp(*this);
      for (int b = 0; b < d_ratio.size(); ++b) {
         tmp.d_ratio[b] = -d_ratio[b];
      }
      return tmp;
   }



   /**
    * @brief Times-equals operator for two integer vectors.
    *
    * @pre getDim() == rhs.getDim()
    */
   MultiIntVector&
   operator *= (
      const MultiIntVector& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      for (int b = 0; b < d_ratio.size(); ++b) {
         d_ratio[b] *= rhs.d_ratio[b];
      }
      return *this;
   }

   /**
    * @brief Times operator for two integer vectors.
    *
    * @pre getDim() == rhs.getDim()
    */
   MultiIntVector
   operator * (
      const MultiIntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      MultiIntVector tmp = *this;
      tmp *= rhs;
      return tmp;
   }

   MultiIntVector
   operator * (
      const int& rhs) const
   {
      MultiIntVector tmp = *this;
      for (int b = 0; b < d_ratio.size(); ++b) {
         tmp.d_ratio[b] *= rhs;
      }
      return tmp;
   }

   void
   ceilingDivide(const MultiIntVector& rhs)
   {
      for (int b = 0; b < d_ratio.size(); ++b) {
         d_ratio[b].ceilingDivide(rhs.d_ratio[b]);
      }
   }

   static MultiIntVector
   ceilingDivide(
      const MultiIntVector& numerator,
      const MultiIntVector& denominator)
   {
      MultiIntVector tmp_num(numerator);
      tmp_num.ceilingDivide(denominator);
      return tmp_num;
   }

   /**
    * @brief Quotient operator for two integer vectors.
    *
    * @pre getDim() == rhs.getDim()
    */
   MultiIntVector
   operator / (
      const MultiIntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      MultiIntVector tmp = *this;
      tmp /= rhs;
      return tmp;
   }

   MultiIntVector&
   operator /= (
      const MultiIntVector rhs)
   {
      for (int b = 0; b < d_ratio.size(); ++b) {
         d_ratio[b] /= rhs.d_ratio[b];
      }
      return *this;
   }

   MultiIntVector
   operator / (
      const int& rhs) const
   {
      MultiIntVector tmp = *this;
      for (int b = 0; b < d_ratio.size(); ++b) {
         tmp.d_ratio[b] /= rhs;
      }
      return tmp;
   }


   /**
    * @brief Returns true if each integer in vector is less or equal to
    *        corresponding integer in comparison vector.
    *
    * @pre getDim() == rhs.getDim()
    */
   bool
   operator <= (
      const MultiIntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      bool result = true;
      for (int b = 0; result && (b < d_ratio.size()); ++b) {
         result = result && (d_ratio[b] <= rhs.d_ratio[b]);
      }
      return result;
   }

   /**
    * @brief Returns true if each integer in vector is less than the
    *        corresponding integer in comparison vector.
    *
    * @pre getDim() == rhs.getDim()
    */
   bool
   operator < (
      const MultiIntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      bool result = true;
      for (int b = 0; result && (b < d_ratio.size()); ++b) {
         result = result && (d_ratio[b] < rhs.d_ratio[b]);
      }
      return result;
   }

   /**
    * @brief Returns true if each integer in vector is greater or equal to
    *        corresponding integer in comparison vector.
    *
    * @pre getDim() == rhs.getDim()
    */
   bool
   operator >= (
      const MultiIntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      bool result = true;
      for (int b = 0; result && (b < d_ratio.size()); ++b) {
         result = result && (d_ratio[b] >= rhs.d_ratio[b]);
      }
      return result;
   }

   bool
   operator > (
      const MultiIntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      bool result = true;
      for (int b = 0; result && (b < d_ratio.size()); ++b) {
         result = result && (d_ratio[b] > rhs.d_ratio[b]);
      }
      return result;
   }

   /**
    * @brief Returns true if two vector objects are equal.  All components
    *        must be the same for equality.
    *
    * @pre getDim() == rhs.getDim()
    */
   bool
   operator == (
      const MultiIntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      bool result = (d_ratio.size() == rhs.d_ratio.size());
      for (int b = 0; result && (b < d_ratio.size()); ++b) {
         result = result && (d_ratio[b] == rhs.d_ratio[b]);
      }
      return result;
   }

   /**
    * @brief Returns true if two vector objects are not equal.  Any of
    *        the components may be different for inequality.
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



   /**
    * @brief Read an integer vector from an input stream.  The format for
    *        the input is (i0,...,in) for an n-dimensional vector.
    */
   friend std::istream&
   operator >> (
      std::istream& s,
      MultiIntVector& rhs);

   /**
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

   static int s_max_blocks;

   std::vector<IntVector> d_ratio;


};


}
}


#endif
