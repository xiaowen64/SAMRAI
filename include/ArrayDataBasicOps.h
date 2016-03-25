//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/mathops/array/ArrayDataBasicOps.h $
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Basic templated operations for array data.
//

#ifndef included_math_ArrayDataBasicOps
#define included_math_ArrayDataBasicOps

#include "SAMRAI_config.h"
#include "tbox/Complex.h"
#include "Box.h"
#include "ArrayData.h"

namespace SAMRAI {
    namespace math {

/**
 * Class ArrayDataBasicOps<DIM> implements a set of basic operations 
 * that apply to numerical data maintained as pdat::ArrayData<DIM> objects.
 * These operations include simple arithmetic operations as well as min
 * and max, etc.  This class provides a single implementation of these 
 * operations that may used to manipulate any of the standard array-based 
 * patch data types defined on a patch.  Note that each member function 
 * accepts a box argument which specifies the portion of the array data 
 * on which the associated operation is performed.   The actual index
 * region on which the operation occurs is the intersection of this box
 * and the boxes of all the pdat::ArrayData<DIM> objects involved. 
 *
 * These operations typically apply only to the numerical standard built-in 
 * types, such as double, float, and int, and the complex type (which may or 
 * may not be a built-in type depending on the C++ compiler).  Thus, this 
 * templated class should only be used to instantiate objects with those
 * types as the template parameter.  Those operations whose implementations
 * depend of the data type are specialized for each numerical type.  To use 
 * this class with other standard types or user-defined types (which may or 
 * may not make sense), the member functions must be specialized so that the 
 * correct operations are performed.
 *
 * @see pdat::ArrayData
 */

template<int DIM, class TYPE>
class ArrayDataBasicOps
{
public:
   /**
    * Empty constructor and destructor.
    */ 
   ArrayDataBasicOps();

   ~ArrayDataBasicOps<DIM,TYPE>();

   /**
    * Set dst = alpha * src, elementwise.
    */
   void scale(pdat::ArrayData<DIM,TYPE>& dst,
              const TYPE& alpha,
              const pdat::ArrayData<DIM,TYPE>& src,
              const hier::Box<DIM>& box) const; 

   /**
    * Set dst = src + alpha, elementwise.
    */
   void addScalar(pdat::ArrayData<DIM,TYPE>& dst,
                  const pdat::ArrayData<DIM,TYPE>& src,
                  const TYPE& alpha,
                  const hier::Box<DIM>& box) const;

   /**
    * Set dst = src1 + src2, elementwise.
    */
   void add(pdat::ArrayData<DIM,TYPE>& dst,
            const pdat::ArrayData<DIM,TYPE>& src1,
            const pdat::ArrayData<DIM,TYPE>& src2,
            const hier::Box<DIM>& box) const;

   /**
    * Set dst = src1 - src2, elementwise.
    */
   void subtract(pdat::ArrayData<DIM,TYPE>& dst,
                 const pdat::ArrayData<DIM,TYPE>& src1,
                 const pdat::ArrayData<DIM,TYPE>& src2,
                 const hier::Box<DIM>& box) const;

   /**
    * Set dst = src1 * src2, elementwise.
    */
   void multiply(pdat::ArrayData<DIM,TYPE>& dst,
                 const pdat::ArrayData<DIM,TYPE>& src1,
                 const pdat::ArrayData<DIM,TYPE>& src2,
                 const hier::Box<DIM>& box) const;

   /**
    * Set dst = src1 / src2, elementwise.  No check for division by zero.
    */
   void divide(pdat::ArrayData<DIM,TYPE>& dst,
               const pdat::ArrayData<DIM,TYPE>& src1,
               const pdat::ArrayData<DIM,TYPE>& src2,
               const hier::Box<DIM>& box) const;

   /**
    * Set dst = 1 / src, elementwise.  No check for division by zero.
    */
   void reciprocal(pdat::ArrayData<DIM,TYPE>& dst,
                   const pdat::ArrayData<DIM,TYPE>& src,
                   const hier::Box<DIM>& box) const;

   /**
    * Set dst = alpha * src1 + beta * src2, elementwise.
    */
   void linearSum(pdat::ArrayData<DIM,TYPE>& dst,
                  const TYPE& alpha,
                  const pdat::ArrayData<DIM,TYPE>& src1,
                  const TYPE& beta,
                  const pdat::ArrayData<DIM,TYPE>& src2,
                  const hier::Box<DIM>& box) const;

   /**
    * Set dst = alpha * src1 + src2, elementwise.
    */
   void axpy(pdat::ArrayData<DIM,TYPE>& dst,
             const TYPE& alpha,
             const pdat::ArrayData<DIM,TYPE>& src1,
             const pdat::ArrayData<DIM,TYPE>& src2,
             const hier::Box<DIM>& box) const;

   /**
    * Set dst = alpha * src1 - src2, elementwise.
    */
   void axmy(pdat::ArrayData<DIM,TYPE>& dst,
             const TYPE& alpha,
             const pdat::ArrayData<DIM,TYPE>& src1,
             const pdat::ArrayData<DIM,TYPE>& src2,
             const hier::Box<DIM>& box) const;

   /**
    * Return the minimum array data entry.  If data is complex, return the
    * array data entry with the minimum norm.
    */
   TYPE min(const pdat::ArrayData<DIM,TYPE>& data,
            const hier::Box<DIM>& box) const;

   /**
    * Return the maximum array data entry.  If data is complex, return the
    * array data entry with the maximum norm.
    */
   TYPE max(const pdat::ArrayData<DIM,TYPE>& data,
            const hier::Box<DIM>& box) const;

   /**
    * Set dst to random values.  If the data is int, each element of dst 
    * is set as dst = mrand48(). If the data is double or float, each
    * element of dst is set as dst = width * drand48() + low.  If the 
    * data is complex, each element of dst is set as dst = dcomplex(rval, ival),
    * where rval = real(width) * drand48() + real(low), and 
    * ival = imag(width) * drand48() + imag(low).
    */
   void setRandomValues(pdat::ArrayData<DIM,TYPE>& dst,
                        const TYPE& width,
                        const TYPE& low,
                        const hier::Box<DIM>& box) const;

private:
   // The following are not implemented:
   ArrayDataBasicOps(const ArrayDataBasicOps<DIM,TYPE>&);
   void operator=(const ArrayDataBasicOps<DIM,TYPE>&);

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "ArrayDataBasicOps.C"
#endif
