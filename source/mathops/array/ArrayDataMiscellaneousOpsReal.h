//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/mathops/array/ArrayDataMiscellaneousOpsReal.h $
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1811 $
// Modified:	$LastChangedDate: 2007-12-20 01:19:26 -0800 (Thu, 20 Dec 2007) $
// Description:	Miscellaneous templated operations for real array data
//

#ifndef included_math_ArrayDataMiscellaneousOpsReal
#define included_math_ArrayDataMiscellaneousOpsReal

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_hier_Box
#include "Box.h"
#endif
#ifndef included_pdat_ArrayData
#include "ArrayData.h"
#endif

namespace SAMRAI {
    namespace math {

/**
 * Class ArrayDataMiscellaneousOpsReal<DIM> provides various operations that 
 * may be applied to arrays of real (double and float) numerical data 
 * values maintained using pdat::ArrayData<DIM> objects.  These operations are 
 * sufficiently different from basic arithmetic and norm operations that we 
 * chose to implement them in a separate class.  However, as in the case of the  * more common operations, the intent of this class is to provide a single 
 * implementation of the operations as they are needed by objects that 
 * manipulate standard array-based patch data types (i.e., cell-centered, 
 * face-centered, node-centered).  Each operation is implemented in two 
 * different ways.  The choice of operation is based on whether control volume 
 * information is to be used to weight the contribution of each data entry 
 * to the calculation.  The use of control volumes is important when these 
 * operations are used in vector kernels where the data resides over multiple 
 * levels of spatial resolution in an AMR hierarchy.  The actual index
 * region on which each operation occurs is the intersection of this box
 * and the boxes of all the pdat::ArrayData<DIM> objects involved.
 *
 * Since these operations are used only by the vector kernels for the KINSOL 
 * and CVODE solver packages at this time, they are intended to be instantiated
 * for the standard built-in types double and float (since those solvers only
 * treat double and float data).  To extend this class to other data types or 
 * to include other operations, the member functions must be specialized or the  * new operations must be added.
 *
 * @see pdat::ArrayData
 */

template<int DIM, class TYPE>
class ArrayDataMiscellaneousOpsReal
{
public:
   /** 
    * Empty constructor and destructor.
    */
   ArrayDataMiscellaneousOpsReal();

   ~ArrayDataMiscellaneousOpsReal<DIM,TYPE>();

   /**
    * Return 1 if \f$\|data2_i\| > 0\f$ and \f$data1_i * data2_i \leq 0\f$, for
    * any \f$i\f$ in the index region, where \f$cvol_i > 0\f$.  Otherwise return 0.
    */
   int 
   computeConstrProdPosWithControlVolume(const pdat::ArrayData<DIM,TYPE>& data1,
                                         const pdat::ArrayData<DIM,TYPE>& data2,
                                         const pdat::ArrayData<DIM,double>& cvol,
                                         const hier::Box<DIM>& box) const;

   /**
    * Return 1 if \f$\|data2_i\| > 0\f$ and \f$data1_i * data2_i \leq 0\f$, for
    * any \f$i\f$ in the index region.  Otherwise return 0.
    */
   int 
   computeConstrProdPos(const pdat::ArrayData<DIM,TYPE>& data1,
                        const pdat::ArrayData<DIM,TYPE>& data2,
                        const hier::Box<DIM>& box) const;

   /**
    * Wherever \f$cvol_i > 0\f$ in the index region, set \f$dst_i = 1\f$
    * if \f$\|src_i\| > \alpha\f$, and \f$dst_i = 0\f$ otherwise.
    */
   void 
   compareToScalarWithControlVolume(pdat::ArrayData<DIM,TYPE>& dst,
                                    const pdat::ArrayData<DIM,TYPE>& src,
                                    const TYPE& alpha,
                                    const pdat::ArrayData<DIM,double>& cvol,
                                    const hier::Box<DIM>& box) const;

   /**
    * Set \f$dst_i = 1\f$ if \f$\|src_i\| > \alpha\f$, and \f$dst_i = 0\f$ otherwise.
    */
   void 
   compareToScalar(pdat::ArrayData<DIM,TYPE>& dst,
                   const pdat::ArrayData<DIM,TYPE>& src,
                   const TYPE& alpha,
                   const hier::Box<DIM>& box) const;

   /**
    * Wherever \f$cvol_i > 0\f$ in the index region, set \f$dst_i = 1/src_i\f$ if
    * \f$src_i \neq 0\f$, and \f$dst_i = 0\f$ otherwise.  If \f$dst_i = 0\f$ anywhere,
    * 0 is the return value.  Otherwise 1 is returned. 
    */
   int 
   testReciprocalWithControlVolume(pdat::ArrayData<DIM,TYPE>& dst,
                                   const pdat::ArrayData<DIM,TYPE>& src,
                                   const pdat::ArrayData<DIM,double>& cvol,
                                   const hier::Box<DIM>& box) const; 

   /**
    * Set \f$dst_i = 1/src_i\f$ if \f$src_i \neq 0\f$, and \f$dst_i = 0\f$ otherwise.  
    * If \f$dst_i = 0\f$ anywhere, 0 is the return value.  Otherwise 1 is returned. 
    */
   int 
   testReciprocal(pdat::ArrayData<DIM,TYPE>& dst,
                  const pdat::ArrayData<DIM,TYPE>& src,
                  const hier::Box<DIM>& box) const;

   /*!
    * @brief Compute max of "conditional" quotients of two arrays.
    *
    * Return the maximum of pointwise "conditional" quotients of the numerator
    * and denominator.
    *
    * The "conditional" quotient is defined as |numerator/denominator|
    * where the denominator is nonzero.  Otherwise, it is defined as
    * |numerator|.
    *
    * @b Note: This method is currently intended to support the
    * PETSc-2.1.6 vector wrapper only.  Please do not use it!
    */
   TYPE 
   maxPointwiseDivide(const pdat::ArrayData<DIM,TYPE>& numer,
		      const pdat::ArrayData<DIM,TYPE>& denom,
		      const hier::Box<DIM>& box) const;

   /*!
    * @brief Compute min of quotients of two arrays.
    *
    * Return the minimum of pointwise quotients of the numerator
    * and denominator.
    *
    * The quotient is defined as (numerator/denominator)
    * where the denominator is nonzero.  When the denominator is zero, the
    * entry is skipped.  If the denominator is always zero, the value of
    * tbox::IEEE::getDBL_MAX() is returned (see @ref SAMRAI::tbox::IEEE).
    *
    * @b Note: This method is currently intended to support the
    * SUNDIALS vector wrapper only.  Please do not use it!
    */
   TYPE 
   minPointwiseDivide(const pdat::ArrayData<DIM,TYPE>& numer,
		      const pdat::ArrayData<DIM,TYPE>& denom,
		      const hier::Box<DIM>& box) const;

private:
   // The following are not implemented:
   ArrayDataMiscellaneousOpsReal(
      const ArrayDataMiscellaneousOpsReal<DIM,TYPE>&);
   void operator=(const ArrayDataMiscellaneousOpsReal<DIM,TYPE>&);

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "ArrayDataMiscellaneousOpsReal.C"
#endif
