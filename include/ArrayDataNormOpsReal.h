//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/mathops/array/ArrayDataNormOpsReal.h $
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Templated norm operations for real data arrays.
//

#ifndef included_math_ArrayDataNormOpsReal
#define included_math_ArrayDataNormOpsReal

#include "SAMRAI_config.h"
#include "Box.h"
#include "ArrayData.h"

namespace SAMRAI {
    namespace math {

/**
 * Class ArrayDataNormOpsReal<DIM> provides a set of common norm operations 
 * that may be applied to arrays of real data values (either float or double)
 * maintained as pdat::ArrayData<DIM> objects.  The intent of this class 
 * is to provide a single implementation of these operations as they are needed 
 * by objects that perform these operations on the standard array-based patch 
 * data types (i.e., cell-centered, face-centered, node-centered).  Each of 
 * the norm operations is implemented in two different ways.  The choice of 
 * operation is based on whether control volume information is to be used to 
 * weight the contribution of each data entry to the norm calculation.  The 
 * use of control volumes is important when these operations are used in 
 * vector kernels where the data resides over multiple levels in an AMR 
 * hierarchy.  Note also that each operation will be performed on the 
 * intersection of the box in the function argument list and the boxes 
 * associated with all pdat::ArrayData<DIM> objects.
 *
 * These operations typically apply only to the numerical standard built-in
 * types, such as double, float, and the complex type (which may or may not 
 * be a built-in type depending on the C++ compiler).  This templated 
 * class should only be used to instantiate objects with double or float as 
 * the template parameter. A similar set of norm operations is implemented 
 * for complex array data in the class ArrayDataNormOpsComplex<DIM>.
 *
 * @see pdat::ArrayData
 */

template<int DIM, class TYPE>
class ArrayDataNormOpsReal
{
public:
   /** 
    * Empty constructor and destructor.
    */
   ArrayDataNormOpsReal();

   ~ArrayDataNormOpsReal<DIM,TYPE>();

   /**
    * Return sum of entries in control volume array.
    */
   double sumControlVolumes(const pdat::ArrayData<DIM,TYPE>& data,
                            const pdat::ArrayData<DIM,double>& cvol,
                            const hier::Box<DIM>& box) const;

   /**
    * Set destination component to absolute value of source component.
    * That is, each destination entry is set to \f$d_i = \| s_i \|\f$.
    */
   void abs(pdat::ArrayData<DIM,TYPE>& dst,
            const pdat::ArrayData<DIM,TYPE>& src,
            const hier::Box<DIM>& box) const;

   /**
    * Return discrete \f$L_1\f$-norm of the data using the control volume to 
    * weight the contribution of each data entry to the sum.  That is, the 
    * return value is the sum \f$\sum_i ( \| data_i \| cvol_i )\f$.  
    */
   double L1NormWithControlVolume(const pdat::ArrayData<DIM,TYPE>& data,
                                  const pdat::ArrayData<DIM,double>& cvol,
                                  const hier::Box<DIM>& box) const;

   /**
    * Return discrete \f$L_1\f$-norm of the data.  That is, the return value is 
    * the sum \f$\sum_i ( \| data_i \| )\f$. 
    */
   double L1Norm(const pdat::ArrayData<DIM,TYPE>& data,
                 const hier::Box<DIM>& box) const;

   /**
    * Return discrete \f$L_2\f$-norm of the data using the control volume to 
    * weight the contribution of each data entry to the sum.  That is, the 
    * return value is the sum \f$\sqrt{ \sum_i ( (data_i)^2 cvol_i ) }\f$. 
    */
   double L2NormWithControlVolume(const pdat::ArrayData<DIM,TYPE>& data,
                                  const pdat::ArrayData<DIM,double>& cvol,
                                  const hier::Box<DIM>& box) const;

   /**
    * Return discrete \f$L_2\f$-norm of the data using the control volume to 
    * weight the contribution of each data entry to the sum.  That is, the 
    * return value is the sum \f$\sqrt{ \sum_i ( (data_i)^2 cvol_i ) }\f$. 
    */
   double L2Norm(const pdat::ArrayData<DIM,TYPE>& data,
                 const hier::Box<DIM>& box) const;

   /**
    * Return discrete weighted \f$L_2\f$-norm of the data using the control 
    * volume to weight the contribution of the data and weight entries to 
    * the sum.  That is, the return value is the sum \f$\sqrt{ \sum_i ( 
    * (data_i * weight_i)^2 cvol_i ) }\f$.
    */
   double weightedL2NormWithControlVolume(const pdat::ArrayData<DIM,TYPE>& data,
                                          const pdat::ArrayData<DIM,TYPE>& weight,
                                          const pdat::ArrayData<DIM,double>& cvol,
                                          const hier::Box<DIM>& box) const;

   /**
    * Return discrete weighted \f$L_2\f$-norm of the data.  That is, the return 
    * value is the sum \f$\sqrt{ \sum_i ( (data_i * weight_i)^2 ) }\f$.
    */
   double weightedL2Norm(const pdat::ArrayData<DIM,TYPE>& data,
                         const pdat::ArrayData<DIM,TYPE>& weight,
                         const hier::Box<DIM>& box) const;

   /**
    * Return the \f$\max\f$-norm of the data using the control volume to weight 
    * the contribution of each data entry to the maximum.  That is, the return 
    * value is \f$\max_i ( \| data_i \| )\f$, where the max is over the data 
    * elements where \f$cvol_i > 0\f$.
    */
   double maxNormWithControlVolume(const pdat::ArrayData<DIM,TYPE>& data,
                                   const pdat::ArrayData<DIM,double>& cvol,
                                   const hier::Box<DIM>& box) const;

   /**
    * Return the \f$\max\f$-norm of the data.  That is, the return value is 
    * \f$\max_i ( \| data_i \| )\f$.
    */
   double maxNorm(const pdat::ArrayData<DIM,TYPE>& data,
                  const hier::Box<DIM>& box) const;

   /**
    * Return the dot product of the two data arrays using the control volume 
    * to weight the contribution of each product to the sum.  That is, the 
    * return value is the sum \f$\sum_i ( data1_i * data2_i * cvol_i )\f$.
    */
   TYPE dotWithControlVolume(const pdat::ArrayData<DIM,TYPE>& data1,
                             const pdat::ArrayData<DIM,TYPE>& data2,
                             const pdat::ArrayData<DIM,double>& cvol,
                             const hier::Box<DIM>& box) const;

   /**
    * Return the dot product of the two data arrays.  That is, the 
    * return value is the sum \f$\sum_i ( data1_i * data2_i )\f$.
    */
   TYPE dot(const pdat::ArrayData<DIM,TYPE>& data1,
            const pdat::ArrayData<DIM,TYPE>& data2,
            const hier::Box<DIM>& box) const;

   /**
    * Return the integral of the function based on the data array.
    * The return value is the sum \f$\sum_i ( data_i * vol_i )\f$.
    */
   TYPE integral(const pdat::ArrayData<DIM,TYPE>& data,
                const pdat::ArrayData<DIM,double>& vol,
                const hier::Box<DIM>& box) const;

private:
   // The following are not implemented:
   ArrayDataNormOpsReal(const ArrayDataNormOpsReal<DIM,TYPE>&);
   void operator=(const ArrayDataNormOpsReal<DIM,TYPE>&);
};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "ArrayDataNormOpsReal.C"
#endif
