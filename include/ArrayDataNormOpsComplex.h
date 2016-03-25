//
// File:	ArrayDataNormOpsComplex.h
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Norm operations for complex data arrays.
//

#ifndef included_tbox_ArrayDataNormOpsComplex
#define included_tbox_ArrayDataNormOpsComplex

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_hier_Box
#include "Box.h"
#endif
#ifndef included_tbox_ArrayData
#include "ArrayData.h"
#endif
#ifndef included_tbox_Complex
#include "tbox/Complex.h"
#endif

namespace SAMRAI {
    namespace math {

/**
 * Class ArrayDataNormOpsComplex<DIM> provides a set of common norm 
 * operations that may be applied to arrays of complex data values 
 * maintained as pdat::ArrayData<DIM> objects.  The intent of this class is to 
 * provide a single implementation of these operations as they are needed 
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
 * Note that a similar set of norm operations is implemented for real array 
 * data (double and float) in the class ArrayDataNormOpsReal<DIM>.
 *
 * @see pdat::ArrayData
 */

template<int DIM> class ArrayDataNormOpsComplex
{
public:
   /** 
    * Empty constructor and destructor.
    */
   ArrayDataNormOpsComplex();

   ~ArrayDataNormOpsComplex<DIM>();

   /**
    * Set destination component to norm of source component.  That is, 
    * each destination entry is set to 
    * \f$d_i = \sqrt{ {real(s_i)}^2 + {imag(s_i)}^2 }\f$.
    */
   void abs(pdat::ArrayData<DIM,double>& dst,
            const pdat::ArrayData<DIM,dcomplex>& src,
            const hier::Box<DIM>& box) const;

   /**
    * Return sum of entries in control volume array.
    */
   double sumControlVolumes(const pdat::ArrayData<DIM,dcomplex>& data,
                            const pdat::ArrayData<DIM,double>& cvol,
                            const hier::Box<DIM>& box) const;

   /**
    * Return discrete \f$L_1\f$-norm of the data using the control volume to
    * weight the contribution of each data entry to the sum.  That is, the
    * return value is the sum \f$\sum_i ( \sqrt{data_i * \bar{data_i}} cvol_i )\f$.
    */
   double L1NormWithControlVolume(const pdat::ArrayData<DIM,dcomplex>& data,
                                  const pdat::ArrayData<DIM,double>& cvol,
                                  const hier::Box<DIM>& box) const;

   /**
    * Return discrete \f$L_1\f$-norm of the data.  That is, the return value is
    * the sum \f$\sum_i ( \sqrt{data_i * \bar{data_i}} )\f$.
    */
   double L1Norm(const pdat::ArrayData<DIM,dcomplex>& data,
                 const hier::Box<DIM>& box) const;

   /**
    * Return discrete \f$L_2\f$-norm of the data using the control volume to
    * weight the contribution of each data entry to the sum.  That is, the
    * return value is the sum \f$\sqrt{ \sum_i ( 
    * data_i * \bar{data_i} cvol_i ) }\f$.
    */
   double L2NormWithControlVolume(const pdat::ArrayData<DIM,dcomplex>& data,
                                  const pdat::ArrayData<DIM,double>& cvol,
                                  const hier::Box<DIM>& box) const;

   /**
    * Return discrete \f$L_2\f$-norm of the data using the control volume to
    * weight the contribution of each data entry to the sum.  That is, the
    * return value is the sum \f$\sqrt{ \sum_i ( data_i * \bar{data_i} ) }\f$.
    */
   double L2Norm(const pdat::ArrayData<DIM,dcomplex>& data,
                 const hier::Box<DIM>& box) const;

   /**
    * Return discrete weighted \f$L_2\f$-norm of the data using the control
    * volume to weight the contribution of the data and weight entries to
    * the sum.  That is, the return value is the sum \f$\sqrt{ \sum_i (
    * (data_i * wgt_i) * \bar{(data_i * wgt_i)} cvol_i ) }\f$.
    */
   double weightedL2NormWithControlVolume(const pdat::ArrayData<DIM,dcomplex>& data,
                                          const pdat::ArrayData<DIM,dcomplex>& wgt,
                                          const pdat::ArrayData<DIM,double>& cvol,
                                          const hier::Box<DIM>& box) const;

   /**
    * Return discrete weighted \f$L_2\f$-norm of the data.  That is, the return
    * value is the sum \f$\sqrt{ \sum_i ( (data_i * wgt_i) * 
    * \bar{(data_i * wgt_i)} cvol_i ) }\f$.
    */
   double weightedL2Norm(const pdat::ArrayData<DIM,dcomplex>& data,
                         const pdat::ArrayData<DIM,dcomplex>& wgt,
                         const hier::Box<DIM>& box) const;

   /**
    * Return the \f$\max\f$-norm of the data using the control volume to weight
    * the contribution of each data entry to the maximum.  That is, the return
    * value is \f$\max_i ( \sqrt{data_i * \bar{data_i}} )\f$, where the max is 
    * over the data elements where \f$cvol_i > 0\f$.
    */ 
   double maxNormWithControlVolume(const pdat::ArrayData<DIM,dcomplex>& data,
                                   const pdat::ArrayData<DIM,double>& cvol,
                                   const hier::Box<DIM>& box) const;

   /**
    * Return the \f$\max\f$-norm of the data.  That is, the return value is
    * \f$\max_i ( \sqrt{data_i * \bar{data_i}} )\f$.
    */
   double maxNorm(const pdat::ArrayData<DIM,dcomplex>& data,
                  const hier::Box<DIM>& box) const;

   /**
    * Return the dot product of the two data arrays using the control volume
    * to weight the contribution of each product to the sum.  That is, the
    * return value is the sum \f$\sum_i ( data1_i * \bar{data2_i} * cvol_i )\f$.
    */
   dcomplex dotWithControlVolume(const pdat::ArrayData<DIM,dcomplex>& data1,
                                 const pdat::ArrayData<DIM,dcomplex>& data2,
                                 const pdat::ArrayData<DIM,double>& cvol,
                                 const hier::Box<DIM>& box) const;

   /**
    * Return the dot product of the two data arrays.  That is, the
    * return value is the sum \f$\sum_i ( data1_i * \bar{data2_i} )\f$.
    */
   dcomplex dot(const pdat::ArrayData<DIM,dcomplex>& data1,
                const pdat::ArrayData<DIM,dcomplex>& data2,
                const hier::Box<DIM>& box) const;

   /**
    * Return the integral of the function based on the data array.
    * The return value is the sum \f$\sum_i ( data_i * vol_i )\f$.
    */
   dcomplex integral(const pdat::ArrayData<DIM,dcomplex>& data,
                const pdat::ArrayData<DIM,double>& vol,
                const hier::Box<DIM>& box) const;

private:
   // The following are not implemented:
   ArrayDataNormOpsComplex(const ArrayDataNormOpsComplex<DIM>&);
   void operator=(const ArrayDataNormOpsComplex<DIM>&);
};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "ArrayDataNormOpsComplex.C"
#endif
