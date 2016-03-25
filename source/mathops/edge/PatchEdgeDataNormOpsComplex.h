//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/mathops/edge/PatchEdgeDataNormOpsComplex.h $
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Norm operations for complex edge-centered data.
//

#ifndef included_math_PatchEdgeDataNormOpsComplex
#define included_math_PatchEdgeDataNormOpsComplex

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif
#ifndef included_tbox_Complex
#include "tbox/Complex.h"
#endif
#ifndef included_math_ArrayDataNormOpsComplex
#include "ArrayDataNormOpsComplex.h"
#endif
#ifndef included_hier_Box
#include "Box.h"
#endif
#ifndef included_pdat_EdgeData
#include "EdgeData.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif

namespace SAMRAI {
    namespace math {

#ifndef NULL
#define NULL (0)
#endif

/**
 * Class PatchEdgeDataNormOpsComplex<DIM> provides a collection of common 
 * norm operations that may be applied to complex edge-centered patch data.  
 * The primary intent of this class is to define part of the interface for the
 * PatchEdgeDataOpsComplex<DIM> class which provides access to operations 
 * that may be used to manipulate complex edge-centered patch data.  Each 
 * member function accepts a box argument indicating the region of index space 
 * on which the operation should be performed.  The operation will be performed
 * on the intersection of this box and those boxes corresponding to the patch
 * data objects.  Also, each operation allows an additional edge-centered 
 * patch data object to be used to represent a control volume that weights the
 * contribution of each data entry in the given norm calculation.  Note that 
 * the control volume patch data must be of type double and have edge-centered 
 * geometry (i.e., the same as the data itself).  The use of control volumes 
 * is important when these operations are used in vector kernels where the 
 * data resides over multiple levels of spatial resolution in an AMR hierarchy.
 * If the control volume is not given in the function call, it will be ignored
 * in the calculation.  Also, note that the depth of the control volume patch
 * data object must be either 1 or be equal to the depth of the other data 
 * objects.
 * 
 * Note that a similar set of norm operations is implemented for real 
 * patch data (double and float) in the class PatchEdgeDataNormOpsReal<DIM>.
 *
 * @see math::ArrayDataNormOpsComplex
 */

template<int DIM>
class PatchEdgeDataNormOpsComplex
{
public:
   /** 
    * Empty constructor and destructor.
    */
   PatchEdgeDataNormOpsComplex();

   virtual ~PatchEdgeDataNormOpsComplex<DIM>();

   /**
    * Return the number of data values for the edge-centered data object
    * in the given box.  Note that it is assumed that the box refers to
    * the cell-centered index space corresponding to the patch hierarchy.
    */
   int numberOfEntries(const tbox::Pointer< pdat::EdgeData<DIM,dcomplex> >& data,
                       const hier::Box<DIM>& box) const;

   /**
    * Return sum of control volume entries for the edge-centered data object.
    */
   double sumControlVolumes(
      const tbox::Pointer< pdat::EdgeData<DIM,dcomplex> >& data, 
      const tbox::Pointer< pdat::EdgeData<DIM,double> >& cvol, 
      const hier::Box<DIM>& box) const;

   /**
    * Set destination component to norm of source component.  That is,
    * each destination entry is set to
    * \f$d_i = \sqrt{ {real(s_i)}^2 + {imag(s_i)}^2 }\f$.
    */
   void abs(tbox::Pointer< pdat::EdgeData<DIM,double> >& dst,
            const tbox::Pointer< pdat::EdgeData<DIM,dcomplex> >& src,
            const hier::Box<DIM>& box) const;

   /**
    * Return discrete \f$L_1\f$-norm of the data using the control volume to
    * weight the contribution of each data entry to the sum.  That is, the
    * return value is the sum \f$\sum_i ( \sqrt{data_i * \bar{data_i}}*cvol_i )\f$.
    * If the control volume is NULL, the return value is 
    * \f$\sum_i ( \sqrt{data_i * \bar{data_i}} )\f$.
    */
   double L1Norm(
      const tbox::Pointer< pdat::EdgeData<DIM,dcomplex> >& data, 
      const hier::Box<DIM>& box, 
      const tbox::Pointer< pdat::EdgeData<DIM,double> > cvol = NULL) const;

   /**
    * Return discrete \f$L_2\f$-norm of the data using the control volume to
    * weight the contribution of each data entry to the sum.  That is, the
    * return value is the sum 
    * \f$\sqrt{ \sum_i ( data_i * \bar{data_i} cvol_i ) }\f$.
    * If the control volume is NULL, the return value is 
    * \f$\sqrt{ \sum_i ( data_i * \bar{data_i} ) }\f$.
    */
   double L2Norm(
      const tbox::Pointer< pdat::EdgeData<DIM,dcomplex> >& data, 
      const hier::Box<DIM>& box, 
      const tbox::Pointer< pdat::EdgeData<DIM,double> > cvol = NULL) const;

   /**
    * Return discrete weighted \f$L_2\f$-norm of the data using the control
    * volume to weight the contribution of the data and weight entries to
    * the sum.  That is, the return value is the sum \f$\sqrt{ \sum_i (
    * (data_i * wgt_i) * \bar{(data_i * wgt_i)} cvol_i ) }\f$.  If the control 
    * volume is NULL, the return value is 
    * \f$\sqrt{ \sum_i ( (data_i * wgt_i) * \bar{(data_i * wgt_i)} cvol_i ) }\f$.
    */
   double weightedL2Norm(
      const tbox::Pointer< pdat::EdgeData<DIM,dcomplex> >& data,
      const tbox::Pointer< pdat::EdgeData<DIM,dcomplex> >& weight,
      const hier::Box<DIM>& box,
      const tbox::Pointer< pdat::EdgeData<DIM,double> > cvol = NULL) const;

   /**
    * Return discrete root mean squared norm of the data.  If the control
    * volume is not NULL, the return value is the \f$L_2\f$-norm divided by
    * the square root of the sum of the control volumes.  Otherwise, the
    * return value is the \f$L_2\f$-norm divided by the square root of the 
    * number of data entries.
    */
   double RMSNorm(
      const tbox::Pointer< pdat::EdgeData<DIM,dcomplex> >& data,
      const hier::Box<DIM>& box,
      const tbox::Pointer< pdat::EdgeData<DIM,double> > cvol = NULL) const;

   /**
    * Return discrete weighted root mean squared norm of the data.  If the 
    * control volume is not NULL, the return value is the weighted \f$L_2\f$-norm 
    * divided by the square root of the sum of the control volumes.  Otherwise, 
    * the return value is the weighted \f$L_2\f$-norm divided by the square root 
    * of the number of data entries.
    */
   double weightedRMSNorm(
      const tbox::Pointer< pdat::EdgeData<DIM,dcomplex> >& data,
      const tbox::Pointer< pdat::EdgeData<DIM,dcomplex> >& weight,
      const hier::Box<DIM>& box,
      const tbox::Pointer< pdat::EdgeData<DIM,double> > cvol = NULL) const;

   /**
    * Return the \f$\max\f$-norm of the data using the control volume to weight
    * the contribution of each data entry to the maximum.  That is, the return
    * value is \f$\max_i ( \sqrt{data_i * \bar{data_i}} )\f$, where the max is 
    * over the data elements where \f$cvol_i > 0\f$.  If the control volume is 
    * NULL, it is ignored during the computation of the maximum.
    */
   double maxNorm(
      const tbox::Pointer< pdat::EdgeData<DIM,dcomplex> >& data,
      const hier::Box<DIM>& box,
      const tbox::Pointer< pdat::EdgeData<DIM,double> > cvol = NULL) const;

   /**
    * Return the dot product of the two data arrays using the control volume
    * to weight the contribution of each product to the sum.  That is, the
    * return value is the sum \f$\sum_i ( data1_i * \bar{data2_i} * cvol_i )\f$.
    * If the control volume is NULL, it is ignored during the summation.
    */
   dcomplex dot(
      const tbox::Pointer< pdat::EdgeData<DIM,dcomplex> >& data1,
      const tbox::Pointer< pdat::EdgeData<DIM,dcomplex> >& data2,
      const hier::Box<DIM>& box,
      const tbox::Pointer< pdat::EdgeData<DIM,double> > cvol = NULL) const;

   /**
    * Return the integral of the function represented by the data array.
    * The return value is the sum \f$\sum_i ( data_i * vol_i )\f$.
    */
   dcomplex integral(
      const tbox::Pointer< pdat::EdgeData<DIM,dcomplex> >& data,
      const hier::Box<DIM>& box,
      const tbox::Pointer< pdat::EdgeData<DIM,double> > vol) const;

private:
   // The following are not implemented:
   PatchEdgeDataNormOpsComplex(const PatchEdgeDataNormOpsComplex<DIM>&);
   void operator=(const PatchEdgeDataNormOpsComplex<DIM>&);

   ArrayDataNormOpsComplex<DIM> d_array_ops;
};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "PatchEdgeDataNormOpsComplex.C"
#endif
