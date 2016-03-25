//
// File:	PatchNodeDataMiscellaneousOpsReal.h
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Templated miscellaneous operations for real node-centered data.
//

#ifndef included_math_PatchNodeDataMiscellaneousOpsReal
#define included_math_PatchNodeDataMiscellaneousOpsReal

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_iostream
#define included_iostream
#include <iostream>
using namespace std;
#endif
#ifndef included_tbox_ArrayDataMiscellaneousOpsReal
#include "ArrayDataMiscellaneousOpsReal.h"
#endif
#ifndef included_hier_Box
#include "Box.h"
#endif
#ifndef included_pdat_NodeData
#include "NodeData.h"
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
 * Class PatchNodeDataMiscellaneousOpsReal<DIM> provides access to a 
 * collection of operations that may be applied to numerical node-centered 
 * patch data of type double and float.  The primary intent of this class is 
 * to provide the interface to these operations for the class 
 * PatchNodeDataOpsReal<DIM> which provides access to a more complete 
 * set of operations that may be used to manipulate node-centered 
 * patch data.  Each member function accepts a box argument
 * indicating the region of index space on which the operation should be
 * performed.  The operation will be performed on the intersection of this
 * box and those boxes corresponding to the patch data objects.  Also, each
 * operation allows an additional node-centered patch data object to be used
 * to represent a control volume that weights the contribution of each data
 * entry in the given norm calculation.  Note that the control volume patch
 * data must be of type double and have node-centered geometry (i.e., the
 * same as the data itself).  The use of control volumes is important when
 * these operations are used in vector kernels where the data resides over
 * multiple levels of spatial resolution in an AMR hierarchy.  If the control
 * volume is not given in the function call, it will be ignored in the
 * calculation.  Also, note that the depth of the control volume patch data
 * object must be either 1 or be equal to the depth of the other data objects.
 *
 * Since these operations are used only by the vector kernels for the KINSOL
 * and PVODE solver packages at this time, they are intended to be instantiated
 * for the standard built-in types double and float (since those solvers only
 * treat double and float data).  To extend this class to other data types or
 * to include other operations, the member functions must be specialized or the
 * new operations must be added.
 *
 * @see math::ArrayDataMiscellaneousOpsReal
 */

template<int DIM, class TYPE>
class PatchNodeDataMiscellaneousOpsReal
{
public:
   /** 
    * Empty constructor and destructor.
    */
   PatchNodeDataMiscellaneousOpsReal();

   virtual ~PatchNodeDataMiscellaneousOpsReal<DIM,TYPE>();

   /**
    * Return 1 if \f$\|data2_i\| > 0\f$ and \f$data1_i * data2_i \leq 0\f$, for
    * any \f$i\f$ in the index region, where \f$cvol_i > 0\f$.  Otherwise return 0.
    * If the control volume is NULL, all values in the index set are used.
    */
   int computeConstrProdPos(const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& data1, 
                            const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& data2,
                            const hier::Box<DIM>& box, 
                            const tbox::Pointer< pdat::NodeData<DIM,double> > cvol =
                                  ((pdat::NodeData<DIM,double>*)NULL)) const;

   /**
    * Wherever \f$cvol_i > 0\f$ in the index region, set \f$dst_i = 1\f$
    * if \f$\|src_i\| > \alpha\f$, and \f$dst_i = 0\f$ otherwise.  If the control 
    * volume is NULL, all values in the index set are considered.
    */
   void compareToScalar(tbox::Pointer< pdat::NodeData<DIM,TYPE> >& dst, 
                        const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& src, 
                        const TYPE& alpha,
                        const hier::Box<DIM>& box, 
                        const tbox::Pointer< pdat::NodeData<DIM,double> > cvol =
                              ((pdat::NodeData<DIM,double>*)NULL)) const;

   /**
    * Wherever \f$cvol_i > 0\f$ in the index region, set \f$dst_i = 1/src_i\f$ if
    * \f$src_i \neq 0\f$, and \f$dst_i = 0\f$ otherwise.  If \f$dst_i = 0\f$ anywhere,
    * 0 is the return value.  Otherwise 1 is returned.  If the control volume
    * all values in the index set are considered.
    */
   int testReciprocal(tbox::Pointer< pdat::NodeData<DIM,TYPE> >& dst, 
                      const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& src,
                      const hier::Box<DIM>& box, 
                      const tbox::Pointer< pdat::NodeData<DIM,double> > cvol =
                            ((pdat::NodeData<DIM,double>*)NULL)) const;

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
   maxPointwiseDivide(const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& numer,
		      const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& denom,
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
    * tbox::IEEE::getFLT_MAX() is returned (see @tbox:IEEE).
    *
    * @b Note: This method is currently intended to support the
    * SUNDIALS vector wrapper only.  Please do not use it!
    */
   TYPE 
   minPointwiseDivide(const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& numer,
		      const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& denom,
		      const hier::Box<DIM>& box) const;

private:
   // The following are not implemented:
   PatchNodeDataMiscellaneousOpsReal(
      const PatchNodeDataMiscellaneousOpsReal<DIM,TYPE>&);
   void operator=(const PatchNodeDataMiscellaneousOpsReal<DIM,TYPE>&);

   ArrayDataMiscellaneousOpsReal<DIM,TYPE> d_array_ops;
};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "PatchNodeDataMiscellaneousOpsReal.C"
#endif
