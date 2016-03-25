//
// File:        solv_NVector.C
// Package:     SAMRAI solvers
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Interface to C++ vector implementation for PVodeTrio package.
//

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_tbox_Utilities
#include "tbox/Utilities.h"
#endif

#if HAVE_KINSOL || HAVE_PVODE

#ifndef included_solv_PVodeTrioAbstractVector
#include "PVodeTrioAbstractVector.h"
#endif

/**
 * External linkage for vector kernel operations needed by the PVodeTrio
 * nonlinear solver package.  This source file must be compiled for 
 * any application that will use PVodeTrio and a vector kernel whose 
 * operations are provided in a subclass of the abstract base class 
 * solv::PVodeTrioAbstractVector.
 *
 * NOTE: These function definitions assume that the vector data will be
 *       of type {\tt double}.  Also, the return type int is used here 
 *       whereas in PVodeTrio "boole" is used (boole is #define'd to be int).
 */

extern "C" {

   typedef SAMRAI::solv::PVodeTrioAbstractVector* N_Vector; 

   /* Duplicate vector structure and allocate data storage for new vector. 
      Note: This function should only be invoked from within the PVodeTrio
            package for producing temporary vectors. */ 
   N_Vector N_VNew(int N, void* vec_in)
   {
      NULL_USE(N);
      N_Vector v = ((N_Vector) vec_in)->makeNewVector();
      return( v );
   }

   /* Free vector structure and associated data. */ 
   void N_VFree(N_Vector v)
   {
      v->freeVector();
   }

   /* Print vector data. */
   void N_VPrint(const N_Vector v)
   {
      v->printVector();
   }

   /* Set vector entries to scalar: v = c  */ 
   void N_VConst(const double c, N_Vector v)
   {
      v->setToScalar(c);
   }

   /* Scale vector entries: z = c * x   */ 
   void N_VScale(const double c, const N_Vector x,
                 N_Vector z)
   {
      z->scaleVector(x, c);
   }

   /* Set z = a * x + b * y */ 
   void N_VLinearSum(const double a, const N_Vector x,
                     const double b, const N_Vector y,
                     N_Vector z)
   {
      z->setLinearSum(a, x, b, y);
   }

   /* Set z_i = x_i * y_i */ 
   void N_VProd(const N_Vector x, const N_Vector y,
                N_Vector z)
   {
      z->pointwiseMultiply(x, y);
   }

   /* Set z_i = x_i / y_i */ 
   void N_VDiv(const N_Vector x, const N_Vector y,
               N_Vector z)
   {
      z->pointwiseDivide(x, y);
   }

   /* Set z_i = | x_i | */ 
   void N_VAbs(const N_Vector x,
               N_Vector z)
   {
      z->setAbs(x);
   }

  /* Set z_i = 1.0 / x_i */
   void N_VInv(const N_Vector x,
               N_Vector z)
   {
      z->pointwiseReciprocal(x);
   }
 
   /* Set z_i = x_i + b */
   void N_VAddConst(const N_Vector x, const double b,
                    N_Vector z)
   {
      z->addScalar(x, b);
   }
 
   /* Return dot product: (x,y) = sum( x_i * y_i ) */
   double N_VDotProd(const N_Vector x, const N_Vector y)
   {
      return( x->dotWith(y) );
   }
 
   /* Return max-norm of vector x */
   double N_VMaxNorm(const N_Vector x)
   {
      return( x->maxNorm() );
   }
 
   /* Return L1-norm of vector x */
   double N_VL1Norm(const N_Vector x)
   {
      return( x->L1Norm() );
   }
 
   /* Return weighted L2-norm of vector x */
   double N_VWL2Norm(const N_Vector x, const N_Vector w)
   {
      return( x->weightedL2Norm(w) );
   }
 
   /* Return weighted RMS-norm of vector x; w is weighting vector. */
   double N_VWrmsNorm(const N_Vector x, const N_Vector w)
   {
      return( x->weightedRMSNorm(w) );
   }
 
   /* Return minimum entry in x. */
   double N_VMin(const N_Vector x)
   {
      return( x->vecMin() );
   }

   /* Return 0 if x_i <> 0.0 and x_i z_i <= 0.0, for some i.  Else, return 1.*/ 
   int N_VConstrProdPos(const N_Vector c, const N_Vector x)
   {
      return( x->constrProdPos(c) );
   }

   /* 
    * Set each entry in vector z based on the vector x as follows:
    * if | x_i | >= c, then z_i = 1.0, else z_i = 0.0.
    */ 
   void N_VCompare(const double c, const N_Vector x,
                   N_Vector z)
   {
      z->compareToScalar(x, c);
   }

   /* 
    * Set each entry of vector z: v_i =  1.0 / x_i, where x_i is an
    * entry in the vector x, unless x_i = 0.0.  If x_i = 0.0, then 
    * return 0.  Otherwise, 1 is returned.
    */ 
   int N_VInvTest(const N_Vector x,
                  N_Vector z)
   {
      return( z->testReciprocal(x) );
   }

}

#endif
