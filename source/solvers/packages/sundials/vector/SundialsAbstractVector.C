//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-1/source/solvers/packages/sundials/vector/SundialsAbstractVector.C $
// Package:     SAMRAI solvers package
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1892 $
// Modified:    $LastChangedDate: 2008-01-23 08:29:36 -0800 (Wed, 23 Jan 2008) $
// Description: Interface to C++ vector implementation for Sundials package.
//

#include "SundialsAbstractVector.h"

#include <cstdlib>

#include <tbox/Utilities.h>

#ifdef HAVE_SUNDIALS

namespace SAMRAI {
   namespace solv {

#define SABSVEC_CAST(v) (static_cast<SAMRAI::solv::SundialsAbstractVector *>(v->content))

SundialsAbstractVector::SundialsAbstractVector()
{
   /* Create N vector */
   d_n_vector = (N_Vector) malloc(sizeof *d_n_vector);
   TBOX_ASSERT(d_n_vector != NULL);

   /* Attach content and ops */
   d_n_vector -> content = this;
   d_n_vector -> ops     = SundialsAbstractVector::createVectorOps();
}

SundialsAbstractVector::~SundialsAbstractVector()
{
#if 0
   if(d_n_vector -> ops)
      free(d_n_vector -> ops);

   if(d_n_vector)
      free(d_n_vector);
#endif
}


N_Vector SundialsAbstractVector::getNVector() {
   return d_n_vector;
}

N_Vector_Ops SundialsAbstractVector::createVectorOps()
{
   /* Create vector operation structure */

   N_Vector_Ops ops;
   
   ops = (N_Vector_Ops) std::calloc(1, sizeof(struct _generic_N_Vector_Ops));

// SGS TODO what about missing fns?
   ops->nvclone           = SundialsAbstractVector::N_VClone_SAMRAI;
//   ops->nvcloneempty      = N_VCloneEmpty_SAMRAI;
   ops->nvdestroy         = N_VDestroy_SAMRAI;
//   ops->nvspace           = N_VSpace_SAMRAI;
//   ops->nvgetarraypointer = N_VGetArrayPointer_SAMRAI;
//   ops->nvsetarraypointer = N_VSetArrayPointer_SAMRAI;
   ops->nvlinearsum       = N_VLinearSum_SAMRAI;
   ops->nvconst           = N_VConst_SAMRAI;
   ops->nvprod            = N_VProd_SAMRAI;
   ops->nvdiv             = N_VDiv_SAMRAI;
   ops->nvscale           = N_VScale_SAMRAI;
   ops->nvabs             = N_VAbs_SAMRAI;
   ops->nvinv             = N_VInv_SAMRAI;
   ops->nvaddconst        = N_VAddConst_SAMRAI;
   ops->nvdotprod         = N_VDotProd_SAMRAI;
   ops->nvmaxnorm         = N_VMaxNorm_SAMRAI;
//   ops->nvwrmsnormmask    = N_VWrmsNormMask_SAMRAI;
   ops->nvwrmsnorm        = N_VWrmsNorm_SAMRAI;
   ops->nvmin             = N_VMin_SAMRAI;
   ops->nvwl2norm         = N_VWL2Norm_SAMRAI;
   ops->nvl1norm          = N_VL1Norm_SAMRAI;
   ops->nvcompare         = N_VCompare_SAMRAI;
   ops->nvinvtest         = N_VInvTest_SAMRAI;
//   ops->nvconstrmask      = N_VConstrMask_SAMRAI;
//   ops->nvminquotient     = N_VMinQuotient_SAMRAI;

   return ops;
}


/* Duplicate vector structure and allocate data storage for new vector. 
   Note: This function should only be invoked from within the Sundials
   package for producing temporary vectors. */ 
N_Vector SundialsAbstractVector::N_VClone_SAMRAI(N_Vector w)
{
   /* Create content, which in this case is the SAMRAI wrapper vector object */
   SAMRAI::solv::SundialsAbstractVector *v
      = SABSVEC_CAST(w)->makeNewVector();

   return(v -> getNVector());
}

/* Free vector structure and associated data. */ 
void SundialsAbstractVector::N_VDestroy_SAMRAI(N_Vector v)
{
   SABSVEC_CAST(v)->freeVector();
   free(v->ops); v->ops = NULL;
   free(v); v = NULL;
}

/* Set z = a * x + b * y */ 
void SundialsAbstractVector::N_VLinearSum_SAMRAI(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z)
{
   SABSVEC_CAST(z)->setLinearSum(a, SABSVEC_CAST(x), b, SABSVEC_CAST(y));
}


/* Set vector entries to scalar: v = c  */ 
void SundialsAbstractVector::N_VConst_SAMRAI(realtype c, N_Vector z)
{
   SABSVEC_CAST(z)->setToScalar(c);
}

/* Set z_i = x_i * y_i */ 
void SundialsAbstractVector::N_VProd_SAMRAI(N_Vector x, N_Vector y, N_Vector z)
{
   SABSVEC_CAST(z)->pointwiseMultiply(SABSVEC_CAST(x), SABSVEC_CAST(y));
}

/* Set z_i = x_i / y_i */ 
void SundialsAbstractVector::N_VDiv_SAMRAI(N_Vector x, N_Vector y, N_Vector z)
{
   SABSVEC_CAST(z)->pointwiseDivide(SABSVEC_CAST(x), SABSVEC_CAST(y));
}

/* Scale vector entries: z = c * x   */ 
void SundialsAbstractVector::N_VScale_SAMRAI(realtype c, N_Vector x, N_Vector z)
{
   SABSVEC_CAST(z)->scaleVector(SABSVEC_CAST(x), c);
}

/* Set z_i = | x_i | */ 
void SundialsAbstractVector::N_VAbs_SAMRAI(N_Vector x, N_Vector z)
{
   SABSVEC_CAST(z)->setAbs(SABSVEC_CAST(x));
}

/* Set z_i = 1.0 / x_i */
void SundialsAbstractVector::N_VInv_SAMRAI(N_Vector x, N_Vector z)
{
   SABSVEC_CAST(z)->pointwiseReciprocal(SABSVEC_CAST(x));
}

/* Set z_i = x_i + b */
void SundialsAbstractVector::N_VAddConst_SAMRAI(N_Vector x, realtype b, N_Vector z)
{
   SABSVEC_CAST(z)->addScalar(SABSVEC_CAST(x), b);
}

/* Return dot product: (x,y) = sum( x_i * y_i ) */
realtype SundialsAbstractVector::N_VDotProd_SAMRAI(N_Vector x, N_Vector y)
{
   return( SABSVEC_CAST(x)->dotWith(SABSVEC_CAST(y)) );
}

 
/* Return max-norm of vector x */
realtype SundialsAbstractVector::N_VMaxNorm_SAMRAI(N_Vector x)
{
   return( SABSVEC_CAST(x)->maxNorm() );
}

/* Return weighted RMS-norm of vector x; w is weighting vector. */
realtype SundialsAbstractVector::N_VWrmsNorm_SAMRAI(N_Vector x, N_Vector w)
{
   return( SABSVEC_CAST(x)->weightedRMSNorm(SABSVEC_CAST(w)) );
}

/* Return minimum entry in x. */
realtype SundialsAbstractVector::N_VMin_SAMRAI(N_Vector x)
{
   return( SABSVEC_CAST(x)->vecMin() );
}

/* Return weighted L2-norm of vector x */
realtype SundialsAbstractVector::N_VWL2Norm_SAMRAI(N_Vector x, N_Vector w)
{
   return( SABSVEC_CAST(x)->weightedL2Norm(SABSVEC_CAST(w)) );
}

/* Return L1-norm of vector x */
realtype SundialsAbstractVector::N_VL1Norm_SAMRAI(N_Vector x)
{
   return( SABSVEC_CAST(x)->L1Norm() );
}

/* 
 * Set each entry in vector z based on the vector x as follows:
 * if | x_i | >= c, then z_i = 1.0, else z_i = 0.0.
 */ 
void SundialsAbstractVector::N_VCompare_SAMRAI(realtype c, N_Vector x, N_Vector z)
{
   SABSVEC_CAST(z)->compareToScalar(SABSVEC_CAST(x), c);
}

/* 
 * Set each entry of vector z: v_i =  1.0 / x_i, where x_i is an
 * entry in the vector x, unless x_i = 0.0.  If x_i = 0.0, then 
 * return 0.  Otherwise, 1 is returned.
 */ 
booleantype SundialsAbstractVector::N_VInvTest_SAMRAI(N_Vector x, N_Vector z)
{
   return( SABSVEC_CAST(z)->testReciprocal(SABSVEC_CAST(x)) );
}

}
}

#endif
