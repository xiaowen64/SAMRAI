/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Interface to C++ vector implementation for Sundials package.
 *
 ************************************************************************/

namespace SAMRAI {
namespace solv {

#define SABSVEC_CAST(v) \
   (static_cast<SundialsAbstractVector *>(v \
                                          -> \
                                          content))

SAMRAI_INLINE_KEYWORD
N_Vector
SundialsAbstractVector::getNVector()
{
   return d_n_vector;
}

/* Duplicate vector structure and allocate data storage for new vector.
 * Note: This function should only be invoked from within the Sundials
 * package for producing temporary vectors. */
SAMRAI_INLINE_KEYWORD
N_Vector
SundialsAbstractVector::N_VClone_SAMRAI(
   N_Vector w)
{
   /* Create content, which in this case is the SAMRAI wrapper vector object */
   SundialsAbstractVector* v = SABSVEC_CAST(w)->makeNewVector();

   return v->getNVector();
}

/* Free vector structure and associated data. */
SAMRAI_INLINE_KEYWORD
void
SundialsAbstractVector::N_VDestroy_SAMRAI(
   N_Vector v)
{
   if (v) {
      SABSVEC_CAST(v)->freeVector();
   }
}

/* Set z = a * x + b * y */
SAMRAI_INLINE_KEYWORD
void
SundialsAbstractVector::N_VLinearSum_SAMRAI(
   realtype a,
   N_Vector x,
   realtype b,
   N_Vector y,
   N_Vector z)
{
   SABSVEC_CAST(z)->setLinearSum(a, SABSVEC_CAST(x), b, SABSVEC_CAST(y));
}

/* Set vector entries to scalar: v = c  */
SAMRAI_INLINE_KEYWORD
void
SundialsAbstractVector::N_VConst_SAMRAI(
   realtype c,
   N_Vector z)
{
   SABSVEC_CAST(z)->setToScalar(c);
}

/* Set z_i = x_i * y_i */
SAMRAI_INLINE_KEYWORD
void
SundialsAbstractVector::N_VProd_SAMRAI(
   N_Vector x,
   N_Vector y,
   N_Vector z)
{
   SABSVEC_CAST(z)->pointwiseMultiply(SABSVEC_CAST(x), SABSVEC_CAST(y));
}

/* Set z_i = x_i / y_i */
SAMRAI_INLINE_KEYWORD
void
SundialsAbstractVector::N_VDiv_SAMRAI(
   N_Vector x,
   N_Vector y,
   N_Vector z)
{
   SABSVEC_CAST(z)->pointwiseDivide(SABSVEC_CAST(x), SABSVEC_CAST(y));
}

/* Scale vector entries: z = c * x   */
SAMRAI_INLINE_KEYWORD
void
SundialsAbstractVector::N_VScale_SAMRAI(
   realtype c,
   N_Vector x,
   N_Vector z)
{
   SABSVEC_CAST(z)->scaleVector(SABSVEC_CAST(x), c);
}

/* Set z_i = | x_i | */
SAMRAI_INLINE_KEYWORD
void
SundialsAbstractVector::N_VAbs_SAMRAI(
   N_Vector x,
   N_Vector z)
{
   SABSVEC_CAST(z)->setAbs(SABSVEC_CAST(x));
}

/* Set z_i = 1.0 / x_i */
SAMRAI_INLINE_KEYWORD
void
SundialsAbstractVector::N_VInv_SAMRAI(
   N_Vector x,
   N_Vector z)
{
   SABSVEC_CAST(z)->pointwiseReciprocal(SABSVEC_CAST(x));
}

/* Set z_i = x_i + b */
SAMRAI_INLINE_KEYWORD
void
SundialsAbstractVector::N_VAddConst_SAMRAI(
   N_Vector x,
   realtype b,
   N_Vector z)
{
   SABSVEC_CAST(z)->addScalar(SABSVEC_CAST(x), b);
}

/* Return dot product: (x,y) = sum( x_i * y_i ) */
SAMRAI_INLINE_KEYWORD
realtype
SundialsAbstractVector::N_VDotProd_SAMRAI(
   N_Vector x,
   N_Vector y)
{
   return SABSVEC_CAST(x)->dotWith(SABSVEC_CAST(y));
}

/* Return max-norm of vector x */
SAMRAI_INLINE_KEYWORD
realtype
SundialsAbstractVector::N_VMaxNorm_SAMRAI(
   N_Vector x)
{
   return SABSVEC_CAST(x)->maxNorm();
}

/* Return weighted RMS-norm of vector x; w is weighting vector. */
SAMRAI_INLINE_KEYWORD
realtype
SundialsAbstractVector::N_VWrmsNorm_SAMRAI(
   N_Vector x,
   N_Vector w)
{
   return SABSVEC_CAST(x)->weightedRMSNorm(SABSVEC_CAST(w));
}

/* Return minimum entry in x. */
SAMRAI_INLINE_KEYWORD
realtype
SundialsAbstractVector::N_VMin_SAMRAI(
   N_Vector x)
{
   return SABSVEC_CAST(x)->vecMin();
}

/* Return weighted L2-norm of vector x */
SAMRAI_INLINE_KEYWORD
realtype
SundialsAbstractVector::N_VWL2Norm_SAMRAI(
   N_Vector x,
   N_Vector w)
{
   return SABSVEC_CAST(x)->weightedL2Norm(SABSVEC_CAST(w));
}

/* Return L1-norm of vector x */
SAMRAI_INLINE_KEYWORD
realtype
SundialsAbstractVector::N_VL1Norm_SAMRAI(
   N_Vector x)
{
   return SABSVEC_CAST(x)->L1Norm();
}

/*
 * Set each entry in vector z based on the vector x as follows:
 * if | x_i | >= c, then z_i = 1.0, else z_i = 0.0.
 */
SAMRAI_INLINE_KEYWORD
void
SundialsAbstractVector::N_VCompare_SAMRAI(
   realtype c,
   N_Vector x,
   N_Vector z)
{
   SABSVEC_CAST(z)->compareToScalar(SABSVEC_CAST(x), c);
}

/*
 * Set each entry of vector z: v_i =  1.0 / x_i, where x_i is an
 * entry in the vector x, unless x_i = 0.0.  If x_i = 0.0, then
 * return 0.  Otherwise, 1 is returned.
 */
SAMRAI_INLINE_KEYWORD
booleantype
SundialsAbstractVector::N_VInvTest_SAMRAI(
   N_Vector x,
   N_Vector z)
{
   return SABSVEC_CAST(z)->testReciprocal(SABSVEC_CAST(x));
}

}
}
