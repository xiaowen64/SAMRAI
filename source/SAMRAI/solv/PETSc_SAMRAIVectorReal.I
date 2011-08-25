/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   "Glue code" between PETSc vector interface and SAMRAI vectors.
 *
 ************************************************************************/

#define C_PSVEC_CAST(x) \
   (dynamic_cast<const PETSc_SAMRAIVectorReal< \
                    TYPE> *>(x))

namespace SAMRAI {
namespace solv {

template<class TYPE>
SAMRAI_INLINE_KEYWORD
double
PETSc_SAMRAIVectorReal<TYPE>::dotWith(
   const PETScAbstractVectorReal<TYPE>* y,
   bool local_only) const
{
   return d_samrai_vector->dot(C_PSVEC_CAST(y)->d_samrai_vector, local_only);
} // dotWith

template<class TYPE>
SAMRAI_INLINE_KEYWORD
double
PETSc_SAMRAIVectorReal<TYPE>::TdotWith(
   const PETScAbstractVectorReal<TYPE>* y,
   bool local_only) const
{
   return d_samrai_vector->dot(C_PSVEC_CAST(y)->d_samrai_vector, local_only);
} // TdotWith

template<class TYPE>
SAMRAI_INLINE_KEYWORD
double
PETSc_SAMRAIVectorReal<TYPE>::L1Norm(
   bool local_only) const
{
   return d_samrai_vector->L1Norm(local_only);
} // L1Norm

template<class TYPE>
SAMRAI_INLINE_KEYWORD
double
PETSc_SAMRAIVectorReal<TYPE>::L2Norm(
   bool local_only) const
{
   return d_samrai_vector->L2Norm(local_only);
} // L2Norm

template<class TYPE>
SAMRAI_INLINE_KEYWORD
double
PETSc_SAMRAIVectorReal<TYPE>::maxNorm(
   bool local_only) const
{
   return d_samrai_vector->maxNorm(local_only);
} // maxNorm

template<class TYPE>
SAMRAI_INLINE_KEYWORD
void
PETSc_SAMRAIVectorReal<TYPE>::scaleVector(
   const TYPE alpha)
{
   d_samrai_vector->scale(alpha, d_samrai_vector);
} // scaleVector

template<class TYPE>
SAMRAI_INLINE_KEYWORD
void
PETSc_SAMRAIVectorReal<TYPE>::copyVector(
   const PETScAbstractVectorReal<TYPE>* v_src)
{
   d_samrai_vector->copyVector(C_PSVEC_CAST(v_src)->d_samrai_vector);
} // copyVector

template<class TYPE>
SAMRAI_INLINE_KEYWORD
void
PETSc_SAMRAIVectorReal<TYPE>::setToScalar(
   const TYPE alpha)
{
   d_samrai_vector->setToScalar(alpha);
} // setToScalar

template<class TYPE>
SAMRAI_INLINE_KEYWORD
void
PETSc_SAMRAIVectorReal<TYPE>::swapWith(
   PETScAbstractVectorReal<TYPE>* v_other)
{
   d_samrai_vector->swapVectors(C_PSVEC_CAST(v_other)->d_samrai_vector);
} // swapWith

template<class TYPE>
SAMRAI_INLINE_KEYWORD
void
PETSc_SAMRAIVectorReal<TYPE>::setAXPY(
   const TYPE alpha,
   const PETScAbstractVectorReal<TYPE>* x)
{
   d_samrai_vector->axpy(alpha, C_PSVEC_CAST(
         x)->d_samrai_vector, d_samrai_vector);
} // setAXPY

template<class TYPE>
SAMRAI_INLINE_KEYWORD
void
PETSc_SAMRAIVectorReal<TYPE>::setAXPBY(
   const TYPE alpha,
   const PETScAbstractVectorReal<TYPE>* x,
   const TYPE beta)
{
   d_samrai_vector->linearSum(alpha, C_PSVEC_CAST(
         x)->d_samrai_vector, beta, d_samrai_vector);
} // setAXPBY

template<class TYPE>
SAMRAI_INLINE_KEYWORD
void
PETSc_SAMRAIVectorReal<TYPE>::setWAXPY(
   const TYPE alpha,
   const PETScAbstractVectorReal<TYPE>* x,
   const PETScAbstractVectorReal<TYPE>* y)
{
   d_samrai_vector->axpy(alpha, C_PSVEC_CAST(x)->d_samrai_vector,
      C_PSVEC_CAST(y)->d_samrai_vector);
} // setWAXPY

template<class TYPE>
SAMRAI_INLINE_KEYWORD
void
PETSc_SAMRAIVectorReal<TYPE>::pointwiseMultiply(
   const PETScAbstractVectorReal<TYPE>* x,
   const PETScAbstractVectorReal<TYPE>* y)
{
   d_samrai_vector->multiply(C_PSVEC_CAST(x)->d_samrai_vector, C_PSVEC_CAST(
         y)->d_samrai_vector);
} // pointwiseMultiply

template<class TYPE>
SAMRAI_INLINE_KEYWORD
void
PETSc_SAMRAIVectorReal<TYPE>::pointwiseDivide(
   const PETScAbstractVectorReal<TYPE>* x,
   const PETScAbstractVectorReal<TYPE>* y)
{
   d_samrai_vector->divide(C_PSVEC_CAST(x)->d_samrai_vector, C_PSVEC_CAST(
         y)->d_samrai_vector);
} // pointwiseDivide

template<class TYPE>
SAMRAI_INLINE_KEYWORD
double
PETSc_SAMRAIVectorReal<TYPE>::maxPointwiseDivide(
   const PETScAbstractVectorReal<TYPE>* y)
{
   return d_samrai_vector->maxPointwiseDivide(C_PSVEC_CAST(y)->d_samrai_vector);
} // maxPointwiseDivide

template<class TYPE>
SAMRAI_INLINE_KEYWORD
void
PETSc_SAMRAIVectorReal<TYPE>::vecMax(
   int& i,
   TYPE& max) const
{
   static const bool interior_only = true;
   max = d_samrai_vector->max(interior_only);
   // Note: This is a bogus index value!
   //       Hopefully, PETSc doesn't use it for anything.
   i = 0;
} // vecMax

template<class TYPE>
SAMRAI_INLINE_KEYWORD
void
PETSc_SAMRAIVectorReal<TYPE>::vecMin(
   int& i,
   TYPE& min) const
{
   static const bool interior_only = true;
   min = d_samrai_vector->min(interior_only);
   // Note: This is a bogus index value!
   //       Hopefully, PETSc doesn't use it for anything.
   i = 0;
} // vecMin

template<class TYPE>
SAMRAI_INLINE_KEYWORD
void
PETSc_SAMRAIVectorReal<TYPE>::setRandomValues(
   const TYPE width,
   const TYPE low)
{
   d_samrai_vector->setRandomValues(width, low);
} // setRandomValues

template<class TYPE>
SAMRAI_INLINE_KEYWORD
void
PETSc_SAMRAIVectorReal<TYPE>::getDataArray(
   TYPE** array)
{
   *array = NULL;
} // getDataArray

template<class TYPE>
SAMRAI_INLINE_KEYWORD
void
PETSc_SAMRAIVectorReal<TYPE>::restoreDataArray(
   TYPE** array)
{
   *array = NULL;
} // restoreDataArray

template<class TYPE>
SAMRAI_INLINE_KEYWORD
int
PETSc_SAMRAIVectorReal<TYPE>::getDataSize() const
{
   // Note: This is a bogus value!
   //       But, PETSc requires some value to be returned.
   //       Hopefully, this will not cause problems.
   return 0;
} // getDataSize

template<class TYPE>
SAMRAI_INLINE_KEYWORD
int
PETSc_SAMRAIVectorReal<TYPE>::getLocalDataSize() const
{
   // Note: This is a bogus value!
   //       But, PETSc requires some value to be returned.
   //       Hopefully, this will not cause problems.
   return 0;
} // getLocalDataSize

}
}
