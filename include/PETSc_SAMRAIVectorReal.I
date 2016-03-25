//
// File:        PETSc_SAMRAIVectorReal.I
// Package:     SAMRAI solvers
// Copyright:   (c) 1998 The Regents of the University of California
// Revision:    $Revision: 7 $
// Modified:    $Date: 2004-11-30 13:18:17 -0800 (Tue, 30 Nov 2004) $
// Description: "Glue code" between PETSc vector interface and SAMRAI vectors.
//

#define SPVEC_CAST(x)  ( ((PETSc_SAMRAIVectorReal<DIM,TYPE>*) x) )

#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace solv {

template<int DIM, class TYPE> inline
tbox::Pointer< SAMRAIVectorReal<DIM,TYPE> > 
PETSc_SAMRAIVectorReal<DIM,TYPE>::getSAMRAIVector()
{
   return(d_samrai_vector); 
}

template<int DIM, class TYPE> inline 
double PETSc_SAMRAIVectorReal<DIM,TYPE>::dotWith(
   const PETScAbstractVectorReal<TYPE>* y,
   bool local_only) const
{
   return( d_samrai_vector->dot(SPVEC_CAST(y)->getSAMRAIVector(),
				local_only) );
}

template<int DIM, class TYPE> inline
double PETSc_SAMRAIVectorReal<DIM,TYPE>::TdotWith(
   const PETScAbstractVectorReal<TYPE>* y,
   bool local_only) const
{
   return( d_samrai_vector->dot(SPVEC_CAST(y)->getSAMRAIVector(),
				local_only) );
}

template<int DIM, class TYPE> inline
double PETSc_SAMRAIVectorReal<DIM,TYPE>::L1Norm(bool local_only) const
{
   return( d_samrai_vector->L1Norm(local_only) );
}

template<int DIM, class TYPE> inline
double PETSc_SAMRAIVectorReal<DIM,TYPE>::L2Norm(bool local_only) const
{
   return( d_samrai_vector->L2Norm(local_only) );
}

template<int DIM, class TYPE> inline
double PETSc_SAMRAIVectorReal<DIM,TYPE>::maxNorm(bool local_only) const
{
   return( d_samrai_vector->maxNorm(local_only) );
}

template<int DIM, class TYPE> inline
void PETSc_SAMRAIVectorReal<DIM,TYPE>::scaleVector(const TYPE alpha)
{
   d_samrai_vector->scale(alpha, d_samrai_vector);
}

template<int DIM, class TYPE> inline
void PETSc_SAMRAIVectorReal<DIM,TYPE>::copyVector(
   const PETScAbstractVectorReal<TYPE>* v_src)
{
   d_samrai_vector->copyVector(SPVEC_CAST(v_src)->getSAMRAIVector());
}

template<int DIM, class TYPE> inline
void PETSc_SAMRAIVectorReal<DIM,TYPE>::setToScalar(const TYPE alpha)
{
   d_samrai_vector->setToScalar(alpha);
}

template<int DIM, class TYPE> inline
void PETSc_SAMRAIVectorReal<DIM,TYPE>::setRandomValues(const TYPE width,
                                                         const TYPE low) 
{
   d_samrai_vector->setRandomValues(width, low);
}

template<int DIM, class TYPE> inline
void PETSc_SAMRAIVectorReal<DIM,TYPE>::swapWith(
   PETScAbstractVectorReal<TYPE>* v_other)
{
   d_samrai_vector->swapVectors(SPVEC_CAST(v_other)->getSAMRAIVector());
}

template<int DIM, class TYPE> inline
void PETSc_SAMRAIVectorReal<DIM,TYPE>::setAXPY(
   const TYPE alpha, 
   const PETScAbstractVectorReal<TYPE>* x)
{
   d_samrai_vector->axpy(alpha, 
                         SPVEC_CAST(x)->getSAMRAIVector(), 
                         d_samrai_vector);
}

template<int DIM, class TYPE> inline
void PETSc_SAMRAIVectorReal<DIM,TYPE>::setAXPBY(
   const TYPE alpha, 
   const PETScAbstractVectorReal<TYPE>* x, 
   const TYPE beta) 
{
   d_samrai_vector->linearSum(alpha, SPVEC_CAST(x)->getSAMRAIVector(),
                              beta, d_samrai_vector);
}

template<int DIM, class TYPE> inline
void PETSc_SAMRAIVectorReal<DIM,TYPE>::setWAXPY(
   const TYPE alpha, 
   const PETScAbstractVectorReal<TYPE>* x,
   const PETScAbstractVectorReal<TYPE>* y)
{
   d_samrai_vector->axpy(alpha, SPVEC_CAST(x)->getSAMRAIVector(), 
                                SPVEC_CAST(y)->getSAMRAIVector());
}

template<int DIM, class TYPE> inline
void PETSc_SAMRAIVectorReal<DIM,TYPE>::pointwiseMultiply(
   const PETScAbstractVectorReal<TYPE>* x,
   const PETScAbstractVectorReal<TYPE>* y)
{
   d_samrai_vector->multiply(SPVEC_CAST(x)->getSAMRAIVector(), 
                             SPVEC_CAST(y)->getSAMRAIVector());
}

template<int DIM, class TYPE> inline
void PETSc_SAMRAIVectorReal<DIM,TYPE>::pointwiseDivide(
   const PETScAbstractVectorReal<TYPE>* x,
   const PETScAbstractVectorReal<TYPE>* y)
{
   d_samrai_vector->divide(SPVEC_CAST(x)->getSAMRAIVector(), 
                           SPVEC_CAST(y)->getSAMRAIVector());
}

template<int DIM, class TYPE> inline
void PETSc_SAMRAIVectorReal<DIM,TYPE>::vecMax(int& i, 
                                                TYPE& max) const
{
   max = d_samrai_vector->max();
   // Note: This is a bogus index value!
   //       Hopefully, PETSc doesn't use it for anything.
   i = 0;
}

template<int DIM, class TYPE> inline
void PETSc_SAMRAIVectorReal<DIM,TYPE>::vecMin(int& i, 
                                                TYPE& min) const
{
   min = d_samrai_vector->min();
   // Note: This is a bogus index value!
   //       Hopefully, PETSc doesn't use it for anything.
   i = 0;
}

template<int DIM, class TYPE> inline
int PETSc_SAMRAIVectorReal<DIM,TYPE>::getDataSize() const
{
   // Note: This is a bogus value!
   //       But, PETSc requires some value to be returned.
   //       Hopefully, this will not cause problems.
   return(0);
}

template<int DIM, class TYPE> inline
int PETSc_SAMRAIVectorReal<DIM,TYPE>::getLocalDataSize() const
{
   // Note: This is a bogus value!
   //       But, PETSc requires some value to be returned.
   //       Hopefully, this will not cause problems.
   return(0);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
