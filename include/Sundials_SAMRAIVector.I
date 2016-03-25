//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/solvers/vectors/Sundials_SAMRAIVector.I $
// Package:     SAMRAI solvers
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1917 $
// Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description: "Glue code" between SAMRAI vector object and Sundials vectors.
//

#define SKVEC_CAST(x)  ( ((Sundials_SAMRAIVector<DIM>*) x) )

#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace solv {

template<int DIM> inline
tbox::Pointer< SAMRAIVectorReal<DIM,double> > 
Sundials_SAMRAIVector<DIM>::getSAMRAIVector()
{
   return(d_samrai_vector);
}

template<int DIM> inline
void Sundials_SAMRAIVector<DIM>::setToScalar(const double c)
{
   d_samrai_vector->setToScalar(c);
}

template<int DIM> inline
void 
Sundials_SAMRAIVector<DIM>::scaleVector(const SundialsAbstractVector* x,
                                       const double c)
{
   d_samrai_vector->scale(c, SKVEC_CAST(x)->getSAMRAIVector());
}

template<int DIM> inline
void
Sundials_SAMRAIVector<DIM>::setLinearSum(const double a, 
                                        const SundialsAbstractVector* x,
                                        const double b, 
                                        const SundialsAbstractVector* y)
{
   d_samrai_vector->linearSum(a, SKVEC_CAST(x)->getSAMRAIVector(),
                              b, SKVEC_CAST(y)->getSAMRAIVector());
}

template<int DIM> inline
void Sundials_SAMRAIVector<DIM>::pointwiseMultiply(
   const SundialsAbstractVector* x, 
   const SundialsAbstractVector* y)
{
   d_samrai_vector->multiply(SKVEC_CAST(x)->getSAMRAIVector(),
                                      SKVEC_CAST(y)->getSAMRAIVector());
}

template<int DIM> inline
void
Sundials_SAMRAIVector<DIM>::pointwiseDivide(const SundialsAbstractVector* x,
                                           const SundialsAbstractVector* y)
{
   d_samrai_vector->divide(SKVEC_CAST(x)->getSAMRAIVector(),
                                    SKVEC_CAST(y)->getSAMRAIVector());
}

template<int DIM> inline
void Sundials_SAMRAIVector<DIM>::setAbs(const SundialsAbstractVector* x)
{
   d_samrai_vector->abs(SKVEC_CAST(x)->getSAMRAIVector());
}

template<int DIM> inline
void Sundials_SAMRAIVector<DIM>::pointwiseReciprocal(
   const SundialsAbstractVector* x)
{
   d_samrai_vector->reciprocal(SKVEC_CAST(x)->getSAMRAIVector());
}

template<int DIM> inline
void Sundials_SAMRAIVector<DIM>::addScalar(const SundialsAbstractVector* x,
                                          const double b)
{
   d_samrai_vector->addScalar(SKVEC_CAST(x)->getSAMRAIVector(), b);
}

template<int DIM> inline
double 
Sundials_SAMRAIVector<DIM>::dotWith(const SundialsAbstractVector* x) const
{
   return( d_samrai_vector->dot(SKVEC_CAST(x)->getSAMRAIVector()) );
}

template<int DIM> inline
double Sundials_SAMRAIVector<DIM>::maxNorm() const
{
   return( d_samrai_vector->maxNorm() );
}

template<int DIM> inline
double Sundials_SAMRAIVector<DIM>::L1Norm() const
{
   return( d_samrai_vector->L1Norm() );
}

template<int DIM> inline
double
Sundials_SAMRAIVector<DIM>::weightedL2Norm(
   const SundialsAbstractVector* x) const
{
   return( d_samrai_vector->weightedL2Norm(SKVEC_CAST(x)->getSAMRAIVector()) );
}

template<int DIM> inline
double
Sundials_SAMRAIVector<DIM>::weightedRMSNorm(const SundialsAbstractVector* x) const
{
   return(d_samrai_vector->weightedRMSNorm(SKVEC_CAST(x)->getSAMRAIVector()));
}

template<int DIM> inline
double Sundials_SAMRAIVector<DIM>::vecMin() const
{
   return( d_samrai_vector->min() );
}

template<int DIM> inline
int Sundials_SAMRAIVector<DIM>::constrProdPos(
   const SundialsAbstractVector* x) const
{
   return( d_samrai_vector->
           computeConstrProdPos(SKVEC_CAST(x)->getSAMRAIVector()) );
}

template<int DIM> inline
void Sundials_SAMRAIVector<DIM>::compareToScalar(
   const SundialsAbstractVector* x, 
   const double c)
{
   d_samrai_vector->compareToScalar(SKVEC_CAST(x)->getSAMRAIVector(), c);
}

template<int DIM> inline
int Sundials_SAMRAIVector<DIM>::testReciprocal(
   const SundialsAbstractVector* x)
{
   return( d_samrai_vector->testReciprocal(SKVEC_CAST(x)->getSAMRAIVector()) );
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
