//
// File:        PVodeTrio_SAMRAIVector.I
// Package:     SAMRAI solvers
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: "Glue code" between SAMRAI vector object and PVodeTrio vectors.
//

#define SKVEC_CAST(x)  ( ((PVodeTrio_SAMRAIVector<DIM>*) x) )

#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace solv {

template<int DIM> inline
tbox::Pointer< SAMRAIVectorReal<DIM,double> > 
PVodeTrio_SAMRAIVector<DIM>::getSAMRAIVector()
{
   return(d_samrai_vector);
}

template<int DIM> inline
void PVodeTrio_SAMRAIVector<DIM>::setToScalar(const double c)
{
   d_samrai_vector->setToScalar(c);
}

template<int DIM> inline
void 
PVodeTrio_SAMRAIVector<DIM>::scaleVector(const PVodeTrioAbstractVector* x,
                                       const double c)
{
   d_samrai_vector->scale(c, SKVEC_CAST(x)->getSAMRAIVector());
}

template<int DIM> inline
void
PVodeTrio_SAMRAIVector<DIM>::setLinearSum(const double a, 
                                        const PVodeTrioAbstractVector* x,
                                        const double b, 
                                        const PVodeTrioAbstractVector* y)
{
   d_samrai_vector->linearSum(a, SKVEC_CAST(x)->getSAMRAIVector(),
                              b, SKVEC_CAST(y)->getSAMRAIVector());
}

template<int DIM> inline
void PVodeTrio_SAMRAIVector<DIM>::pointwiseMultiply(
   const PVodeTrioAbstractVector* x, 
   const PVodeTrioAbstractVector* y)
{
   d_samrai_vector->multiply(SKVEC_CAST(x)->getSAMRAIVector(),
                                      SKVEC_CAST(y)->getSAMRAIVector());
}

template<int DIM> inline
void
PVodeTrio_SAMRAIVector<DIM>::pointwiseDivide(const PVodeTrioAbstractVector* x,
                                           const PVodeTrioAbstractVector* y)
{
   d_samrai_vector->divide(SKVEC_CAST(x)->getSAMRAIVector(),
                                    SKVEC_CAST(y)->getSAMRAIVector());
}

template<int DIM> inline
void PVodeTrio_SAMRAIVector<DIM>::setAbs(const PVodeTrioAbstractVector* x)
{
   d_samrai_vector->abs(SKVEC_CAST(x)->getSAMRAIVector());
}

template<int DIM> inline
void PVodeTrio_SAMRAIVector<DIM>::pointwiseReciprocal(
   const PVodeTrioAbstractVector* x)
{
   d_samrai_vector->reciprocal(SKVEC_CAST(x)->getSAMRAIVector());
}

template<int DIM> inline
void PVodeTrio_SAMRAIVector<DIM>::addScalar(const PVodeTrioAbstractVector* x,
                                          const double b)
{
   d_samrai_vector->addScalar(SKVEC_CAST(x)->getSAMRAIVector(), b);
}

template<int DIM> inline
double 
PVodeTrio_SAMRAIVector<DIM>::dotWith(const PVodeTrioAbstractVector* x) const
{
   return( d_samrai_vector->dot(SKVEC_CAST(x)->getSAMRAIVector()) );
}

template<int DIM> inline
double PVodeTrio_SAMRAIVector<DIM>::maxNorm() const
{
   return( d_samrai_vector->maxNorm() );
}

template<int DIM> inline
double PVodeTrio_SAMRAIVector<DIM>::L1Norm() const
{
   return( d_samrai_vector->L1Norm() );
}

template<int DIM> inline
double
PVodeTrio_SAMRAIVector<DIM>::weightedL2Norm(
   const PVodeTrioAbstractVector* x) const
{
   return( d_samrai_vector->weightedL2Norm(SKVEC_CAST(x)->getSAMRAIVector()) );
}

template<int DIM> inline
double
PVodeTrio_SAMRAIVector<DIM>::weightedRMSNorm(const PVodeTrioAbstractVector* x) const
{
   return(d_samrai_vector->weightedRMSNorm(SKVEC_CAST(x)->getSAMRAIVector()));
}

template<int DIM> inline
double PVodeTrio_SAMRAIVector<DIM>::vecMin() const
{
   return( d_samrai_vector->min() );
}

template<int DIM> inline
int PVodeTrio_SAMRAIVector<DIM>::constrProdPos(
   const PVodeTrioAbstractVector* x) const
{
   return( d_samrai_vector->
           computeConstrProdPos(SKVEC_CAST(x)->getSAMRAIVector()) );
}

template<int DIM> inline
void PVodeTrio_SAMRAIVector<DIM>::compareToScalar(
   const PVodeTrioAbstractVector* x, 
   const double c)
{
   d_samrai_vector->compareToScalar(SKVEC_CAST(x)->getSAMRAIVector(), c);
}

template<int DIM> inline
int PVodeTrio_SAMRAIVector<DIM>::testReciprocal(
   const PVodeTrioAbstractVector* x)
{
   return( d_samrai_vector->testReciprocal(SKVEC_CAST(x)->getSAMRAIVector()) );
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
