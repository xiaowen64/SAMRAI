/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   "Glue code" between SAMRAI vector object and Sundials vectors. 
 *
 ************************************************************************/

#define SKVEC_CAST(x) (((Sundials_SAMRAIVector *)x))

namespace SAMRAI {
namespace solv {

SAMRAI_INLINE_KEYWORD
tbox::Pointer<SAMRAIVectorReal<double> >
Sundials_SAMRAIVector::getSAMRAIVector()
{
   return d_samrai_vector;
}

SAMRAI_INLINE_KEYWORD
void Sundials_SAMRAIVector::setToScalar(
   const double c)
{
   d_samrai_vector->setToScalar(c);
}

SAMRAI_INLINE_KEYWORD
void
Sundials_SAMRAIVector::scaleVector(
   const SundialsAbstractVector* x,
   const double c)
{
   d_samrai_vector->scale(c, SKVEC_CAST(x)->getSAMRAIVector());
}

SAMRAI_INLINE_KEYWORD
void
Sundials_SAMRAIVector::setLinearSum(
   const double a,
   const SundialsAbstractVector* x,
   const double b,
   const SundialsAbstractVector* y)
{
   d_samrai_vector->linearSum(a, SKVEC_CAST(x)->getSAMRAIVector(),
      b, SKVEC_CAST(y)->getSAMRAIVector());
}

SAMRAI_INLINE_KEYWORD
void Sundials_SAMRAIVector::pointwiseMultiply(
   const SundialsAbstractVector* x,
   const SundialsAbstractVector* y)
{
   d_samrai_vector->multiply(SKVEC_CAST(x)->getSAMRAIVector(),
      SKVEC_CAST(y)->getSAMRAIVector());
}

SAMRAI_INLINE_KEYWORD
void
Sundials_SAMRAIVector::pointwiseDivide(
   const SundialsAbstractVector* x,
   const SundialsAbstractVector* y)
{
   d_samrai_vector->divide(SKVEC_CAST(x)->getSAMRAIVector(),
      SKVEC_CAST(y)->getSAMRAIVector());
}

SAMRAI_INLINE_KEYWORD
void Sundials_SAMRAIVector::setAbs(
   const SundialsAbstractVector* x)
{
   d_samrai_vector->abs(SKVEC_CAST(x)->getSAMRAIVector());
}

SAMRAI_INLINE_KEYWORD
void Sundials_SAMRAIVector::pointwiseReciprocal(
   const SundialsAbstractVector* x)
{
   d_samrai_vector->reciprocal(SKVEC_CAST(x)->getSAMRAIVector());
}

SAMRAI_INLINE_KEYWORD
void Sundials_SAMRAIVector::addScalar(
   const SundialsAbstractVector* x,
   const double b)
{
   d_samrai_vector->addScalar(SKVEC_CAST(x)->getSAMRAIVector(), b);
}

SAMRAI_INLINE_KEYWORD
double
Sundials_SAMRAIVector::dotWith(
   const SundialsAbstractVector* x) const
{
   return d_samrai_vector->dot(SKVEC_CAST(x)->getSAMRAIVector());
}

SAMRAI_INLINE_KEYWORD
double Sundials_SAMRAIVector::maxNorm() const
{
   return d_samrai_vector->maxNorm();
}

SAMRAI_INLINE_KEYWORD
double Sundials_SAMRAIVector::L1Norm() const
{
   return d_samrai_vector->L1Norm();
}

SAMRAI_INLINE_KEYWORD
double
Sundials_SAMRAIVector::weightedL2Norm(
   const SundialsAbstractVector* x) const
{
   return d_samrai_vector->weightedL2Norm(SKVEC_CAST(x)->getSAMRAIVector());
}

SAMRAI_INLINE_KEYWORD
double
Sundials_SAMRAIVector::weightedRMSNorm(
   const SundialsAbstractVector* x) const
{
   return d_samrai_vector->weightedRMSNorm(SKVEC_CAST(x)->getSAMRAIVector());
}

SAMRAI_INLINE_KEYWORD
double Sundials_SAMRAIVector::vecMin() const
{
   return d_samrai_vector->min();
}

SAMRAI_INLINE_KEYWORD
int Sundials_SAMRAIVector::constrProdPos(
   const SundialsAbstractVector* x) const
{
   return d_samrai_vector->
          computeConstrProdPos(SKVEC_CAST(x)->getSAMRAIVector());
}

SAMRAI_INLINE_KEYWORD
void Sundials_SAMRAIVector::compareToScalar(
   const SundialsAbstractVector* x,
   const double c)
{
   d_samrai_vector->compareToScalar(SKVEC_CAST(x)->getSAMRAIVector(), c);
}

SAMRAI_INLINE_KEYWORD
int Sundials_SAMRAIVector::testReciprocal(
   const SundialsAbstractVector* x)
{
   return d_samrai_vector->testReciprocal(SKVEC_CAST(x)->getSAMRAIVector());
}

}
}
