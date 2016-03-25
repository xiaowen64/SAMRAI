//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/solvers/vectors/PETSc_SAMRAIVectorReal.I $
// Package:     SAMRAI solvers
// Copyright:   (c) 1998 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1732 $
// Modified:    $LastChangedDate: 2007-12-04 09:01:48 -0800 (Tue, 04 Dec 2007) $
// Description: "Glue code" between PETSc vector interface and SAMRAI vectors.
//

#define C_PSVEC_CAST(x) (dynamic_cast<const PETSc_SAMRAIVectorReal<DIM,TYPE>*>(x))

#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace solv {


template <int DIM, class TYPE>
inline
double
PETSc_SAMRAIVectorReal<DIM, TYPE>::dotWith(
    const PETScAbstractVectorReal<TYPE>* y,
    bool local_only) const
{
    return d_samrai_vector->dot(C_PSVEC_CAST(y)->d_samrai_vector, local_only);
}// dotWith

template <int DIM, class TYPE>
inline
double
PETSc_SAMRAIVectorReal<DIM,TYPE>::TdotWith(
    const PETScAbstractVectorReal<TYPE>* y,
    bool local_only) const
{
    return d_samrai_vector->dot(C_PSVEC_CAST(y)->d_samrai_vector, local_only);
}// TdotWith

template <int DIM, class TYPE>
inline
double
PETSc_SAMRAIVectorReal<DIM,TYPE>::L1Norm(
    bool local_only) const
{
    return d_samrai_vector->L1Norm(local_only);
}// L1Norm

template <int DIM, class TYPE>
inline
double
PETSc_SAMRAIVectorReal<DIM,TYPE>::L2Norm(
    bool local_only) const
{
    return d_samrai_vector->L2Norm(local_only);
}// L2Norm

template <int DIM, class TYPE>
inline
double
PETSc_SAMRAIVectorReal<DIM,TYPE>::maxNorm(
    bool local_only) const
{
    return d_samrai_vector->maxNorm(local_only);
}// maxNorm

template <int DIM, class TYPE>
inline
void
PETSc_SAMRAIVectorReal<DIM,TYPE>::scaleVector(
    const TYPE alpha)
{
    d_samrai_vector->scale(alpha, d_samrai_vector);
    return;
}// scaleVector

template <int DIM, class TYPE>
inline
void
PETSc_SAMRAIVectorReal<DIM,TYPE>::copyVector(
    const PETScAbstractVectorReal<TYPE>* v_src)
{
    d_samrai_vector->copyVector(C_PSVEC_CAST(v_src)->d_samrai_vector);
    return;
}// copyVector

template <int DIM, class TYPE>
inline
void
PETSc_SAMRAIVectorReal<DIM,TYPE>::setToScalar(
    const TYPE alpha)
{
    d_samrai_vector->setToScalar(alpha);
    return;
}// setToScalar

template <int DIM, class TYPE>
inline
void
PETSc_SAMRAIVectorReal<DIM,TYPE>::swapWith(
    PETScAbstractVectorReal<TYPE>* v_other)
{
    d_samrai_vector->swapVectors(C_PSVEC_CAST(v_other)->d_samrai_vector);
    return;
}// swapWith

template <int DIM, class TYPE>
inline
void
PETSc_SAMRAIVectorReal<DIM,TYPE>::setAXPY(
    const TYPE alpha,
    const PETScAbstractVectorReal<TYPE>* x)
{
    d_samrai_vector->axpy(alpha, C_PSVEC_CAST(x)->d_samrai_vector, d_samrai_vector);
    return;
}// setAXPY

template <int DIM, class TYPE>
inline
void
PETSc_SAMRAIVectorReal<DIM,TYPE>::setAXPBY(
    const TYPE alpha,
    const PETScAbstractVectorReal<TYPE>* x,
    const TYPE beta)
{
    d_samrai_vector->linearSum(alpha, C_PSVEC_CAST(x)->d_samrai_vector, beta, d_samrai_vector);
    return;
}// setAXPBY

template <int DIM, class TYPE>
inline
void
PETSc_SAMRAIVectorReal<DIM,TYPE>::setWAXPY(
    const TYPE alpha,
    const PETScAbstractVectorReal<TYPE>* x,
    const PETScAbstractVectorReal<TYPE>* y)
{
    d_samrai_vector->axpy(alpha, C_PSVEC_CAST(x)->d_samrai_vector, C_PSVEC_CAST(y)->d_samrai_vector);
    return;
}// setWAXPY

template <int DIM, class TYPE>
inline
void
PETSc_SAMRAIVectorReal<DIM,TYPE>::pointwiseMultiply(
    const PETScAbstractVectorReal<TYPE>* x,
    const PETScAbstractVectorReal<TYPE>* y)
{
    d_samrai_vector->multiply(C_PSVEC_CAST(x)->d_samrai_vector, C_PSVEC_CAST(y)->d_samrai_vector);
    return;
}// pointwiseMultiply

template <int DIM, class TYPE>
inline
void
PETSc_SAMRAIVectorReal<DIM,TYPE>::pointwiseDivide(
    const PETScAbstractVectorReal<TYPE>* x,
    const PETScAbstractVectorReal<TYPE>* y)
{
    d_samrai_vector->divide(C_PSVEC_CAST(x)->d_samrai_vector, C_PSVEC_CAST(y)->d_samrai_vector);
    return;
}// pointwiseDivide

template <int DIM, class TYPE>
inline
double
PETSc_SAMRAIVectorReal<DIM,TYPE>::maxPointwiseDivide(
    const PETScAbstractVectorReal<TYPE>* y)
{
    return d_samrai_vector->maxPointwiseDivide(C_PSVEC_CAST(y)->d_samrai_vector);
}// maxPointwiseDivide

template <int DIM, class TYPE>
inline
void
PETSc_SAMRAIVectorReal<DIM,TYPE>::vecMax(
    int& i,
    TYPE& max) const
{
    static const bool interior_only = true;
    max = d_samrai_vector->max(interior_only);
    // Note: This is a bogus index value!
    //       Hopefully, PETSc doesn't use it for anything.
    i = 0;
    return;
}// vecMax

template <int DIM, class TYPE>
inline
void
PETSc_SAMRAIVectorReal<DIM,TYPE>::vecMin(
    int& i,
    TYPE& min) const
{
    static const bool interior_only = true;
    min = d_samrai_vector->min(interior_only);
    // Note: This is a bogus index value!
    //       Hopefully, PETSc doesn't use it for anything.
    i = 0;
    return;
}// vecMin

template <int DIM, class TYPE>
inline
void
PETSc_SAMRAIVectorReal<DIM,TYPE>::setRandomValues(
    const TYPE width,
    const TYPE low)
{
    d_samrai_vector->setRandomValues(width, low);
    return;
}// setRandomValues

template <int DIM, class TYPE>
inline
void
PETSc_SAMRAIVectorReal<DIM,TYPE>::getDataArray(
    TYPE** array)
{
    *array = NULL;
    return;
}// getDataArray

template <int DIM, class TYPE>
inline
void
PETSc_SAMRAIVectorReal<DIM,TYPE>::restoreDataArray(
    TYPE** array)
{
    *array = NULL;
    return;
}// restoreDataArray

template <int DIM, class TYPE>
inline
int
PETSc_SAMRAIVectorReal<DIM,TYPE>::getDataSize() const
{
    // Note: This is a bogus value!
    //       But, PETSc requires some value to be returned.
    //       Hopefully, this will not cause problems.
    return 0;
}// getDataSize

template <int DIM, class TYPE>
inline
int
PETSc_SAMRAIVectorReal<DIM,TYPE>::getLocalDataSize() const
{
    // Note: This is a bogus value!
    //       But, PETSc requires some value to be returned.
    //       Hopefully, this will not cause problems.
    return 0;
}// getLocalDataSize


#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
