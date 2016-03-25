//
// File:        PETScAbstractVectorReal.C
// Package:     SAMRAI solvers
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Interface to C++ vector implementation for PETSc package.
//

#include "PETScAbstractVectorReal.h"

#ifdef HAVE_PETSC


#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

extern "C" {
#include "vecimpl.h"
}
#include "tbox/IOStream.h"
#include "tbox/PIO.h"
#include "tbox/Utilities.h"

#ifndef NULL
#define NULL (0)
#endif

namespace SAMRAI {
    namespace solv {

#define ABSPVEC_CAST(v)     ( ((PETScAbstractVectorReal<TYPE>*)(v->data)) )

/*
*************************************************************************
*                                                                       *
* Static member functions to link with PETSc package.                   *
*                                                                       *
*************************************************************************
*/

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::destroyVec(Vec v) 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(v == (Vec)NULL));
#endif
   ABSPVEC_CAST(v)->freeVector();
   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::viewVec(Vec v, PetscViewer view)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(v == (Vec)NULL));
#endif
   (void) view;
   ABSPVEC_CAST(v)->viewVector();
   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::duplicateVec(Vec v_in, Vec* v_new)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(v_in == (Vec)NULL));
#endif
   PETScAbstractVectorReal<TYPE>* v = ABSPVEC_CAST(v_in)->makeNewVector();  
   *v_new = v->getPETScVector();
   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::duplicateVecs(Vec v_in, 
                                                      int n, Vec** varr_new)
{
   int ierr = 0;
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(v_in == (Vec)NULL));
#endif
   ierr = PetscMalloc( n * sizeof(Vec *), varr_new);
   PETSC_SAMRAI_ERROR(ierr);
   // Get rid of KCC warning
   if(ierr);
   for (int i = 0; i < n; i++) {
      duplicateVec(v_in, *varr_new+i);
   }
   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::destroyVecs(const Vec* v_arr, int n)
{
   int i;
#ifdef DEBUG_CHECK_ASSERTIONS
   for (i = 0; i < n; i++) {
      assert(!(v_arr[i] == (Vec)NULL));
   }
#endif
   for (i = 0; i < n; i++) {
      ABSPVEC_CAST(v_arr[i])->freeVector();
   }

   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::dotProduct(Vec x, Vec y, TYPE* dp)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(x == (Vec)NULL));
   assert(!(y == (Vec)NULL));
#endif
   *dp = ABSPVEC_CAST(x)->dotWith(ABSPVEC_CAST(y));  
   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::dotProductT(Vec x, Vec y, TYPE* dp)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(x == (Vec)NULL));
   assert(!(y == (Vec)NULL));
#endif
   *dp = ABSPVEC_CAST(x)->TdotWith(ABSPVEC_CAST(y));
   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::dotProductM(int n, Vec x, const Vec* y,
                                                    TYPE* dp)
{
   int i;
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(x == (Vec)NULL));
   for (i = 0; i < n; i++) {
      assert(!(y[i] == (Vec)NULL));
   }
#endif
   for (i = 0; i < n; i++) {
      dp[i] = ABSPVEC_CAST(x)->dotWith(ABSPVEC_CAST(y[i]));
   }
   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::dotProductMT(int n, Vec x, const Vec* y,
                                                     TYPE* dp)
{
   int i;
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(x == (Vec)NULL));
   for (i = 0; i < n; i++) {
      assert(!(y[i] == (Vec)NULL));
   }
#endif
   for (i = 0; i < n; i++) {
      dp[i] = ABSPVEC_CAST(x)->TdotWith(ABSPVEC_CAST(y[i]));
   }
   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::vecNorm(Vec x, 
                                                NormType n_type, double* norm)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(x == (Vec)NULL));
#endif
   if (n_type == NORM_1) {
      *norm = ABSPVEC_CAST(x)->L1Norm();
   } else if (n_type == NORM_2) {
      *norm = ABSPVEC_CAST(x)->L2Norm();
   } else if (n_type == NORM_INFINITY) {
      *norm =  ABSPVEC_CAST(x)->maxNorm();
   } else if (n_type == NORM_1_AND_2) {
      norm[0] = ABSPVEC_CAST(x)->L1Norm();
      norm[1] = ABSPVEC_CAST(x)->L2Norm(); 
   } else {
      TBOX_ERROR("vector norm type undefined!");
   }
   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::scaleVec(const TYPE* alpha, Vec x) 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(x == (Vec)NULL));
#endif
   ABSPVEC_CAST(x)->scaleVector(*alpha);
   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::copyVec(Vec v_src, Vec v_dst)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(v_src == (Vec)NULL));
   assert(!(v_dst == (Vec)NULL));
#endif
   ABSPVEC_CAST(v_dst)->copyVector(ABSPVEC_CAST(v_src));
   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::setVec(const TYPE* alpha, Vec x)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(x == (Vec)NULL));
#endif
   ABSPVEC_CAST(x)->setToScalar(*alpha);
   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::swapVecs(Vec x, Vec y)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(x == (Vec)NULL));
   assert(!(y == (Vec)NULL));
#endif
   ABSPVEC_CAST(x)->swapWith(ABSPVEC_CAST(y));
   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::computeAXPY(const TYPE* alpha, 
                                                    Vec x, Vec y)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(x == (Vec)NULL));
   assert(!(y == (Vec)NULL));
#endif
   ABSPVEC_CAST(y)->setAXPY(*alpha, ABSPVEC_CAST(x));
   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::computeAXPBY(const TYPE* alpha,
                                                     const TYPE* beta, 
                                                     Vec x, Vec y)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(x == (Vec)NULL));
   assert(!(y == (Vec)NULL));
#endif
   ABSPVEC_CAST(y)->setAXPBY(*alpha, ABSPVEC_CAST(x), *beta);
   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::computeMAXPY(int n, const TYPE* alpha, 
                                                     Vec x, Vec* y)
{
   int i;
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(x == (Vec)NULL));
   for (i = 0; i < n; i++) {
      assert(!(y[i] == (Vec)NULL));
   }
#endif
   for (i = 0; i < n; i++) {
      ABSPVEC_CAST(x)->setAXPY(alpha[i], ABSPVEC_CAST(y[i]));
   }
   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::computeAYPX(const TYPE* alpha,
                                                    Vec x, Vec y)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(x == (Vec)NULL));
   assert(!(y == (Vec)NULL));
#endif
   ABSPVEC_CAST(y)->setAXPBY(1.0, ABSPVEC_CAST(x), *alpha);
   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::computeWAXPY(const TYPE* alpha,
                                                     Vec x, Vec y, 
                                                     Vec w)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(x == (Vec)NULL));
   assert(!(y == (Vec)NULL));
   assert(!(w == (Vec)NULL));
#endif
   ABSPVEC_CAST(w)->setWAXPY(*alpha, ABSPVEC_CAST(x), 
                                     ABSPVEC_CAST(y));
   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::pointwiseMultVecs(Vec x, Vec y, Vec w)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(x == (Vec)NULL));
   assert(!(y == (Vec)NULL));
   assert(!(w == (Vec)NULL));
#endif
   ABSPVEC_CAST(w)->pointwiseMultiply(ABSPVEC_CAST(x), ABSPVEC_CAST(y));
   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::pointwiseDivideVecs(Vec x, Vec y, Vec w)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(x == (Vec)NULL));
   assert(!(y == (Vec)NULL));
   assert(!(w == (Vec)NULL));
#endif
   ABSPVEC_CAST(w)->pointwiseDivide(ABSPVEC_CAST(x), ABSPVEC_CAST(y));
   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::maxPointwiseDivideVecs(Vec x, Vec y, double *max)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(x == (Vec)NULL));
   assert(!(y == (Vec)NULL));
#endif
   *max = ABSPVEC_CAST(x)->maxPointwiseDivide(ABSPVEC_CAST(y));
   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::vectorMax(Vec x, int* i, double* max)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(x == (Vec)NULL));
#endif
   ABSPVEC_CAST(x)->vecMax(*i, *max);
   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::vectorMin(Vec x, int* i, double* min)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(x == (Vec)NULL));
#endif
   ABSPVEC_CAST(x)->vecMin(*i, *min);
   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::setVecRandom(PetscRandom r, Vec x)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(x == (Vec)NULL));
#endif
   (void) r;
   TYPE low = 0.0;
   TYPE width = 1.0;
// if (r->iset == PETSC_TRUE) {
//    width = r->width;
//    low   = r->low;
// }
   ABSPVEC_CAST(x)->setRandomValues(width, low);
   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::getVecArray(Vec x, TYPE** array)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(x == (Vec)NULL));
#endif
   ABSPVEC_CAST(x)->getDataArray(*array);
   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::getVecSize(Vec x, int* size)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(x == (Vec)NULL));
#endif
   *size = ABSPVEC_CAST(x)->getDataSize();
   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::getLocalVecSize(Vec x, int* size)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(x == (Vec)NULL));
#endif
   *size = ABSPVEC_CAST(x)->getLocalDataSize();
   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::setVecValues(
   Vec x, int ni, const int* indices, const TYPE* vals, InsertMode mode)
{
   (void) x;
   (void) ni;
   (void) indices;
   (void) vals;
   (void) mode;
   TBOX_ERROR("PETScAbstractVectorReal<TYPE>::setVecValues undefined!");   
   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::beginVecAssembly(Vec x)
{
   (void) x;
   TBOX_ERROR("PETScAbstractVectorReal<TYPE>::beginVecAssembly undefined!");
   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::endVecAssembly(Vec x)
{
   (void) x;
   TBOX_ERROR("PETScAbstractVectorReal<TYPE>::endVecAssembly undefined!");
   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::restoreVecArray(Vec x, TYPE** array)
{
   (void) x;
   (void) array;
   TBOX_ERROR("PETScAbstractVectorReal<TYPE>::restoreVecArray undefined!");
   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::setVecOption(Vec x, VecOption option)
{
   (void) x;
   (void) option;
   TBOX_ERROR("PETScAbstractVectorReal<TYPE>::setVecOption undefined!");
   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::setVecValuesBlocked(
   Vec x, int n, const int* nb, const TYPE* vals, InsertMode mode)
{
   (void) x;
   (void) n;
   (void) nb;
   (void) vals;
   (void) mode;
   TBOX_ERROR("PETScAbstractVectorReal<TYPE>::setVecValuesBlocked undefined!");
   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::placeArray(Vec x, const TYPE *vals) 
{
   TBOX_ERROR("PETScAbstractVectorReal<TYPE>::placeArray undefined!");
   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::replaceArray(Vec x, const TYPE *vals) 
{
   TBOX_ERROR("PETScAbstractVectorReal<TYPE>::replaceArray undefined!");
   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::dotProductLocal(Vec x, Vec y, TYPE* dp) 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(x == (Vec)NULL));
   assert(!(y == (Vec)NULL));
#endif
   *dp = ABSPVEC_CAST(x)->dotWith(ABSPVEC_CAST(y),true);
   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::dotProductLocalT(
   Vec x, Vec y, TYPE* dp) 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(x == (Vec)NULL));
   assert(!(y == (Vec)NULL));
#endif
   *dp = ABSPVEC_CAST(x)->TdotWith(ABSPVEC_CAST(y),true);
   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::vecNormLocal(
   Vec x, NormType n_type, double* norm) 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(x == (Vec)NULL));
#endif
   if (n_type == NORM_1) {
      *norm = ABSPVEC_CAST(x)->L1Norm(true);
   } else if (n_type == NORM_2) {
      *norm = ABSPVEC_CAST(x)->L2Norm(true);
   } else if (n_type == NORM_INFINITY) {
      *norm =  ABSPVEC_CAST(x)->maxNorm(true);
   } else if (n_type == NORM_1_AND_2) {
      norm[0] = ABSPVEC_CAST(x)->L1Norm(true);
      norm[1] = ABSPVEC_CAST(x)->L2Norm(true); 
   } else {
      TBOX_ERROR("vector norm type undefined!");
   }
   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::loadIntoVector(PetscViewer view, Vec x) 
{
   TBOX_ERROR("PETScAbstractVectorReal<TYPE>::loadIntoVector undefined!");
   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::reciprocal(Vec x) 
{
   TBOX_ERROR("PETScAbstractVectorReal<TYPE>::reciprocal undefined!");
   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::viewNative(Vec x, PetscViewer view) 
{
   TBOX_ERROR("PETScAbstractVectorReal<TYPE>::viewNative undefined!");
   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::conjugate(Vec x) 
{
   TBOX_ERROR("PETScAbstractVectorReal<TYPE>::conjugate undefined!");
   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::setLocalToGlobalMapping(
   Vec x, ISLocalToGlobalMapping mapping) 
{
   TBOX_ERROR("PETScAbstractVectorReal<TYPE>::setLocalToGlobalMapping undefined!");
   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::setValuesLocal(
   Vec x, int count, const int *indices, const TYPE *vals, InsertMode mode) 
{
   TBOX_ERROR("PETScAbstractVectorReal<TYPE>::setValuesLocal undefined!");
   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::resetArray(Vec x) 
{
   TBOX_ERROR("PETScAbstractVectorReal<TYPE>::resetArray undefined!");
   return(0);
}

template <class TYPE>
int PETScAbstractVectorReal<TYPE>::setFromOptions(Vec x) 
{
   TBOX_ERROR("PETScAbstractVectorReal<TYPE>::setFromOptions undefined!");
   return(0);
}


/*
*************************************************************************
*                                                                       *
* PETScAbstractVectorReal constructor and destructor.              *
*                                                                       *
*************************************************************************
*/

template <class TYPE>
PETScAbstractVectorReal<TYPE>::PETScAbstractVectorReal(
   bool vector_created_via_duplicate)
: d_vector_created_via_duplicate(vector_created_via_duplicate)
{
   int ierr = 0;

   // To get rid of KCC warning
   if(ierr);

// ierr = VecCreate(MPI_COMM_SELF, &d_petsc_vector); PETSC_SAMRAI_ERROR(ierr);
   ierr = VecCreate(PETSC_COMM_SELF, &d_petsc_vector); PETSC_SAMRAI_ERROR(ierr);

   string my_name = "PETScAbstractVectorReal";
   if (d_petsc_vector->type_name) {
      ierr = PetscFree(d_petsc_vector->type_name);  PETSC_SAMRAI_ERROR(ierr);
   }
   ierr = PetscMalloc(
      (my_name.size()+1)*sizeof(char), &(d_petsc_vector->type_name));
   (void) strcpy(d_petsc_vector->type_name, my_name.c_str());


   // Set PETSc vector data to this abstract vector object
   d_petsc_vector->data                   = this; 

   d_petsc_vector->n                      = 0;
   d_petsc_vector->N                      = 0;

   // The following are vector data members defined in the PETSc "_p_Vec" 
   // structure, but are not used here. These are all set to
   // "-1," "PETSC_FALSE," or "PETSC_NULL" in ${PETSC_DIR}/src/vec/interface,
   // so we won't mess with them here.  Besides, some of these are
   // structs, so if we did mess with them we might cause a memory leak.
   #if 0
   d_petsc_vector->bs                     = 0;
   d_petsc_vector->map                    = 0;
   d_petsc_vector->mapping                = 0;
   d_petsc_vector->bmapping               = 0;
   d_petsc_vector->array_gotton           = PETSC_FALSE;
   d_petsc_vector->stash                  = 0;
   d_petsc_vector->bstash                 = 0;
   d_petsc_vector->petscnative            = PETSC_FALSE;
   d_petsc_vector->esivec                 = 0;
   #endif

   // Assign vector operations to PETSc vector object.
   d_petsc_vector->ops->destroy           = PETScAbstractVectorReal<TYPE>::
                                            destroyVec;
   d_petsc_vector->ops->view              = PETScAbstractVectorReal<TYPE>::
                                            viewVec;
   d_petsc_vector->ops->duplicate         = PETScAbstractVectorReal<TYPE>::
                                            duplicateVec;
   d_petsc_vector->ops->duplicatevecs     = PETScAbstractVectorReal<TYPE>::
                                            duplicateVecs;
   d_petsc_vector->ops->destroyvecs       = PETScAbstractVectorReal<TYPE>::
                                            destroyVecs;
   d_petsc_vector->ops->dot               = PETScAbstractVectorReal<TYPE>::
                                            dotProduct;
   d_petsc_vector->ops->tdot              = PETScAbstractVectorReal<TYPE>::
                                            dotProductT;
   d_petsc_vector->ops->mdot              = PETScAbstractVectorReal<TYPE>::
                                            dotProductM;
   d_petsc_vector->ops->mtdot             = PETScAbstractVectorReal<TYPE>::
                                            dotProductMT;
   d_petsc_vector->ops->norm              = PETScAbstractVectorReal<TYPE>::
                                            vecNorm;
   d_petsc_vector->ops->scale             = PETScAbstractVectorReal<TYPE>::
                                            scaleVec;
   d_petsc_vector->ops->copy              = PETScAbstractVectorReal<TYPE>::
                                            copyVec;
   d_petsc_vector->ops->set               = PETScAbstractVectorReal<TYPE>::
                                            setVec;
   d_petsc_vector->ops->swap              = PETScAbstractVectorReal<TYPE>::
                                            swapVecs;
   d_petsc_vector->ops->axpy              = PETScAbstractVectorReal<TYPE>::
                                            computeAXPY;
   d_petsc_vector->ops->axpby             = PETScAbstractVectorReal<TYPE>::
                                            computeAXPBY;
   d_petsc_vector->ops->maxpy             = PETScAbstractVectorReal<TYPE>::
                                            computeMAXPY;
   d_petsc_vector->ops->aypx              = PETScAbstractVectorReal<TYPE>::
                                            computeAYPX;
   d_petsc_vector->ops->waxpy             = PETScAbstractVectorReal<TYPE>::
                                            computeWAXPY;
   d_petsc_vector->ops->pointwisemult     = PETScAbstractVectorReal<TYPE>::
                                            pointwiseMultVecs;
   d_petsc_vector->ops->pointwisedivide   = PETScAbstractVectorReal<TYPE>::
                                            pointwiseDivideVecs;
#if (PETSC_VERSION_MAJOR>=2)&&(PETSC_VERSION_MINOR>=1)&&(PETSC_VERSION_SUBMINOR>=5)
   d_petsc_vector->ops->maxpointwisedivide= PETScAbstractVectorReal<TYPE>::
                                            maxPointwiseDivideVecs;
#endif
   d_petsc_vector->ops->max               = PETScAbstractVectorReal<TYPE>::
                                            vectorMax;
   d_petsc_vector->ops->min               = PETScAbstractVectorReal<TYPE>::
                                            vectorMin;
   d_petsc_vector->ops->setrandom         = PETScAbstractVectorReal<TYPE>::
                                            setVecRandom;
   d_petsc_vector->ops->getarray          = PETScAbstractVectorReal<TYPE>::
                                            getVecArray;
   d_petsc_vector->ops->getsize           = PETScAbstractVectorReal<TYPE>::
                                            getVecSize;
   d_petsc_vector->ops->getlocalsize      = PETScAbstractVectorReal<TYPE>::
                                            getLocalVecSize;

   // The remaining functions will result in program abort.
   d_petsc_vector->ops->setvalues         = PETScAbstractVectorReal<TYPE>::
                                            setVecValues;
   d_petsc_vector->ops->assemblybegin     = PETScAbstractVectorReal<TYPE>::
                                            beginVecAssembly;
   d_petsc_vector->ops->assemblyend       = PETScAbstractVectorReal<TYPE>::
                                            endVecAssembly;
   d_petsc_vector->ops->restorearray      = PETScAbstractVectorReal<TYPE>::
                                            restoreVecArray;
   d_petsc_vector->ops->setoption         = PETScAbstractVectorReal<TYPE>::
                                            setVecOption;
   d_petsc_vector->ops->setvaluesblocked  = PETScAbstractVectorReal<TYPE>::
                                            setVecValuesBlocked;

   d_petsc_vector->ops->placearray        = PETScAbstractVectorReal<TYPE>::
                                            placeArray;
   d_petsc_vector->ops->replacearray      = PETScAbstractVectorReal<TYPE>::
                                            replaceArray;
   d_petsc_vector->ops->dot_local         = PETScAbstractVectorReal<TYPE>::
                                            dotProductLocal;
   d_petsc_vector->ops->tdot_local        = PETScAbstractVectorReal<TYPE>::
                                            dotProductLocalT;
   d_petsc_vector->ops->norm_local        = PETScAbstractVectorReal<TYPE>::
                                            vecNormLocal;
   d_petsc_vector->ops->loadintovector    = PETScAbstractVectorReal<TYPE>::
                                            loadIntoVector;
   d_petsc_vector->ops->reciprocal        = PETScAbstractVectorReal<TYPE>::
                                            reciprocal;
   d_petsc_vector->ops->viewnative        = PETScAbstractVectorReal<TYPE>::
                                            viewNative;
   d_petsc_vector->ops->conjugate         = PETScAbstractVectorReal<TYPE>::
                                            conjugate;
   d_petsc_vector->ops->setlocaltoglobalmapping 
                                          = PETScAbstractVectorReal<TYPE>::
                                            setLocalToGlobalMapping;
   d_petsc_vector->ops->setvalueslocal    = PETScAbstractVectorReal<TYPE>::
                                            setValuesLocal;
   d_petsc_vector->ops->resetarray        = PETScAbstractVectorReal<TYPE>::
                                            resetArray;
   d_petsc_vector->ops->setfromoptions    = PETScAbstractVectorReal<TYPE>::
                                            setFromOptions;
}

template <class TYPE>
PETScAbstractVectorReal<TYPE>::~PETScAbstractVectorReal()
{
   int ierr = 0;
   // Get rid of KCC warning
   if(ierr);
   if (!d_vector_created_via_duplicate) {
      d_petsc_vector->ops->destroy = 0;
      ierr = VecDestroy(d_petsc_vector); PETSC_SAMRAI_ERROR(ierr);
   }
   d_petsc_vector = 0;
}

template <class TYPE>
Vec PETScAbstractVectorReal<TYPE>::getPETScVector()
{
   return( d_petsc_vector );    
}

}
}

#endif
