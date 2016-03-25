//
// File:        PETSc_SAMRAIVectorReal.C
// Package:     SAMRAI solvers
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 179 $
// Modified:    $Date: 2005-01-20 14:50:51 -0800 (Thu, 20 Jan 2005) $
// Description: "Glue code" between PETSc vector interface and SAMRAI vectors.
//

#ifndef included_solv_PETSc_SAMRAIVectorReal_C
#define included_solv_PETSc_SAMRAIVectorReal_C

#include "PETSc_SAMRAIVectorReal.h"

#ifdef HAVE_PETSC

#include <stdlib.h>


#include "tbox/Utilities.h"
#include "tbox/IOStream.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

#include "tbox/PIO.h"
extern "C" {
#include "vecimpl.h"
}


#ifndef NULL
#define NULL (0)
#endif

#ifdef DEBUG_NO_INLINE
#include "PETSc_SAMRAIVectorReal.I"
#endif
namespace SAMRAI {
    namespace solv {

/*
*************************************************************************
*                                                                       *
* Static public member functions.                                       *
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
Vec PETSc_SAMRAIVectorReal<DIM,TYPE>::createPETScVector(
   tbox::Pointer< SAMRAIVectorReal<DIM,TYPE> > samrai_vec)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!samrai_vec.isNull());
#endif
   bool vector_created_via_duplicate = false;
   PETSc_SAMRAIVectorReal<DIM,TYPE>* spv = 
      new PETSc_SAMRAIVectorReal<DIM,TYPE>(samrai_vec, 
                                            vector_created_via_duplicate);
   return( spv->getPETScVector() );
} 

template<int DIM, class TYPE>
void PETSc_SAMRAIVectorReal<DIM,TYPE>::destroyPETScVector(Vec petsc_vec)
{
   if (petsc_vec) {
       delete ((PETSc_SAMRAIVectorReal<DIM,TYPE>*)(petsc_vec->data));
   }
}

template<int DIM, class TYPE>
tbox::Pointer< SAMRAIVectorReal<DIM,TYPE> > 
PETSc_SAMRAIVectorReal<DIM,TYPE>::getSAMRAIVector(Vec petsc_vec)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(petsc_vec == (Vec)NULL));
#endif
   return( ((PETSc_SAMRAIVectorReal<DIM,TYPE>*)(petsc_vec->data))->
                                                 getSAMRAIVector() );
}

/*
*************************************************************************
*                                                                       *
* Protected constructor and destructor for PETSc_SAMRAIVectorReal<DIM>.*
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
PETSc_SAMRAIVectorReal<DIM,TYPE>::PETSc_SAMRAIVectorReal(
   tbox::Pointer< SAMRAIVectorReal<DIM,TYPE> > samrai_vector,
   bool vector_created_via_duplicate)
: PETScAbstractVectorReal<TYPE>(vector_created_via_duplicate)
{
   d_vector_created_via_duplicate = vector_created_via_duplicate;
   d_samrai_vector = samrai_vector;
}

template<int DIM, class TYPE>
PETSc_SAMRAIVectorReal<DIM,TYPE>::~PETSc_SAMRAIVectorReal()
{
}

/*
*************************************************************************
*                                                                       *
* Other member functions						*
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
PETScAbstractVectorReal<TYPE>* 
PETSc_SAMRAIVectorReal<DIM,TYPE>::makeNewVector()
{
   tbox::Pointer< SAMRAIVectorReal<DIM,TYPE> > sam_vec =
      d_samrai_vector->cloneVector(d_samrai_vector->getName());
   sam_vec->allocateVectorData();
   bool vector_created_via_duplicate = true;
   PETSc_SAMRAIVectorReal<DIM,TYPE>* out_vec =
      new PETSc_SAMRAIVectorReal<DIM,TYPE>(sam_vec,
                                             vector_created_via_duplicate);
   return( out_vec );
}

template<int DIM, class TYPE>
void PETSc_SAMRAIVectorReal<DIM,TYPE>::freeVector()
{
   if (d_vector_created_via_duplicate) {
      d_samrai_vector->freeVectorComponents();
      d_samrai_vector.setNull();
      Vec petsc_vec = this -> getPETScVector(); 
      if (petsc_vec) { 
         delete ((PETSc_SAMRAIVectorReal<DIM,TYPE>*)(petsc_vec->data));
      }
   }
}

template<int DIM, class TYPE>
void PETSc_SAMRAIVectorReal<DIM,TYPE>::viewVector() const
{
   ostream& s = d_samrai_vector->getOutputStream();
   s << "\nPrinting PETSc_SAMRAIVectorReal<DIM>..." 
     << "\nSAMRAI vector structure and data: " << endl;
   d_samrai_vector->print(s);
   s << "\n" << endl;
}

template<int DIM, class TYPE>
void PETSc_SAMRAIVectorReal<DIM,TYPE>::getDataArray(TYPE* array)
{
   (void) array;
   TBOX_ERROR("PETSc_SAMRAIVectorReal<DIM,TYPE>::getDataArray undefined!");
}

template<int DIM, class TYPE>
double PETSc_SAMRAIVectorReal<DIM,TYPE>::maxPointwiseDivide(
   const PETScAbstractVectorReal<TYPE>* y)
{
   tbox::Pointer<SAMRAIVectorReal<DIM,TYPE> > y_vector
      = SPVEC_CAST(y)->getSAMRAIVector();
   TYPE max = d_samrai_vector->maxPointwiseDivide(y_vector);
   return max;
}

}
}
#endif
#endif
