//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/solvers/vectors/PETSc_SAMRAIVectorReal.C $
// Package:     SAMRAI solvers
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1917 $
// Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description: "Glue code" between PETSc vector interface and SAMRAI vectors.
//

#ifndef included_solv_PETSc_SAMRAIVectorReal_C
#define included_solv_PETSc_SAMRAIVectorReal_C

#include "PETSc_SAMRAIVectorReal.h"

#ifdef HAVE_PETSC

#include <stdlib.h>


#include "tbox/Utilities.h"
#include "tbox/IOStream.h"


#include "tbox/PIO.h"

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

template <int DIM, class TYPE>
Vec
PETSc_SAMRAIVectorReal<DIM,TYPE>::createPETScVector(
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<DIM,TYPE> > samrai_vec,
    MPI_Comm comm)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!samrai_vec.isNull());
#endif

    static const bool vector_created_via_duplicate = false;

    PETSc_SAMRAIVectorReal<DIM,TYPE>* psv = new PETSc_SAMRAIVectorReal<DIM,TYPE>(
        samrai_vec, vector_created_via_duplicate, comm);

    return psv->getPETScVector();
}

template <int DIM, class TYPE>
void
PETSc_SAMRAIVectorReal<DIM,TYPE>::destroyPETScVector(
    Vec petsc_vec)
{
    if (petsc_vec != static_cast<Vec>(NULL))
    {
        PETSc_SAMRAIVectorReal<DIM,TYPE>* psv = static_cast<PETSc_SAMRAIVectorReal<DIM,TYPE>*>(petsc_vec->data);

#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(psv != NULL);
#endif

        delete psv;
    }
    return;
}

template <int DIM, class TYPE>
SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<DIM,TYPE> >
PETSc_SAMRAIVectorReal<DIM,TYPE>::getSAMRAIVector(
    Vec petsc_vec)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(petsc_vec != static_cast<Vec>(NULL));
#endif

    PETSc_SAMRAIVectorReal<DIM,TYPE>* psv = static_cast<PETSc_SAMRAIVectorReal<DIM,TYPE>*>(petsc_vec->data);

#ifdef DEBUG_CHECK_TBOX_ASSERTIONS
    TBOX_ASSERT(psv != NULL);
#endif

    return psv->d_samrai_vector;
}

/*
*************************************************************************
*                                                                       *
* Protected constructor and destructor for PETSc_SAMRAIVectorReal<DIM>.*
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
PETSc_SAMRAIVectorReal<DIM,TYPE>::PETSc_SAMRAIVectorReal(
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<DIM,TYPE> > samrai_vector,
    bool vector_created_via_duplicate,
    MPI_Comm comm)
    : PETScAbstractVectorReal<TYPE>(vector_created_via_duplicate,comm),
      d_samrai_vector(samrai_vector),
      d_vector_created_via_duplicate(vector_created_via_duplicate)
{
    // intentionally blank
   return;
}

template <int DIM, class TYPE>
PETSc_SAMRAIVectorReal<DIM,TYPE>::~PETSc_SAMRAIVectorReal()
{
   // intentionally blank
   return;
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

   Vec petsc_vec = PETSc_SAMRAIVectorReal<DIM,TYPE>::getPETScVector();
   MPI_Comm comm;
   int ierr = PetscObjectGetComm(reinterpret_cast<PetscObject>(petsc_vec),&comm); PETSC_SAMRAI_ERROR(ierr);

   tbox::Pointer< SAMRAIVectorReal<DIM,TYPE> > sam_vec =
      d_samrai_vector->cloneVector(d_samrai_vector->getName());
   sam_vec->allocateVectorData();
   const bool vector_created_via_duplicate = true;
   PETSc_SAMRAIVectorReal<DIM,TYPE>* out_vec =
      new PETSc_SAMRAIVectorReal<DIM,TYPE>(sam_vec,
					   vector_created_via_duplicate, 
					   comm);
   return( out_vec );
}

template<int DIM, class TYPE>
void PETSc_SAMRAIVectorReal<DIM,TYPE>::freeVector()
{

   if (d_vector_created_via_duplicate) {
      d_samrai_vector->freeVectorComponents();
      d_samrai_vector.setNull();
      Vec petsc_vec = this -> getPETScVector(); 

#ifdef DEBUG_CHECK_TBOX_ASSERTIONS
      TBOX_ASSERT(petsc_vec != static_cast<Vec>(NULL));
#endif
      delete ((PETSc_SAMRAIVectorReal<DIM,TYPE>*)(petsc_vec->data));
   }
}

template<int DIM, class TYPE>
void PETSc_SAMRAIVectorReal<DIM,TYPE>::viewVector() const
{
   std::ostream& s = d_samrai_vector->getOutputStream();
   s << "\nPrinting PETSc_SAMRAIVectorReal<DIM>..." 
     << "\nSAMRAI vector structure and data: " << std::endl;
   d_samrai_vector->print(s);
   s << "\n" << std::endl;
}

}
}
#endif
#endif
