/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   "Glue code" between PETSc vector interface and SAMRAI vectors.
 *
 ************************************************************************/

#ifndef included_solv_PETSc_SAMRAIVectorReal_C
#define included_solv_PETSc_SAMRAIVectorReal_C

#include "SAMRAI/solv/PETSc_SAMRAIVectorReal.h"

#ifdef HAVE_PETSC

#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/IOStream.h"
#include "SAMRAI/tbox/PIO.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/solv/PETSc_SAMRAIVectorReal.I"
#endif

#include <cstdlib>

namespace SAMRAI {
namespace solv {

/*
 *************************************************************************
 *
 * Static public member functions.
 *
 *************************************************************************
 */

template<class TYPE>
Vec
PETSc_SAMRAIVectorReal<TYPE>::createPETScVector(
   boost::shared_ptr<SAMRAIVectorReal<TYPE> > samrai_vec,
   MPI_Comm comm)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(samrai_vec);
#endif

   static const bool vector_created_via_duplicate = false;

   PETSc_SAMRAIVectorReal<TYPE>* psv = new PETSc_SAMRAIVectorReal<TYPE>(
         samrai_vec, vector_created_via_duplicate, comm);

   return psv->getPETScVector();
}

template<class TYPE>
void
PETSc_SAMRAIVectorReal<TYPE>::destroyPETScVector(
   Vec petsc_vec)
{
   if (petsc_vec != static_cast<Vec>(NULL)) {
      PETSc_SAMRAIVectorReal<TYPE>* psv =
         static_cast<PETSc_SAMRAIVectorReal<TYPE> *>(petsc_vec->data);

#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(psv != NULL);
#endif

      delete psv;
   }
}

template<class TYPE>
boost::shared_ptr<SAMRAIVectorReal<TYPE> >
PETSc_SAMRAIVectorReal<TYPE>::getSAMRAIVector(
   Vec petsc_vec)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(petsc_vec != static_cast<Vec>(NULL));
#endif

   PETSc_SAMRAIVectorReal<TYPE>* psv =
      static_cast<PETSc_SAMRAIVectorReal<TYPE> *>(petsc_vec->data);

#ifdef DEBUG_CHECK_TBOX_ASSERTIONS
   TBOX_ASSERT(psv != NULL);
#endif

   return psv->d_samrai_vector;
}

/*
 *************************************************************************
 *
 * Protected constructor and destructor for PETSc_SAMRAIVectorReal.
 *
 *************************************************************************
 */

template<class TYPE>
PETSc_SAMRAIVectorReal<TYPE>::PETSc_SAMRAIVectorReal(
   boost::shared_ptr<SAMRAIVectorReal<TYPE> > samrai_vector,
   bool vector_created_via_duplicate,
   MPI_Comm comm):
   PETScAbstractVectorReal<TYPE>(vector_created_via_duplicate, comm),
   d_samrai_vector(samrai_vector),
   d_vector_created_via_duplicate(vector_created_via_duplicate)
{
   // intentionally blank
}

template<class TYPE>
PETSc_SAMRAIVectorReal<TYPE>::~PETSc_SAMRAIVectorReal()
{
   // intentionally blank
}

/*
 *************************************************************************
 *
 * Other member functions
 *
 *************************************************************************
 */

template<class TYPE>
PETScAbstractVectorReal<TYPE> *
PETSc_SAMRAIVectorReal<TYPE>::makeNewVector()
{

   Vec petsc_vec = PETSc_SAMRAIVectorReal<TYPE>::getPETScVector();
   MPI_Comm comm;
   int ierr = PetscObjectGetComm(reinterpret_cast<PetscObject>(petsc_vec),
         &comm);
   PETSC_SAMRAI_ERROR(ierr);

   boost::shared_ptr<SAMRAIVectorReal<TYPE> > sam_vec(
      d_samrai_vector->cloneVector(d_samrai_vector->getName()));
   sam_vec->allocateVectorData();
   const bool vector_created_via_duplicate = true;
   PETSc_SAMRAIVectorReal<TYPE>* out_vec =
      new PETSc_SAMRAIVectorReal<TYPE>(sam_vec,
                                       vector_created_via_duplicate,
                                       comm);
   return out_vec;
}

template<class TYPE>
void PETSc_SAMRAIVectorReal<TYPE>::freeVector()
{

   if (d_vector_created_via_duplicate) {
      d_samrai_vector->freeVectorComponents();
      d_samrai_vector.reset();
      Vec petsc_vec = this->getPETScVector();

#ifdef DEBUG_CHECK_TBOX_ASSERTIONS
      TBOX_ASSERT(petsc_vec != static_cast<Vec>(NULL));
#endif
      delete ((PETSc_SAMRAIVectorReal<TYPE> *)(petsc_vec->data));
   }
}

template<class TYPE>
void PETSc_SAMRAIVectorReal<TYPE>::viewVector() const
{
   std::ostream& s = d_samrai_vector->getOutputStream();
   s << "\nPrinting PETSc_SAMRAIVectorReal..."
     << "\nSAMRAI vector structure and data: " << std::endl;
   d_samrai_vector->print(s);
   s << "\n" << std::endl;
}

}
}
#endif
#endif
