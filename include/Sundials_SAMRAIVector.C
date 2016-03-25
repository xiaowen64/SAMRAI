//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/solvers/vectors/Sundials_SAMRAIVector.C $
// Package:     SAMRAI solvers
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1917 $
// Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description: "Glue code" between SAMRAI vector object and Sundials vector.
//

#ifndef included_solv_Sundials_SAMRAIVector_C
#define included_solv_Sundials_SAMRAIVector_C

#include "Sundials_SAMRAIVector.h"

#ifdef HAVE_SUNDIALS


#ifndef NULL
#define NULL (0)
#endif

#ifdef DEBUG_NO_INLINE
#include "Sundials_SAMRAIVector.I"
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

template<int DIM> SundialsAbstractVector *Sundials_SAMRAIVector<DIM>::createSundialsVector(
   tbox::Pointer< SAMRAIVectorReal<DIM,double> > samrai_vec)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!samrai_vec.isNull());
#endif
   SundialsAbstractVector* skv = new Sundials_SAMRAIVector<DIM>(samrai_vec);
   
   return( skv );
}

template<int DIM> void Sundials_SAMRAIVector<DIM>::destroySundialsVector(
   SundialsAbstractVector *sundials_vec)
{
   if (sundials_vec) {
      delete ( dynamic_cast<Sundials_SAMRAIVector<DIM>*>(sundials_vec) );
   }
}

template<int DIM> tbox::Pointer< SAMRAIVectorReal<DIM,double> >
Sundials_SAMRAIVector<DIM>::getSAMRAIVector(
   SundialsAbstractVector* sundials_vec)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!(sundials_vec == (SundialsAbstractVector*)NULL));
#endif
   return( ( dynamic_cast<Sundials_SAMRAIVector<DIM>*>(sundials_vec))->getSAMRAIVector() );
}

template<int DIM> tbox::Pointer< SAMRAIVectorReal<DIM,double> >
Sundials_SAMRAIVector<DIM>::getSAMRAIVector(
   N_Vector sundials_vec)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!(sundials_vec == NULL));
#endif
// sgs
   return ( static_cast<Sundials_SAMRAIVector<DIM>*>(sundials_vec -> content) -> getSAMRAIVector() );
}

/*
*************************************************************************
*                                                                       *
* Constructor and destructor for Sundials_SAMRAIVector<DIM>.          *
*                                                                       *
*************************************************************************
*/

template<int DIM>  Sundials_SAMRAIVector<DIM>::Sundials_SAMRAIVector(
   tbox::Pointer< SAMRAIVectorReal<DIM,double> > samrai_vector)
:
   SundialsAbstractVector(),
   d_samrai_vector(samrai_vector)
{
}

template<int DIM>  Sundials_SAMRAIVector<DIM>::~Sundials_SAMRAIVector()
{
}

/*
*************************************************************************
*                                                                       *
* Other miscellaneous member functions					*
*                                                                       *
*************************************************************************
*/

template<int DIM> SundialsAbstractVector* Sundials_SAMRAIVector<DIM>::makeNewVector() 
{
   Sundials_SAMRAIVector<DIM>* out_vec =
      new Sundials_SAMRAIVector<DIM>(d_samrai_vector->cloneVector("out_vec"));
   out_vec->getSAMRAIVector()->allocateVectorData();
   return( out_vec );
}

template<int DIM> void Sundials_SAMRAIVector<DIM>::freeVector()
{
   d_samrai_vector->freeVectorComponents();
   d_samrai_vector.setNull();
   delete this;
}

template<int DIM> void Sundials_SAMRAIVector<DIM>::printVector() const
{
   std::ostream& s = d_samrai_vector->getOutputStream();
   s << "\nPrinting Sundials_SAMRAIVector<DIM>..."
     << "\nthis = " << (Sundials_SAMRAIVector<DIM>*)this
     << "\nSAMRAI vector structure and data: " << std::endl;
   d_samrai_vector->print(s);
   s << "\n" << std::endl;
}

}
}

#endif 
#endif
