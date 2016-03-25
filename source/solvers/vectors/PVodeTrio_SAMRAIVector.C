//
// File:        PVodeTrio_SAMRAIVector.C
// Package:     SAMRAI solvers
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: "Glue code" between SAMRAI vector object and PVodeTrio vector.
//

#ifndef included_solv_PVodeTrio_SAMRAIVector_C
#define included_solv_PVodeTrio_SAMRAIVector_C

#include "PVodeTrio_SAMRAIVector.h"

#if HAVE_KINSOL || HAVE_PVODE

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

#ifndef NULL
#define NULL (0)
#endif

#ifdef DEBUG_NO_INLINE
#include "PVodeTrio_SAMRAIVector.I"
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

template<int DIM> PVodeTrioAbstractVector* PVodeTrio_SAMRAIVector<DIM>::createPVodeTrioVector(
   tbox::Pointer< SAMRAIVectorReal<DIM,double> > samrai_vec)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!samrai_vec.isNull());
#endif
   PVodeTrioAbstractVector* skv = new PVodeTrio_SAMRAIVector<DIM>(samrai_vec);
   return( skv );
}

template<int DIM> void PVodeTrio_SAMRAIVector<DIM>::destroyPVodeTrioVector(
   PVodeTrioAbstractVector* pvode_trio_vec)
{
   if (pvode_trio_vec) {
      delete ( ((PVodeTrio_SAMRAIVector<DIM>*)(pvode_trio_vec)) );
   }
}

template<int DIM> tbox::Pointer< SAMRAIVectorReal<DIM,double> >
PVodeTrio_SAMRAIVector<DIM>::getSAMRAIVector(
   PVodeTrioAbstractVector* pvode_trio_vec)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(pvode_trio_vec == (PVodeTrioAbstractVector*)NULL));
#endif
   return( ((PVodeTrio_SAMRAIVector<DIM>*)(pvode_trio_vec))->getSAMRAIVector() );
}

/*
*************************************************************************
*                                                                       *
* Constructor and destructor for PVodeTrio_SAMRAIVector<DIM>.          *
*                                                                       *
*************************************************************************
*/

template<int DIM>  PVodeTrio_SAMRAIVector<DIM>::PVodeTrio_SAMRAIVector(
   tbox::Pointer< SAMRAIVectorReal<DIM,double> > samrai_vector)
:
   PVodeTrioAbstractVector(),
   d_samrai_vector(samrai_vector)
{
}

template<int DIM>  PVodeTrio_SAMRAIVector<DIM>::~PVodeTrio_SAMRAIVector()
{
}

/*
*************************************************************************
*                                                                       *
* Other miscellaneous member functions					*
*                                                                       *
*************************************************************************
*/

template<int DIM> PVodeTrioAbstractVector* PVodeTrio_SAMRAIVector<DIM>::makeNewVector() 
{
   PVodeTrio_SAMRAIVector<DIM>* out_vec =
      new PVodeTrio_SAMRAIVector<DIM>(d_samrai_vector->cloneVector("out_vec"));
   out_vec->getSAMRAIVector()->allocateVectorData();
   return( out_vec );
}

template<int DIM> void PVodeTrio_SAMRAIVector<DIM>::freeVector()
{
   d_samrai_vector->freeVectorComponents();
   d_samrai_vector.setNull();
   delete this;
}

template<int DIM> void PVodeTrio_SAMRAIVector<DIM>::printVector() const
{
   ostream& s = d_samrai_vector->getOutputStream();
   s << "\nPrinting PVodeTrio_SAMRAIVector<DIM>..."
     << "\nthis = " << (PVodeTrio_SAMRAIVector<DIM>*)this
     << "\nSAMRAI vector structure and data: " << endl;
   d_samrai_vector->print(s);
   s << "\n" << endl;
}

}
}

#endif 
#endif
