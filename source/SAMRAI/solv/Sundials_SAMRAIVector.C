/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   "Glue code" between SAMRAI vector object and Sundials vector. 
 *
 ************************************************************************/

#ifndef included_solv_Sundials_SAMRAIVector_C
#define included_solv_Sundials_SAMRAIVector_C

#include "SAMRAI/solv/Sundials_SAMRAIVector.h"

#ifdef HAVE_SUNDIALS

#ifndef SAMRAI_INLINE
#include "SAMRAI/solv/Sundials_SAMRAIVector.I"
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

SundialsAbstractVector *Sundials_SAMRAIVector::createSundialsVector(
   tbox::Pointer<SAMRAIVectorReal<double> > samrai_vec)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!samrai_vec.isNull());
#endif
   SundialsAbstractVector* skv = new Sundials_SAMRAIVector(samrai_vec);

   return skv;
}

void Sundials_SAMRAIVector::destroySundialsVector(
   SundialsAbstractVector* sundials_vec)
{
   if (sundials_vec) {
      delete (dynamic_cast<Sundials_SAMRAIVector *>(sundials_vec));
   }
}

tbox::Pointer<SAMRAIVectorReal<double> >
Sundials_SAMRAIVector::getSAMRAIVector(
   SundialsAbstractVector* sundials_vec)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!(sundials_vec == (SundialsAbstractVector *)NULL));
#endif
   return (dynamic_cast<Sundials_SAMRAIVector *>(sundials_vec))->
          getSAMRAIVector();
}

tbox::Pointer<SAMRAIVectorReal<double> >
Sundials_SAMRAIVector::getSAMRAIVector(
   N_Vector sundials_vec)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!(sundials_vec == NULL));
#endif
// sgs
   return static_cast<Sundials_SAMRAIVector *>(sundials_vec->content)->
          getSAMRAIVector();
}

/*
 *************************************************************************
 *                                                                       *
 * Constructor and destructor for Sundials_SAMRAIVector.          *
 *                                                                       *
 *************************************************************************
 */

Sundials_SAMRAIVector::Sundials_SAMRAIVector(
   tbox::Pointer<SAMRAIVectorReal<double> > samrai_vector):
   SundialsAbstractVector(),
   d_samrai_vector(samrai_vector)
{
}

Sundials_SAMRAIVector::~Sundials_SAMRAIVector()
{
}

/*
 *************************************************************************
 *                                                                       *
 * Other miscellaneous member functions					*
 *                                                                       *
 *************************************************************************
 */

SundialsAbstractVector *Sundials_SAMRAIVector::makeNewVector()
{
   Sundials_SAMRAIVector* out_vec =
      new Sundials_SAMRAIVector(d_samrai_vector->cloneVector("out_vec"));
   out_vec->getSAMRAIVector()->allocateVectorData();
   return out_vec;
}

void Sundials_SAMRAIVector::freeVector()
{
   d_samrai_vector->freeVectorComponents();
   d_samrai_vector.setNull();
   delete this;
}

void Sundials_SAMRAIVector::printVector() const
{
   std::ostream& s = d_samrai_vector->getOutputStream();
   s << "\nPrinting Sundials_SAMRAIVector..."
     << "\nthis = " << (Sundials_SAMRAIVector *)this
     << "\nSAMRAI vector structure and data: " << std::endl;
   d_samrai_vector->print(s);
   s << "\n" << std::endl;
}

}
}

#endif
#endif
