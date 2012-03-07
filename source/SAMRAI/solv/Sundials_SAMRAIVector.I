/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   "Glue code" between SAMRAI vector object and Sundials vectors.
 *
 ************************************************************************/

#define SKVEC_CAST(x) (((Sundials_SAMRAIVector *)x))

namespace SAMRAI {
namespace solv {

SAMRAI_INLINE_KEYWORD
boost::shared_ptr<SAMRAIVectorReal<double> >
Sundials_SAMRAIVector::getSAMRAIVector()
{
   return d_samrai_vector;
}

}
}
