//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/toolbox/base/IEEE.I $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	IEEE routines to set up handlers and get signaling NaNs
//

#include "tbox/MathUtilities.h"

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
   namespace tbox {

/*
 * Routines to get signaling NaN value.
 */

inline
float IEEE::getSignalingFloatNaN()
{
   return(MathUtilities<float>::getSignalingNaN());
}

inline
double IEEE::getSignalingNaN()
{
   return(MathUtilities<double>::getSignalingNaN());
}

inline
dcomplex IEEE::getSignalingComplexNaN()
{
   return(MathUtilities<dcomplex>::getSignalingNaN());
}


/*
 * Routines to set signaling NaN value.
 */

inline
void IEEE::setNaN(float& f)
{
   f = MathUtilities<float>::getSignalingNaN();
}

inline
void IEEE::setNaN(double& d)
{
   d = MathUtilities<double>::getSignalingNaN();
}

inline
void IEEE::setNaN(dcomplex& dc)
{
   dc = MathUtilities<dcomplex>::getSignalingNaN();
}

/*
 * Routines to check against signaling NaN value.
 */

inline
bool IEEE::isNaN(const float& f)
{
   return(MathUtilities<float>::isNaN(f));
}

inline
bool IEEE::isNaN(const double& d)
{
   return(MathUtilities<double>::isNaN(d));
}

inline
bool IEEE::isNaN(const dcomplex& dc)
{
   return(MathUtilities<dcomplex>::isNaN(dc));
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}
