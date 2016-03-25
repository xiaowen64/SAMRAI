//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-0/source/toolbox/templates/special/MathUtilitiesInst.C $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2195 $
// Modified:	$LastChangedDate: 2008-05-14 11:33:30 -0700 (Wed, 14 May 2008) $
// Description:	MathUtilities routines to set up handlers and get signaling NaNs
//

#include "tbox/MathUtilities.h"

#ifdef HAVE_CMATH_ISNAN
#include <cmath>
#include <math.h>
#else
#include <math.h>
#endif

#include <float.h>
#include <limits.h>
#include <stdlib.h>

#include "tbox/Complex.h"

namespace SAMRAI {
   namespace tbox {

/*
 * TODO this should be updated to correct syntax for instantiating.
 */

/*
 * WARNING: This method should not be used.
 *
 * Force library templates to be created.
 *
 */
void DontUseThisFunction_MathInstantiations(void) {

   int    i = 1;
   double d = 1.0;
   float  f = 1.0;

   /*
    * These were needed by gcc 3.x compilers to force
    * instantiations of a pow helper template.
    */
   d = std::pow(d,d);
   d = std::pow(d,i);
   f = std::pow(f,f);
}

}

}

