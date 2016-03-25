#ifndef included_solv_RobinBcCoefStrategy_C
#define included_solv_RobinBcCoefStrategy_C

/*
 * File:        RobinBcCoefStrategy.C
 * Package:     SAMRAI solvers
 * Copyright:	(c) 1997-2005 The Regents of the University of California
 * Revision:	$Revision: 453 $
 * Modified:	$Date: 2005-06-16 10:19:28 -0700 (Thu, 16 Jun 2005) $
 * Description:	Operator class for solving scalar Poisson using FAC
 */

#include "RobinBcCoefStrategy.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

namespace SAMRAI {
    namespace solv {

/*
********************************************************************
* Default constructor.                                             *
********************************************************************
*/

template<int DIM>  RobinBcCoefStrategy<DIM>::RobinBcCoefStrategy() {
   return;
}

/*
********************************************************************
*   Destructor.                                                    *
********************************************************************
*/

template<int DIM>  RobinBcCoefStrategy<DIM>::~RobinBcCoefStrategy() {
   return;
}



}
}
#endif
