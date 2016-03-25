#ifndef included_solv_RobinBcCoefStrategy_C
#define included_solv_RobinBcCoefStrategy_C

/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/solvers/poisson/RobinBcCoefStrategy.C $
 * Package:     SAMRAI solvers
 * Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
 * Revision:	$LastChangedRevision: 1917 $
 * Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
 * Description:	Operator class for solving scalar Poisson using FAC
 */

#include "RobinBcCoefStrategy.h"


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
