#ifndef included_solv_RobinBcCoefStrategy_C
#define included_solv_RobinBcCoefStrategy_C

/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/solvers/poisson/RobinBcCoefStrategy.C $
 * Package:     SAMRAI solvers
 * Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
 * Revision:	$LastChangedRevision: 1704 $
 * Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
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
