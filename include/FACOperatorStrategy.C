//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/solvers/FAC/FACOperatorStrategy.C $
// Package:	SAMRAI solvers
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	Interface to user-defined operations used in FAC solve.
//

#ifndef included_solv_FACOperatorStrategy_C
#define included_solv_FACOperatorStrategy_C

#include "FACOperatorStrategy.h"

namespace SAMRAI {
    namespace solv {

template<int DIM>  FACOperatorStrategy<DIM>::FACOperatorStrategy() {
}

template<int DIM>  FACOperatorStrategy<DIM>::~FACOperatorStrategy() {
}

template<int DIM> void FACOperatorStrategy<DIM>::postprocessOneCycle(
   int fac_cycle_num ,
   const SAMRAIVectorReal<DIM,double> &current_soln ,
   const SAMRAIVectorReal<DIM,double> &residual )
{
  NULL_USE(fac_cycle_num);
  NULL_USE(current_soln);
  NULL_USE(residual);
  return;
}

template<int DIM> void FACOperatorStrategy<DIM>::initializeOperatorState(
   const SAMRAIVectorReal<DIM,double> &solution ,
   const SAMRAIVectorReal<DIM,double> &rhs )
{
   NULL_USE(solution);
   NULL_USE(rhs);
   return;
}

template<int DIM> void FACOperatorStrategy<DIM>::deallocateOperatorState()
{
   return;
}


}
}
#endif
