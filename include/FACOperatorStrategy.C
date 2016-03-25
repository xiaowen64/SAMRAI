//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/solvers/FAC/FACOperatorStrategy.C $
// Package:	SAMRAI solvers
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
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
  return;
}

template<int DIM> void FACOperatorStrategy<DIM>::initializeOperatorState(
   const SAMRAIVectorReal<DIM,double> &solution ,
   const SAMRAIVectorReal<DIM,double> &rhs )
{
   return;
}

template<int DIM> void FACOperatorStrategy<DIM>::deallocateOperatorState()
{
   return;
}


}
}
#endif
