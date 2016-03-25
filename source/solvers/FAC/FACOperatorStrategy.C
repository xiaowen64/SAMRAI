//
// File:	FACOperatorStrategy.C
// Package:	SAMRAI solvers
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
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
