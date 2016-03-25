//
// File:        NonlinearSolverStrategy.C
// Package:     SAMRAI solvers
// Copyright:   (c) 1997-2000 The Regents of the University of California
// Revision:    $Revision: 47 $
// Modified:    $Date: 2004-12-09 16:08:57 -0800 (Thu, 09 Dec 2004) $
// Description: Interface between implicit integrator and nonlinear solver.
//
 
#ifndef included_solv_NonlinearSolverStrategy_C
#define included_solv_NonlinearSolverStrategy_C

#include "NonlinearSolverStrategy.h"

namespace SAMRAI {
    namespace solv {

template<int DIM>  NonlinearSolverStrategy<DIM>::NonlinearSolverStrategy()
{
}
 
template<int DIM>  NonlinearSolverStrategy<DIM>::~NonlinearSolverStrategy()
{
}

}
}
#endif
