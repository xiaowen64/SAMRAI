//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/solvers/nonlinear/NonlinearSolverStrategy.C $
// Package:     SAMRAI solvers
// Copyright:   (c) 1997-2000 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
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
