 //
// File:        ImplicitEquationStrategy.C
// Package:     SAMRAI algorithms
// Copyright:   (c) 1997-2000 The Regents of the University of California
// Revision:    $Revision: 47 $
// Modified:    $Date: 2004-12-09 16:08:57 -0800 (Thu, 09 Dec 2004) $
// Description: Interface between implicit integrator and equations to solve.
//
 
#ifndef included_algs_ImplicitEquationStrategy_C
#define included_algs_ImplicitEquationStrategy_C

#include "ImplicitEquationStrategy.h"

namespace SAMRAI {
    namespace algs {

template<int DIM>  ImplicitEquationStrategy<DIM>::ImplicitEquationStrategy()
{
}
 
template<int DIM>  ImplicitEquationStrategy<DIM>::~ImplicitEquationStrategy()
{
}

}
}
#endif
