 //
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/algorithm/implicit/ImplicitEquationStrategy.C $
// Package:     SAMRAI algorithms
// Copyright:   (c) 1997-2000 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
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
