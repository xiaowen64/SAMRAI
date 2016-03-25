//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/algorithm/method_of_lines/MethodOfLinesPatchStrategy.C $
// Package:     SAMRAI algorithms
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1917 $
// Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description: Interface to application-specific patch functions to support
//              MethodOfLines integration algorithm
//
 
#ifndef included_algs_MethodOfLinesPatchStrategy_C
#define included_algs_MethodOfLinesPatchStrategy_C

#include "MethodOfLinesPatchStrategy.h"
#include "VariableDatabase.h"

namespace SAMRAI {
    namespace algs {


/*
*************************************************************************
*                                                                       *
* Note: hier::Variable contexts should be consistent with those in            *
*       MethodOfLinesIntegrator<DIM> class.                            *
*                                                                       *
*************************************************************************
*/ 

template<int DIM>  MethodOfLinesPatchStrategy<DIM>::MethodOfLinesPatchStrategy()
{
   d_interior_with_ghosts.setNull();
   d_interior.setNull();
}
 
template<int DIM>  MethodOfLinesPatchStrategy<DIM>::~MethodOfLinesPatchStrategy()
{
}
  
}
}
#endif
