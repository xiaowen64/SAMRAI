//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/algorithm/method_of_lines/MethodOfLinesPatchStrategy.C $
// Package:     SAMRAI algorithms
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
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
