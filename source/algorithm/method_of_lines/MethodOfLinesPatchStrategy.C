//
// File:        MethodOfLinesPatchStrategy.C
// Package:     SAMRAI algorithms
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
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
