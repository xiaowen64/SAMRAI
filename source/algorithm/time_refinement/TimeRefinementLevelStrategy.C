//
// File:        TimeRefinementLevelStrategy.C
// Package:     SAMRAI algorithms
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Interface to level routines for time-refinement integrator.
//

#ifndef included_algs_TimeRefinementLevelStrategy_C
#define included_algs_TimeRefinementLevelStrategy_C

#include "TimeRefinementLevelStrategy.h"

namespace SAMRAI {
    namespace algs {

/*
*************************************************************************
*                                                                       *
* Constructor and destructor for TimeRefinementLevelStrategy<DIM>.     *
*                                                                       *
*************************************************************************
*/

template<int DIM>  TimeRefinementLevelStrategy<DIM>::TimeRefinementLevelStrategy() 
{
}

template<int DIM>  TimeRefinementLevelStrategy<DIM>::~TimeRefinementLevelStrategy()
{
}

}
}
#endif
