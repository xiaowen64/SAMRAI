//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/algorithm/time_refinement/TimeRefinementLevelStrategy.C $
// Package:     SAMRAI algorithms
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1917 $
// Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
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
