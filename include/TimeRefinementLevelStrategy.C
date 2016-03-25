//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/algorithm/time_refinement/TimeRefinementLevelStrategy.C $
// Package:     SAMRAI algorithms
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
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
