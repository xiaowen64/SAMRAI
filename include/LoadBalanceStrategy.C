//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/mesh/load_balance/LoadBalanceStrategy.C $
// Package:     SAMRAI mesh generation
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1917 $
// Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description: Strategy interface for box load balancing routines.
//

#ifndef included_mesh_LoadBalanceStrategy_C
#define included_mesh_LoadBalanceStrategy_C

#include "LoadBalanceStrategy.h"

namespace SAMRAI {
    namespace mesh {

/*
*************************************************************************
*									*
* The constructor and destructor for LoadBalanceStrategy do        *
* nothing that could be considered even remotely interesting.		*
*									*
*************************************************************************
*/

template<int DIM>  LoadBalanceStrategy<DIM>::LoadBalanceStrategy()
{
}

template<int DIM>  LoadBalanceStrategy<DIM>::~LoadBalanceStrategy()
{
}

}
}

#endif
