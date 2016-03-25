//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/mesh/load_balance/LoadBalanceStrategy.C $
// Package:     SAMRAI mesh generation
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
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
