//
// File:        LoadBalanceStrategy.C
// Package:     SAMRAI mesh generation
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
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
