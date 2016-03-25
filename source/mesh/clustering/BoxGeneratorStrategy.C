//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/mesh/clustering/BoxGeneratorStrategy.C $
// Package:     SAMRAI mesh generation
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1917 $
// Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description: Strategy interface for box generation routines.
//

#ifndef included_mesh_BoxGeneratorStrategy_C
#define included_mesh_BoxGeneratorStrategy_C

#include "BoxGeneratorStrategy.h"

namespace SAMRAI {
    namespace mesh {


/*
*************************************************************************
*                                                                       *
* Default constructor and destructor for BoxGeneratorStrategy<DIM>.    *
*                                                                       *
*************************************************************************
*/

template<int DIM>  BoxGeneratorStrategy<DIM>::BoxGeneratorStrategy()
{
}

template<int DIM>  BoxGeneratorStrategy<DIM>::~BoxGeneratorStrategy()
{
}

}
}

#endif
