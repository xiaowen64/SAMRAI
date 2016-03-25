//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/mesh/clustering/BoxGeneratorStrategy.C $
// Package:     SAMRAI mesh generation
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
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
