//
// File:        BoxGeneratorStrategy.C
// Package:     SAMRAI mesh generation
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
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
