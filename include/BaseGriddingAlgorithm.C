//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/mesh/gridding/BaseGriddingAlgorithm.C $
// Package:     SAMRAI mesh
// Copyright:   (c) 1997-2003 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: AMR hierarchy generation and regridding routines.
//

#ifndef included_mesh_BaseGriddingAlgorithm_C
#define included_mesh_BaseGriddingAlgorithm_C

#include "BaseGriddingAlgorithm.h"

namespace SAMRAI {
   namespace mesh {

/*
*************************************************************************
*                                                                       *
* Constructor and destructor for BaseGriddingAlgorithm<DIM>.            *
*                                                                       *
*************************************************************************
*/
template<int DIM> BaseGriddingAlgorithm<DIM>::BaseGriddingAlgorithm()
{
}

/*
*************************************************************************
*                                                                       *
* Destructor tells the tbox_RestartManager to remove this object from   *
* the list of restart items.                                            *
*                                                                       *
*************************************************************************
*/
template<int DIM> BaseGriddingAlgorithm<DIM>::~BaseGriddingAlgorithm()
{
}

       
}
}
#endif
