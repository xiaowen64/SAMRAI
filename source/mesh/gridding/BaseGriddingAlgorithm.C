//
// File:        BaseGriddingAlgorithm.C
// Package:     SAMRAI mesh
// Copyright:   (c) 1997-2003 The Regents of the University of California
// Revision:    $Revision: 47 $
// Modified:    $Date: 2004-12-09 16:08:57 -0800 (Thu, 09 Dec 2004) $
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
