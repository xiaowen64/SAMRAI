//
// File:   $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/hierarchy/patches/BasePatchLevel.C $
// Package:   SAMRAI hierarchy
// Copyright:   (c) 1997-2003 Lawrence Livermore National Security, LLC
// Revision:   $LastChangedRevision: 1704 $
// Modified:   $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:   An abstract base class for a level of the AMR hierarchy
//

#include "BasePatchLevel.h"

namespace SAMRAI {
   namespace hier {


/*
*************************************************************************
*                                                                       *
* Constructor and destructor for abstract base class                    *
*                                                                       *
*************************************************************************
*/

template<int DIM> BasePatchLevel<DIM>::BasePatchLevel()
{
}

template<int DIM> BasePatchLevel<DIM>::~BasePatchLevel()
{
}

}
}
