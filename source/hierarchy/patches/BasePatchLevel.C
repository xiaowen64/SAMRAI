//
// File:   $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-0/source/hierarchy/patches/BasePatchLevel.C $
// Package:   SAMRAI hierarchy
// Copyright:   (c) 1997-2003 Lawrence Livermore National Security, LLC
// Revision:   $LastChangedRevision: 2195 $
// Modified:   $LastChangedDate: 2008-05-14 11:33:30 -0700 (Wed, 14 May 2008) $
// Description:   An abstract base class for a level of the AMR hierarchy
//

#ifndef included_hier_BasePatchLevel_C
#define included_hier_BasePatchLevel_C

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

#endif
