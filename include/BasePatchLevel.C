//
// File:   BasePatchLevel.C
// Package:   SAMRAI hierarchy
// Copyright:   (c) 1997-2003 The Regents of the University of California
// Revision:   $Revision: 7 $
// Modified:   $Date: 2004-11-30 13:18:17 -0800 (Tue, 30 Nov 2004) $
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
