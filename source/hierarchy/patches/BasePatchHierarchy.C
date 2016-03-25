//
// File:	BasePatchHierarchy.C
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2003 The Regents of the University of California
// Revision:	$Revision: 47 $
// Modified:	$Date: 2004-12-09 16:08:57 -0800 (Thu, 09 Dec 2004) $
// Description:	An abstract base class for hierarchies
//

#ifndef included_hier_BasePatchHierarchy_C
#define included_hier_BasePatchHierarchy_C

#include "BasePatchHierarchy.h"


namespace SAMRAI {
   namespace hier {

/*
*************************************************************************
*									*
* Constructor and destructor for abstract base class			*
*									*
*************************************************************************
*/

template<int DIM> BasePatchHierarchy<DIM>::BasePatchHierarchy()
{
}

template<int DIM> BasePatchHierarchy<DIM>::~BasePatchHierarchy()
{
}

}
}

#endif
