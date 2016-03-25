//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/toolbox/database/Serializable.C $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	An abstract base class for objects than can be serialized
//

#include "tbox/Serializable.h"

#ifdef DEBUG_NO_INLINE
#include "tbox/Serializable.I"
#endif

namespace SAMRAI {
   namespace tbox {


Serializable::~Serializable()
{
}

}
}
