//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/database/Serializable.C $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
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
