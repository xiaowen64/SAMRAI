//
// File:	Serializable.C
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
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
