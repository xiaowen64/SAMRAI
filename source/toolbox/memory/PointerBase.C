//
// File:	PointerBase.C
// Package:	SAMRAI toolbox for memory management
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: A smart pointer base class with RTTI
//

#include "tbox/PointerBase.h"

#ifdef DEBUG_NO_INLINE
#include "tbox/PointerBase.I"
#endif

namespace SAMRAI {
   namespace tbox {


PointerBase::~PointerBase()
{
}

}
}
