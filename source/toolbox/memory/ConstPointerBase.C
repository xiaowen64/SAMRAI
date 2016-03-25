//
// File:	ConstPointerBase.C
// Package:	SAMRAI toolbox for memory management
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: A smart const pointer base class with RTTI
//

#include "tbox/ConstPointerBase.h"

#ifdef DEBUG_NO_INLINE
#include "tbox/ConstPointerBase.I"
#endif

namespace SAMRAI {
   namespace tbox {


ConstPointerBase::~ConstPointerBase()
{
}

}
}
