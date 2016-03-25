//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/toolbox/memory/PointerBase.C $
// Package:	SAMRAI toolbox for memory management
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
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
