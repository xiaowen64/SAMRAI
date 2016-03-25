//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/memory/PointerBase.C $
// Package:	SAMRAI toolbox for memory management
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
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
