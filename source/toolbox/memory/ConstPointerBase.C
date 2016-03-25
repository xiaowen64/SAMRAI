//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/toolbox/memory/ConstPointerBase.C $
// Package:	SAMRAI toolbox for memory management
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
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
