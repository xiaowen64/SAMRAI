//
// File:	ScratchArena.C
// Package:	SAMRAI toolbox for memory management
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Arena memory manager for scratch space
//

#include "tbox/ScratchArena.h"
#include "tbox/Utilities.h"

#ifdef DEBUG_NO_INLINE
#include "tbox/ScratchArena.I"
#endif

namespace SAMRAI {
   namespace tbox {

ScratchArena::~ScratchArena()
{
}

void *ScratchArena::alloc(const size_t bytes)
{
   void *p = new char[bytes];
   if (p == ((void *) NULL)) {
      TBOX_ERROR("ScratchArena::alloc(size_t) error ...\n"
                 << "Out of memory: size = " << bytes);
   }
   return(p);
}

void ScratchArena::free(void *p)
{
   delete [] ((char *) p);
}

}
}
