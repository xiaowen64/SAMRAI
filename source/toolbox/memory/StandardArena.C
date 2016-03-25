//
// File:	StandardArena.C
// Package:	SAMRAI toolbox for memory management
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Arena memory manager for standard allocation requests
//

#include "tbox/StandardArena.h"
#include "tbox/Utilities.h"

#ifdef DEBUG_NO_INLINE
#include "tbox/StandardArena.I"
#endif

namespace SAMRAI {
   namespace tbox {


StandardArena::~StandardArena()
{
}

void *StandardArena::alloc(const size_t bytes)
{
   void *p = new char[bytes];
   if (p == ((void *) NULL)) {
      TBOX_ERROR("StandardArena::alloc(size_t) error ...\n"
                 << "Out of memory: size = " << bytes);
   }
   return(p);
}

void StandardArena::free(void *p)
{
   delete [] ((char *) p);
}

}
}
