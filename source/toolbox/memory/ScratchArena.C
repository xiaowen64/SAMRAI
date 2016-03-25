//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/toolbox/memory/ScratchArena.C $
// Package:	SAMRAI toolbox for memory management
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
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
