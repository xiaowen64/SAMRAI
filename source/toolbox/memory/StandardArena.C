//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/memory/StandardArena.C $
// Package:	SAMRAI toolbox for memory management
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
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
