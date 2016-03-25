//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/memory/FixedArena.C $
// Package:	SAMRAI toolbox for memory management
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Fixed-size arena for efficient memory management
//

#include "tbox/FixedArena.h"

#include "tbox/ArenaManager.h"
#include "tbox/Utilities.h"

namespace SAMRAI {
   namespace tbox {

FixedArena::FixedArena(const size_t bytes)
{
   if (bytes > 0) {
      d_arena = ArenaManager::getManager()->getStandardAllocator();
      d_pool  = d_arena->alloc(bytes);
      d_size  = bytes;
      d_used  = 0;
   } else {
      d_pool = NULL;
      d_size = 0;
      d_used = 0;
   }
}

FixedArena::~FixedArena()
{
   if (d_pool) d_arena->free(d_pool);
}

void *FixedArena::alloc(const size_t bytes)
{
   const size_t allocate = Arena::align(bytes);
   if (allocate + d_used > d_size) {
      TBOX_ERROR("FixedArena::alloc(size) error ...\n"
                 << "Out of memory: size = " << allocate);
   }
   void *ptr = ((char *) d_pool) + d_used;
   d_used += allocate;
   return(ptr);
}

void FixedArena::free(void *p)
{
   NULL_USE(p);
}

Pointer<Arena> FixedArena::allocateArena(const size_t bytes)
{
   return(Pointer<Arena>(new FixedArena(bytes)));
}

void FixedArena::printClassData(std::ostream& os) const
{
   os << "Pointer to fixed memory arena = " << d_pool << std::endl;
   os << "Number of bytes allocated from arena = " << d_used << std::endl;
   os << "Maximum size of fixed memory arena = " << d_size << std::endl;
}

}
}
