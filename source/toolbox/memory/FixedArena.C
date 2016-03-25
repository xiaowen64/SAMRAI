//
// File:	FixedArena.C
// Package:	SAMRAI toolbox for memory management
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
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

void FixedArena::printClassData(ostream& os) const
{
   os << "Pointer to fixed memory arena = " << d_pool << endl;
   os << "Number of bytes allocated from arena = " << d_used << endl;
   os << "Maximum size of fixed memory arena = " << d_size << endl;
}

}
}
