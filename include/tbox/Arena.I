//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/memory/Arena.I $
// Package:	SAMRAI toolbox for memory management
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Virtual base class for Arena memory management
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
   namespace tbox {


inline
Arena::Arena()
{
}

inline
size_t Arena::align(const size_t bytes)
{
   size_t aligned = bytes + ArenaAllocationAlignment - 1;
   aligned -= aligned % ArenaAllocationAlignment;
   return(aligned);
}

}
}

inline
void *operator new(size_t bytes, SAMRAI::tbox::Arena *arena)
{
   return(arena->alloc(bytes));
}

#ifdef LACKS_NEW_PLACEMENT_OPERATOR
inline
void *operator new(size_t bytes, void *ptr)
{
   (void) bytes;
   return(ptr);
}
#endif

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
