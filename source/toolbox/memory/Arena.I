//
// File:	Arena.I
// Package:	SAMRAI toolbox for memory management
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
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
