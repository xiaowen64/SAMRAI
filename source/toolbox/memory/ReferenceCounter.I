//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/toolbox/memory/ReferenceCounter.I $
// Package:	SAMRAI toolbox for memory management
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	Reference counting class for Array and Pointer
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
   namespace tbox {


#ifndef NULL
#define NULL (0)
#endif

inline
ReferenceCounter::ReferenceCounter()
{
   d_references = 1;
   d_arena      = (Arena *) NULL;
   d_counter    = (ReferenceCounter *) NULL;
}

inline
Arena *ReferenceCounter::getArena()
{
   return(d_arena);
}

inline
bool ReferenceCounter::deleteReference()
{
   return(--d_references == 0);
}

inline
void ReferenceCounter::addReference()
{
   d_references++;
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}

