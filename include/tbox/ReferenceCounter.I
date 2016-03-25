//
// File:	ReferenceCounter.I
// Package:	SAMRAI toolbox for memory management
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
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

