//
// File:	ScratchArena.I
// Package:	SAMRAI toolbox for memory management
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Arena memory manager for scratch space
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
   namespace tbox {


inline
ScratchArena::ScratchArena()
{
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}
