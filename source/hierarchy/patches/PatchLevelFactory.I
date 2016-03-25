//
// File:	PatchLevelFactory.I
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Abstract factory class for creating patch level objects
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace hier {

template<int DIM> inline
PatchLevelFactory<DIM>::PatchLevelFactory()
{
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
