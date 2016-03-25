//
// File:	Serializable.I
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	An abstract base class for objects than can be serialized
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
   namespace tbox {

inline
Serializable::Serializable()
{
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}
