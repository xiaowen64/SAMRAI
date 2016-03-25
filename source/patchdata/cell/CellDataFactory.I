//
// File:	CellDataFactory.I
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Factory class for creating cell data objects
//

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace pdat {

template<int DIM, class TYPE> inline
int CellDataFactory<DIM,TYPE>::getDefaultDepth() const
{
   return(d_depth);
}

template<int DIM, class TYPE> inline
void CellDataFactory<DIM,TYPE>::setDefaultDepth(const int depth)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(depth > 0);
#endif
   d_depth = depth;
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
