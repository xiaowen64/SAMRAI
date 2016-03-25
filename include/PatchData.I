//
// File:	PatchData.I
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 301 $
// Modified:	$Date: 2005-04-25 10:32:08 -0700 (Mon, 25 Apr 2005) $
// Description:	Abstract base class for patch data objects
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
    namespace hier {

template<int DIM> inline
const Box<DIM>& PatchData<DIM>::getBox() const
{
   return(d_box);
}

template<int DIM> inline
const Box<DIM>& PatchData<DIM>::getGhostBox() const
{
   return(d_ghost_box);
}

template<int DIM> inline
void PatchData<DIM>::setGhostBox(const Box<DIM>& ghost_box)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(ghost_box * d_box == d_box);
#endif
   d_ghost_box = ghost_box;
}

template<int DIM> inline
const IntVector<DIM>& PatchData<DIM>::getGhostCellWidth() const
{
   return(d_ghosts);
}

template<int DIM> inline
void PatchData<DIM>::setTime(const double timestamp)
{
   d_timestamp = timestamp;
}

template<int DIM> inline
double PatchData<DIM>::getTime() const
{
   return(d_timestamp);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
