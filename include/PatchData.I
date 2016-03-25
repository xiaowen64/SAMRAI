//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/hierarchy/variables/PatchData.I $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Abstract base class for patch data objects
//


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
   TBOX_ASSERT(ghost_box * d_box == d_box);
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
