//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/hierarchy/patches/BoundaryBox.I $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	BoundaryBox representing a portion of the physical boundary 
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
   namespace hier {


template<int DIM> inline
BoundaryBox<DIM>::BoundaryBox() : d_bdry_type(-1), d_location_index(-1)
{
}

template<int DIM> inline
BoundaryBox<DIM>::BoundaryBox(const BoundaryBox<DIM>& boundary_box) 
: d_box(boundary_box.d_box),
  d_bdry_type(boundary_box.d_bdry_type),
  d_location_index(boundary_box.d_location_index),
  d_is_mblk_singularity(boundary_box.d_is_mblk_singularity)
{
}

template<int DIM> inline
const Box<DIM>& BoundaryBox<DIM>::getBox() const
{
   return (d_box);
}

template<int DIM> inline
int BoundaryBox<DIM>::getBoundaryType() const
{
   return (d_bdry_type);
}

template<int DIM> inline
int BoundaryBox<DIM>::getLocationIndex() const
{
   return (d_location_index);
}

template<int DIM> inline
void BoundaryBox<DIM>::setIsMultiblockSingularity(bool is_mblk_singularity)
{
   d_is_mblk_singularity = is_mblk_singularity;
}
                                                                                
template<int DIM> inline
bool BoundaryBox<DIM>::getIsMultiblockSingularity() const
{
   return (d_is_mblk_singularity);
}

template<int DIM> inline
BoundaryBox<DIM>& BoundaryBox<DIM>::operator=(
   const BoundaryBox<DIM>& boundary_box)
{
   d_box = boundary_box.d_box;
   d_bdry_type = boundary_box.d_bdry_type;
   d_location_index = boundary_box.d_location_index;
   d_is_mblk_singularity = boundary_box.d_is_mblk_singularity;

   return(*this);
}


#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}
