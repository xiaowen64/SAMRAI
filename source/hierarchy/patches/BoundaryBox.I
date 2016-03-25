//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/patches/BoundaryBox.I $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2050 $
// Modified:	$LastChangedDate: 2008-03-12 13:47:16 -0700 (Wed, 12 Mar 2008) $
// Description:	BoundaryBox representing a portion of the physical boundary 
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif

#include "BoundaryLookupTable.h"

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
typename BoundaryBox<DIM>::BoundaryDirection
BoundaryBox<DIM>::getBoundaryDirection(const int dir) const
{
   BoundaryLookupTable<DIM>* blut =
      BoundaryLookupTable<DIM>::getLookupTable();

   int bdry_dir = blut->getBoundaryDirections(d_bdry_type)[d_location_index](dir);

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(bdry_dir == -1 || bdry_dir == 0 || bdry_dir == 1);  
#endif

   BoundaryDirection retval;

   if (bdry_dir == -1) {
      retval = LOWER;
   } else if (bdry_dir == 0) {
      retval = MIDDLE;
   } else {
      retval = UPPER;
   }

   return (retval);
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
