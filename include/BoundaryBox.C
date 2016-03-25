//
// File:	BoundaryBox.C
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	BoundaryBox representing a portion of the physical boundary
//

#ifndef included_hier_BoundaryBox_C
#define included_hier_BoundaryBox_C

#include "BoundaryBox.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

#include "BoundaryLookupTable.h"

#ifdef DEBUG_NO_INLINE
#include "BoundaryBox.I"
#endif


namespace SAMRAI {
   namespace hier {


template<int DIM>  BoundaryBox<DIM>::BoundaryBox(const Box<DIM>& box,
                                     const int bdry_type,
                                     const int location_index)
:  d_box(box)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   BoundaryLookupTable<DIM>* blut =
      BoundaryLookupTable<DIM>::getLookupTable();


   const tbox::Array<int>& location_index_max = blut->getMaxLocationIndices();

   assert ( (bdry_type >= 1) && (bdry_type <= DIM) );

   assert (location_index >= 0);
   assert (location_index < location_index_max[bdry_type-1]);
#endif

   d_bdry_type = bdry_type;

   d_location_index = location_index;

   d_is_mblk_singularity = false;
}

template<int DIM>  BoundaryBox<DIM>::~BoundaryBox()
{
}

}
}

#endif
