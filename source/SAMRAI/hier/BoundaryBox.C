/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   BoundaryBox representing a portion of the physical boundary
 *
 ************************************************************************/

#ifndef included_hier_BoundaryBox_C
#define included_hier_BoundaryBox_C

#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/BoundaryLookupTable.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/hier/BoundaryBox.I"
#endif

namespace SAMRAI {
namespace hier {

BoundaryBox::BoundaryBox(
   const Box& box,
   const int bdry_type,
   const int location_index):
   d_dim(box.getDim()),
   d_box(box)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   BoundaryLookupTable* blut =
      BoundaryLookupTable::getLookupTable(d_dim);
   const tbox::Array<int>& location_index_max = blut->getMaxLocationIndices();

   TBOX_ASSERT((bdry_type >= 1) && (bdry_type <= d_dim.getValue()));
   TBOX_ASSERT(location_index >= 0);
   TBOX_ASSERT(location_index < location_index_max[bdry_type - 1]);
#endif

   d_bdry_type = bdry_type;

   d_location_index = location_index;

   d_is_mblk_singularity = false;
}

BoundaryBox::~BoundaryBox()
{
}

}
}

#endif
