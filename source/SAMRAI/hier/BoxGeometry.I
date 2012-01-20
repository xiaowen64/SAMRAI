/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Box geometry description for overlap computations
 *
 ************************************************************************/

#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace hier {

SAMRAI_INLINE_KEYWORD
BoxGeometry::BoxGeometry()
{
}

SAMRAI_INLINE_KEYWORD
boost::shared_ptr<BoxOverlap>
BoxGeometry::calculateOverlap(
   const BoxGeometry& src_geometry,
   const Box& src_mask,
   const Box& fill_box,
   const bool overwrite_interior,
   const Transformation& transformation,
   const BoxContainer& dst_restrict_boxes) const
{
   return this->calculateOverlap(
      *this, src_geometry, src_mask, fill_box,
      overwrite_interior, transformation, true,
      dst_restrict_boxes);
}

}
}
