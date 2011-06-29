/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   hier 
 *
 ************************************************************************/

#ifndef included_pdat_EdgeOverlap_C
#define included_pdat_EdgeOverlap_C

#include "SAMRAI/pdat/EdgeOverlap.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/pdat/EdgeOverlap.I"
#endif

namespace SAMRAI {
namespace pdat {

EdgeOverlap::EdgeOverlap(
   const tbox::Array<hier::BoxList>& boxes,
   const hier::Transformation& transformation):
   d_is_overlap_empty(true),
   d_transformation(transformation)
{

   d_dst_boxes.resizeArray(boxes.getSize());

   for (int d = 0; d < boxes.getSize(); d++) {
      d_dst_boxes[d] = boxes[d];
      if (!d_dst_boxes[d].isEmpty()) d_is_overlap_empty = false;
   }
}

EdgeOverlap::~EdgeOverlap()
{
}

bool EdgeOverlap::isOverlapEmpty() const
{
   return d_is_overlap_empty;
}

}
}
#endif
