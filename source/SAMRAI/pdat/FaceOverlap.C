/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_FaceOverlap_C
#define included_pdat_FaceOverlap_C

#include "SAMRAI/pdat/FaceOverlap.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/pdat/FaceOverlap.I"
#endif

namespace SAMRAI {
namespace pdat {

FaceOverlap::FaceOverlap(
   const tbox::Array<hier::BoxList>& boxes,
   const hier::Transformation& transformation):
   d_is_overlap_empty(true),
   d_transformation(transformation)
{
   const tbox::Dimension& dim = d_transformation.getOffset().getDim();
   d_dst_boxes.resizeArray(boxes.getSize(), hier::BoxList(dim));

   for (int d = 0; d < dim.getValue(); d++) {
      d_dst_boxes[d] = boxes[d];
      if (!d_dst_boxes[d].isEmpty()) d_is_overlap_empty = false;
   }

}

FaceOverlap::~FaceOverlap()
{
}

bool FaceOverlap::isOverlapEmpty() const
{
   return d_is_overlap_empty;
}

}
}
#endif
