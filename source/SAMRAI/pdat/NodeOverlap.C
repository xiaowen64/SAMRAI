/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_NodeOverlap_C
#define included_pdat_NodeOverlap_C

#include "SAMRAI/pdat/NodeOverlap.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/pdat/NodeOverlap.I"
#endif

namespace SAMRAI {
namespace pdat {

NodeOverlap::NodeOverlap(
   const hier::BoxContainer& boxes,
   const hier::Transformation& transformation):
   d_is_overlap_empty(boxes.isEmpty()),
   d_transformation(transformation),
   d_dst_boxes(boxes)
{
}

NodeOverlap::~NodeOverlap()
{
}

bool NodeOverlap::isOverlapEmpty() const
{
   return d_is_overlap_empty;
}

}
}
#endif
