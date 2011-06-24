/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   hier 
 *
 ************************************************************************/

#ifndef included_pdat_CellOverlap_C
#define included_pdat_CellOverlap_C

#include "SAMRAI/pdat/CellOverlap.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/pdat/CellOverlap.I"
#endif

namespace SAMRAI {
namespace pdat {

CellOverlap::CellOverlap(
   const hier::BoxList& boxes,
   const hier::Transformation& transformation):
   d_is_overlap_empty(boxes.isEmpty()),
   d_transformation(transformation),
   d_dst_boxes(boxes)
{
}

CellOverlap::~CellOverlap()
{
}

bool CellOverlap::isOverlapEmpty() const
{
   return d_is_overlap_empty;
}

void CellOverlap::print(
   std::ostream& os) const
{
   os << "CellOverlap boxes:";
   for (hier::BoxList::Iterator b(d_dst_boxes); b; b++) {
      const hier::Box& box = b();
      os << "  " << box << std::endl;
   }
}

}
}
#endif
