/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Generic utilities for boundary box calculus.
 *
 ************************************************************************/

namespace SAMRAI {
namespace hier {

SAMRAI_INLINE_KEYWORD
void
BoundaryBoxUtils::setBoundaryBox(
   const BoundaryBox& bbox)
{
   d_bbox = bbox;
   computeOutwardShift();
}

SAMRAI_INLINE_KEYWORD
const BoundaryBox&
BoundaryBoxUtils::getBoundaryBox() const
{
   return d_bbox;
}

SAMRAI_INLINE_KEYWORD
const IntVector&
BoundaryBoxUtils::getOutwardShift() const
{
   return d_outward;
}

SAMRAI_INLINE_KEYWORD
int
BoundaryBoxUtils::normalDir() const
{
   return d_bbox.getLocationIndex() / 2;
}

}
}
