/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

namespace SAMRAI {
namespace pdat {

SAMRAI_INLINE_KEYWORD
const hier::Box& NodeGeometry::getBox() const
{
   return d_box;
}

SAMRAI_INLINE_KEYWORD
const hier::IntVector& NodeGeometry::getGhosts() const
{
   return d_ghosts;
}

SAMRAI_INLINE_KEYWORD
hier::Box NodeGeometry::toNodeBox(
   const hier::Box& box)
{
   return box.empty() ? box : hier::Box(box.lower(), box.upper() + 1);
}

}
}
