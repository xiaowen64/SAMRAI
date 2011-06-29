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
const hier::Box& SideGeometry::getBox() const
{
   return d_box;
}

SAMRAI_INLINE_KEYWORD
const hier::IntVector& SideGeometry::getGhosts() const
{
   return d_ghosts;
}

SAMRAI_INLINE_KEYWORD
const hier::IntVector& SideGeometry::getDirectionVector() const
{
   return d_directions;
}

}
}
