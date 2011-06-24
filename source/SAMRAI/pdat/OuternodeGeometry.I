/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   hier 
 *
 ************************************************************************/

namespace SAMRAI {
namespace pdat {

SAMRAI_INLINE_KEYWORD
const hier::Box& OuternodeGeometry::getBox() const
{
   return d_box;
}

SAMRAI_INLINE_KEYWORD
const hier::IntVector& OuternodeGeometry::getGhosts() const
{
   return d_ghosts;
}

}
}
