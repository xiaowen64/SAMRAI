/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Simple Cartesian grid geometry for an AMR hierarchy. 
 *
 ************************************************************************/

namespace SAMRAI {
namespace geom {

SAMRAI_INLINE_KEYWORD
const double *CartesianPatchGeometry::getDx() const
{
   return d_dx;
}

SAMRAI_INLINE_KEYWORD
const double *CartesianPatchGeometry::getXLower() const
{
   return d_x_lo;
}

SAMRAI_INLINE_KEYWORD
const double *CartesianPatchGeometry::getXUpper() const
{
   return d_x_up;
}

}
}
