/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Extension of a std
 *
 ************************************************************************/
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace hier {

SAMRAI_INLINE_KEYWORD
NeighborhoodSet::NeighborhoodSet(const tbox::Dimension& dim):
   d_dim(dim),
   d_map()
{
}

SAMRAI_INLINE_KEYWORD
NeighborhoodSet::NeighborhoodSet(
   const NeighborhoodSet& r):
   d_dim(r.d_dim),
   d_map(r.d_map)
{
}

SAMRAI_INLINE_KEYWORD
NeighborhoodSet::~NeighborhoodSet()
{
}

SAMRAI_INLINE_KEYWORD
NeighborhoodSet& NeighborhoodSet::operator = (
   const NeighborhoodSet& r)
{
   d_map = r.d_map;
   return *this;
}

SAMRAI_INLINE_KEYWORD
bool NeighborhoodSet::operator == (
   const NeighborhoodSet& r) const
{
   return (d_map == r.d_map && d_dim == r.d_dim);
}

SAMRAI_INLINE_KEYWORD
bool NeighborhoodSet::operator != (
   const NeighborhoodSet& r) const
{
   return (d_map != r.d_map || d_dim != r.d_dim);
}

SAMRAI_INLINE_KEYWORD
void NeighborhoodSet::swap(
   NeighborhoodSet& r)
{
   d_map.swap(r.d_map);
}

}
}
