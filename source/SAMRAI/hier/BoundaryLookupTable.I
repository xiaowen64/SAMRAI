/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:    
 *
 ************************************************************************/

#include "SAMRAI/tbox/StartupShutdownManager.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace hier {

/*
 *************************************************************************
 *                                                                       *
 * Static lookup table member functions.                                 *
 *                                                                       *
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
BoundaryLookupTable *BoundaryLookupTable::getLookupTable(
   const tbox::Dimension& dim)
{
   if (!s_lookup_table_instance[dim.getValue() - 1]) {
      s_lookup_table_instance[dim.getValue() - 1] = new BoundaryLookupTable(dim);
   }
   return s_lookup_table_instance[dim.getValue() - 1];
}

/*
 *************************************************************************
 *                                                                       *
 * Inline accessory functions to retrieve data members.                  *
 *                                                                       *
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
const tbox::Array<int>& BoundaryLookupTable::getMaxLocationIndices() const
{
   return d_max_li;
}

SAMRAI_INLINE_KEYWORD
const tbox::Array<int>& BoundaryLookupTable::getDirections(
   const int loc,
   const int codim) const
{
   TBOX_ASSERT((codim > 0) && (codim <= d_dim.getValue()));
   TBOX_ASSERT((loc >= 0) && (loc < d_max_li[codim - 1]));

   int iloc = loc / (1 << codim);
   return d_table[codim - 1][iloc];
}

SAMRAI_INLINE_KEYWORD
bool BoundaryLookupTable::isUpper(
   const int loc,
   const int codim,
   const int index) const
{
   TBOX_ASSERT((codim > 0) && (codim <= d_dim.getValue()));
   TBOX_ASSERT((loc >= 0) && (loc < d_max_li[codim - 1]));
   TBOX_ASSERT((index >= 0) && (index < codim));

   return (loc % (1 << codim)) & (1 << (index));
}

SAMRAI_INLINE_KEYWORD
bool BoundaryLookupTable::isLower(
   const int loc,
   const int codim,
   const int index) const
{
   TBOX_ASSERT((codim > 0) && (codim <= d_dim.getValue()));
   TBOX_ASSERT((loc >= 0) && (loc < d_max_li[codim - 1]));
   TBOX_ASSERT((index >= 0) && (index < codim));

   return !isUpper(loc, codim, index);
}

/*
 *************************************************************************
 *                                                                       *
 * Get array of IntVectors indication location of a boundary relative    *
 * to a patch.                                                           *
 *                                                                       *
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
const tbox::Array<IntVector>&
BoundaryLookupTable::getBoundaryDirections(
   int codim) const
{
   TBOX_ASSERT((codim > 0) && (codim <= d_dim.getValue()));
   return d_bdry_dirs[codim - 1];
}

}
}
