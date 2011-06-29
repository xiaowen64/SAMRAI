/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   BoundaryBox representing a portion of the physical boundary 
 *
 ************************************************************************/

#include "SAMRAI/hier/BoundaryLookupTable.h"

namespace SAMRAI {
namespace hier {

SAMRAI_INLINE_KEYWORD
const tbox::Dimension& BoundaryBox::getDim() const
{
   return d_dim;
}

SAMRAI_INLINE_KEYWORD
BoundaryBox::BoundaryBox(
   const tbox::Dimension& dim):
   d_dim(dim),
   d_box(dim),
   d_bdry_type(-1),
   d_location_index(-1)
{
   TBOX_DIM_ASSERT_CHECK_DIM(d_dim);
}

SAMRAI_INLINE_KEYWORD
BoundaryBox::BoundaryBox(
   const BoundaryBox& boundary_box):
   d_dim(boundary_box.getDim()),
   d_box(boundary_box.d_box),
   d_bdry_type(boundary_box.d_bdry_type),
   d_location_index(boundary_box.d_location_index),
   d_is_mblk_singularity(boundary_box.d_is_mblk_singularity)
{
}

SAMRAI_INLINE_KEYWORD
const Box& BoundaryBox::getBox() const
{
   return d_box;
}

SAMRAI_INLINE_KEYWORD
int BoundaryBox::getBoundaryType() const
{
   return d_bdry_type;
}

SAMRAI_INLINE_KEYWORD
int BoundaryBox::getLocationIndex() const
{
   return d_location_index;
}

SAMRAI_INLINE_KEYWORD
void BoundaryBox::setIsMultiblockSingularity(
   bool is_mblk_singularity)
{
   d_is_mblk_singularity = is_mblk_singularity;
}

SAMRAI_INLINE_KEYWORD
bool BoundaryBox::getIsMultiblockSingularity() const
{
   return d_is_mblk_singularity;
}

SAMRAI_INLINE_KEYWORD
BoundaryBox::BoundaryOrientation BoundaryBox::getBoundaryOrientation(
   const int dir) const
{
   TBOX_ASSERT(dir < d_dim.getValue());

   BoundaryLookupTable* blut =
      BoundaryLookupTable::getLookupTable(d_dim);

   int bdry_dir =
      blut->getBoundaryDirections(d_bdry_type)[d_location_index](dir);

   TBOX_ASSERT(bdry_dir == -1 || bdry_dir == 0 || bdry_dir == 1);

   BoundaryOrientation retval;

   if (bdry_dir == -1) {
      retval = LOWER;
   } else if (bdry_dir == 0) {
      retval = MIDDLE;
   } else {
      retval = UPPER;
   }

   return retval;
}

SAMRAI_INLINE_KEYWORD
BoundaryBox& BoundaryBox::operator = (
   const BoundaryBox& boundary_box)
{
   d_box = boundary_box.d_box;
   d_bdry_type = boundary_box.d_bdry_type;
   d_location_index = boundary_box.d_location_index;
   d_is_mblk_singularity = boundary_box.d_is_mblk_singularity;

   return *this;
}

}
}
