/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Abstract base class for patch data objects 
 *
 ************************************************************************/

namespace SAMRAI {
namespace hier {

SAMRAI_INLINE_KEYWORD
const tbox::Dimension& PatchData::getDim() const
{
   return d_box.getDim();
}

SAMRAI_INLINE_KEYWORD
const Box& PatchData::getBox() const
{
   return d_box;
}

SAMRAI_INLINE_KEYWORD
const Box& PatchData::getGhostBox() const
{
   return d_ghost_box;
}

SAMRAI_INLINE_KEYWORD
void PatchData::setGhostBox(
   const Box& ghost_box)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(d_box, ghost_box);
   TBOX_ASSERT((ghost_box * d_box).isSpatiallyEqual(d_box));

   d_ghost_box = ghost_box;
}

SAMRAI_INLINE_KEYWORD
const IntVector& PatchData::getGhostCellWidth() const
{
   return d_ghosts;
}

SAMRAI_INLINE_KEYWORD
void PatchData::setTime(
   const double timestamp)
{
   d_timestamp = timestamp;
}

SAMRAI_INLINE_KEYWORD
double PatchData::getTime() const
{
   return d_timestamp;
}

}
}
