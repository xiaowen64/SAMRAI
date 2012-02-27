/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   For describing coarse-fine boundary interfaces
 *
 ************************************************************************/

namespace SAMRAI {
namespace hier {

SAMRAI_INLINE_KEYWORD
void
CoarseFineBoundary::clear()
{
   d_boundary_boxes.clear();
}

SAMRAI_INLINE_KEYWORD
const tbox::Array<BoundaryBox>&
CoarseFineBoundary::getNodeBoundaries(
   const GlobalId& global_id,
   const BlockId& block_id) const
{
   return getBoundaries(global_id, d_dim.getValue(), block_id);
}

SAMRAI_INLINE_KEYWORD
const tbox::Array<BoundaryBox>&
CoarseFineBoundary::getEdgeBoundaries(
   const GlobalId& global_id,
   const BlockId& block_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_dim.getValue() < 2) {
      TBOX_ERROR("CoarseFineBoundary::getEdgeBoundaries():  There are\n"
         << "no edge boundaries in " << d_dim << "d.\n");
   }
#endif
   return getBoundaries(global_id, d_dim.getValue() - 1, block_id);
}

SAMRAI_INLINE_KEYWORD
const tbox::Array<BoundaryBox>&
CoarseFineBoundary::getFaceBoundaries(
   const GlobalId& global_id,
   const BlockId& block_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_dim.getValue() < 3) {
      TBOX_ERROR("CoarseFineBoundary::getFaceBoundaries():  There are\n"
         << "no face boundaries in " << d_dim << "d.\n");
   }
#endif
   return getBoundaries(global_id, d_dim.getValue() - 2, block_id);
}

SAMRAI_INLINE_KEYWORD
CoarseFineBoundary&
CoarseFineBoundary::operator = (
   const CoarseFineBoundary& rhs)
{
   d_initialized = rhs.d_initialized;
   d_boundary_boxes = rhs.d_boundary_boxes;
   return *this;
}

}
}
