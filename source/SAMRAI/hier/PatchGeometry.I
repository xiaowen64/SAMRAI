/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Base class for geometry management on patches
 *
 ************************************************************************/

#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)

#endif

namespace hier {

SAMRAI_INLINE_KEYWORD
bool& PatchGeometry::TwoDimBool::operator () (
   int dim,
   int side)
{
   TBOX_ASSERT(dim >= 0 && dim < d_dim.getValue());
   TBOX_ASSERT(side == 0 || side == 1);

   return d_data[2 * dim + side];
}

SAMRAI_INLINE_KEYWORD
const bool& PatchGeometry::TwoDimBool::operator () (
   int dim,
   int side) const
{
   TBOX_ASSERT(dim >= 0 && dim < d_dim.getValue());
   TBOX_ASSERT(side == 0 || side == 1);

   return d_data[2 * dim + side];
}

SAMRAI_INLINE_KEYWORD
const tbox::Array<tbox::Array<BoundaryBox> > PatchGeometry::getPatchBoundaries()
const
{
   return d_patch_boundaries.getArrays();
}

SAMRAI_INLINE_KEYWORD
const IntVector& PatchGeometry::getRatio() const
{
   return d_ratio_to_level_zero;
}

SAMRAI_INLINE_KEYWORD
bool PatchGeometry::intersectsPhysicalBoundary() const
{
   return d_has_regular_boundary;
}

SAMRAI_INLINE_KEYWORD
bool PatchGeometry::getTouchesRegularBoundary() const
{
   return d_has_regular_boundary;
}

SAMRAI_INLINE_KEYWORD
bool PatchGeometry::getTouchesPeriodicBoundary() const
{
   return d_has_periodic_boundary;
}

SAMRAI_INLINE_KEYWORD
bool
PatchGeometry::getTouchesRegularBoundary(
   int axis,
   int upperlower) const
{
   TBOX_ASSERT(axis >= 0 && axis < d_dim.getValue());
   TBOX_ASSERT(upperlower == 0 || upperlower == 1);

   return d_touches_regular_bdry(axis, upperlower);
}

SAMRAI_INLINE_KEYWORD
bool
PatchGeometry::getTouchesPeriodicBoundary(
   int axis,
   int upperlower) const
{
   TBOX_ASSERT(axis >= 0 && axis < d_dim.getValue());
   TBOX_ASSERT(upperlower == 0 || upperlower == 1);

   return d_touches_periodic_bdry(axis, upperlower);
}

SAMRAI_INLINE_KEYWORD
const tbox::Array<BoundaryBox>&
PatchGeometry::getNodeBoundaries() const
{
   return d_patch_boundaries[d_dim.getValue() - 1];
}

SAMRAI_INLINE_KEYWORD
const tbox::Array<BoundaryBox>&
PatchGeometry::getEdgeBoundaries() const
{
   if (d_dim.getValue() < 2) {
      TBOX_ERROR("PatchGeometry error in getEdgeBoundary...\n"
         << "DIM < 2 not supported." << std::endl);
   }

   // The "funny" indexing prevents a warning when compiling for
   // DIM < 2.  This code is only reached if DIM >= 2 when
   // executing.
   return d_patch_boundaries[d_dim.getValue() < 2 ? 0 : d_dim.getValue() - 2];
}

SAMRAI_INLINE_KEYWORD
const tbox::Array<BoundaryBox>&
PatchGeometry::getFaceBoundaries() const
{
   if (d_dim.getValue() < 3) {
      TBOX_ERROR("PatchGeometry error in getFaceBoundary...\n"
         << "DIM < 3 not supported." << std::endl);
   }

   // The "funny" indexing prevents a warning when compiling for
   // DIM < 3.  This code is only reached if DIM >= 3 when
   // executing.
   return d_patch_boundaries[d_dim.getValue() < 3 ? 0 : d_dim.getValue() - 3];
}

SAMRAI_INLINE_KEYWORD
const tbox::Array<BoundaryBox>&
PatchGeometry::getCodimensionBoundaries(
   const int codim) const
{
   TBOX_ASSERT((codim > 0) && (codim <= d_dim.getValue()));

   return d_patch_boundaries[codim - 1];
}

SAMRAI_INLINE_KEYWORD
void PatchGeometry::TwoDimBool::setAll(
   bool v)
{
   for (int i = 0; i < 2 * d_dim.getValue(); ++i) {
      d_data[i] = v;
   }
}

SAMRAI_INLINE_KEYWORD
PatchGeometry::TwoDimBool::TwoDimBool():
   d_dim(tbox::Dimension::getInvalidDimension())
{
}

SAMRAI_INLINE_KEYWORD
PatchGeometry::TwoDimBool::TwoDimBool(
   const tbox::Dimension& dim):
   d_dim(dim)
{
   TBOX_DIM_ASSERT_CHECK_DIM(dim);
   setAll(false);
}

SAMRAI_INLINE_KEYWORD
PatchGeometry::TwoDimBool::TwoDimBool(
   const tbox::Dimension& dim,
   bool v):
   d_dim(dim)
{
   for (int i = 0; i < 2 * d_dim.getValue(); ++i) {
      d_data[i] = v;
   }
}

SAMRAI_INLINE_KEYWORD
const tbox::Dimension& PatchGeometry::TwoDimBool::getDim() const {
   return d_dim;
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif
