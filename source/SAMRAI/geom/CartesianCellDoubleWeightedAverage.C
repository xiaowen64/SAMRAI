/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Weighted averaging operator for cell-centered double data on
 *                a Cartesian mesh. 
 *
 ************************************************************************/

#ifndef included_geom_CartesianCellDoubleWeightedAverage_C
#define included_geom_CartesianCellDoubleWeightedAverage_C

#include "SAMRAI/geom/CartesianCellDoubleWeightedAverage.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/Utilities.h"

#include <float.h>
#include <math.h>

/*
 *************************************************************************
 *                                                                       *
 * External declarations for FORTRAN  routines.                          *
 *                                                                       *
 *************************************************************************
 */
extern "C" {

#ifdef __INTEL_COMPILER
#pragma warning (disable:1419)
#endif

// in cartcoarsen1d.f:
void F77_FUNC(cartwgtavgcelldoub1d, CARTWGTAVGCELLDOUB1D) (const int &,
   const int &,
   const int &, const int &,
   const int &, const int &,
   const int *, const double *, const double *,
   const double *, double *);
// in cartcoarsen2d.f:
void F77_FUNC(cartwgtavgcelldoub2d, CARTWGTAVGCELLDOUB2D) (const int &,
   const int &, const int &, const int &,
   const int &, const int &, const int &, const int &,
   const int &, const int &, const int &, const int &,
   const int *, const double *, const double *,
   const double *, double *);
// in cartcoarsen3d.f:
void F77_FUNC(cartwgtavgcelldoub3d, CARTWGTAVGCELLDOUB3D) (const int &,
   const int &, const int &,
   const int &, const int &, const int &,
   const int &, const int &, const int &,
   const int &, const int &, const int &,
   const int &, const int &, const int &,
   const int &, const int &, const int &,
   const int *, const double *, const double *,
   const double *, double *);
}

namespace SAMRAI {
namespace geom {

CartesianCellDoubleWeightedAverage::CartesianCellDoubleWeightedAverage(
   const tbox::Dimension& dim):
   hier::CoarsenOperator(dim, "CONSERVATIVE_COARSEN")
{
}

CartesianCellDoubleWeightedAverage::~CartesianCellDoubleWeightedAverage()
{
}

bool CartesianCellDoubleWeightedAverage::findCoarsenOperator(
   const tbox::Pointer<hier::Variable>& var,
   const std::string& op_name) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *var);

   const tbox::Pointer<pdat::CellVariable<double> > cast_var(var);
   if (!cast_var.isNull() && (op_name == getOperatorName())) {
      return true;
   } else {
      return false;
   }
}

int CartesianCellDoubleWeightedAverage::getOperatorPriority() const
{
   return 0;
}

hier::IntVector
CartesianCellDoubleWeightedAverage::getStencilWidth() const {
   return hier::IntVector::getZero(getDim());
}

void CartesianCellDoubleWeightedAverage::coarsen(
   hier::Patch& coarse,
   const hier::Patch& fine,
   const int dst_component,
   const int src_component,
   const hier::Box& coarse_box,
   const hier::IntVector& ratio) const
{
   const tbox::Dimension& dim(getDim());

   TBOX_DIM_ASSERT_CHECK_DIM_ARGS4(dim, coarse, fine, coarse_box, ratio);

   tbox::Pointer<pdat::CellData<double> >
   fdata = fine.getPatchData(src_component);
   tbox::Pointer<pdat::CellData<double> >
   cdata = coarse.getPatchData(dst_component);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!fdata.isNull());
   TBOX_ASSERT(!cdata.isNull());
   TBOX_ASSERT(cdata->getDepth() == fdata->getDepth());
#endif

   const hier::Index filo = fdata->getGhostBox().lower();
   const hier::Index fihi = fdata->getGhostBox().upper();
   const hier::Index cilo = cdata->getGhostBox().lower();
   const hier::Index cihi = cdata->getGhostBox().upper();

   const tbox::Pointer<CartesianPatchGeometry> fgeom =
      fine.getPatchGeometry();
   const tbox::Pointer<CartesianPatchGeometry> cgeom =
      coarse.getPatchGeometry();

   const hier::Index ifirstc = coarse_box.lower();
   const hier::Index ilastc = coarse_box.upper();

   for (int d = 0; d < cdata->getDepth(); d++) {
      if ((dim == tbox::Dimension(1))) {
         F77_FUNC(cartwgtavgcelldoub1d, CARTWGTAVGCELLDOUB1D) (ifirstc(0),
            ilastc(0),
            filo(0), fihi(0),
            cilo(0), cihi(0),
            &ratio[0],
            fgeom->getDx(),
            cgeom->getDx(),
            fdata->getPointer(d),
            cdata->getPointer(d));
      } else if ((dim == tbox::Dimension(2))) {
         F77_FUNC(cartwgtavgcelldoub2d, CARTWGTAVGCELLDOUB2D) (ifirstc(0),
            ifirstc(1), ilastc(0), ilastc(1),
            filo(0), filo(1), fihi(0), fihi(1),
            cilo(0), cilo(1), cihi(0), cihi(1),
            &ratio[0],
            fgeom->getDx(),
            cgeom->getDx(),
            fdata->getPointer(d),
            cdata->getPointer(d));
      } else if ((dim == tbox::Dimension(3))) {
         F77_FUNC(cartwgtavgcelldoub3d, CARTWGTAVGCELLDOUB3D) (ifirstc(0),
            ifirstc(1), ifirstc(2),
            ilastc(0), ilastc(1), ilastc(2),
            filo(0), filo(1), filo(2),
            fihi(0), fihi(1), fihi(2),
            cilo(0), cilo(1), cilo(2),
            cihi(0), cihi(1), cihi(2),
            &ratio[0],
            fgeom->getDx(),
            cgeom->getDx(),
            fdata->getPointer(d),
            cdata->getPointer(d));
      } else {
         TBOX_ERROR("CartesianCellDoubleWeightedAverage error...\n"
            << "dim > 3 not supported." << std::endl);

      }
   }
}

}
}
#endif
