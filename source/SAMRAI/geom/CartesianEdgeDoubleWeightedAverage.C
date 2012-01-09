/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Weighted averaging operator for edge-centered double data on
 *                a Cartesian mesh.
 *
 ************************************************************************/

#ifndef included_geom_CartesianEdgeDoubleWeightedAverage_C
#define included_geom_CartesianEdgeDoubleWeightedAverage_C

#include "SAMRAI/geom/CartesianEdgeDoubleWeightedAverage.h"

#include <float.h>
#include <math.h>
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/EdgeData.h"
#include "SAMRAI/pdat/EdgeVariable.h"
#include "SAMRAI/tbox/Utilities.h"

/*
 *************************************************************************
 *
 * External declarations for FORTRAN  routines.
 *
 *************************************************************************
 */
extern "C" {

#ifdef __INTEL_COMPILER
#pragma warning (disable:1419)
#endif

// in cartcoarsen1d.f:
void F77_FUNC(cartwgtavgedgedoub1d, CARTWGTAVGEDGEDOUB1D) (const int&,
   const int&,
   const int&, const int&,
   const int&, const int&,
   const int *, const double *, const double *,
   const double *, double *);
// in cartcoarsen2d.f:
void F77_FUNC(cartwgtavgedgedoub2d0, CARTWGTAVGEDGEDOUB2D0) (const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *, const double *, const double *,
   const double *, double *);

void F77_FUNC(cartwgtavgedgedoub2d1, CARTWGTAVGEDGEDOUB2D1) (const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *, const double *, const double *,
   const double *, double *);
// in cartcoarsen3d.f:
void F77_FUNC(cartwgtavgedgedoub3d0, CARTWGTAVGEDGEDOUB3D0) (const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *, const double *, const double *,
   const double *, double *);
void F77_FUNC(cartwgtavgedgedoub3d1, CARTWGTAVGEDGEDOUB3D1) (const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *, const double *, const double *,
   const double *, double *);
void F77_FUNC(cartwgtavgedgedoub3d2, CARTWGTAVGEDGEDOUB3D2) (const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *, const double *, const double *,
   const double *, double *);
}

namespace SAMRAI {
namespace geom {

// using namespace std;

CartesianEdgeDoubleWeightedAverage::CartesianEdgeDoubleWeightedAverage(
   const tbox::Dimension& dim):
   hier::CoarsenOperator(dim, "CONSERVATIVE_COARSEN")
{
}

CartesianEdgeDoubleWeightedAverage::~CartesianEdgeDoubleWeightedAverage()
{
}

bool CartesianEdgeDoubleWeightedAverage::findCoarsenOperator(
   const tbox::Pointer<hier::Variable>& var,
   const std::string& op_name) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *var);

   const tbox::Pointer<pdat::EdgeVariable<double> > cast_var(
      var,
      tbox::__dynamic_cast_tag());
   if (cast_var && (op_name == getOperatorName())) {
      return true;
   } else {
      return false;
   }
}

int CartesianEdgeDoubleWeightedAverage::getOperatorPriority() const
{
   return 0;
}

hier::IntVector
CartesianEdgeDoubleWeightedAverage::getStencilWidth() const {
   return hier::IntVector::getZero(getDim());
}

void CartesianEdgeDoubleWeightedAverage::coarsen(
   hier::Patch& coarse,
   const hier::Patch& fine,
   const int dst_component,
   const int src_component,
   const hier::Box& coarse_box,
   const hier::IntVector& ratio) const
{
   const tbox::Dimension& dim(getDim());

   TBOX_DIM_ASSERT_CHECK_DIM_ARGS4(dim, coarse, fine, coarse_box, ratio);

   tbox::Pointer<pdat::EdgeData<double> > fdata(
       fine.getPatchData(src_component),
      tbox::__dynamic_cast_tag());
   tbox::Pointer<pdat::EdgeData<double> > cdata(
      coarse.getPatchData(dst_component),
      tbox::__dynamic_cast_tag());
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(fdata);
   TBOX_ASSERT(cdata);
   TBOX_ASSERT(cdata->getDepth() == fdata->getDepth());
#endif

   const hier::Index filo = fdata->getGhostBox().lower();
   const hier::Index fihi = fdata->getGhostBox().upper();
   const hier::Index cilo = cdata->getGhostBox().lower();
   const hier::Index cihi = cdata->getGhostBox().upper();

   const tbox::Pointer<CartesianPatchGeometry> fgeom(
      fine.getPatchGeometry(),
      tbox::__dynamic_cast_tag());
   const tbox::Pointer<CartesianPatchGeometry> cgeom(
      coarse.getPatchGeometry(),
      tbox::__dynamic_cast_tag());

   const hier::Index ifirstc = coarse_box.lower();
   const hier::Index ilastc = coarse_box.upper();

   for (int d = 0; d < cdata->getDepth(); d++) {
      if ((dim == tbox::Dimension(1))) {
         F77_FUNC(cartwgtavgedgedoub1d, CARTWGTAVGEDGEDOUB1D) (ifirstc(0),
            ilastc(0),
            filo(0), fihi(0),
            cilo(0), cihi(0),
            &ratio[0],
            fgeom->getDx(),
            cgeom->getDx(),
            fdata->getPointer(0, d),
            cdata->getPointer(0, d));
      } else if ((dim == tbox::Dimension(2))) {
         F77_FUNC(cartwgtavgedgedoub2d0, CARTWGTAVGEDGEDOUB2D0) (ifirstc(0),
            ifirstc(1), ilastc(0), ilastc(1),
            filo(0), filo(1), fihi(0), fihi(1),
            cilo(0), cilo(1), cihi(0), cihi(1),
            &ratio[0],
            fgeom->getDx(),
            cgeom->getDx(),
            fdata->getPointer(0, d),
            cdata->getPointer(0, d));
         F77_FUNC(cartwgtavgedgedoub2d1, CARTWGTAVGEDGEDOUB2D1) (ifirstc(0),
            ifirstc(1), ilastc(0), ilastc(1),
            filo(0), filo(1), fihi(0), fihi(1),
            cilo(0), cilo(1), cihi(0), cihi(1),
            &ratio[0],
            fgeom->getDx(),
            cgeom->getDx(),
            fdata->getPointer(1, d),
            cdata->getPointer(1, d));
      } else if ((dim == tbox::Dimension(3))) {
         F77_FUNC(cartwgtavgedgedoub3d0, CARTWGTAVGEDGEDOUB3D0) (ifirstc(0),
            ifirstc(1), ifirstc(2),
            ilastc(0), ilastc(1), ilastc(2),
            filo(0), filo(1), filo(2),
            fihi(0), fihi(1), fihi(2),
            cilo(0), cilo(1), cilo(2),
            cihi(0), cihi(1), cihi(2),
            &ratio[0],
            fgeom->getDx(),
            cgeom->getDx(),
            fdata->getPointer(0, d),
            cdata->getPointer(0, d));
         F77_FUNC(cartwgtavgedgedoub3d1, CARTWGTAVGEDGEDOUB3D1) (ifirstc(0),
            ifirstc(1), ifirstc(2),
            ilastc(0), ilastc(1), ilastc(2),
            filo(0), filo(1), filo(2),
            fihi(0), fihi(1), fihi(2),
            cilo(0), cilo(1), cilo(2),
            cihi(0), cihi(1), cihi(2),
            &ratio[0],
            fgeom->getDx(),
            cgeom->getDx(),
            fdata->getPointer(1, d),
            cdata->getPointer(1, d));
         F77_FUNC(cartwgtavgedgedoub3d2, CARTWGTAVGEDGEDOUB3D2) (ifirstc(0),
            ifirstc(1), ifirstc(2),
            ilastc(0), ilastc(1), ilastc(2),
            filo(0), filo(1), filo(2),
            fihi(0), fihi(1), fihi(2),
            cilo(0), cilo(1), cilo(2),
            cihi(0), cihi(1), cihi(2),
            &ratio[0],
            fgeom->getDx(),
            cgeom->getDx(),
            fdata->getPointer(2, d),
            cdata->getPointer(2, d));
      } else {
         TBOX_ERROR("CartesianEdgeDoubleWeightedAverage error...\n"
            << "dim > 3 not supported." << std::endl);
      }

   }
}

}
}
#endif
