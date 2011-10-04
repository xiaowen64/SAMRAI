/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Linear refine operator for node-centered complex data on
 *                a Cartesian mesh.
 *
 ************************************************************************/

#ifndef included_geom_CartesianNodeComplexLinearRefine_C
#define included_geom_CartesianNodeComplexLinearRefine_C

#include "SAMRAI/geom/CartesianNodeComplexLinearRefine.h"
#include "SAMRAI/tbox/Complex.h"

#include <float.h>
#include <math.h>
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/pdat/NodeVariable.h"
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

// in cartrefine1d.f:
void F77_FUNC(cartlinrefnodecplx1d, CARTLINREFNODECPLX1D) (const int&,
   const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *, const double *, const double *,
   const dcomplex *, dcomplex *);
// in cartrefine2d.f:
void F77_FUNC(cartlinrefnodecplx2d, CARTLINREFNODECPLX2D) (const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *, const double *, const double *,
   const dcomplex *, dcomplex *);
// in cartrefine3d.f:
void F77_FUNC(cartlinrefnodecplx3d, CARTLINREFNODECPLX3D) (const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *, const double *, const double *,
   const dcomplex *, dcomplex *);
}

namespace SAMRAI {
namespace geom {

// using namespace std;

CartesianNodeComplexLinearRefine::CartesianNodeComplexLinearRefine(
   const tbox::Dimension& dim):
   hier::RefineOperator(dim, "LINEAR_REFINE")
{
}

CartesianNodeComplexLinearRefine::~CartesianNodeComplexLinearRefine()
{
}

bool CartesianNodeComplexLinearRefine::findRefineOperator(
   const tbox::Pointer<hier::Variable>& var,
   const std::string& op_name) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *var);

   const tbox::Pointer<pdat::NodeVariable<dcomplex> > cast_var(var);
   if (!cast_var.isNull() && (op_name == getOperatorName())) {
      return true;
   } else {
      return false;
   }
}

int CartesianNodeComplexLinearRefine::getOperatorPriority() const
{
   return 0;
}

hier::IntVector
CartesianNodeComplexLinearRefine::getStencilWidth() const {
   return hier::IntVector::getZero(getDim());
}

void CartesianNodeComplexLinearRefine::refine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const int dst_component,
   const int src_component,
   const hier::BoxOverlap& fine_overlap,
   const hier::IntVector& ratio) const
{
   const pdat::NodeOverlap* t_overlap =
      dynamic_cast<const pdat::NodeOverlap *>(&fine_overlap);

   TBOX_ASSERT(t_overlap != NULL);

   const hier::BoxList& boxes = t_overlap->getDestinationBoxList();
   for (hier::BoxList::ConstIterator b(boxes); b != boxes.end(); ++b) {
      hier::Box fine_box(b());
      fine_box.growUpper(hier::IntVector(ratio.getDim(), -1));
      refine(fine,
         coarse,
         dst_component,
         src_component,
         fine_box,
         ratio);
   }
}

void CartesianNodeComplexLinearRefine::refine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const int dst_component,
   const int src_component,
   const hier::Box& fine_box,
   const hier::IntVector& ratio) const
{
   const tbox::Dimension& dim(getDim());
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS4(dim, fine, coarse, fine_box, ratio);

   tbox::Pointer<pdat::NodeData<dcomplex> >
   cdata = coarse.getPatchData(src_component);
   tbox::Pointer<pdat::NodeData<dcomplex> >
   fdata = fine.getPatchData(dst_component);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!cdata.isNull());
   TBOX_ASSERT(!fdata.isNull());
   TBOX_ASSERT(cdata->getDepth() == fdata->getDepth());
#endif

   const hier::Box cgbox(cdata->getGhostBox());

   const hier::Index cilo = cgbox.lower();
   const hier::Index cihi = cgbox.upper();
   const hier::Index filo = fdata->getGhostBox().lower();
   const hier::Index fihi = fdata->getGhostBox().upper();

   const tbox::Pointer<CartesianPatchGeometry> cgeom =
      coarse.getPatchGeometry();
   const tbox::Pointer<CartesianPatchGeometry> fgeom =
      fine.getPatchGeometry();

   const hier::Box coarse_box = hier::Box::coarsen(fine_box, ratio);
   const hier::Index ifirstc = coarse_box.lower();
   const hier::Index ilastc = coarse_box.upper();
   const hier::Index ifirstf = fine_box.lower();
   const hier::Index ilastf = fine_box.upper();

   for (int d = 0; d < fdata->getDepth(); d++) {
      if ((dim == tbox::Dimension(1))) {
         F77_FUNC(cartlinrefnodecplx1d, CARTLINREFNODECPLX1D) (ifirstc(0),
            ilastc(0),
            ifirstf(0), ilastf(0),
            cilo(0), cihi(0),
            filo(0), fihi(0),
            &ratio[0],
            cgeom->getDx(),
            fgeom->getDx(),
            cdata->getPointer(d),
            fdata->getPointer(d));
      } else if ((dim == tbox::Dimension(2))) {
         F77_FUNC(cartlinrefnodecplx2d, CARTLINREFNODECPLX2D) (ifirstc(0),
            ifirstc(1), ilastc(0), ilastc(1),
            ifirstf(0), ifirstf(1), ilastf(0), ilastf(1),
            cilo(0), cilo(1), cihi(0), cihi(1),
            filo(0), filo(1), fihi(0), fihi(1),
            &ratio[0],
            cgeom->getDx(),
            fgeom->getDx(),
            cdata->getPointer(d),
            fdata->getPointer(d));
      } else if ((dim == tbox::Dimension(3))) {
         F77_FUNC(cartlinrefnodecplx3d, CARTLINREFNODECPLX3D) (ifirstc(0),
            ifirstc(1), ifirstc(2),
            ilastc(0), ilastc(1), ilastc(2),
            ifirstf(0), ifirstf(1), ifirstf(2),
            ilastf(0), ilastf(1), ilastf(2),
            cilo(0), cilo(1), cilo(2),
            cihi(0), cihi(1), cihi(2),
            filo(0), filo(1), filo(2),
            fihi(0), fihi(1), fihi(2),
            &ratio[0],
            cgeom->getDx(),
            fgeom->getDx(),
            cdata->getPointer(d),
            fdata->getPointer(d));
      } else {
         TBOX_ERROR("CartesianNodeComplexLinearRefine error...\n"
            << "dim > 3 not supported." << std::endl);
      }
   }
}

}
}
#endif
