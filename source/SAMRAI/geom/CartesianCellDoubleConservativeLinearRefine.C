/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2019 Lawrence Livermore National Security, LLC
 * Description:   Conservative linear refine operator for cell-centered
 *                double data on a Cartesian mesh.
 *
 ************************************************************************/
#include "SAMRAI/geom/CartesianCellDoubleConservativeLinearRefine.h"
#include <float.h>
#include <math.h>
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/ForAll.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/Utilities.h"

#include "SAMRAI/tbox/AllocatorDatabase.h"

#define MAX(a, b) (((b) > (a)) ? (b) : (a))
#define MIN(a, b) (((b) < (a)) ? (b) : (a))
#define ABS abs
#define SQRT sqrt
#define COPYSIGN copysign

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
void SAMRAI_F77_FUNC(cartclinrefcelldoub1d, CARTCLINREFCELLDOUB1D) (const int&,
   const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *, const double *, const double *,
   const double *, double *,
   double *, double *);
// in cartrefine2d.f:
void SAMRAI_F77_FUNC(cartclinrefcelldoub2d, CARTCLINREFCELLDOUB2D) (const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *, const double *, const double *,
   const double *, double *,
   double *, double *, double *, double *);
// in cartrefine3d.f:
void SAMRAI_F77_FUNC(cartclinrefcelldoub3d, CARTCLINREFCELLDOUB3D) (const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *, const double *, const double *,
   const double *, double *,
   double *, double *, double *,
   double *, double *, double *);
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace geom {


CartesianCellDoubleConservativeLinearRefine::
CartesianCellDoubleConservativeLinearRefine():
   hier::RefineOperator("CONSERVATIVE_LINEAR_REFINE")
{
}

CartesianCellDoubleConservativeLinearRefine::~
CartesianCellDoubleConservativeLinearRefine()
{
}

int
CartesianCellDoubleConservativeLinearRefine::getOperatorPriority() const
{
   return 0;
}

hier::IntVector
CartesianCellDoubleConservativeLinearRefine::getStencilWidth(const tbox::Dimension& dim) const
{
   return hier::IntVector::getOne(dim);
}

void
CartesianCellDoubleConservativeLinearRefine::refine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const int dst_component,
   const int src_component,
   const hier::BoxOverlap& fine_overlap,
   const hier::IntVector& ratio) const
{
   const pdat::CellOverlap* t_overlap =
      CPP_CAST<const pdat::CellOverlap *>(&fine_overlap);

   TBOX_ASSERT(t_overlap != 0);

   const hier::BoxContainer& boxes = t_overlap->getDestinationBoxContainer();
   for (hier::BoxContainer::const_iterator b = boxes.begin();
        b != boxes.end(); ++b) {
      refine(fine,
         coarse,
         dst_component,
         src_component,
         *b,
         ratio);
   }
}

void
CartesianCellDoubleConservativeLinearRefine::refine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const int dst_component,
   const int src_component,
   const hier::Box& fine_box,
   const hier::IntVector& ratio) const
{
   RANGE_PUSH("ConservativeLinearRefine::refine", 3);

   const tbox::Dimension& dim(fine.getDim());
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY3(dim, coarse, fine_box, ratio);

   std::shared_ptr<pdat::CellData<double> > cdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         coarse.getPatchData(src_component)));
   std::shared_ptr<pdat::CellData<double> > fdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         fine.getPatchData(dst_component)));
   TBOX_ASSERT(cdata);
   TBOX_ASSERT(fdata);
   TBOX_ASSERT(cdata->getDepth() == fdata->getDepth());

   const hier::Box cgbox(cdata->getGhostBox());

   const hier::Index& cilo = cgbox.lower();
   const hier::Index& cihi = cgbox.upper();
   const hier::Index& filo = fdata->getGhostBox().lower();
   const hier::Index& fihi = fdata->getGhostBox().upper();

   const std::shared_ptr<CartesianPatchGeometry> cgeom(
      SAMRAI_SHARED_PTR_CAST<CartesianPatchGeometry, hier::PatchGeometry>(
         coarse.getPatchGeometry()));
   const std::shared_ptr<CartesianPatchGeometry> fgeom(
      SAMRAI_SHARED_PTR_CAST<CartesianPatchGeometry, hier::PatchGeometry>(
         fine.getPatchGeometry()));

   TBOX_ASSERT(cgeom);
   TBOX_ASSERT(fgeom);

   const hier::Box coarse_box = hier::Box::coarsen(fine_box, ratio);
   const hier::Index& ifirstc = coarse_box.lower();
   const hier::Index& ilastc = coarse_box.upper();
   const hier::Index& ifirstf = fine_box.lower();
   const hier::Index& ilastf = fine_box.upper();

   const hier::IntVector tmp_ghosts(dim, 0);


   SAMRAI::hier::Box diff_box = coarse_box;
   diff_box.growUpper(SAMRAI::hier::IntVector::getOne(dim));

   tbox::AllocatorDatabase* alloc_db = tbox::AllocatorDatabase::getDatabase();

   pdat::CellData<double> diff(diff_box, dim.getValue(), tmp_ghosts, alloc_db->getDevicePool());
   pdat::CellData<double> slope(cgbox, dim.getValue(), tmp_ghosts, alloc_db->getDevicePool());

   for (int d = 0; d < fdata->getDepth(); ++d) {
      if ((dim == tbox::Dimension(1))) {
         std::vector<double> diff0(cgbox.numberCells(0) + 1);

         SAMRAI_F77_FUNC(cartclinrefcelldoub1d, CARTCLINREFCELLDOUB1D) (ifirstc(0),
            ilastc(0),
            ifirstf(0), ilastf(0),
            cilo(0), cihi(0),
            filo(0), fihi(0),
            &ratio[0],
            cgeom->getDx(),
            fgeom->getDx(),
            cdata->getPointer(d),
            fdata->getPointer(d),
            &diff0[0], slope.getPointer());
      } else if ((dim == tbox::Dimension(2))) {

#if defined(HAVE_RAJA)
      auto fine_array = fdata->getView<2>(d);
      auto coarse_array = cdata->getView<2>(d);

      auto diff0 = diff.getView<2>(0);
      auto diff1 = diff.getView<2>(1);

      auto slope0 = slope.getView<2>(0);
      auto slope1 = slope.getView<2>(1);

      const double* fdx = fgeom->getDx();
      const double* cdx = cgeom->getDx();

      const int fdx0 = fdx[0];
      const int fdx1 = fdx[1];
      const int cdx0 = cdx[0];
      const int cdx1 = cdx[1];

      const int r0 = ratio[0];
      const int r1 = ratio[1];

      pdat::parallel_for_all(coarse_box, [=] SAMRAI_HOST_DEVICE (int k, int j) {
         diff0(j,k) = coarse_array(j,k) - coarse_array(j-1,k);
         diff1(j,k) = coarse_array(j,k) - coarse_array(j,k-1);
      });

      pdat::parallel_for_all(coarse_box, [=] SAMRAI_HOST_DEVICE (int k, int j) {
         const double coef2j = 0.5*(diff0(j+1,k)+diff0(j,k));
         const double boundj = 2.0*MIN(ABS(diff0(j+1,k)),ABS(diff0(j,k)));

         if (diff0(j,k)*diff0(j+1,k) > 0.0) {
            slope0(j,k) = COPYSIGN(MIN(ABS(coef2j),boundj),coef2j)/cdx0;
         } else {
            slope0(j,k) = 0.0;
         }

         const double coef2k = 0.5*(diff1(j,k+1)+diff1(j,k));
         const double boundk = 2.0*MIN(ABS(diff1(j,k+1)),ABS(diff1(j,k)));

         if (diff1(j,k)*diff1(j,k+1) > 0.0) {
            slope1(j,k) = COPYSIGN(MIN(ABS(coef2k),boundk),coef2k)/cdx1;
         } else {
            slope1(j,k) = 0.0;
         }
      });

      pdat::parallel_for_all(fine_box, [=] SAMRAI_HOST_DEVICE (int k, int j) {
         const int ic1 = (k < 0) ? (k+1)/r1-1 : k/r1;
         const int ic0 = (j < 0) ? (j+1)/r0-1 : j/r0;

         const int ir0 = j - ic0*r0;
         const int ir1 = k - ic1*r1;

         const double deltax1 = (static_cast<double>(ir1)+0.5)*fdx1-cdx1*0.5;
         const double deltax0 = (static_cast<double>(ir0)+0.5)*fdx0-cdx0*0.5;

         fine_array(j,k) = coarse_array(ic0,ic1) + slope0(ic0, ic1)*deltax0 + slope1(ic0,ic1)*deltax1;
      });
#else
         std::vector<double> diff1(cgbox.numberCells(1) + 1);
         pdat::CellData<double> slope1(cgbox, 1, tmp_ghosts);

         SAMRAI_F77_FUNC(cartclinrefcelldoub2d, CARTCLINREFCELLDOUB2D) (ifirstc(0),
            ifirstc(1), ilastc(0), ilastc(1),
            ifirstf(0), ifirstf(1), ilastf(0), ilastf(1),
            cilo(0), cilo(1), cihi(0), cihi(1),
            filo(0), filo(1), fihi(0), fihi(1),
            &ratio[0],
            cgeom->getDx(),
            fgeom->getDx(),
            cdata->getPointer(d),
            fdata->getPointer(d),
            &diff0[0], slope.getPointer(0),
            &diff1[0], slope.getPointer(1));
#endif
      } else if ((dim == tbox::Dimension(3))) {

         std::vector<double> diff0(cgbox.numberCells(0) + 1);
         std::vector<double> diff1(cgbox.numberCells(1) + 1);
         pdat::CellData<double> slope1(cgbox, 1, tmp_ghosts);

         std::vector<double> diff2(cgbox.numberCells(2) + 1);
         pdat::CellData<double> slope2(cgbox, 1, tmp_ghosts);

         SAMRAI_F77_FUNC(cartclinrefcelldoub3d, CARTCLINREFCELLDOUB3D) (ifirstc(0),
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
            fdata->getPointer(d),
            &diff0[0], slope.getPointer(0),
            &diff1[0], slope.getPointer(1),
            &diff2[0], slope.getPointer(2));
      } else {
         TBOX_ERROR("CartesianCellDoubleConservativeLinearRefine error...\n"
            << "dim > 3 not supported." << std::endl);

      }
   }

   RANGE_POP;
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
