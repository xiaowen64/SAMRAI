/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2018 Lawrence Livermore National Security, LLC
 * Description:   Linear refine operator for node-centered double data on
 *                a Cartesian mesh.
 *
 ************************************************************************/
#include "SAMRAI/geom/CartesianNodeDoubleLinearRefine.h"

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
void SAMRAI_F77_FUNC(cartlinrefnodedoub1d, CARTLINREFNODEDOUB1D) (const int&,
   const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *, const double *, const double *,
   const double *, double *);
// in cartrefine2d.f:
void SAMRAI_F77_FUNC(cartlinrefnodedoub2d, CARTLINREFNODEDOUB2D) (const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *, const double *, const double *,
   const double *, double *);
// in cartrefine3d.f:
void SAMRAI_F77_FUNC(cartlinrefnodedoub3d, CARTLINREFNODEDOUB3D) (const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
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


CartesianNodeDoubleLinearRefine::CartesianNodeDoubleLinearRefine():
   hier::RefineOperator("LINEAR_REFINE")
{
}

CartesianNodeDoubleLinearRefine::~CartesianNodeDoubleLinearRefine()
{
}

int
CartesianNodeDoubleLinearRefine::getOperatorPriority() const
{
   return 0;
}

hier::IntVector
CartesianNodeDoubleLinearRefine::getStencilWidth(const tbox::Dimension& dim) const
{
   return hier::IntVector::getZero(dim);
}

void
CartesianNodeDoubleLinearRefine::refine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const int dst_component,
   const int src_component,
   const hier::BoxOverlap& fine_overlap,
   const hier::IntVector& ratio) const
{
   const pdat::NodeOverlap* t_overlap =
      CPP_CAST<const pdat::NodeOverlap *>(&fine_overlap);

   TBOX_ASSERT(t_overlap != 0);

   const hier::BoxContainer& boxes = t_overlap->getDestinationBoxContainer();
   for (hier::BoxContainer::const_iterator b = boxes.begin();
        b != boxes.end(); ++b) {
      hier::Box fine_box(*b);
      fine_box.growUpper(hier::IntVector(ratio.getDim(), -1));
      refine(fine,
         coarse,
         dst_component,
         src_component,
         fine_box,
         ratio);
   }
}

void
CartesianNodeDoubleLinearRefine::refine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const int dst_component,
   const int src_component,
   const hier::Box& fine_box,
   const hier::IntVector& ratio) const
{
   const tbox::Dimension& dim(fine.getDim());
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY3(dim, coarse, fine_box, ratio);

   std::shared_ptr<pdat::NodeData<double> > cdata(
      SAMRAI_SHARED_PTR_CAST<pdat::NodeData<double>, hier::PatchData>(
         coarse.getPatchData(src_component)));
   std::shared_ptr<pdat::NodeData<double> > fdata(
      SAMRAI_SHARED_PTR_CAST<pdat::NodeData<double>, hier::PatchData>(
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

   for (int d = 0; d < fdata->getDepth(); ++d) {
      if ((dim == tbox::Dimension(1))) {
         SAMRAI_F77_FUNC(cartlinrefnodedoub1d, CARTLINREFNODEDOUB1D) (ifirstc(0),
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
#if defined(HAVE_RAJA)
      // this kernel has iterations writing to the same address
      // we're going to use atomics until kernel undergoes rewrite or we can verify thread-safe  
      auto fine_array = fdata->getView<2>(d);
      auto coarse_array = cdata->getView<2>(d);

      const int r0 = ratio[0];
      const int r1 = ratio[1];

      const int filo0 = filo(0);
      const int filo1 = filo(1);
      const int fihi0 = fihi(0);
      const int fihi1 = fihi(1);

      pdat::parallel_for_all_x(coarse_box, [=] SAMRAI_HOST_DEVICE (int j /*fast*/, int k) {
         const int ic0 = j;
         const int ic1 = k;
         const int if0 = ic0*r0;
         const int if1 = ic1*r1;
         
         const double realrat0 = 1.0/static_cast<double>(r0);
         const double realrat1 = 1.0/static_cast<double>(r1);

         for(int ir1 = 0; ir1 <= r1; ++ir1) {
           int ie1 = if1+ir1;
           //fprintf(stderr,"ir1=%d if1=%d ie1=%d filo1=%d fihi1=%d\n",ir1,if1,ie1,filo1,fihi1);
           if(ie1 >= filo1 && ie1 <= fihi1+1) {
             for(int ir0 = 0; ir0 <= r0; ++ir0) {
               int ie0 = if0+ir0;

               if(ie0 >= filo0 && ie0 <= fihi0+1) {
                 double x = double(ir0) * realrat0;
                 double y = double(ir1) * realrat1;
                 //fprintf(stderr,"ie0=%d ie1=%d filo0=%d fihi0=%d filo1=%d fihi1=%d r0=%d r1=%d\n",
                 //        ie0,ie1,filo0,fihi0,filo1,fihi1,r0,r1); 
                 double fineValue =
                   (coarse_array(ic0,ic1)*(1.0-x) +
                    coarse_array(ic0+1,ic1)*x)*(1.0-y) +
                   (coarse_array(ic0,ic1+1)*(1.0-x) +
                    coarse_array(ic0+1,ic1+1)*x)*y;
                 //fine_array(ie0,ie1) = fineValue;
                 RAJA::atomic::atomicExchange<RAJA::atomic::auto_atomic>(&fine_array(ie0,ie1),fineValue);
                 //fprintf(stderr,"fine_array(%d,%d) = %f @ %p\n",ie0,ie1,fine_array(ie0,ie1),&fine_array(ie0,ie1));
               }
             }
           }
         }
      });
#else // Fortran Dimension 2

         SAMRAI_F77_FUNC(cartlinrefnodedoub2d, CARTLINREFNODEDOUB2D) (ifirstc(0),
            ifirstc(1), ilastc(0), ilastc(1),
            ifirstf(0), ifirstf(1), ilastf(0), ilastf(1),
            cilo(0), cilo(1), cihi(0), cihi(1),
            filo(0), filo(1), fihi(0), fihi(1),
            &ratio[0],
            cgeom->getDx(),
            fgeom->getDx(),
            cdata->getPointer(d),
            fdata->getPointer(d));

         //exit(-1);
#endif // test for RAJA

      } else if ((dim == tbox::Dimension(3))) {
#if defined(HAVE_RAJA)
      // this kernel has iterations writing to the same address
      // we're going to use atomics until kernel undergoes rewrite or we can verify thread-safe  
      auto fine_array = fdata->getView<3>(d);
      auto coarse_array = cdata->getView<3>(d);

      const int r0 = ratio[0];
      const int r1 = ratio[1];
      const int r2 = ratio[2];

      const int filo0 = filo(0);
      const int filo1 = filo(1);
      const int filo2 = filo(2);

      const int fihi0 = fihi(0);
      const int fihi1 = fihi(1);
      const int fihi2 = fihi(2);

      pdat::parallel_for_all_x(coarse_box, [=] SAMRAI_HOST_DEVICE (int i /*fast*/, int j, int k) {
         const int ic0 = i;
         const int ic1 = j;
         const int ic2 = k;

         const int if0 = ic0*r0;
         const int if1 = ic1*r1;
         const int if2 = ic2*r2;
         
         const double realrat0 = 1.0/static_cast<double>(r0);
         const double realrat1 = 1.0/static_cast<double>(r1);
         const double realrat2 = 1.0/static_cast<double>(r2);

         for(int ir2 = 0; ir2 <= r2; ++ir2) {
           int ie2 = if2+ir2;
           if(ie2 >= filo2 && ie2 <= (fihi2+1)) {
             for(int ir1 = 0; ir1 <= r1; ++ir1) {
               int ie1 = if1+ir1;
               if(ie1 >= filo1 && ie1 <= fihi1+1) {
                 for(int ir0 = 0; ir0 <= r0; ++ir0) {
                   int ie0 = if0+ir0;
                   if(ie0 >= filo0 && ie0 <= fihi0+1) {
                     double x = double(ir0) * realrat0;
                     double y = double(ir1) * realrat1;
                     double z = double(ir2) * realrat2;
                     double fineValue =
                       ((coarse_array(ic0,ic1,ic2)*(1.0-x) +
                         coarse_array(ic0+1,ic1,ic2)*x)*(1.0-y)
                       + (coarse_array(ic0,ic1+1,ic2)*(1.0-x) +
                         coarse_array(ic0+1,ic1+1,ic2)*x)*y) * (1.0-z) +
                       ((coarse_array(ic0,ic1,ic2+1)*(1.0-x) +
                         coarse_array(ic0+1,ic1,ic2+1)*x)*(1.0-y)
                       + (coarse_array(ic0,ic1+1,ic2+1)*(1.0-x) +
                          coarse_array(ic0+1,ic1+1,ic2+1)*x)*y)*z; 
                     RAJA::atomic::atomicExchange<RAJA::atomic::auto_atomic>(&fine_array(ie0,ie1,ie2),fineValue);
                   }
                 }
               }
             }
           }
         }
      });
#else // Fortran Dimension 3
         SAMRAI_F77_FUNC(cartlinrefnodedoub3d, CARTLINREFNODEDOUB3D) (ifirstc(0),
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
         //exit(-1);
#endif // test for RAJA

      } else {
         TBOX_ERROR("CartesianNodeDoubleLinearRefine error...\n"
            << "dim > 3 not supported." << std::endl);
      }
   }
}

}
}
