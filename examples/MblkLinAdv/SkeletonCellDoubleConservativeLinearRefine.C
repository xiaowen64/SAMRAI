//
// File:	SkeletonCellDoubleConservativeLinearRefine.C
// Package:	SAMRAI geometry
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Conservative linear refine operator for cell-centered 
//              double data on a Skeleton mesh.
//


#include "SkeletonCellDoubleConservativeLinearRefine.h"

#include<float.h>
#include<math.h>
#include "BlockPatchGeometry.h"
#include "Index.h"
#include "CellData.h"
#include "CellVariable.h"
#include "tbox/ArenaManager.h"
#include "tbox/Array.h"
#include "tbox/Utilities.h"

/*
*************************************************************************
*                                                                       *
* External declarations for FORTRAN  routines.                          *
*                                                                       *
*************************************************************************
*/

extern "C" {
// in cartrefine1d.f:
   void cartclinrefcelldoub1d_( const int&, const int&,
                                const int&, const int&,
                                const int&, const int&,
                                const int&, const int&,
                                const int*, const double*, const double*,
                                const double*, double*,
                                double*, double* );
// in cartrefine2d.f:
   void cartclinrefcelldoub2d_( const int&, const int&, const int&, const int&,
                                const int&, const int&, const int&, const int&,
                                const int&, const int&, const int&, const int&,
                                const int&, const int&, const int&, const int&,
                                const int*, const double*, const double*,
                                const double*, double*,
                                double*, double*, double*, double* );
// in cartrefine3d.f:
   void cartclinrefcelldoub3d_( const int&, const int&, const int&,
                                const int&, const int&, const int&,
                                const int&, const int&, const int&,
                                const int&, const int&, const int&,
                                const int&, const int&, const int&,
                                const int&, const int&, const int&,
                                const int&, const int&, const int&,
                                const int&, const int&, const int&,
                                const int*, const double*, const double*,
                                const double*, double*,
                                double*, double*, double*,
                                double*, double*, double* );
}

using namespace SAMRAI;


SkeletonCellDoubleConservativeLinearRefine::SkeletonCellDoubleConservativeLinearRefine()
: xfer::RefineOperator<NDIM>()
{
   d_name_id = "CONSERVATIVE_LINEAR_REFINE";

   const int max_levels = 10;
   d_dx.resizeArray(max_levels);
   for (int n = 0; n < max_levels; n++) {
      d_dx[n].resizeArray(NDIM);
      for (int i = 0; i < NDIM; i++) {
         d_dx[n][i] = 1.;
      }
   }
   
}

SkeletonCellDoubleConservativeLinearRefine::~SkeletonCellDoubleConservativeLinearRefine()
{
}

bool SkeletonCellDoubleConservativeLinearRefine::findRefineOperator(
   const tbox::Pointer< hier::Variable<NDIM> >& var,
   const string &op_name) const
{
   const tbox::Pointer< pdat::CellVariable<NDIM,double> > cast_var(var);
   if ( !cast_var.isNull() && (op_name == d_name_id) ) {
      return(true);
   } else {
      return(false);
   }
}

const string&
SkeletonCellDoubleConservativeLinearRefine::getOperatorName() const
{
   return(d_name_id);
}

int 
SkeletonCellDoubleConservativeLinearRefine::getOperatorPriority() const
{
   return(0);
}

hier::IntVector<NDIM> 
SkeletonCellDoubleConservativeLinearRefine::getStencilWidth() const {
   return(hier::IntVector<NDIM>(1));
}

void SkeletonCellDoubleConservativeLinearRefine::refine(
   hier::Patch<NDIM>& fine, 
   const hier::Patch<NDIM>& coarse, 
   const int dst_component, 
   const int src_component, 
   const hier::Box<NDIM>& fine_box, 
   const hier::IntVector<NDIM>& ratio) const
{
   tbox::Pointer< pdat::CellData<NDIM,double> >
      cdata = coarse.getPatchData(src_component);
   tbox::Pointer< pdat::CellData<NDIM,double> >
      fdata = fine.getPatchData(dst_component);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!cdata.isNull());
   TBOX_ASSERT(!fdata.isNull());
   TBOX_ASSERT(cdata->getDepth() == fdata->getDepth());
#endif

   const hier::Box<NDIM> cgbox(cdata->getGhostBox());

   const hier::Index<NDIM> cilo = cgbox.lower();
   const hier::Index<NDIM> cihi = cgbox.upper();
   const hier::Index<NDIM> filo = fdata->getGhostBox().lower();
   const hier::Index<NDIM> fihi = fdata->getGhostBox().upper();

   const tbox::Pointer<geom::BlockPatchGeometry<NDIM> > cgeom =
                                                    coarse.getPatchGeometry();
   const tbox::Pointer<geom::BlockPatchGeometry<NDIM> > fgeom =
                                                    fine.getPatchGeometry();

   const hier::Box<NDIM> coarse_box = hier::Box<NDIM>::coarsen(fine_box, ratio);
   const hier::Index<NDIM> ifirstc = coarse_box.lower();
   const hier::Index<NDIM> ilastc = coarse_box.upper();
   const hier::Index<NDIM> ifirstf = fine_box.lower();
   const hier::Index<NDIM> ilastf = fine_box.upper();    

   tbox::Pointer<tbox::Arena> 
   memory_arena = tbox::ArenaManager::getManager()->getScratchAllocator();

   const hier::IntVector<NDIM> tmp_ghosts(0);
   tbox::Array<double> diff0(cgbox.numberCells(0)+1, memory_arena);
   pdat::CellData<NDIM,double> slope0(cgbox, 1, tmp_ghosts, memory_arena);

   int flev_num = fine.getPatchLevelNumber();
   int clev_num = coarse.getPatchLevelNumber();

   // deal with levels not in hierarchy
   if (flev_num < 0) flev_num = clev_num + 1;
   if (clev_num < 0) clev_num = flev_num - 1;

   double cdx[NDIM];
   double fdx[NDIM];
   getDx(clev_num, cdx);
   getDx(flev_num, fdx);

   for (int d = 0; d < fdata->getDepth(); d++) {
      if (NDIM == 1) {
	 cartclinrefcelldoub1d_(ifirstc(0),ilastc(0),
				ifirstf(0),ilastf(0),
				cilo(0),cihi(0),
				filo(0),fihi(0),
				ratio,
				cdx,
				fdx,
				cdata->getPointer(d),
				fdata->getPointer(d),
				diff0.getPointer(),slope0.getPointer());
      } else if (NDIM == 2) {
	 
	 tbox::Array<double> diff1(cgbox.numberCells(1)+1, memory_arena);
	 pdat::CellData<NDIM,double> slope1(cgbox, 1, tmp_ghosts, memory_arena);

	 cartclinrefcelldoub2d_(ifirstc(0),ifirstc(1),ilastc(0),ilastc(1),
				ifirstf(0),ifirstf(1),ilastf(0),ilastf(1),
				cilo(0),cilo(1),cihi(0),cihi(1),
				filo(0),filo(1),fihi(0),fihi(1),
				ratio,
				cdx,
				fdx,
                                cdata->getPointer(d),
				fdata->getPointer(d),
				diff0.getPointer(),slope0.getPointer(),
				diff1.getPointer(),slope1.getPointer());
      } else if (NDIM == 3) {

	 tbox::Array<double> diff1(cgbox.numberCells(1)+1, memory_arena);
	 pdat::CellData<NDIM,double> slope1(cgbox, 1, tmp_ghosts, memory_arena);

	 tbox::Array<double> diff2(cgbox.numberCells(2)+1, memory_arena);
	 pdat::CellData<NDIM,double> slope2(cgbox, 1, tmp_ghosts, memory_arena);

	 cartclinrefcelldoub3d_(ifirstc(0),ifirstc(1),ifirstc(2),
				ilastc(0),ilastc(1),ilastc(2),
				ifirstf(0),ifirstf(1),ifirstf(2),
				ilastf(0),ilastf(1),ilastf(2),
				cilo(0),cilo(1),cilo(2),
				cihi(0),cihi(1),cihi(2),
				filo(0),filo(1),filo(2),
				fihi(0),fihi(1),fihi(2),
				ratio,
				cdx,
				fdx,
				cdata->getPointer(d),
				fdata->getPointer(d),
				diff0.getPointer(),slope0.getPointer(),
				diff1.getPointer(),slope1.getPointer(),
				diff2.getPointer(),slope2.getPointer());
      } else {
	 TBOX_ERROR("SkeletonCellDoubleConservativeLinearRefine error...\n"
		    << "NDIM > 3 not supported." << endl);

      }
   }
}

void SkeletonCellDoubleConservativeLinearRefine::setDx(
   const int level_number,
   const double* dx)
{
   if (level_number >= d_dx.getSize()) {
      d_dx.resizeArray(level_number+1);
      d_dx[level_number].resizeArray(NDIM);
      for (int i = 0; i < NDIM; i++) {
         d_dx[level_number][i] = dx[i];
      }
   }
}

void SkeletonCellDoubleConservativeLinearRefine::getDx(
   const int level_number,
   double* dx) const
{
   for (int i = 0; i < NDIM; i++) {
      dx[i] = d_dx[level_number][i];
   }
}





