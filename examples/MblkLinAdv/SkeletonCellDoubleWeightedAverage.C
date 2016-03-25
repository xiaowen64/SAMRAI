//
// File:	SkeletonCellDoubleWeightedAverage.C
// Package:	SAMRAI geometry
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description: Weighted averaging operator for cell-centered double data on 
//              a Skeleton mesh.
//

#include "SkeletonCellDoubleWeightedAverage.h"

#include<float.h>
#include<math.h>
#include "BlockPatchGeometry.h"
#include "Index.h"
#include "CellData.h"
#include "CellVariable.h"
#include "tbox/Utilities.h"

/*
*************************************************************************
*                                                                       *
* External declarations for FORTRAN  routines.                          *
*                                                                       *
*************************************************************************
*/
extern "C" {
// in cartcoarsen1d.f:
   void cartwgtavgcelldoub1d_( const int&, const int&,
                               const int&, const int&,
                               const int&, const int&,
                               const int*, const double*, const double*,
                               const double*, double* );
// in cartcoarsen2d.f:
   void cartwgtavgcelldoub2d_( const int&, const int&, const int&, const int&,
                               const int&, const int&, const int&, const int&,
                               const int&, const int&, const int&, const int&,
                               const int*, const double*, const double*,
                               const double*, double* );
// in cartcoarsen3d.f:
   void cartwgtavgcelldoub3d_( const int&, const int&, const int&,
                               const int&, const int&, const int&,
                               const int&, const int&, const int&,
                               const int&, const int&, const int&,
                               const int&, const int&, const int&,
                               const int&, const int&, const int&,
                               const int*, const double*, const double*,
                               const double*, double* );
}

using namespace SAMRAI;

SkeletonCellDoubleWeightedAverage::SkeletonCellDoubleWeightedAverage()
: xfer::CoarsenOperator<NDIM>()
{
   d_name_id = "CONSERVATIVE_COARSEN";
}

SkeletonCellDoubleWeightedAverage::~SkeletonCellDoubleWeightedAverage()
{
}

bool SkeletonCellDoubleWeightedAverage::findCoarsenOperator(
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
SkeletonCellDoubleWeightedAverage::getOperatorName() const
{
   return(d_name_id);
}

int SkeletonCellDoubleWeightedAverage::getOperatorPriority() const
{
   return(0);
}

hier::IntVector<NDIM> 
SkeletonCellDoubleWeightedAverage::getStencilWidth() const {
   return(hier::IntVector<NDIM>(0));
}

void SkeletonCellDoubleWeightedAverage::coarsen(
   hier::Patch<NDIM>& coarse, 
   const hier::Patch<NDIM>& fine, 
   const int dst_component, 
   const int src_component, 
   const hier::Box<NDIM>& coarse_box, 
   const hier::IntVector<NDIM>& ratio) const 
{
   tbox::Pointer< pdat::CellData<NDIM,double> > 
      fdata = fine.getPatchData(src_component);
   tbox::Pointer< pdat::CellData<NDIM,double> > 
      cdata = coarse.getPatchData(dst_component);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!fdata.isNull());
   TBOX_ASSERT(!cdata.isNull());
   TBOX_ASSERT(cdata->getDepth() == fdata->getDepth());
#endif

   const hier::Index<NDIM> filo = fdata->getGhostBox().lower();
   const hier::Index<NDIM> fihi = fdata->getGhostBox().upper();
   const hier::Index<NDIM> cilo = cdata->getGhostBox().lower();
   const hier::Index<NDIM> cihi = cdata->getGhostBox().upper();

   const tbox::Pointer<geom::BlockPatchGeometry<NDIM> > fgeom =
                                                    fine.getPatchGeometry();
   const tbox::Pointer<geom::BlockPatchGeometry<NDIM> > cgeom =
                                                    coarse.getPatchGeometry();

   const hier::Index<NDIM> ifirstc = coarse_box.lower();
   const hier::Index<NDIM> ilastc = coarse_box.upper();

   int flev_num = fine.getPatchLevelNumber();
   int clev_num = coarse.getPatchLevelNumber();

   // deal with levels not in hierarchy
   if (flev_num < 0) flev_num = clev_num + 1;
   if (clev_num < 0) clev_num = flev_num - 1;

   double cdx[NDIM];
   double fdx[NDIM];
   getDx(clev_num, cdx);
   getDx(flev_num, fdx);

   for (int d = 0; d < cdata->getDepth(); d++) {
      if (NDIM == 1) {
	 cartwgtavgcelldoub1d_(ifirstc(0),ilastc(0),
			       filo(0),fihi(0),
			       cilo(0),cihi(0),
			       ratio,
			       fdx,
			       cdx,
			       fdata->getPointer(d),
			       cdata->getPointer(d));
      } else if (NDIM == 2) {
	 cartwgtavgcelldoub2d_(ifirstc(0),ifirstc(1),ilastc(0),ilastc(1),
			       filo(0),filo(1),fihi(0),fihi(1),
			       cilo(0),cilo(1),cihi(0),cihi(1),
			       ratio,
			       fdx,
			       cdx,
			       fdata->getPointer(d),
			       cdata->getPointer(d));
      } else if (NDIM == 3) {
	 cartwgtavgcelldoub3d_(ifirstc(0),ifirstc(1),ifirstc(2),
			       ilastc(0),ilastc(1),ilastc(2),
			       filo(0),filo(1),filo(2),
			       fihi(0),fihi(1),fihi(2),
			       cilo(0),cilo(1),cilo(2),
			       cihi(0),cihi(1),cihi(2),
			       ratio,
			       fdx,
			       cdx,
			       fdata->getPointer(d),
			       cdata->getPointer(d));
      } else {
	 TBOX_ERROR("SkeletonCellDoubleWeightedAverage error...\n"
		    << "NDIM > 3 not supported." << endl);

      }
   }
}

void SkeletonCellDoubleWeightedAverage::setDx(
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

void SkeletonCellDoubleWeightedAverage::getDx(
   const int level_number,
   double* dx) const
{
   for (int i = 0; i < NDIM; i++) {
      dx[i] = d_dx[level_number][i];
   }
}


