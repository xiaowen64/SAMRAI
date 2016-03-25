//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/geometry/cartesian/operators/outerface/CartesianOuterfaceDoubleWeightedAverage.C $
// Package:	SAMRAI geometry
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description: Weighted averaging operator for outerface double data on 
//              a Cartesian mesh.
//

#ifndef included_geom_CartesianOuterfaceDoubleWeightedAverage_C
#define included_geom_CartesianOuterfaceDoubleWeightedAverage_C

#include "CartesianOuterfaceDoubleWeightedAverage.h"

#include<float.h>
#include<math.h>
#include "CartesianPatchGeometry.h"
#include "Index.h"
#include "OuterfaceData.h"
#include "OuterfaceVariable.h"
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
   void cartwgtavgoutfacedoub1d_( const int&, const int&,
                                  const int&, const int&,
                                  const int&, const int&,
                                  const int*, const double*, const double*,
                                  const double*, double* );
// in cartcoarsen2d.f:
   void cartwgtavgoutfacedoub2d0_( const int&, const int&, 
                                   const int&, const int&,
                                   const int&, const int&, 
                                   const int&, const int&,
                                   const int&, const int&, 
                                   const int&, const int&,
                                   const int*, const double*, const double*,
                                   const double*, double* );

   void cartwgtavgoutfacedoub2d1_( const int&, const int&,
                                   const int&, const int&,
                                   const int&, const int&, 
                                   const int&, const int&,
                                   const int&, const int&, 
                                   const int&, const int&,
                                   const int*, const double*, const double*,
                                   const double*, double* );
// in cartcoarsen3d.f:
   void cartwgtavgoutfacedoub3d0_( const int&, const int&, const int&,
                                   const int&, const int&, const int&,
                                   const int&, const int&, const int&,
                                   const int&, const int&, const int&,
                                   const int&, const int&, const int&,
                                   const int&, const int&, const int&,
                                   const int*, const double*, const double*,
                                   const double*, double* );
   void cartwgtavgoutfacedoub3d1_( const int&, const int&, const int&,
                                   const int&, const int&, const int&,
                                   const int&, const int&, const int&,
                                   const int&, const int&, const int&,
                                   const int&, const int&, const int&,
                                   const int&, const int&, const int&,
                                   const int*, const double*, const double*,
                                   const double*, double* );
   void cartwgtavgoutfacedoub3d2_( const int&, const int&, const int&,
                                   const int&, const int&, const int&,
                                   const int&, const int&, const int&,
                                   const int&, const int&, const int&,
                                   const int&, const int&, const int&,
                                   const int&, const int&, const int&,
                                   const int*, const double*, const double*,
                                   const double*, double* );
}

namespace SAMRAI {
    namespace geom {

template<int DIM> CartesianOuterfaceDoubleWeightedAverage<DIM>::CartesianOuterfaceDoubleWeightedAverage()
: xfer::CoarsenOperator<DIM>()
{
   d_name_id = "CONSERVATIVE_COARSEN";
}

template<int DIM> CartesianOuterfaceDoubleWeightedAverage<DIM>::~CartesianOuterfaceDoubleWeightedAverage()
{
}

template<int DIM> bool CartesianOuterfaceDoubleWeightedAverage<DIM>::findCoarsenOperator(
   const tbox::Pointer< hier::Variable<DIM> >& var,
   const std::string &op_name) const
{
   const tbox::Pointer< pdat::OuterfaceVariable<DIM,double> > cast_var(var);
   if ( !cast_var.isNull() && (op_name == d_name_id) ) {
      return(true);
   } else {
      return(false);
   }
}

template<int DIM> const std::string&
CartesianOuterfaceDoubleWeightedAverage<DIM>::getOperatorName() const
{
   return(d_name_id);
}

template<int DIM> int CartesianOuterfaceDoubleWeightedAverage<DIM>::getOperatorPriority() const
{
   return(0);
}

template<int DIM> hier::IntVector<DIM> 
CartesianOuterfaceDoubleWeightedAverage<DIM>::getStencilWidth() const {
   return(hier::IntVector<DIM>(0));
}

template<int DIM> void CartesianOuterfaceDoubleWeightedAverage<DIM>::coarsen(
   hier::Patch<DIM>& coarse, 
   const hier::Patch<DIM>& fine, 
   const int dst_component, 
   const int src_component, 
   const hier::Box<DIM>& coarse_box, 
   const hier::IntVector<DIM>& ratio) const 
{
   tbox::Pointer< pdat::OuterfaceData<DIM,double> > 
      fdata = fine.getPatchData(src_component);
   tbox::Pointer< pdat::OuterfaceData<DIM,double> > 
      cdata = coarse.getPatchData(dst_component);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!fdata.isNull());
   TBOX_ASSERT(!cdata.isNull());
   TBOX_ASSERT(cdata->getDepth() == fdata->getDepth());
#endif

   const hier::Index<DIM> filo = fdata->getGhostBox().lower();
   const hier::Index<DIM> fihi = fdata->getGhostBox().upper();
   const hier::Index<DIM> cilo = cdata->getGhostBox().lower();
   const hier::Index<DIM> cihi = cdata->getGhostBox().upper();

   const tbox::Pointer<CartesianPatchGeometry<DIM> > fgeom =
                                                    fine.getPatchGeometry();
   const tbox::Pointer<CartesianPatchGeometry<DIM> > cgeom =
                                                    coarse.getPatchGeometry();

   const hier::Index<DIM> ifirstc = coarse_box.lower();
   const hier::Index<DIM> ilastc = coarse_box.upper();

   for (int d = 0; d < cdata->getDepth(); d++) {
      // loop over lower and upper outerface arrays
      for (int i = 0; i < 2; i++) {
	 if (DIM == 1) {
	    cartwgtavgoutfacedoub1d_(ifirstc(0),ilastc(0),
				     filo(0),fihi(0),
				     cilo(0),cihi(0),
				     ratio,
				     fgeom->getDx(),
				     cgeom->getDx(),
				     fdata->getPointer(0,i,d),
				     cdata->getPointer(0,i,d));
	 } else if (DIM == 2) {
	    cartwgtavgoutfacedoub2d0_(ifirstc(0),ifirstc(1),ilastc(0),ilastc(1),
				      filo(0),filo(1),fihi(0),fihi(1),
				      cilo(0),cilo(1),cihi(0),cihi(1),
				      ratio,
				      fgeom->getDx(),
				      cgeom->getDx(),
				      fdata->getPointer(0,i,d),
				      cdata->getPointer(0,i,d));
	    cartwgtavgoutfacedoub2d1_(ifirstc(0),ifirstc(1),ilastc(0),ilastc(1),
				      filo(0),filo(1),fihi(0),fihi(1),
				      cilo(0),cilo(1),cihi(0),cihi(1),
				      ratio,
				      fgeom->getDx(),
				      cgeom->getDx(),
				      fdata->getPointer(1,i,d),
				      cdata->getPointer(1,i,d));
	 } else if (DIM == 3) {
	    cartwgtavgoutfacedoub3d0_(ifirstc(0),ifirstc(1),ifirstc(2),
				      ilastc(0),ilastc(1),ilastc(2),
				      filo(0),filo(1),filo(2),
				      fihi(0),fihi(1),fihi(2),
				      cilo(0),cilo(1),cilo(2),
				      cihi(0),cihi(1),cihi(2),
				      ratio,
				      fgeom->getDx(),
				      cgeom->getDx(),
				      fdata->getPointer(0,i,d),
				      cdata->getPointer(0,i,d));
	    cartwgtavgoutfacedoub3d1_(ifirstc(0),ifirstc(1),ifirstc(2),
				      ilastc(0),ilastc(1),ilastc(2),
				      filo(0),filo(1),filo(2),
				      fihi(0),fihi(1),fihi(2),
				      cilo(0),cilo(1),cilo(2),
				      cihi(0),cihi(1),cihi(2),
				      ratio,
				      fgeom->getDx(),
				      cgeom->getDx(),
				      fdata->getPointer(1,i,d),
				      cdata->getPointer(1,i,d));
	    cartwgtavgoutfacedoub3d2_(ifirstc(0),ifirstc(1),ifirstc(2),
				      ilastc(0),ilastc(1),ilastc(2),
				      filo(0),filo(1),filo(2),
				      fihi(0),fihi(1),fihi(2),
				      cilo(0),cilo(1),cilo(2),
				      cihi(0),cihi(1),cihi(2),
				      ratio,
				      fgeom->getDx(),
				      cgeom->getDx(),
				      fdata->getPointer(2,i,d),
				      cdata->getPointer(2,i,d));
	 } else {
	    TBOX_ERROR("CartesianOuterfaceDoubleWeightedAverage error...\n"
		       << "DIM > 3 not supported." << std::endl);
	 }
      }
   }
}

}
}
#endif
