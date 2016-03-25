//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/geometry/cartesian/operators/side/CartesianSideFloatWeightedAverage.C $
// Package:	SAMRAI geometry
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description: Weighted averaging operator for side-centered float data on 
//              a Cartesian mesh.
//

#ifndef included_geom_CartesianSideFloatWeightedAverage_C
#define included_geom_CartesianSideFloatWeightedAverage_C

#include "CartesianSideFloatWeightedAverage.h"

#include<float.h>
#include<math.h>
#include "CartesianPatchGeometry.h"
#include "Index.h"
#include "SideData.h"
#include "SideVariable.h"
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
   void cartwgtavgsideflot1d_( const int&, const int&,
                               const int&, const int&,
                               const int&, const int&,
                               const int*, const double*, const double*,
                               const float*, float* );
// in cartcoarsen2d.f:
   void cartwgtavgsideflot2d0_( const int&, const int&, const int&, const int&,
                                const int&, const int&, const int&, const int&,
                                const int&, const int&, const int&, const int&,
                                const int*, const double*, const double*,
                                const float*, float* );

   void cartwgtavgsideflot2d1_( const int&, const int&, const int&, const int&,
                                const int&, const int&, const int&, const int&,
                                const int&, const int&, const int&, const int&,
                                const int*, const double*, const double*,
                                const float*, float* );
// in cartcoarsen3d.f:
   void cartwgtavgsideflot3d0_( const int&, const int&, const int&,
                                const int&, const int&, const int&,
                                const int&, const int&, const int&,
                                const int&, const int&, const int&,
                                const int&, const int&, const int&,
                                const int&, const int&, const int&,
                                const int*, const double*, const double*,
                                const float*, float* );
   void cartwgtavgsideflot3d1_( const int&, const int&, const int&,
                                const int&, const int&, const int&,
                                const int&, const int&, const int&,
                                const int&, const int&, const int&,
                                const int&, const int&, const int&,
                                const int&, const int&, const int&,
                                const int*, const double*, const double*,
                                const float*, float* );
   void cartwgtavgsideflot3d2_( const int&, const int&, const int&,
                                const int&, const int&, const int&,
                                const int&, const int&, const int&,
                                const int&, const int&, const int&,
                                const int&, const int&, const int&,
                                const int&, const int&, const int&,
                                const int*, const double*, const double*,
                                const float*, float* );
}

namespace SAMRAI {
    namespace geom {

template<int DIM> CartesianSideFloatWeightedAverage<DIM>::CartesianSideFloatWeightedAverage()
: xfer::CoarsenOperator<DIM>()
{
   d_name_id = "CONSERVATIVE_COARSEN";
}

template<int DIM> CartesianSideFloatWeightedAverage<DIM>::~CartesianSideFloatWeightedAverage()
{
}

template<int DIM> bool CartesianSideFloatWeightedAverage<DIM>::findCoarsenOperator(
   const tbox::Pointer< hier::Variable<DIM> >& var,
   const std::string &op_name) const
{
   const tbox::Pointer< pdat::SideVariable<DIM,float> > cast_var(var);
   if ( !cast_var.isNull() && (op_name == d_name_id) ) {
      return(true);
   } else {
      return(false);
   }
}

template<int DIM> const std::string&
CartesianSideFloatWeightedAverage<DIM>::getOperatorName() const
{
   return(d_name_id);
}

template<int DIM> int CartesianSideFloatWeightedAverage<DIM>::getOperatorPriority() const
{
   return(0);
}

template<int DIM> hier::IntVector<DIM> 
CartesianSideFloatWeightedAverage<DIM>::getStencilWidth() const {
   return(hier::IntVector<DIM>(0));
}

template<int DIM> void CartesianSideFloatWeightedAverage<DIM>::coarsen(
   hier::Patch<DIM>& coarse, 
   const hier::Patch<DIM>& fine, 
   const int dst_component, 
   const int src_component, 
   const hier::Box<DIM>& coarse_box, 
   const hier::IntVector<DIM>& ratio) const 
{
   tbox::Pointer< pdat::SideData<DIM,float> > 
      fdata = fine.getPatchData(src_component);
   tbox::Pointer< pdat::SideData<DIM,float> > 
      cdata = coarse.getPatchData(dst_component);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!fdata.isNull());
   TBOX_ASSERT(!cdata.isNull());
   TBOX_ASSERT(cdata->getDepth() == fdata->getDepth());
#endif
   const hier::IntVector<DIM>& directions = cdata->getDirectionVector();
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(directions ==
          hier::IntVector<DIM>::min(directions, fdata->getDirectionVector())); 
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
      if (DIM == 1) {
	 if (directions(0)) {
	    cartwgtavgsideflot1d_(ifirstc(0),ilastc(0),
				  filo(0),fihi(0),
				  cilo(0),cihi(0),
				  ratio,
				  fgeom->getDx(),
				  cgeom->getDx(),
                               fdata->getPointer(0,d),
				  cdata->getPointer(0,d));
	 }
      } else if (DIM == 2) {
	 if (directions(0)) {
	    cartwgtavgsideflot2d0_(ifirstc(0),ifirstc(1),ilastc(0),ilastc(1),
				   filo(0),filo(1),fihi(0),fihi(1),
				   cilo(0),cilo(1),cihi(0),cihi(1),
				   ratio,
				   fgeom->getDx(),
				   cgeom->getDx(),
				   fdata->getPointer(0,d),
                                cdata->getPointer(0,d));
	 }
	 if (directions(1)) {
	    cartwgtavgsideflot2d1_(ifirstc(0),ifirstc(1),ilastc(0),ilastc(1),
				   filo(0),filo(1),fihi(0),fihi(1),
				   cilo(0),cilo(1),cihi(0),cihi(1),
				   ratio,
				   fgeom->getDx(),
				   cgeom->getDx(),
				   fdata->getPointer(1,d),
				   cdata->getPointer(1,d));
	 }
      } else if (DIM == 3) {
	 if (directions(0)) {
	    cartwgtavgsideflot3d0_(ifirstc(0),ifirstc(1),ifirstc(2),
				   ilastc(0),ilastc(1),ilastc(2),
				   filo(0),filo(1),filo(2),
				   fihi(0),fihi(1),fihi(2),
				   cilo(0),cilo(1),cilo(2),
				   cihi(0),cihi(1),cihi(2),
				   ratio,
				   fgeom->getDx(),
				   cgeom->getDx(),
				   fdata->getPointer(0,d),
				   cdata->getPointer(0,d));
	 }
	 if (directions(1)) {
	    cartwgtavgsideflot3d1_(ifirstc(0),ifirstc(1),ifirstc(2),
				   ilastc(0),ilastc(1),ilastc(2),
				   filo(0),filo(1),filo(2),
				   fihi(0),fihi(1),fihi(2),
				   cilo(0),cilo(1),cilo(2),
				   cihi(0),cihi(1),cihi(2),
				   ratio,
				   fgeom->getDx(),
				   cgeom->getDx(),
				   fdata->getPointer(1,d),
				   cdata->getPointer(1,d));
	 }
	 if (directions(2)) {
	    cartwgtavgsideflot3d2_(ifirstc(0),ifirstc(1),ifirstc(2),
				   ilastc(0),ilastc(1),ilastc(2),
				   filo(0),filo(1),filo(2),
				   fihi(0),fihi(1),fihi(2),
				   cilo(0),cilo(1),cilo(2),
				   cihi(0),cihi(1),cihi(2),
				   ratio,
				   fgeom->getDx(),
				   cgeom->getDx(),
				   fdata->getPointer(2,d),
				   cdata->getPointer(2,d));
	 }
      } else {
	 TBOX_ERROR("CartesianSideFloatWeightedAverage error...\n"
		    << "DIM > 3 not supported." << std::endl);
      }
   }
}

}
}
#endif
