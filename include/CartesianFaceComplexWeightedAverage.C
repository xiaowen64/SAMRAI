//
// File:	CartesianFaceComplexWeightedAverage.C
// Package:	SAMRAI geometry
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Weighted averaging operator for face-centered complex data on 
//              a Cartesian mesh.
//

#ifndef included_geom_CartesianFaceComplexWeightedAverage_C
#define included_geom_CartesianFaceComplexWeightedAverage_C

#include "CartesianFaceComplexWeightedAverage.h"
#include "tbox/Complex.h"

#include<float.h>
#include<math.h>
#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif
#include "CartesianPatchGeometry.h"
#include "Index.h"
#include "FaceData.h"
#include "FaceVariable.h"
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
   void cartwgtavgfacecplx1d_( const int&, const int&,
                               const int&, const int&,
                               const int&, const int&,
                               const int*, const double*, const double*,
                               const dcomplex*, dcomplex* );
// in cartcoarsen2d.f:
   void cartwgtavgfacecplx2d0_( const int&, const int&, const int&, const int&,
                                const int&, const int&, const int&, const int&,
                                const int&, const int&, const int&, const int&,
                                const int*, const double*, const double*,
                                const dcomplex*, dcomplex* );

   void cartwgtavgfacecplx2d1_( const int&, const int&, const int&, const int&,
                                const int&, const int&, const int&, const int&,
                                const int&, const int&, const int&, const int&,
                                const int*, const double*, const double*,
                                const dcomplex*, dcomplex* );
// in cartcoarsen3d.f:
   void cartwgtavgfacecplx3d0_( const int&, const int&, const int&,
                                const int&, const int&, const int&,
                                const int&, const int&, const int&,
                                const int&, const int&, const int&,
                                const int&, const int&, const int&,
                                const int&, const int&, const int&,
                                const int*, const double*, const double*,
                                const dcomplex*, dcomplex* );
   void cartwgtavgfacecplx3d1_( const int&, const int&, const int&,
                                const int&, const int&, const int&,
                                const int&, const int&, const int&,
                                const int&, const int&, const int&,
                                const int&, const int&, const int&,
                                const int&, const int&, const int&,
                                const int*, const double*, const double*,
                                const dcomplex*, dcomplex* );
   void cartwgtavgfacecplx3d2_( const int&, const int&, const int&,
                                const int&, const int&, const int&,
                                const int&, const int&, const int&,
                                const int&, const int&, const int&,
                                const int&, const int&, const int&,
                                const int&, const int&, const int&,
                                const int*, const double*, const double*,
                                const dcomplex*, dcomplex* );
}

namespace SAMRAI {
    namespace geom {

template<int DIM> CartesianFaceComplexWeightedAverage<DIM>::CartesianFaceComplexWeightedAverage()
: xfer::CoarsenOperator<DIM>()
{
   d_name_id = "CONSERVATIVE_COARSEN";
}

template<int DIM> CartesianFaceComplexWeightedAverage<DIM>::~CartesianFaceComplexWeightedAverage()
{
}

template<int DIM> bool CartesianFaceComplexWeightedAverage<DIM>::findCoarsenOperator(
   const tbox::Pointer< hier::Variable<DIM> >& var,
   const string &op_name) const
{
   const tbox::Pointer< pdat::FaceVariable<DIM,dcomplex> > cast_var(var);
   if ( !cast_var.isNull() && (op_name == d_name_id) ) {
      return(true);
   } else {
      return(false);
   }
}

template<int DIM> const string&
CartesianFaceComplexWeightedAverage<DIM>::getOperatorName() const
{
   return(d_name_id);
}

template<int DIM> int CartesianFaceComplexWeightedAverage<DIM>::getOperatorPriority() const
{
   return(0);
}

template<int DIM> hier::IntVector<DIM> 
CartesianFaceComplexWeightedAverage<DIM>::getStencilWidth() const {
   return(hier::IntVector<DIM>(0));
}

template<int DIM> void CartesianFaceComplexWeightedAverage<DIM>::coarsen(
   hier::Patch<DIM>& coarse, 
   const hier::Patch<DIM>& fine, 
   const int dst_component, 
   const int src_component, 
   const hier::Box<DIM>& coarse_box, 
   const hier::IntVector<DIM>& ratio) const 
{
   tbox::Pointer< pdat::FaceData<DIM,dcomplex> > 
      fdata = fine.getPatchData(src_component);
   tbox::Pointer< pdat::FaceData<DIM,dcomplex> > 
      cdata = coarse.getPatchData(dst_component);
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!fdata.isNull());
   assert(!cdata.isNull());
   assert(cdata->getDepth() == fdata->getDepth());
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
	 cartwgtavgfacecplx1d_(ifirstc(0),ilastc(0),
			       filo(0),fihi(0),
			       cilo(0),cihi(0),
			       ratio,
			       fgeom->getDx(),
			       cgeom->getDx(),
			       fdata->getPointer(0,d),
			       cdata->getPointer(0,d));
      } else if (DIM == 2) {
	 cartwgtavgfacecplx2d0_(ifirstc(0),ifirstc(1),ilastc(0),ilastc(1),
				filo(0),filo(1),fihi(0),fihi(1),
				cilo(0),cilo(1),cihi(0),cihi(1),
				ratio,
				fgeom->getDx(),
				cgeom->getDx(),
				fdata->getPointer(0,d),
				cdata->getPointer(0,d));
	 cartwgtavgfacecplx2d1_(ifirstc(0),ifirstc(1),ilastc(0),ilastc(1),
				filo(0),filo(1),fihi(0),fihi(1),
				cilo(0),cilo(1),cihi(0),cihi(1),
				ratio,
				fgeom->getDx(),
				cgeom->getDx(),
				fdata->getPointer(1,d),
				cdata->getPointer(1,d));
      } else if (DIM == 3) {
	 cartwgtavgfacecplx3d0_(ifirstc(0),ifirstc(1),ifirstc(2),
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
	 cartwgtavgfacecplx3d1_(ifirstc(0),ifirstc(1),ifirstc(2),
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
	 cartwgtavgfacecplx3d2_(ifirstc(0),ifirstc(1),ifirstc(2),
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
      } else {
	 TBOX_ERROR("CartesianFaceComplexWeightedAverage error...\n"
		    << "DIM > 3 not supported." << endl);
      }

   }
}

}
}
#endif
