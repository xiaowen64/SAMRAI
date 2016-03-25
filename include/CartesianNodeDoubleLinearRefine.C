//
// File:	CartesianNodeDoubleLinearRefine.C
// Package:	SAMRAI geometry
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Linear refine operator for node-centered double data on 
//              a Cartesian mesh.
//

#ifndef included_geom_CartesianNodeDoubleLinearRefine_C
#define included_geom_CartesianNodeDoubleLinearRefine_C

#include "CartesianNodeDoubleLinearRefine.h"

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
#include "NodeData.h"
#include "NodeVariable.h"
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
   void cartlinrefnodedoub1d_( const int&, const int&,
                               const int&, const int&,
                               const int&, const int&,
                               const int&, const int&,
                               const int*, const double*, const double*,
                               const double*, double* );
// in cartrefine2d.f:
   void cartlinrefnodedoub2d_( const int&, const int&, const int&, const int&,
                               const int&, const int&, const int&, const int&,
                               const int&, const int&, const int&, const int&,
                               const int&, const int&, const int&, const int&,
                               const int*, const double*, const double*,
                               const double*, double* );
// in cartrefine3d.f:
   void cartlinrefnodedoub3d_( const int&, const int&, const int&,
                               const int&, const int&, const int&,
                               const int&, const int&, const int&,
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

template<int DIM> CartesianNodeDoubleLinearRefine<DIM>::CartesianNodeDoubleLinearRefine()
: xfer::RefineOperator<DIM>()
{
   d_name_id = "LINEAR_REFINE";
}

template<int DIM> CartesianNodeDoubleLinearRefine<DIM>::~CartesianNodeDoubleLinearRefine()
{
}

template<int DIM> bool CartesianNodeDoubleLinearRefine<DIM>::findRefineOperator(
   const tbox::Pointer< hier::Variable<DIM> >& var,
   const string &op_name) const
{
   const tbox::Pointer< pdat::NodeVariable<DIM,double> > cast_var(var);
   if ( !cast_var.isNull() && (op_name == d_name_id) ) {
      return(true);
   } else {
      return(false);
   }
}

template<int DIM> const string&
CartesianNodeDoubleLinearRefine<DIM>::getOperatorName() const
{
   return(d_name_id);
}

template<int DIM> int CartesianNodeDoubleLinearRefine<DIM>::getOperatorPriority() const
{
   return(0);
}

template<int DIM> hier::IntVector<DIM> 
CartesianNodeDoubleLinearRefine<DIM>::getStencilWidth() const {
   return(hier::IntVector<DIM>(0));
}

template<int DIM> void CartesianNodeDoubleLinearRefine<DIM>::refine(
   hier::Patch<DIM>& fine, 
   const hier::Patch<DIM>& coarse, 
   const int dst_component, 
   const int src_component, 
   const hier::Box<DIM>& fine_box, 
   const hier::IntVector<DIM>& ratio) const
{
   tbox::Pointer< pdat::NodeData<DIM,double> >
      cdata = coarse.getPatchData(src_component);
   tbox::Pointer< pdat::NodeData<DIM,double> >
      fdata = fine.getPatchData(dst_component);
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!cdata.isNull());
   assert(!fdata.isNull());
   assert(cdata->getDepth() == fdata->getDepth());
#endif

   const hier::Box<DIM> cgbox(cdata->getGhostBox());

   const hier::Index<DIM> cilo = cgbox.lower();
   const hier::Index<DIM> cihi = cgbox.upper();
   const hier::Index<DIM> filo = fdata->getGhostBox().lower();
   const hier::Index<DIM> fihi = fdata->getGhostBox().upper();

   const tbox::Pointer<CartesianPatchGeometry<DIM> > cgeom =
                                                    coarse.getPatchGeometry();
   const tbox::Pointer<CartesianPatchGeometry<DIM> > fgeom =
                                                    fine.getPatchGeometry();

   const hier::Box<DIM> coarse_box = hier::Box<DIM>::coarsen(fine_box, ratio);
   const hier::Index<DIM> ifirstc = coarse_box.lower();
   const hier::Index<DIM> ilastc = coarse_box.upper();
   const hier::Index<DIM> ifirstf = fine_box.lower();
   const hier::Index<DIM> ilastf = fine_box.upper();

   for (int d = 0; d < fdata->getDepth(); d++) {
      if (DIM == 1) {
	 cartlinrefnodedoub1d_(ifirstc(0),ilastc(0),
			       ifirstf(0),ilastf(0),
			       cilo(0),cihi(0),
			       filo(0),fihi(0),
			       ratio,
			       cgeom->getDx(),
			       fgeom->getDx(),
			       cdata->getPointer(d),
			       fdata->getPointer(d));
      } else if (DIM == 2) {
	 cartlinrefnodedoub2d_(ifirstc(0),ifirstc(1),ilastc(0),ilastc(1),
			       ifirstf(0),ifirstf(1),ilastf(0),ilastf(1),
			       cilo(0),cilo(1),cihi(0),cihi(1),
			       filo(0),filo(1),fihi(0),fihi(1),
			       ratio,
			       cgeom->getDx(),
			       fgeom->getDx(),
			       cdata->getPointer(d),
			       fdata->getPointer(d));
      } else if (DIM == 3) {
	 cartlinrefnodedoub3d_(ifirstc(0),ifirstc(1),ifirstc(2),
			       ilastc(0),ilastc(1),ilastc(2),
			       ifirstf(0),ifirstf(1),ifirstf(2),
			       ilastf(0),ilastf(1),ilastf(2),
			       cilo(0),cilo(1),cilo(2),
			       cihi(0),cihi(1),cihi(2),
			       filo(0),filo(1),filo(2),
			       fihi(0),fihi(1),fihi(2),
			       ratio,
			       cgeom->getDx(),
			       fgeom->getDx(),
			       cdata->getPointer(d),
			       fdata->getPointer(d));
      } else {
	 TBOX_ERROR("CartesianNodeDoubleLinearRefine error...\n"
		    << "DIM > 3 not supported." << endl);
      }
   }
}

}
}
#endif
