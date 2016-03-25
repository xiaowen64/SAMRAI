//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/geometry/cartesian/operators/cell/CartesianCellDoubleLinearRefine.C $
// Package:	SAMRAI geometry
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Linear refine operator for cell-centered double data on 
//              a Cartesian mesh.
//

#ifndef included_geom_CartesianCellDoubleLinearRefine_C
#define included_geom_CartesianCellDoubleLinearRefine_C

#include "CartesianCellDoubleLinearRefine.h"

#include<float.h>
#include<math.h>
#include "CartesianPatchGeometry.h"
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
// in cartrefine1d.f:
   void cartlinrefcelldoub1d_( const int&, const int&,
                               const int&, const int&,
                               const int&, const int&,
                               const int&, const int&,
                               const int*, const double*, const double*,
                               const double*, double* );
// in cartrefine2d.f:
   void cartlinrefcelldoub2d_( const int&, const int&, const int&, const int&,
                               const int&, const int&, const int&, const int&,
                               const int&, const int&, const int&, const int&,
                               const int&, const int&, const int&, const int&,
                               const int*, const double*, const double*,
                               const double*, double* );
// in cartrefine3d.f:
   void cartlinrefcelldoub3d_( const int&, const int&, const int&,
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


template<int DIM> CartesianCellDoubleLinearRefine<DIM>::CartesianCellDoubleLinearRefine()
: xfer::RefineOperator<DIM>()
{
   d_name_id = "LINEAR_REFINE";
}

template<int DIM> CartesianCellDoubleLinearRefine<DIM>::~CartesianCellDoubleLinearRefine()
{
}

template<int DIM> bool CartesianCellDoubleLinearRefine<DIM>::findRefineOperator(
   const tbox::Pointer< hier::Variable<DIM> >& var,
   const std::string &op_name) const
{
   const tbox::Pointer< pdat::CellVariable<DIM,double> > cast_var(var);
   if ( !cast_var.isNull() && (op_name == d_name_id) ) {
      return(true);
   } else {
      return(false);
   }
}

template<int DIM> const std::string&
CartesianCellDoubleLinearRefine<DIM>::getOperatorName() const
{
   return(d_name_id);
}

template<int DIM> int CartesianCellDoubleLinearRefine<DIM>::getOperatorPriority() const
{
   return(0);
}

template<int DIM> hier::IntVector<DIM> 
CartesianCellDoubleLinearRefine<DIM>::getStencilWidth() const {
   return(hier::IntVector<DIM>(1));
}

template<int DIM> void CartesianCellDoubleLinearRefine<DIM>::refine(
   hier::Patch<DIM>& fine, 
   const hier::Patch<DIM>& coarse, 
   const int dst_component, 
   const int src_component, 
   const hier::Box<DIM>& fine_box, 
   const hier::IntVector<DIM>& ratio) const
{
   tbox::Pointer< pdat::CellData<DIM,double> >
      cdata = coarse.getPatchData(src_component);
   tbox::Pointer< pdat::CellData<DIM,double> >
      fdata = fine.getPatchData(dst_component);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!cdata.isNull());
   TBOX_ASSERT(!fdata.isNull());
   TBOX_ASSERT(cdata->getDepth() == fdata->getDepth());
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
	 cartlinrefcelldoub1d_(ifirstc(0),ilastc(0),
			       ifirstf(0),ilastf(0),
			       cilo(0),cihi(0),
			       filo(0),fihi(0),
			       ratio,
			       cgeom->getDx(),
			       fgeom->getDx(),
			       cdata->getPointer(d),
			       fdata->getPointer(d));
      } else if (DIM == 2) {
	 cartlinrefcelldoub2d_(ifirstc(0),ifirstc(1),ilastc(0),ilastc(1),
			       ifirstf(0),ifirstf(1),ilastf(0),ilastf(1),
			       cilo(0),cilo(1),cihi(0),cihi(1),
			       filo(0),filo(1),fihi(0),fihi(1),
			       ratio,
			       cgeom->getDx(),
			       fgeom->getDx(),
			       cdata->getPointer(d),
			       fdata->getPointer(d));
      } else if (DIM == 3) {
	 cartlinrefcelldoub3d_(ifirstc(0),ifirstc(1),ifirstc(2),
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
	 TBOX_ERROR("CartesianCellDoubleLinearRefine error...\n"
		    << "DIM > 3 not supported." << std::endl);

      }
   }
}

}
}
#endif
