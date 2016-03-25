//
// File:	CartesianEdgeDoubleConservativeLinearRefine.C
// Package:	SAMRAI geometry
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Conservative linear refine operator for edge-centered 
//              double data on a Cartesian mesh.
//

#ifndef included_geom_CartesianEdgeDoubleConservativeLinearRefine_C
#define included_geom_CartesianEdgeDoubleConservativeLinearRefine_C

#include "CartesianEdgeDoubleConservativeLinearRefine.h"
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
#include "EdgeData.h"
#include "EdgeVariable.h"
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
   void cartclinrefedgedoub1d_( const int&, const int&,
                                const int&, const int&,
                                const int&, const int&,
                                const int&, const int&,
                                const int*, const double*, const double*,
                                const double*, double*,
                                double*, double* );
// in cartrefine2d.f:
   void cartclinrefedgedoub2d0_( const int&, const int&, const int&, const int&,
                                 const int&, const int&, const int&, const int&,
                                 const int&, const int&, const int&, const int&,
                                 const int&, const int&, const int&, const int&,
                                 const int*, const double*, const double*,
                                 const double*, double*,
                                 double*, double*, double*, double* );
   void cartclinrefedgedoub2d1_( const int&, const int&, const int&, const int&,
                                 const int&, const int&, const int&, const int&,
                                 const int&, const int&, const int&, const int&,
                                 const int&, const int&, const int&, const int&,
                                 const int*, const double*, const double*,
                                 const double*, double*,
                                 double*, double*, double*, double* );
// in cartrefine3d.f:
   void cartclinrefedgedoub3d0_( const int&, const int&, const int&,
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
   void cartclinrefedgedoub3d1_( const int&, const int&, const int&,
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
   void cartclinrefedgedoub3d2_( const int&, const int&, const int&,
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

namespace SAMRAI {
    namespace geom {

template<int DIM> CartesianEdgeDoubleConservativeLinearRefine<DIM>::CartesianEdgeDoubleConservativeLinearRefine()
: xfer::RefineOperator<DIM>()
{
   d_name_id = "CONSERVATIVE_LINEAR_REFINE";
}

template<int DIM> CartesianEdgeDoubleConservativeLinearRefine<DIM>::~CartesianEdgeDoubleConservativeLinearRefine()
{
}

template<int DIM> bool CartesianEdgeDoubleConservativeLinearRefine<DIM>::findRefineOperator(
   const tbox::Pointer< hier::Variable<DIM> >& var,
   const string &op_name) const
{
   const tbox::Pointer< pdat::EdgeVariable<DIM,double> > cast_var(var);
   if ( !cast_var.isNull() && (op_name == d_name_id) ) {
      return(true);
   } else {
      return(false);
   }
}

template<int DIM> const string&
CartesianEdgeDoubleConservativeLinearRefine<DIM>::getOperatorName() const
{
   return(d_name_id);
}

template<int DIM> int 
CartesianEdgeDoubleConservativeLinearRefine<DIM>::getOperatorPriority() const
{
   return(0);
}

template<int DIM> hier::IntVector<DIM> 
CartesianEdgeDoubleConservativeLinearRefine<DIM>::getStencilWidth() const {
   return(hier::IntVector<DIM>(1));
}

template<int DIM> void CartesianEdgeDoubleConservativeLinearRefine<DIM>::refine(
   hier::Patch<DIM>& fine, 
   const hier::Patch<DIM>& coarse, 
   const int dst_component, 
   const int src_component, 
   const hier::Box<DIM>& fine_box, 
   const hier::IntVector<DIM>& ratio) const
{
   tbox::Pointer< pdat::EdgeData<DIM,double> >
      cdata = coarse.getPatchData(src_component);
   tbox::Pointer< pdat::EdgeData<DIM,double> >
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

   tbox::Pointer<tbox::Arena> 
   memory_arena = tbox::ArenaManager::getManager()->getScratchAllocator();

   const hier::IntVector<DIM> tmp_ghosts(0);
   tbox::Array<double> diff0(cgbox.numberCells(0)+2, memory_arena);
   pdat::EdgeData<DIM,double> slope0(cgbox, 1, tmp_ghosts, memory_arena);

   for (int d = 0; d < fdata->getDepth(); d++) {
      if (DIM == 1) {
	 cartclinrefedgedoub1d_(ifirstc(0),ilastc(0),
				ifirstf(0),ilastf(0),
				cilo(0),cihi(0),
				filo(0),fihi(0),
				ratio,
				cgeom->getDx(),
				fgeom->getDx(),
				cdata->getPointer(0,d),
				fdata->getPointer(0,d),
				diff0.getPointer(),slope0.getPointer(0));
      } else if (DIM == 2) {
	 tbox::Array<double> diff1(cgbox.numberCells(1)+2, memory_arena);
	 pdat::EdgeData<DIM,double> slope1(cgbox, 1, tmp_ghosts, memory_arena);

	 cartclinrefedgedoub2d0_(ifirstc(0),ifirstc(1),ilastc(0),ilastc(1),
				 ifirstf(0),ifirstf(1),ilastf(0),ilastf(1),
				 cilo(0),cilo(1),cihi(0),cihi(1),
				 filo(0),filo(1),fihi(0),fihi(1),
				 ratio,
				 cgeom->getDx(),
				 fgeom->getDx(),
				 cdata->getPointer(0,d),
				 fdata->getPointer(0,d),
				 diff0.getPointer(),slope0.getPointer(0),
				 diff1.getPointer(),slope1.getPointer(0));
	 cartclinrefedgedoub2d1_(ifirstc(0),ifirstc(1),ilastc(0),ilastc(1),
				 ifirstf(0),ifirstf(1),ilastf(0),ilastf(1),
				 cilo(0),cilo(1),cihi(0),cihi(1),
				 filo(0),filo(1),fihi(0),fihi(1),
				 ratio,
				 cgeom->getDx(),
				 fgeom->getDx(),
				 cdata->getPointer(1,d),
				 fdata->getPointer(1,d),
				 diff1.getPointer(),slope1.getPointer(1),
				 diff0.getPointer(),slope0.getPointer(1));
      } else if (DIM == 3) {
	 tbox::Array<double> diff1(cgbox.numberCells(1)+2, memory_arena);
	 pdat::EdgeData<DIM,double> slope1(cgbox, 1, tmp_ghosts, memory_arena);

	 tbox::Array<double> diff2(cgbox.numberCells(2)+2, memory_arena);
	 pdat::EdgeData<DIM,double> slope2(cgbox, 1, tmp_ghosts, memory_arena);


	 cartclinrefedgedoub3d0_(ifirstc(0),ifirstc(1),ifirstc(2),
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
				 cdata->getPointer(0,d),
				 fdata->getPointer(0,d),
				 diff0.getPointer(),slope0.getPointer(0),
				 diff1.getPointer(),slope1.getPointer(0),
				 diff2.getPointer(),slope2.getPointer(0));
	 cartclinrefedgedoub3d1_(ifirstc(0),ifirstc(1),ifirstc(2),
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
				 cdata->getPointer(1,d),
				 fdata->getPointer(1,d),
				 diff1.getPointer(),slope1.getPointer(1),
				 diff2.getPointer(),slope2.getPointer(1),
				 diff0.getPointer(),slope0.getPointer(1));
	 cartclinrefedgedoub3d2_(ifirstc(0),ifirstc(1),ifirstc(2),
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
				 cdata->getPointer(2,d),
				 fdata->getPointer(2,d),
				 diff2.getPointer(),slope2.getPointer(2),
				 diff0.getPointer(),slope0.getPointer(2),
				 diff1.getPointer(),slope1.getPointer(2));
      } else {
	 TBOX_ERROR("CartesianEdgeDoubleConservativeLinearRefine error...\n"
		    << "DIM > 3 not supported." << endl);
      }

   }
}

}
}
#endif
