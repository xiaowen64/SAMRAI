//
// File:	CartesianCellFloatConservativeLinearRefine.C
// Package:	SAMRAI geometry
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Conservative linear refine operator for cell-centered float 
//              data on a Cartesian mesh.
//

#ifndef included_geom_CartesianCellFloatConservativeLinearRefine_C
#define included_geom_CartesianCellFloatConservativeLinearRefine_C

#include "CartesianCellFloatConservativeLinearRefine.h"
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
   void cartclinrefcellflot1d_( const int&, const int&,
                                const int&, const int&,
                                const int&, const int&,
                                const int&, const int&,
                                const int*, const double*, const double*,
                                const float*, float*,
                                float*, float* );
// in cartrefine2d.f:
   void cartclinrefcellflot2d_( const int&, const int&, const int&, const int&,
                                const int&, const int&, const int&, const int&,
                                const int&, const int&, const int&, const int&,
                                const int&, const int&, const int&, const int&,
                                const int*, const double*, const double*,
                                const float*, float*,
                                float*, float*, float*, float* );
// in cartrefine3d.f:
   void cartclinrefcellflot3d_( const int&, const int&, const int&,
                                const int&, const int&, const int&,
                                const int&, const int&, const int&,
                                const int&, const int&, const int&,
                                const int&, const int&, const int&,
                                const int&, const int&, const int&,
                                const int&, const int&, const int&,
                                const int&, const int&, const int&,
                                const int*, const double*, const double*,
                                const float*, float*,
                                float*, float*, float*,
                                float*, float*, float* );
}

namespace SAMRAI {
    namespace geom {

template<int DIM> CartesianCellFloatConservativeLinearRefine<DIM>::CartesianCellFloatConservativeLinearRefine()
: xfer::RefineOperator<DIM>()
{
   d_name_id = "CONSERVATIVE_LINEAR_REFINE";
}

template<int DIM> CartesianCellFloatConservativeLinearRefine<DIM>::~CartesianCellFloatConservativeLinearRefine()
{
}

template<int DIM> bool CartesianCellFloatConservativeLinearRefine<DIM>::findRefineOperator(
   const tbox::Pointer< hier::Variable<DIM> >& var,
   const string &op_name) const
{
   const tbox::Pointer< pdat::CellVariable<DIM,float> > cast_var(var);
   if ( !cast_var.isNull() && (op_name == d_name_id) ) {
      return(true);
   } else {
      return(false);
   }
}

template<int DIM> const string&
CartesianCellFloatConservativeLinearRefine<DIM>::getOperatorName() const
{
   return(d_name_id);
}

template<int DIM> int 
CartesianCellFloatConservativeLinearRefine<DIM>::getOperatorPriority() const
{
   return(0);
}

template<int DIM> hier::IntVector<DIM> 
CartesianCellFloatConservativeLinearRefine<DIM>::getStencilWidth() const {
   return(hier::IntVector<DIM>(1));
}

template<int DIM> void CartesianCellFloatConservativeLinearRefine<DIM>::refine(
   hier::Patch<DIM>& fine, 
   const hier::Patch<DIM>& coarse, 
   const int dst_component, 
   const int src_component, 
   const hier::Box<DIM>& fine_box, 
   const hier::IntVector<DIM>& ratio) const
{
   tbox::Pointer< pdat::CellData<DIM,float> >
      cdata = coarse.getPatchData(src_component);
   tbox::Pointer< pdat::CellData<DIM,float> >
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
   tbox::Array<float> diff0(cgbox.numberCells(0)+1, memory_arena);
   pdat::CellData<DIM,float> slope0(cgbox, 1, tmp_ghosts, memory_arena);

   for (int d = 0; d < fdata->getDepth(); d++) {
      if (DIM == 1) {
	 cartclinrefcellflot1d_(ifirstc(0),ilastc(0),
				ifirstf(0),ilastf(0),
				cilo(0),cihi(0),
				filo(0),fihi(0),
				ratio,
				cgeom->getDx(),
				fgeom->getDx(),
				cdata->getPointer(d),
				fdata->getPointer(d),
				diff0.getPointer(),slope0.getPointer()); 
      } else if (DIM == 2) {
	 tbox::Array<float> diff1(cgbox.numberCells(1)+1, memory_arena);
	 pdat::CellData<DIM,float> slope1(cgbox, 1, tmp_ghosts, memory_arena);

	 cartclinrefcellflot2d_(ifirstc(0),ifirstc(1),ilastc(0),ilastc(1),
				ifirstf(0),ifirstf(1),ilastf(0),ilastf(1),
				cilo(0),cilo(1),cihi(0),cihi(1),
				filo(0),filo(1),fihi(0),fihi(1),
				ratio,
				cgeom->getDx(),
				fgeom->getDx(),
				cdata->getPointer(d),
				fdata->getPointer(d),
				diff0.getPointer(),slope0.getPointer(),
				diff1.getPointer(),slope1.getPointer()); 
      } else if (DIM == 3) {
	 tbox::Array<float> diff1(cgbox.numberCells(1)+1, memory_arena);
	 pdat::CellData<DIM,float> slope1(cgbox, 1, tmp_ghosts, memory_arena);
	 
	 tbox::Array<float> diff2(cgbox.numberCells(2)+1, memory_arena);
	 pdat::CellData<DIM,float> slope2(cgbox, 1, tmp_ghosts, memory_arena);

	 cartclinrefcellflot3d_(ifirstc(0),ifirstc(1),ifirstc(2),
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
				fdata->getPointer(d),
				diff0.getPointer(),slope0.getPointer(),
				diff1.getPointer(),slope1.getPointer(),
				diff2.getPointer(),slope2.getPointer());
      } else {
	 TBOX_ERROR("CartesianFloatConservativeLinearRefine error...\n"
		    << "DIM > 3 not supported." << endl);

      }
   }
}

}
}
#endif
