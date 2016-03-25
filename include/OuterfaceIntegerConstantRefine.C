//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/operators/constant/outerface/OuterfaceIntegerConstantRefine.C $
// Package:	SAMRAI patchdata
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Constant refine operator for outerface integer data on 
//              a  mesh.
//

#ifndef included_pdat_OuterfaceIntegerConstantRefine_C
#define included_pdat_OuterfaceIntegerConstantRefine_C

#include "OuterfaceIntegerConstantRefine.h"

#include<float.h>
#include<math.h>
#include "tbox/Utilities.h"
#include "Index.h"
#include "OuterfaceData.h"
#include "OuterfaceVariable.h"

/*
*************************************************************************
*                                                                       *
* External declarations for FORTRAN  routines.                          *
*                                                                       *
*************************************************************************
*/
extern "C" {
// in conrefine1d.f:
   void conrefoutfaceintg1d_( const int&, const int&,
                              const int&, const int&,
                              const int&, const int&,
                              const int&, const int&,
                              const int*,
                              const int*, int* );
// in conrefine2d.f:
   void conrefoutfaceintg2d0_( const int&, const int&, 
                               const int&, const int&,
                               const int&, const int&,
                               const int&, const int&,
                               const int&, const int&,
                               const int&, const int&,
                               const int&, const int&,
                               const int&, const int&,
                               const int*,
                               const int*, int* );
   void conrefoutfaceintg2d1_( const int&, const int&,
                               const int&, const int&,
                               const int&, const int&,
                               const int&, const int&,
                               const int&, const int&,
                               const int&, const int&,
                               const int&, const int&,
                               const int&, const int&,
                               const int*,
                               const int*, int* );
// in conrefine3d.f:
   void conrefoutfaceintg3d0_( const int&, const int&, const int&,
                               const int&, const int&, const int&,
                               const int&, const int&, const int&,
                               const int&, const int&, const int&,
                               const int&, const int&, const int&,
                               const int&, const int&, const int&,
                               const int&, const int&, const int&,
                               const int&, const int&, const int&,
                               const int*,
                               const int*, int* );
   void conrefoutfaceintg3d1_( const int&, const int&, const int&,
                               const int&, const int&, const int&,
                               const int&, const int&, const int&,
                               const int&, const int&, const int&,
                               const int&, const int&, const int&,
                               const int&, const int&, const int&,
                               const int&, const int&, const int&,
                               const int&, const int&, const int&,
                               const int*,
                               const int*, int* );
   void conrefoutfaceintg3d2_( const int&, const int&, const int&,
                               const int&, const int&, const int&,
                               const int&, const int&, const int&,
                               const int&, const int&, const int&,
                               const int&, const int&, const int&,
                               const int&, const int&, const int&,
                               const int&, const int&, const int&,
                               const int&, const int&, const int&,
                               const int*,
                               const int*, int* );
}

namespace SAMRAI {
    namespace pdat {

template<int DIM> OuterfaceIntegerConstantRefine<DIM>::OuterfaceIntegerConstantRefine()
: xfer::RefineOperator<DIM>()
{
   d_name_id = "CONSTANT_REFINE";
}

template<int DIM> OuterfaceIntegerConstantRefine<DIM>::~OuterfaceIntegerConstantRefine()
{
}

template<int DIM> bool OuterfaceIntegerConstantRefine<DIM>::findRefineOperator(
   const tbox::Pointer< hier::Variable<DIM> >& var,
   const std::string &op_name) const
{
   const tbox::Pointer< OuterfaceVariable<DIM,int> > cast_var(var);
   if ( !cast_var.isNull() && (op_name == d_name_id) ) {
      return(true);
   } else {
      return(false);
   }
}

template<int DIM> const std::string&
OuterfaceIntegerConstantRefine<DIM>::getOperatorName() const
{
   return(d_name_id);
}

template<int DIM> int OuterfaceIntegerConstantRefine<DIM>::getOperatorPriority() const
{
   return(0);
}

template<int DIM> hier::IntVector<DIM> 
OuterfaceIntegerConstantRefine<DIM>::getStencilWidth() const {
   return(hier::IntVector<DIM>(0));
}

template<int DIM> void OuterfaceIntegerConstantRefine<DIM>::refine(
   hier::Patch<DIM>& fine, 
   const hier::Patch<DIM>& coarse, 
   const int dst_component, 
   const int src_component, 
   const hier::Box<DIM>& fine_box, 
   const hier::IntVector<DIM>& ratio) const
{
   tbox::Pointer< OuterfaceData<DIM,int> >
      cdata = coarse.getPatchData(src_component);
   tbox::Pointer< OuterfaceData<DIM,int> >
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

   const hier::Box<DIM> coarse_box = hier::Box<DIM>::coarsen(fine_box, ratio);
   const hier::Index<DIM> ifirstc = coarse_box.lower();
   const hier::Index<DIM> ilastc = coarse_box.upper();
   const hier::Index<DIM> ifirstf = fine_box.lower();
   const hier::Index<DIM> ilastf = fine_box.upper();

   for (int d = 0; d < fdata->getDepth(); d++) {
      // loop over lower and upper outerface arrays
      for (int i = 0; i < 2; i++) {
	 if (DIM == 1) {
	    conrefoutfaceintg1d_(ifirstc(0),ilastc(0),
				 ifirstf(0),ilastf(0),
				 cilo(0),cihi(0),
				 filo(0),fihi(0),
				 ratio,
				 cdata->getPointer(0,i,d),
				 fdata->getPointer(0,i,d));
	 } else if (DIM == 2) {
	    conrefoutfaceintg2d0_(ifirstc(0),ifirstc(1),ilastc(0),ilastc(1),
				  ifirstf(0),ifirstf(1),ilastf(0),ilastf(1),
				  cilo(0),cilo(1),cihi(0),cihi(1),
				  filo(0),filo(1),fihi(0),fihi(1),
				  ratio,
				  cdata->getPointer(0,i,d),
				  fdata->getPointer(0,i,d));
	    conrefoutfaceintg2d1_(ifirstc(0),ifirstc(1),ilastc(0),ilastc(1),
				  ifirstf(0),ifirstf(1),ilastf(0),ilastf(1),
				  cilo(0),cilo(1),cihi(0),cihi(1),
				  filo(0),filo(1),fihi(0),fihi(1),
				  ratio,
				  cdata->getPointer(1,i,d),
				  fdata->getPointer(1,i,d));
	 } else if (DIM == 3) {
	    conrefoutfaceintg3d0_(ifirstc(0),ifirstc(1),ifirstc(2),
				  ilastc(0),ilastc(1),ilastc(2),
				  ifirstf(0),ifirstf(1),ifirstf(2),
				  ilastf(0),ilastf(1),ilastf(2),
				  cilo(0),cilo(1),cilo(2),
				  cihi(0),cihi(1),cihi(2),
				  filo(0),filo(1),filo(2),
				  fihi(0),fihi(1),fihi(2),
				  ratio,
				  cdata->getPointer(0,i,d),
				  fdata->getPointer(0,i,d));
	    conrefoutfaceintg3d1_(ifirstc(0),ifirstc(1),ifirstc(2),
				  ilastc(0),ilastc(1),ilastc(2),
				  ifirstf(0),ifirstf(1),ifirstf(2),
				  ilastf(0),ilastf(1),ilastf(2),
				  cilo(0),cilo(1),cilo(2),
				  cihi(0),cihi(1),cihi(2),
				  filo(0),filo(1),filo(2),
				  fihi(0),fihi(1),fihi(2),
				  ratio,
				  cdata->getPointer(1,i,d),
				  fdata->getPointer(1,i,d));
	    conrefoutfaceintg3d2_(ifirstc(0),ifirstc(1),ifirstc(2),
				  ilastc(0),ilastc(1),ilastc(2),
				  ifirstf(0),ifirstf(1),ifirstf(2),
				  ilastf(0),ilastf(1),ilastf(2),
				  cilo(0),cilo(1),cilo(2),
				  cihi(0),cihi(1),cihi(2),
				  filo(0),filo(1),filo(2),
				  fihi(0),fihi(1),fihi(2),
				  ratio,
				  cdata->getPointer(2,i,d),
				  fdata->getPointer(2,i,d));
	 } else {
	    TBOX_ERROR("OuterfaceIntegerConstantRefine::refine DIM > 3 not supported" << std::endl);
	 }
      }
   }
}

}
}
#endif
