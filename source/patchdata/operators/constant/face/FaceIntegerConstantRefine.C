//
// File:	FaceIntegerConstantRefine.C
// Package:	SAMRAI patchdata
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Constant refine operator for face-centered integer data on 
//              a  mesh.
//

#ifndef included_pdat_FaceIntegerConstantRefine_C
#define included_pdat_FaceIntegerConstantRefine_C

#include "FaceIntegerConstantRefine.h"

#include<float.h>
#include<math.h>
#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif
#include "tbox/Utilities.h"
#include "Index.h"
#include "FaceData.h"
#include "FaceVariable.h"

/*
*************************************************************************
*                                                                       *
* External declarations for FORTRAN  routines.                          *
*                                                                       *
*************************************************************************
*/
extern "C" {
// in cartrefine1d.f:
   void conreffaceintg1d_( const int&, const int&,
                           const int&, const int&,
                           const int&, const int&,
                           const int&, const int&,
                           const int*,
                           const int*, int* );
// in cartrefine2d.f:
   void conreffaceintg2d0_( const int&, const int&, const int&, const int&,
                            const int&, const int&, const int&, const int&,
                            const int&, const int&, const int&, const int&,
                            const int&, const int&, const int&, const int&,
                            const int*,
                            const int*, int* );
   void conreffaceintg2d1_( const int&, const int&, const int&, const int&,
                            const int&, const int&, const int&, const int&,
                            const int&, const int&, const int&, const int&,
                            const int&, const int&, const int&, const int&,
                            const int*,
                            const int*, int* );
// in cartrefine3d.f:
   void conreffaceintg3d0_( const int&, const int&, const int&,
                            const int&, const int&, const int&,
                            const int&, const int&, const int&,
                            const int&, const int&, const int&,
                            const int&, const int&, const int&,
                            const int&, const int&, const int&,
                            const int&, const int&, const int&,
                            const int&, const int&, const int&,
                            const int*,
                            const int*, int* );
   void conreffaceintg3d1_( const int&, const int&, const int&,
                            const int&, const int&, const int&,
                            const int&, const int&, const int&,
                            const int&, const int&, const int&,
                            const int&, const int&, const int&,
                            const int&, const int&, const int&,
                            const int&, const int&, const int&,
                            const int&, const int&, const int&,
                            const int*,
                            const int*, int* );
   void conreffaceintg3d2_( const int&, const int&, const int&,
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

template<int DIM> FaceIntegerConstantRefine<DIM>::FaceIntegerConstantRefine()
: xfer::RefineOperator<DIM>()
{
   d_name_id = "CONSTANT_REFINE";
}

template<int DIM> FaceIntegerConstantRefine<DIM>::~FaceIntegerConstantRefine()
{
}

template<int DIM> bool FaceIntegerConstantRefine<DIM>::findRefineOperator(
   const tbox::Pointer< hier::Variable<DIM> >& var,
   const string &op_name) const
{
   const tbox::Pointer< FaceVariable<DIM,int> > cast_var(var);
   if ( !cast_var.isNull() && (op_name == d_name_id) ) {
      return(true);
   } else {
      return(false);
   }
}

template<int DIM> const string&
FaceIntegerConstantRefine<DIM>::getOperatorName() const
{
   return(d_name_id);
}

template<int DIM> int FaceIntegerConstantRefine<DIM>::getOperatorPriority() const
{
   return(0);
}

template<int DIM> hier::IntVector<DIM> 
FaceIntegerConstantRefine<DIM>::getStencilWidth() const {
   return(hier::IntVector<DIM>(0));
}

template<int DIM> void FaceIntegerConstantRefine<DIM>::refine(
   hier::Patch<DIM>& fine, 
   const hier::Patch<DIM>& coarse, 
   const int dst_component, 
   const int src_component, 
   const hier::Box<DIM>& fine_box, 
   const hier::IntVector<DIM>& ratio) const
{
   tbox::Pointer< FaceData<DIM,int> >
      cdata = coarse.getPatchData(src_component);
   tbox::Pointer< FaceData<DIM,int> >
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

   const hier::Box<DIM> coarse_box = hier::Box<DIM>::coarsen(fine_box, ratio);
   const hier::Index<DIM> ifirstc = coarse_box.lower();
   const hier::Index<DIM> ilastc = coarse_box.upper();
   const hier::Index<DIM> ifirstf = fine_box.lower();
   const hier::Index<DIM> ilastf = fine_box.upper();

   for (int d = 0; d < fdata->getDepth(); d++) {
      if (DIM == 1) {
	 conreffaceintg1d_(ifirstc(0),ilastc(0),
			   ifirstf(0),ilastf(0),
			   cilo(0),cihi(0),
			   filo(0),fihi(0),
			   ratio,
			   cdata->getPointer(0,d),
			   fdata->getPointer(0,d));
      } else if (DIM == 2) {
	 conreffaceintg2d0_(ifirstc(0),ifirstc(1),ilastc(0),ilastc(1),
			    ifirstf(0),ifirstf(1),ilastf(0),ilastf(1),
			    cilo(0),cilo(1),cihi(0),cihi(1),
			    filo(0),filo(1),fihi(0),fihi(1),
			    ratio,
			    cdata->getPointer(0,d),
			    fdata->getPointer(0,d));
	 conreffaceintg2d1_(ifirstc(0),ifirstc(1),ilastc(0),ilastc(1),
			    ifirstf(0),ifirstf(1),ilastf(0),ilastf(1),
			    cilo(0),cilo(1),cihi(0),cihi(1),
			    filo(0),filo(1),fihi(0),fihi(1),
			    ratio,
			    cdata->getPointer(1,d),
			    fdata->getPointer(1,d));
      } else if (DIM == 3) {
	 conreffaceintg3d0_(ifirstc(0),ifirstc(1),ifirstc(2),
			    ilastc(0),ilastc(1),ilastc(2),
			    ifirstf(0),ifirstf(1),ifirstf(2),
			    ilastf(0),ilastf(1),ilastf(2),
			    cilo(0),cilo(1),cilo(2),
			    cihi(0),cihi(1),cihi(2),
			    filo(0),filo(1),filo(2),
			    fihi(0),fihi(1),fihi(2),
			    ratio,
			    cdata->getPointer(0,d),
			    fdata->getPointer(0,d));
	 conreffaceintg3d1_(ifirstc(0),ifirstc(1),ifirstc(2),
			    ilastc(0),ilastc(1),ilastc(2),
			    ifirstf(0),ifirstf(1),ifirstf(2),
			    ilastf(0),ilastf(1),ilastf(2),
			    cilo(0),cilo(1),cilo(2),
			    cihi(0),cihi(1),cihi(2),
			    filo(0),filo(1),filo(2),
			    fihi(0),fihi(1),fihi(2),
			    ratio,
			    cdata->getPointer(1,d),
			    fdata->getPointer(1,d));
	 conreffaceintg3d2_(ifirstc(0),ifirstc(1),ifirstc(2),
			    ilastc(0),ilastc(1),ilastc(2),
			    ifirstf(0),ifirstf(1),ifirstf(2),
			    ilastf(0),ilastf(1),ilastf(2),
			    cilo(0),cilo(1),cilo(2),
			    cihi(0),cihi(1),cihi(2),
			    filo(0),filo(1),filo(2),
			    fihi(0),fihi(1),fihi(2),
			    ratio,
			    cdata->getPointer(2,d),
			    fdata->getPointer(2,d));
      } else {
	 TBOX_ERROR("FaceIntegerConstanRefine::refine DIM > 3 not supported" << endl);
      }
   }
}

}
}
#endif
