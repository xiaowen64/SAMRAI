//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/operators/constant/face/FaceFloatConstantRefine.C $
// Package:	SAMRAI patchdata
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description: Constant refine operator for face-centered float data on 
//              a  mesh.
//

#ifndef included_pdat_FaceFloatConstantRefine_C
#define included_pdat_FaceFloatConstantRefine_C

#include "FaceFloatConstantRefine.h"

#include<float.h>
#include<math.h>
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
// in refine1d.f:
   void conreffaceflot1d_( const int&, const int&,
                           const int&, const int&,
                           const int&, const int&,
                           const int&, const int&,
                           const int*,
                           const float*, float* );
// in refine2d.f:
   void conreffaceflot2d0_( const int&, const int&, const int&, const int&,
                            const int&, const int&, const int&, const int&,
                            const int&, const int&, const int&, const int&,
                            const int&, const int&, const int&, const int&,
                            const int*,
                            const float*, float* );
   void conreffaceflot2d1_( const int&, const int&, const int&, const int&,
                            const int&, const int&, const int&, const int&,
                            const int&, const int&, const int&, const int&,
                            const int&, const int&, const int&, const int&,
                            const int*,
                            const float*, float* );
// in refine3d.f:
   void conreffaceflot3d0_( const int&, const int&, const int&,
                            const int&, const int&, const int&,
                            const int&, const int&, const int&,
                            const int&, const int&, const int&,
                            const int&, const int&, const int&,
                            const int&, const int&, const int&,
                            const int&, const int&, const int&,
                            const int&, const int&, const int&,
                            const int*,
                            const float*, float* );
   void conreffaceflot3d1_( const int&, const int&, const int&,
                            const int&, const int&, const int&,
                            const int&, const int&, const int&,
                            const int&, const int&, const int&,
                            const int&, const int&, const int&,
                            const int&, const int&, const int&,
                            const int&, const int&, const int&,
                            const int&, const int&, const int&,
                            const int*,
                            const float*, float* );
   void conreffaceflot3d2_( const int&, const int&, const int&,
                            const int&, const int&, const int&,
                            const int&, const int&, const int&,
                            const int&, const int&, const int&,
                            const int&, const int&, const int&,
                            const int&, const int&, const int&,
                            const int&, const int&, const int&,
                            const int&, const int&, const int&,
                            const int*,
                            const float*, float* );
}

namespace SAMRAI {
    namespace pdat {

template<int DIM> FaceFloatConstantRefine<DIM>::FaceFloatConstantRefine()
: xfer::RefineOperator<DIM>()
{
   d_name_id = "CONSTANT_REFINE";
}

template<int DIM> FaceFloatConstantRefine<DIM>::~FaceFloatConstantRefine()
{
}

template<int DIM> bool FaceFloatConstantRefine<DIM>::findRefineOperator(
   const tbox::Pointer< hier::Variable<DIM> >& var,
   const std::string &op_name) const
{
   const tbox::Pointer< FaceVariable<DIM,float> > cast_var(var);
   if ( !cast_var.isNull() && (op_name == d_name_id) ) {
      return(true);
   } else {
      return(false);
   }
}

template<int DIM> const std::string&
FaceFloatConstantRefine<DIM>::getOperatorName() const
{
   return(d_name_id);
}

template<int DIM> int FaceFloatConstantRefine<DIM>::getOperatorPriority() const
{
   return(0);
}

template<int DIM> hier::IntVector<DIM> 
FaceFloatConstantRefine<DIM>::getStencilWidth() const {
   return(hier::IntVector<DIM>(0));
}

template<int DIM> void FaceFloatConstantRefine<DIM>::refine(
   hier::Patch<DIM>& fine, 
   const hier::Patch<DIM>& coarse, 
   const int dst_component, 
   const int src_component, 
   const hier::Box<DIM>& fine_box, 
   const hier::IntVector<DIM>& ratio) const
{
   tbox::Pointer< FaceData<DIM,float> >
      cdata = coarse.getPatchData(src_component);
   tbox::Pointer< FaceData<DIM,float> >
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
      if (DIM == 1) {
	 conreffaceflot1d_(ifirstc(0),ilastc(0),
			   ifirstf(0),ilastf(0),
			   cilo(0),cihi(0),
			   filo(0),fihi(0),
			   ratio,
			   cdata->getPointer(0,d),
			   fdata->getPointer(0,d));
      } else if (DIM == 2) {
	 conreffaceflot2d0_(ifirstc(0),ifirstc(1),ilastc(0),ilastc(1),
			    ifirstf(0),ifirstf(1),ilastf(0),ilastf(1),
			    cilo(0),cilo(1),cihi(0),cihi(1),
			    filo(0),filo(1),fihi(0),fihi(1),
			    ratio,
			    cdata->getPointer(0,d),
			    fdata->getPointer(0,d));
	 conreffaceflot2d1_(ifirstc(0),ifirstc(1),ilastc(0),ilastc(1),
			    ifirstf(0),ifirstf(1),ilastf(0),ilastf(1),
			    cilo(0),cilo(1),cihi(0),cihi(1),
			    filo(0),filo(1),fihi(0),fihi(1),
			    ratio,
			    cdata->getPointer(1,d),
			    fdata->getPointer(1,d));
      } else if (DIM == 3) {
	 conreffaceflot3d0_(ifirstc(0),ifirstc(1),ifirstc(2),
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
	 conreffaceflot3d1_(ifirstc(0),ifirstc(1),ifirstc(2),
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
	 conreffaceflot3d2_(ifirstc(0),ifirstc(1),ifirstc(2),
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
	 TBOX_ERROR("FaceFloatConstantRefine::refine DIM > 3 not supported" << std::endl);
      }
   }
}

}
}
#endif
