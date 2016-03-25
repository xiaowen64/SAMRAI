//
// File:	NodeIntegerConstantAverage.C
// Package:	SAMRAI patchdata
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Constant averaging operator for node-centered integer data on 
//              a  mesh.
//

#ifndef included_pdat_NodeIntegerConstantAverage_C
#define included_pdat_NodeIntegerConstantAverage_C

#include "NodeIntegerConstantAverage.h"

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
#include "NodeData.h"
#include "NodeVariable.h"

/*
*************************************************************************
*                                                                       *
* External declarations for FORTRAN  routines.                          *
*                                                                       *
*************************************************************************
*/
extern "C" {
// in concoarsen1d.f:
   void conavgnodeintg1d_( const int&, const int&,
                           const int&, const int&,
                           const int&, const int&,
                           const int*,
                           const int*, int* );
// in concoarsen2d.f:
   void conavgnodeintg2d_( const int&, const int&, const int&, const int&,
                           const int&, const int&, const int&, const int&,
                           const int&, const int&, const int&, const int&,
                           const int*,
                           const int*, int* );
// in concoarsen3d.f:
   void conavgnodeintg3d_( const int&, const int&, const int&,
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

template<int DIM> NodeIntegerConstantAverage<DIM>::NodeIntegerConstantAverage()
: xfer::CoarsenOperator<DIM>()
{
   d_name_id = "CONSTANT_COARSEN";
}

template<int DIM> NodeIntegerConstantAverage<DIM>::~NodeIntegerConstantAverage()
{
}

template<int DIM> bool NodeIntegerConstantAverage<DIM>::findCoarsenOperator(
   const tbox::Pointer< hier::Variable<DIM> >& var,
   const string &op_name) const
{
   const tbox::Pointer< NodeVariable<DIM,int> > cast_var(var);
   if ( !cast_var.isNull() && (op_name == d_name_id) ) {
      return(true);
   } else {
      return(false);
   }
}

template<int DIM> const string&
NodeIntegerConstantAverage<DIM>::getOperatorName() const
{
   return(d_name_id);
}

template<int DIM> int NodeIntegerConstantAverage<DIM>::getOperatorPriority() const
{
   return(0);
}

template<int DIM> hier::IntVector<DIM> 
NodeIntegerConstantAverage<DIM>::getStencilWidth() const {
   return(hier::IntVector<DIM>(0));
}

template<int DIM> void NodeIntegerConstantAverage<DIM>::coarsen(
   hier::Patch<DIM>& coarse, 
   const hier::Patch<DIM>& fine, 
   const int dst_component, 
   const int src_component, 
   const hier::Box<DIM>& coarse_box, 
   const hier::IntVector<DIM>& ratio) const 
{
   tbox::Pointer< NodeData<DIM,int> > 
      fdata = fine.getPatchData(src_component);
   tbox::Pointer< NodeData<DIM,int> > 
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

   const hier::Index<DIM> ifirstc = coarse_box.lower();
   const hier::Index<DIM> ilastc = coarse_box.upper();

   for (int d = 0; d < cdata->getDepth(); d++) {
      if (DIM == 1) {
	 conavgnodeintg1d_(ifirstc(0),ilastc(0),
			   filo(0),fihi(0),
			   cilo(0),cihi(0),
			   ratio,
			   fdata->getPointer(d),
			   cdata->getPointer(d));
      } else if (DIM == 2) {
	 conavgnodeintg2d_(ifirstc(0),ifirstc(1),ilastc(0),ilastc(1),
			   filo(0),filo(1),fihi(0),fihi(1),
			   cilo(0),cilo(1),cihi(0),cihi(1),
			   ratio,
			   fdata->getPointer(d),
			   cdata->getPointer(d));
      } else if (DIM == 3) {
	 conavgnodeintg3d_(ifirstc(0),ifirstc(1),ifirstc(2),
			   ilastc(0),ilastc(1),ilastc(2),
			   filo(0),filo(1),filo(2),
			   fihi(0),fihi(1),fihi(2),
			   cilo(0),cilo(1),cilo(2),
			   cihi(0),cihi(1),cihi(2),
			   ratio,
			   fdata->getPointer(d),
			   cdata->getPointer(d));
      } else {
	 TBOX_ERROR("NodeIntegerConstantRefine::coarsen DIM > 3 not supported" << endl);
      }
   }
}

}
}
#endif
