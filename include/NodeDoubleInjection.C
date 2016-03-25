//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/operators/constant/node/NodeDoubleInjection.C $
// Package:	SAMRAI patchdata
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Constant averaging operator for node-centered double data on
//              a  mesh.
//

#ifndef included_pdat_NodeDoubleInjection_C
#define included_pdat_NodeDoubleInjection_C

#include "NodeDoubleInjection.h"

#include<float.h>
#include<math.h>
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
   void conavgnodedoub1d_( const int&, const int&,
                           const int&, const int&,
                           const int&, const int&,
                           const int*,
                           const double*, double* );
// in concoarsen2d.f:
   void conavgnodedoub2d_( const int&, const int&, const int&, const int&,
                           const int&, const int&, const int&, const int&,
                           const int&, const int&, const int&, const int&,
                           const int*,
                           const double*, double* );
// in concoarsen3d.f:
   void conavgnodedoub3d_( const int&, const int&, const int&,
                           const int&, const int&, const int&,
                           const int&, const int&, const int&,
                           const int&, const int&, const int&,
                           const int&, const int&, const int&,
                           const int&, const int&, const int&,
                           const int*,
                           const double*, double* );
}

namespace SAMRAI {
    namespace pdat {

template<int DIM> NodeDoubleInjection<DIM>::NodeDoubleInjection()
: xfer::CoarsenOperator<DIM>()
{
   d_name_id = "CONSTANT_COARSEN";
}

template<int DIM> NodeDoubleInjection<DIM>::~NodeDoubleInjection()
{
}

template<int DIM> bool NodeDoubleInjection<DIM>::findCoarsenOperator(
   const tbox::Pointer< hier::Variable<DIM> >& var,
   const std::string &op_name) const
{
   const tbox::Pointer< NodeVariable<DIM,double> > cast_var(var);
   if ( !cast_var.isNull() && (op_name == d_name_id) ) {
      return(true);
   } else {
      return(false);
   }
}

template<int DIM> const std::string&
NodeDoubleInjection<DIM>::getOperatorName() const
{
   return(d_name_id);
}

template<int DIM> int NodeDoubleInjection<DIM>::getOperatorPriority() const
{
   return(0);
}

template<int DIM> hier::IntVector<DIM>
NodeDoubleInjection<DIM>::getStencilWidth() const {
   return(hier::IntVector<DIM>(0));
}

template<int DIM> void NodeDoubleInjection<DIM>::coarsen(
   hier::Patch<DIM>& coarse,
   const hier::Patch<DIM>& fine,
   const int dst_component,
   const int src_component,
   const hier::Box<DIM>& coarse_box,
   const hier::IntVector<DIM>& ratio) const
{
   tbox::Pointer< NodeData<DIM,double> >
      fdata = fine.getPatchData(src_component);
   tbox::Pointer< NodeData<DIM,double> >
      cdata = coarse.getPatchData(dst_component);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!fdata.isNull());
   TBOX_ASSERT(!cdata.isNull());
   TBOX_ASSERT(cdata->getDepth() == fdata->getDepth());
#endif

   const hier::Index<DIM> filo = fdata->getGhostBox().lower();
   const hier::Index<DIM> fihi = fdata->getGhostBox().upper();
   const hier::Index<DIM> cilo = cdata->getGhostBox().lower();
   const hier::Index<DIM> cihi = cdata->getGhostBox().upper();

   const hier::Index<DIM> ifirstc = coarse_box.lower();
   const hier::Index<DIM> ilastc = coarse_box.upper();

   for (int d = 0; d < cdata->getDepth(); d++) {
      if (DIM == 1) {
	 conavgnodedoub1d_(ifirstc(0),ilastc(0),
			   filo(0),fihi(0),
			   cilo(0),cihi(0),
			   ratio,
			   fdata->getPointer(d),
			   cdata->getPointer(d));
      } else if (DIM == 2) {
	 conavgnodedoub2d_(ifirstc(0),ifirstc(1),ilastc(0),ilastc(1),
			   filo(0),filo(1),fihi(0),fihi(1),
			   cilo(0),cilo(1),cihi(0),cihi(1),
			   ratio,
			   fdata->getPointer(d),
			   cdata->getPointer(d));
      } else if (DIM == 3) {
	 conavgnodedoub3d_(ifirstc(0),ifirstc(1),ifirstc(2),
			   ilastc(0),ilastc(1),ilastc(2),
			   filo(0),filo(1),filo(2),
			   fihi(0),fihi(1),fihi(2),
			   cilo(0),cilo(1),cilo(2),
			   cihi(0),cihi(1),cihi(2),
			   ratio,
			   fdata->getPointer(d),
			   cdata->getPointer(d));
      } else {
	 TBOX_ERROR("NodeDoubleConstantRefine::coarsen DIM > 3 not supported" << std::endl);
      }
   }
}

}
}
#endif
