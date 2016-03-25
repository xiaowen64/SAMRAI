//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/operators/constant/outernode/OuternodeDoubleConstantCoarsen.C $
// Package:	SAMRAI patchdata
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: ConstantCoarsen averaging operator for outernode-centered 
//              double data on a mesh.
//

#include "OuternodeDoubleConstantCoarsen.h"

#include <float.h>
#include <math.h>

#ifdef DEBUG_CHECK_ASSERTIONS
#include "tbox/Utilities.h"
#endif
#include "Index.h"
#include "OuternodeData.h"
#include "OuternodeVariable.h"

/*
*************************************************************************
*                                                                       *
* External declarations for FORTRAN  routines.                          *
*                                                                       *
*************************************************************************
*/
extern "C" {
// in pdat_concoarsen1d.f:
void conavgouternodedoub1d_(  const int&, const int&, 
                              const int&, const int&,
                              const int&, const int&,
                              const int*, 
                              const double*, double*);
// in pdat_concoarsen2d.f:
void conavgouternodedoub2d0_(  const int&, const int&, 
                               const int&, const int&,
                               const int&, const int&, 
                               const int&, const int&,
                               const int&, const int&,
                               const int&, const int&,
                               const int*, 
                               const double*, double*);
void conavgouternodedoub2d1_(  const int&, const int&, 
                               const int&, const int&,
                               const int&, const int&, 
                               const int&, const int&,
                               const int&, const int&,
                               const int&, const int&,
                               const int*, 
                               const double*, double*);
   

// in pdat_concoarsen3d.f:
void conavgouternodedoub3d0_(  const int&, const int&, const int&, 
                               const int&, const int&, const int&,
                               const int&, const int&, const int&,
                               const int&, const int&, const int&,
                               const int&, const int&, const int&,
                               const int&, const int&, const int&,
                               const int*, 
                               const double*, double*);
void conavgouternodedoub3d1_(  const int&, const int&, const int&, 
                               const int&, const int&, const int&,
                               const int&, const int&, const int&,
                               const int&, const int&, const int&,
                               const int&, const int&, const int&,
                               const int&, const int&, const int&,
                               const int*, 
                               const double*, double*);
void conavgouternodedoub3d2_(  const int&, const int&, const int&, 
                               const int&, const int&, const int&,
                               const int&, const int&, const int&,
                               const int&, const int&, const int&,
                               const int&, const int&, const int&,
                               const int&, const int&, const int&,
                               const int*, 
                               const double*, double*);
}


namespace SAMRAI {
   namespace pdat {


template<int DIM> OuternodeDoubleConstantCoarsen<DIM>::OuternodeDoubleConstantCoarsen()
: xfer::CoarsenOperator<DIM>()
{
   d_name_id = "CONSTANT_COARSEN";
   return;
}

template<int DIM> OuternodeDoubleConstantCoarsen<DIM>::~OuternodeDoubleConstantCoarsen()
{
   return;
}

template<int DIM> bool OuternodeDoubleConstantCoarsen<DIM>::findCoarsenOperator(
   const tbox::Pointer<hier::Variable<DIM> >& var,
   const std::string &op_name) const
{
   const tbox::Pointer< pdat::OuternodeVariable<DIM, double> > cast_var(var);
   if ( !cast_var.isNull() && (op_name == d_name_id) ) {
      return(true);
   } else {
      return(false);
   }
}

template<int DIM> const std::string&
OuternodeDoubleConstantCoarsen<DIM>::getOperatorName() const
{
   return(d_name_id);
}

template<int DIM> int OuternodeDoubleConstantCoarsen<DIM>::getOperatorPriority() const
{
   return(0);
}

template<int DIM> hier::IntVector<DIM>
OuternodeDoubleConstantCoarsen<DIM>::getStencilWidth() const {
   return(hier::IntVector<DIM>(0));
}



template<int DIM> void OuternodeDoubleConstantCoarsen<DIM>::coarsen(
   hier::Patch<DIM>& coarse, 
   const hier::Patch<DIM>& fine, 
   const int dst_component, 
   const int src_component, 
   const hier::Box<DIM>& coarse_box, 
   const hier::IntVector<DIM>& ratio) const 
{
   tbox::Pointer< pdat::OuternodeData<DIM,double> > 
      fdata = fine.getPatchData(src_component);
   tbox::Pointer< pdat::OuternodeData<DIM,double> > 
      cdata = coarse.getPatchData(dst_component);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!fdata.isNull());
   TBOX_ASSERT(!cdata.isNull());
   TBOX_ASSERT(cdata->getDepth() == fdata->getDepth());
#endif

   const hier::Index<DIM> filo = fine.getBox().lower();
   const hier::Index<DIM> fihi = fine.getBox().upper();
   const hier::Index<DIM> cilo = coarse.getBox().lower();
   const hier::Index<DIM> cihi = coarse.getBox().upper();

   const hier::Index<DIM> ifirstc = coarse_box.lower();
   const hier::Index<DIM> ilastc = coarse_box.upper();

   for ( int i = 0; i < 2; i++ ) {     

      for ( int axis = 0; axis < DIM; axis++ ) {     

         if ( cdata->dataExists(axis) ) {

            for (int d = 0; d < cdata->getDepth(); d++) {

	       if (DIM == 1) {
		  conavgouternodedoub1d_(ifirstc(0),ilastc(0),
					 filo(0), fihi(0),
					 cilo(0),cihi(0),
					 ratio,
					 fdata->getPointer(axis,i,d),
					 cdata->getPointer(axis,i,d));
	       } else if (DIM == 2) {
		  if (axis == 0) {
		     conavgouternodedoub2d0_(ifirstc(0), ifirstc(1), 
					     ilastc(0), ilastc(1),
					     filo(0), filo(1), 
					     fihi(0), fihi(1),
					     cilo(0), cilo(1), 
					     cihi(0), cihi(1),
					     ratio,
					     fdata->getPointer(axis,i,d),
					     cdata->getPointer(axis,i,d));
		  }
		  
		  if (axis == 1) {
		     conavgouternodedoub2d1_(ifirstc(0), ifirstc(1), 
					     ilastc(0), ilastc(1),
					     filo(0), filo(1), 
					     fihi(0), fihi(1),
					     cilo(0), cilo(1), 
					     cihi(0), cihi(1),
					     ratio,
					     fdata->getPointer(axis,i,d),
					     cdata->getPointer(axis,i,d));
		  }
	       } else if (DIM == 3) {
		  if (axis == 0) {
		     conavgouternodedoub3d0_(ifirstc(0), ifirstc(1), ifirstc(2), 
					     ilastc(0), ilastc(1), ilastc(2),
					     filo(0), filo(1), filo(2),
					     fihi(0), fihi(1), fihi(2),
					     cilo(0), cilo(1), cilo(2),
					     cihi(0), cihi(1), cihi(2),
					     ratio,
					     fdata->getPointer(axis,i,d),
					     cdata->getPointer(axis,i,d));
		  }
		  if (axis == 1) {
		     conavgouternodedoub3d1_(ifirstc(0), ifirstc(1), ifirstc(2), 
					     ilastc(0), ilastc(1), ilastc(2),
					     filo(0), filo(1), filo(2),
					     fihi(0), fihi(1), fihi(2),
					     cilo(0), cilo(1), cilo(2),
					     cihi(0), cihi(1), cihi(2),
					     ratio,
					     fdata->getPointer(axis,i,d),
					     cdata->getPointer(axis,i,d));
		  }
		  if (axis == 2) {
		     conavgouternodedoub3d2_(ifirstc(0), ifirstc(1), ifirstc(2), 
					     ilastc(0), ilastc(1), ilastc(2),
					     filo(0), filo(1), filo(2),
					     fihi(0), fihi(1), fihi(2),
					     cilo(0), cilo(1), cilo(2),
					     cihi(0), cihi(1), cihi(2),
					     ratio,
					     fdata->getPointer(axis,i,d),
					     cdata->getPointer(axis,i,d));
		  }
	       }
		  
	    }
         }
      }
   }
}


}
}
