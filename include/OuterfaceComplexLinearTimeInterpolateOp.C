//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/operators/time_interpolate/outerface/OuterfaceComplexLinearTimeInterpolateOp.C $
// Package:	SAMRAI patchdata
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Linear time interp operator for complex outerface data.
//




#ifndef included_pdat_OuterfaceComplexLinearTimeInterpolateOp_C
#define included_pdat_OuterfaceComplexLinearTimeInterpolateOp_C

#include "OuterfaceComplexLinearTimeInterpolateOp.h"
#include "tbox/Complex.h"

#include "Box.h"
#include "Index.h"
#include "OuterfaceData.h"
#include "OuterfaceVariable.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"

/*
*************************************************************************
*                                                                       *
* External declarations for FORTRAN  routines.                          *
*                                                                       *
*************************************************************************
*/
extern "C" {
// in lintimint1d.f:
   void lintimeintoutfacecmplx1d_( const int&, const int&,
                                   const int&, const int&,
                                   const int&, const int&,
                                   const int&, const int&,
                                   const double&,
                                   const dcomplex*, const dcomplex*,
                                   dcomplex* );
// in lintimint2d.f:
   void lintimeintoutfacecmplx2d0_( const int&, const int&,
                                    const int&, const int&,
                                    const int&, const int&,
                                    const int&, const int&,
                                    const int&, const int&,
                                    const int&, const int&,
                                    const int&, const int&,
                                    const int&, const int&,
                                    const double&,
                                    const dcomplex*, const dcomplex*,
                                    dcomplex* );
   void lintimeintoutfacecmplx2d1_( const int&, const int&,
                                    const int&, const int&,
                                    const int&, const int&,
                                    const int&, const int&,
                                    const int&, const int&,
                                    const int&, const int&,
                                    const int&, const int&,
                                    const int&, const int&,
                                    const double&,
                                    const dcomplex*, const dcomplex*,
                                    dcomplex* );
// in lintimint3d.f:
   void lintimeintoutfacecmplx3d0_( const int&, const int&, const int&,
                                    const int&, const int&, const int&,
                                    const int&, const int&, const int&,
                                    const int&, const int&, const int&,
                                    const int&, const int&, const int&,
                                    const int&, const int&, const int&,
                                    const int&, const int&, const int&,
                                    const int&, const int&, const int&,
                                    const double&,
                                    const dcomplex*, const dcomplex*,
                                    dcomplex* );
   void lintimeintoutfacecmplx3d1_( const int&, const int&, const int&,
                                    const int&, const int&, const int&,
                                    const int&, const int&, const int&,
                                    const int&, const int&, const int&,
                                    const int&, const int&, const int&,
                                    const int&, const int&, const int&,
                                    const int&, const int&, const int&,
                                    const int&, const int&, const int&,
                                    const double&,
                                    const dcomplex*, const dcomplex*,
                                    dcomplex* );
   void lintimeintoutfacecmplx3d2_( const int&, const int&, const int&,
                                    const int&, const int&, const int&,
                                    const int&, const int&, const int&,
                                    const int&, const int&, const int&,
                                    const int&, const int&, const int&,
                                    const int&, const int&, const int&,
                                    const int&, const int&, const int&,
                                    const int&, const int&, const int&,
                                    const double&,
                                    const dcomplex*, const dcomplex*,
                                    dcomplex* );
}

namespace SAMRAI {
    namespace pdat {

template<int DIM> OuterfaceComplexLinearTimeInterpolateOp<DIM>::OuterfaceComplexLinearTimeInterpolateOp()
: xfer::TimeInterpolateOperator<DIM>()
{
}

template<int DIM> OuterfaceComplexLinearTimeInterpolateOp<DIM>::~OuterfaceComplexLinearTimeInterpolateOp()
{
}

template<int DIM> bool OuterfaceComplexLinearTimeInterpolateOp<DIM>::findTimeInterpolateOperator(
   const tbox::Pointer< hier::Variable<DIM> >& var,
   const std::string &op_name) const
{
   const tbox::Pointer< OuterfaceVariable<DIM,dcomplex> > cast_var(var);
   if ( !cast_var.isNull() && (op_name == "STD_LINEAR_TIME_INTERPOLATE") ) {
      return(true);
   } else {
      return(false);
   }
}

template<int DIM> void OuterfaceComplexLinearTimeInterpolateOp<DIM>::timeInterpolate(
   hier::PatchData<DIM>& dst_data, 
   const hier::Box<DIM>& where, 
   const hier::PatchData<DIM>& src_data_old, 
   const hier::PatchData<DIM>& src_data_new) const
{

   const OuterfaceData<DIM,dcomplex> *old_dat =
      dynamic_cast<const OuterfaceData<DIM,dcomplex> *>(&src_data_old);
   const OuterfaceData<DIM,dcomplex> *new_dat =
      dynamic_cast<const OuterfaceData<DIM,dcomplex> *>(&src_data_new);
   OuterfaceData<DIM,dcomplex> *dst_dat = 
      dynamic_cast<OuterfaceData<DIM,dcomplex> *>(&dst_data);

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( old_dat != NULL );
   TBOX_ASSERT( new_dat != NULL );
   TBOX_ASSERT( dst_dat != NULL );
   TBOX_ASSERT( where*old_dat->getGhostBox() == where );
   TBOX_ASSERT( where*new_dat->getGhostBox() == where );
   TBOX_ASSERT( where*dst_dat->getGhostBox() == where );
#endif

   const hier::Index<DIM> old_ilo = old_dat->getGhostBox().lower();
   const hier::Index<DIM> old_ihi = old_dat->getGhostBox().upper();
   const hier::Index<DIM> new_ilo = new_dat->getGhostBox().lower();
   const hier::Index<DIM> new_ihi = new_dat->getGhostBox().upper();

   const hier::Index<DIM> dst_ilo = dst_dat->getGhostBox().lower();
   const hier::Index<DIM> dst_ihi = dst_dat->getGhostBox().upper();

   const hier::Index<DIM> ifirst = where.lower();
   const hier::Index<DIM> ilast = where.upper();

   const double old_time = old_dat->getTime();
   const double new_time = new_dat->getTime();
   const double dst_time = dst_dat->getTime();
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((old_time < dst_time || tbox::MathUtilities<double>::equalEps(old_time,dst_time)) && 
          (dst_time < new_time || tbox::MathUtilities<double>::equalEps(dst_time,new_time)));
#endif

   double tfrac = dst_time - old_time;
   double denom = new_time - old_time;
   if ( denom > tbox::MathUtilities<double>::getMin() ) {
      tfrac /= denom;
   } else {
      tfrac = 0.0;
   }

   for (int d = 0; d < dst_dat->getDepth(); d++) {
      // loop over lower and upper outerface arrays
      for (int i = 0; i < 2; i++) {
	 if (DIM == 1) {
	    lintimeintoutfacecmplx1d_(ifirst(0),ilast(0),
				      old_ilo(0),old_ihi(0),
				      new_ilo(0),new_ihi(0),
				      dst_ilo(0),dst_ihi(0),
				      tfrac,
				      old_dat->getPointer(0,i,d),
				      new_dat->getPointer(0,i,d),
				      dst_dat->getPointer(0,i,d));
	 } else if (DIM == 2) {
	    lintimeintoutfacecmplx2d0_(ifirst(0),ifirst(1),ilast(0),ilast(1),
				       old_ilo(0),old_ilo(1),old_ihi(0),old_ihi(1),
				       new_ilo(0),new_ilo(1),new_ihi(0),new_ihi(1),
				       dst_ilo(0),dst_ilo(1),dst_ihi(0),dst_ihi(1),
				       tfrac,
				       old_dat->getPointer(0,i,d),
				       new_dat->getPointer(0,i,d),
				       dst_dat->getPointer(0,i,d));
	    lintimeintoutfacecmplx2d1_(ifirst(0),ifirst(1),ilast(0),ilast(1),
				       old_ilo(0),old_ilo(1),old_ihi(0),old_ihi(1),
				       new_ilo(0),new_ilo(1),new_ihi(0),new_ihi(1),
				       dst_ilo(0),dst_ilo(1),dst_ihi(0),dst_ihi(1),
				       tfrac,
				       old_dat->getPointer(1,i,d),
				       new_dat->getPointer(1,i,d),
				       dst_dat->getPointer(1,i,d));
	 } else if (DIM == 3) {
	    lintimeintoutfacecmplx3d0_(ifirst(0),ifirst(1),ifirst(2),
				       ilast(0),ilast(1),ilast(2),
				       old_ilo(0),old_ilo(1),old_ilo(2),
				       old_ihi(0),old_ihi(1),old_ihi(2),
				       new_ilo(0),new_ilo(1),new_ilo(2),
				       new_ihi(0),new_ihi(1),new_ihi(2),
				       dst_ilo(0),dst_ilo(1),dst_ilo(2),
				       dst_ihi(0),dst_ihi(1),dst_ihi(2),
				       tfrac,
				       old_dat->getPointer(0,i,d),
				       new_dat->getPointer(0,i,d),
				       dst_dat->getPointer(0,i,d));
	    lintimeintoutfacecmplx3d1_(ifirst(0),ifirst(1),ifirst(2),
				       ilast(0),ilast(1),ilast(2),
				       old_ilo(0),old_ilo(1),old_ilo(2),
				       old_ihi(0),old_ihi(1),old_ihi(2),
				       new_ilo(0),new_ilo(1),new_ilo(2),
				       new_ihi(0),new_ihi(1),new_ihi(2),
				       dst_ilo(0),dst_ilo(1),dst_ilo(2),
				       dst_ihi(0),dst_ihi(1),dst_ihi(2),
				       tfrac,
				       old_dat->getPointer(1,i,d),
				       new_dat->getPointer(1,i,d),
				       dst_dat->getPointer(1,i,d));
	    lintimeintoutfacecmplx3d2_(ifirst(0),ifirst(1),ifirst(2),
				       ilast(0),ilast(1),ilast(2),
				       old_ilo(0),old_ilo(1),old_ilo(2),
				       old_ihi(0),old_ihi(1),old_ihi(2),
				       new_ilo(0),new_ilo(1),new_ilo(2),
				       new_ihi(0),new_ihi(1),new_ihi(2),
				       dst_ilo(0),dst_ilo(1),dst_ilo(2),
				       dst_ihi(0),dst_ihi(1),dst_ihi(2),
				       tfrac,
				       old_dat->getPointer(2,i,d),
				       new_dat->getPointer(2,i,d),
				       dst_dat->getPointer(2,i,d));
	 } else {
	    TBOX_ERROR("OuterfaceComplexLinearTimeInterpolateOp::TimeInterpolate DIM > 3 not supported" << std::endl);
	 }
      }
   }
}

}
}
#endif
