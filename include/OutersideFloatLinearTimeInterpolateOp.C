//
// File:	OutersideFloatLinearTimeInterpolateOp.C
// Package:	SAMRAI patchdata
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Linear time interp operator for float outerside patch data.
//




#ifndef included_pdat_OutersideFloatLinearTimeInterpolateOp_C
#define included_pdat_OutersideFloatLinearTimeInterpolateOp_C

#include "OutersideFloatLinearTimeInterpolateOp.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif
#include "Box.h"
#include "Index.h"
#include "OutersideData.h"
#include "OutersideVariable.h"
#include "tbox/IEEE.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

/*
*************************************************************************
*                                                                       *
* External declarations for FORTRAN  routines.                          *
*                                                                       *
*************************************************************************
*/
extern "C" {
// in lintimint1d.f:
   void lintimeintoutsidefloat1d_( const int&, const int&,
                                   const int&, const int&,
                                   const int&, const int&,
                                   const int&, const int&,
                                   const double&,
                                   const float*, const float*,
                                   float* );
// in lintimint2d.f:
   void lintimeintoutsidefloat2d0_( const int&, const int&,
                                    const int&, const int&,
                                    const int&, const int&,
                                    const int&, const int&,
                                    const int&, const int&,
                                    const int&, const int&,
                                    const int&, const int&,
                                    const int&, const int&,
                                    const double&,
                                    const float*, const float*,
                                    float* );
   void lintimeintoutsidefloat2d1_( const int&, const int&,
                                    const int&, const int&,
                                    const int&, const int&,
                                    const int&, const int&,
                                    const int&, const int&,
                                    const int&, const int&,
                                    const int&, const int&,
                                    const int&, const int&,
                                    const double&,
                                    const float*, const float*,
                                    float* );
// in lintimint3d.f:
   void lintimeintoutsidefloat3d0_( const int&, const int&, const int&,
                                    const int&, const int&, const int&,
                                    const int&, const int&, const int&,
                                    const int&, const int&, const int&,
                                    const int&, const int&, const int&,
                                    const int&, const int&, const int&,
                                    const int&, const int&, const int&,
                                    const int&, const int&, const int&,
                                    const double&,
                                    const float*, const float*,
                                    float* );
   void lintimeintoutsidefloat3d1_( const int&, const int&, const int&,
                                    const int&, const int&, const int&,
                                    const int&, const int&, const int&,
                                    const int&, const int&, const int&,
                                    const int&, const int&, const int&,
                                    const int&, const int&, const int&,
                                    const int&, const int&, const int&,
                                    const int&, const int&, const int&,
                                    const double&,
                                    const float*, const float*,
                                    float* );
   void lintimeintoutsidefloat3d2_( const int&, const int&, const int&,
                                    const int&, const int&, const int&,
                                    const int&, const int&, const int&,
                                    const int&, const int&, const int&,
                                    const int&, const int&, const int&,
                                    const int&, const int&, const int&,
                                    const int&, const int&, const int&,
                                    const int&, const int&, const int&,
                                    const double&,
                                    const float*, const float*,
                                    float* );
}

namespace SAMRAI {
    namespace pdat {

template<int DIM> OutersideFloatLinearTimeInterpolateOp<DIM>::OutersideFloatLinearTimeInterpolateOp()
: xfer::TimeInterpolateOperator<DIM>()
{
}

template<int DIM> OutersideFloatLinearTimeInterpolateOp<DIM>::~OutersideFloatLinearTimeInterpolateOp()
{
}

template<int DIM> bool OutersideFloatLinearTimeInterpolateOp<DIM>::findTimeInterpolateOperator(
   const tbox::Pointer< hier::Variable<DIM> >& var,
   const string &op_name) const
{
   const tbox::Pointer< OutersideVariable<DIM,float> > cast_var(var);
   if ( !cast_var.isNull() && (op_name == "STD_LINEAR_TIME_INTERPOLATE") ) {
      return(true);
   } else {
      return(false);
   }
}

template<int DIM> void OutersideFloatLinearTimeInterpolateOp<DIM>::timeInterpolate(
   hier::PatchData<DIM>& dst_data, 
   const hier::Box<DIM>& where, 
   const hier::PatchData<DIM>& src_data_old, 
   const hier::PatchData<DIM>& src_data_new) const
{

   const OutersideData<DIM,float> *old_dat =
      dynamic_cast<const OutersideData<DIM,float> *>(&src_data_old);
   const OutersideData<DIM,float> *new_dat =
      dynamic_cast<const OutersideData<DIM,float> *>(&src_data_new);
   OutersideData<DIM,float> *dst_dat = 
      dynamic_cast<OutersideData<DIM,float> *>(&dst_data);

#ifdef DEBUG_CHECK_ASSERTIONS
   assert( old_dat != NULL );
   assert( new_dat != NULL );
   assert( dst_dat != NULL );
   assert( where*old_dat->getGhostBox() == where );
   assert( where*new_dat->getGhostBox() == where );
   assert( where*dst_dat->getGhostBox() == where );
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
   assert((old_time < dst_time || tbox::Utilities::deq(old_time,dst_time)) && 
          (dst_time < new_time || tbox::Utilities::deq(dst_time,new_time)));
#endif

   double tfrac = dst_time - old_time;
   double denom = new_time - old_time;
   if ( denom > tbox::IEEE::getDBL_MIN() ) {
      tfrac /= denom;
   } else {
      tfrac = 0.0;
   }

   for (int d = 0; d < dst_dat->getDepth(); d++) {
      // loop over lower and upper outerside arrays
      for (int i = 0; i < 2; i++) {
	 if (DIM == 1) {
	    lintimeintoutsidefloat1d_(ifirst(0),ilast(0),
				      old_ilo(0),old_ihi(0),
				      new_ilo(0),new_ihi(0),
				      dst_ilo(0),dst_ihi(0),
				      tfrac,
				      old_dat->getPointer(0,i,d),
				      new_dat->getPointer(0,i,d),
				      dst_dat->getPointer(0,i,d));
	 } else if (DIM == 2) {
	    lintimeintoutsidefloat2d0_(ifirst(0),ifirst(1),ilast(0),ilast(1),
				       old_ilo(0),old_ilo(1),old_ihi(0),old_ihi(1),
				       new_ilo(0),new_ilo(1),new_ihi(0),new_ihi(1),
				       dst_ilo(0),dst_ilo(1),dst_ihi(0),dst_ihi(1),
				       tfrac,
				       old_dat->getPointer(0,i,d),
				       new_dat->getPointer(0,i,d),
				       dst_dat->getPointer(0,i,d));
	    lintimeintoutsidefloat2d1_(ifirst(0),ifirst(1),ilast(0),ilast(1),
				       old_ilo(0),old_ilo(1),old_ihi(0),old_ihi(1),
				       new_ilo(0),new_ilo(1),new_ihi(0),new_ihi(1),
				       dst_ilo(0),dst_ilo(1),dst_ihi(0),dst_ihi(1),
				       tfrac,
				       old_dat->getPointer(1,i,d),
				       new_dat->getPointer(1,i,d),
				       dst_dat->getPointer(1,i,d));
	 } else if (DIM == 3) {
	    lintimeintoutsidefloat3d0_(ifirst(0),ifirst(1),ifirst(2),
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
	    lintimeintoutsidefloat3d1_(ifirst(0),ifirst(1),ifirst(2),
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
	    lintimeintoutsidefloat3d2_(ifirst(0),ifirst(1),ifirst(2),
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
	    TBOX_ERROR("OutersideFloatLinearTimeInterpolateOp::TimeInterpolate DIM > 3 not supported" << endl);
	 }
      }
   }
}

}
}
#endif
