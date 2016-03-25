//
// File:	FaceComplexLinearTimeInterpolateOp.C
// Package:	SAMRAI patchdata
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Linear time interp operator for face-centered complex data.
//

#ifndef included_pdat_FaceComplexLinearTimeInterpolateOp_C
#define included_pdat_FaceComplexLinearTimeInterpolateOp_C

#include "FaceComplexLinearTimeInterpolateOp.h"
#include "tbox/Complex.h"


#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif
#include "Box.h"
#include "Index.h"
#include "FaceData.h"
#include "FaceVariable.h"
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
   void lintimeintfacecmplx1d_( const int&, const int&,
                                const int&, const int&,
                                const int&, const int&,
                                const int&, const int&,
                                const double&,
                                const dcomplex*, const dcomplex*,
                                dcomplex* );
// in lintimint2d.f:
   void lintimeintfacecmplx2d0_( const int&, const int&,
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
   void lintimeintfacecmplx2d1_( const int&, const int&,
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
   void lintimeintfacecmplx3d0_( const int&, const int&, const int&,
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
   void lintimeintfacecmplx3d1_( const int&, const int&, const int&,
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
   void lintimeintfacecmplx3d2_( const int&, const int&, const int&,
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

template<int DIM> FaceComplexLinearTimeInterpolateOp<DIM>::FaceComplexLinearTimeInterpolateOp()
: xfer::TimeInterpolateOperator<DIM>()
{
}

template<int DIM> FaceComplexLinearTimeInterpolateOp<DIM>::~FaceComplexLinearTimeInterpolateOp()
{
}

template<int DIM> bool FaceComplexLinearTimeInterpolateOp<DIM>::findTimeInterpolateOperator(
   const tbox::Pointer< hier::Variable<DIM> >& var,
   const string &op_name) const
{
   const tbox::Pointer< FaceVariable<DIM,dcomplex> > cast_var(var);
   if ( !cast_var.isNull() && (op_name == "STD_LINEAR_TIME_INTERPOLATE") ) {
      return(true);
   } else {
      return(false);
   }
}

template<int DIM> void FaceComplexLinearTimeInterpolateOp<DIM>::timeInterpolate(
   hier::PatchData<DIM>& dst_data, 
   const hier::Box<DIM>& where, 
   const hier::PatchData<DIM>& src_data_old, 
   const hier::PatchData<DIM>& src_data_new) const
{

   const FaceData<DIM,dcomplex> *old_dat =
      dynamic_cast<const FaceData<DIM,dcomplex> *>(&src_data_old);
   const FaceData<DIM,dcomplex> *new_dat =
      dynamic_cast<const FaceData<DIM,dcomplex> *>(&src_data_new);
   FaceData<DIM,dcomplex> *dst_dat = 
      dynamic_cast<FaceData<DIM,dcomplex> *>(&dst_data);

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
      if (DIM == 1) {
	 lintimeintfacecmplx1d_(ifirst(0),ilast(0),
				old_ilo(0),old_ihi(0),
				new_ilo(0),new_ihi(0),
				dst_ilo(0),dst_ihi(0),
				tfrac,
				old_dat->getPointer(0,d),
				new_dat->getPointer(0,d),
				dst_dat->getPointer(0,d));
      } else if (DIM == 2) {
	 lintimeintfacecmplx2d0_(ifirst(0),ifirst(1),ilast(0),ilast(1),
				 old_ilo(0),old_ilo(1),old_ihi(0),old_ihi(1),
				 new_ilo(0),new_ilo(1),new_ihi(0),new_ihi(1),
				 dst_ilo(0),dst_ilo(1),dst_ihi(0),dst_ihi(1),
				 tfrac,
				 old_dat->getPointer(0,d),
				 new_dat->getPointer(0,d),
				 dst_dat->getPointer(0,d));
	 lintimeintfacecmplx2d1_(ifirst(0),ifirst(1),ilast(0),ilast(1),
				 old_ilo(0),old_ilo(1),old_ihi(0),old_ihi(1),
				 new_ilo(0),new_ilo(1),new_ihi(0),new_ihi(1),
				 dst_ilo(0),dst_ilo(1),dst_ihi(0),dst_ihi(1),
				 tfrac,
				 old_dat->getPointer(1,d),
				 new_dat->getPointer(1,d),
				 dst_dat->getPointer(1,d));
      } else if (DIM == 3) {
	 lintimeintfacecmplx3d0_(ifirst(0),ifirst(1),ifirst(2),
				 ilast(0),ilast(1),ilast(2),
				 old_ilo(0),old_ilo(1),old_ilo(2),
				 old_ihi(0),old_ihi(1),old_ihi(2),
				 new_ilo(0),new_ilo(1),new_ilo(2),
				 new_ihi(0),new_ihi(1),new_ihi(2),
				 dst_ilo(0),dst_ilo(1),dst_ilo(2),
				 dst_ihi(0),dst_ihi(1),dst_ihi(2),
				 tfrac,
				 old_dat->getPointer(0,d),
				 new_dat->getPointer(0,d),
				 dst_dat->getPointer(0,d));
	 lintimeintfacecmplx3d1_(ifirst(0),ifirst(1),ifirst(2),
				 ilast(0),ilast(1),ilast(2),
				 old_ilo(0),old_ilo(1),old_ilo(2),
				 old_ihi(0),old_ihi(1),old_ihi(2),
				 new_ilo(0),new_ilo(1),new_ilo(2),
				 new_ihi(0),new_ihi(1),new_ihi(2),
				 dst_ilo(0),dst_ilo(1),dst_ilo(2),
				 dst_ihi(0),dst_ihi(1),dst_ihi(2),
				 tfrac,
				 old_dat->getPointer(1,d),
				 new_dat->getPointer(1,d),
				 dst_dat->getPointer(1,d));
	 lintimeintfacecmplx3d2_(ifirst(0),ifirst(1),ifirst(2),
				 ilast(0),ilast(1),ilast(2),
				 old_ilo(0),old_ilo(1),old_ilo(2),
				 old_ihi(0),old_ihi(1),old_ihi(2),
				 new_ilo(0),new_ilo(1),new_ilo(2),
				 new_ihi(0),new_ihi(1),new_ihi(2),
				 dst_ilo(0),dst_ilo(1),dst_ilo(2),
				 dst_ihi(0),dst_ihi(1),dst_ihi(2),
				 tfrac,
				 old_dat->getPointer(2,d),
				 new_dat->getPointer(2,d),
				 dst_dat->getPointer(2,d));
      } else {
	 TBOX_ERROR("FaceComplexLinearTimeInterpolateOp::TimeInterpolate DIM > 3 not supported" << endl);
      }
   }
}

}
}
#endif
