//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/transfer/operators/Geometry.C $
// Package:	SAMRAI transfer
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Base class for interface between transfer ops and geometry.
//

#ifndef included_xfer_Geometry_C
#define included_xfer_Geometry_C

#include "Geometry.h"

#include "tbox/Utilities.h"

#ifndef NULL
#define NULL (0)
#endif

namespace SAMRAI {
    namespace xfer {

/*
*************************************************************************
*                                                                       *
* Constructor and destructor for Geometry objects.                 *
*                                                                       *
*************************************************************************
*/

template<int DIM>  Geometry<DIM>::Geometry(const std::string &object_name) 
: hier::GridGeometry<DIM>(object_name)
{
}

template<int DIM>  Geometry<DIM>::~Geometry()
{
}

/*
*************************************************************************
*                                                                       *
* Add operator to appropriate lookup list.                              * 
*                                                                       *
*************************************************************************
*/

template<int DIM> void Geometry<DIM>::addSpatialCoarsenOperator(
   tbox::Pointer< CoarsenOperator<DIM> > coarsen_op)
{
   d_coarsen_operators.addItem(coarsen_op);
}

template<int DIM> void Geometry<DIM>::addSpatialRefineOperator(
   tbox::Pointer< RefineOperator<DIM> > refine_op)
{
   d_refine_operators.addItem(refine_op);
}

template<int DIM> void Geometry<DIM>::addTimeInterpolateOperator(
   tbox::Pointer< TimeInterpolateOperator<DIM> > time_op)
{
   d_time_operators.addItem(time_op);
}

/*
*************************************************************************
*                                                                       *
* Search operator lists for operator matching request.                  *
*                                                                       *
*************************************************************************
*/

template<int DIM> tbox::Pointer< CoarsenOperator<DIM> >
Geometry<DIM>::lookupCoarsenOperator(
   const tbox::Pointer< hier::Variable<DIM> >& var,
   const std::string& op_name) const
{
   tbox::Pointer< CoarsenOperator<DIM> > coarsen_op = NULL;
   bool found_op = false;

   if ( (op_name == "NO_COARSEN") || 
        (op_name == "USER_DEFINED_COARSEN") ||
        (op_name.empty()) ) { 
      found_op = true;
   } else {

      typename tbox::List< tbox::Pointer< CoarsenOperator<DIM> > >::Iterator
         lop = d_coarsen_operators.listStart();

      while( coarsen_op.isNull() && lop ) {
         if ( lop()->findCoarsenOperator(var, op_name) ) {
            found_op = true;
            coarsen_op = lop();
         }
         lop++;
      }
   }

   if ( !found_op ) {
      TBOX_ERROR("Coarsen operator '" << op_name <<
                 "' not found for variable '" << var->getName());
   }

   return(coarsen_op);
}

template<int DIM> tbox::Pointer< RefineOperator<DIM> >
Geometry<DIM>::lookupRefineOperator(
   const tbox::Pointer< hier::Variable<DIM> >& var,
   const std::string& op_name) const
{
   tbox::Pointer< RefineOperator<DIM> > refine_op = NULL;
   bool found_op = false;

   if ( (op_name == "NO_REFINE") ||
        (op_name == "USER_DEFINED_REFINE") ||
        (op_name.empty()) ) {
      found_op = true;
   } else {

      typename tbox::List< tbox::Pointer< RefineOperator<DIM> > >::Iterator
         lop = d_refine_operators.listStart();

      while( refine_op.isNull() && lop ) {
         if ( lop()->findRefineOperator(var, op_name) ) {
            found_op = true;
            refine_op = lop();
         }
         lop++;
      }
   }

   if ( !found_op ) {
      TBOX_ERROR("Refine operator '" << op_name <<
                 "' not found for variable '" << var->getName());
   }

   return(refine_op);
}

template<int DIM> tbox::Pointer< TimeInterpolateOperator<DIM> >
Geometry<DIM>::lookupTimeInterpolateOperator(
   const tbox::Pointer< hier::Variable<DIM> >& var,
   const std::string& op_name) const
{
   tbox::Pointer< TimeInterpolateOperator<DIM> > time_op = NULL;
   bool found_op = false;

   if ( (op_name == "NO_TIME_INTERPOLATE") ||
        (op_name.empty()) ) {
      found_op = true;
   } else {

      typename tbox::List< tbox::Pointer< TimeInterpolateOperator<DIM> > >::Iterator
         lop = d_time_operators.listStart();

      while( time_op.isNull() && lop ) {
         if ( lop()->findTimeInterpolateOperator(var, op_name) ) {
            found_op = true;
            time_op = lop();
         }
         lop++;
      }
   }

   if ( !found_op ) {
      TBOX_ERROR("Time interpolation operator '" << op_name <<
                 "' not found for variable '" << var->getName());
   }

   return(time_op);
}

/*
*************************************************************************
*                                                                       *
* Print CartesianGridGeometry class data.                               *
*                                                                       *
*************************************************************************
*/

template<int DIM> void Geometry<DIM>::printClassData(std::ostream& os) const
{
   os << "printing Geometry<DIM> data..." << std::endl;
   os << "Geometry<DIM>: this = " << (Geometry<DIM>*)this << std::endl;

   os << "Coarsen operator list: " << std::endl;
   typename tbox::List< tbox::Pointer< CoarsenOperator<DIM> > >::Iterator
      cop = d_coarsen_operators.listStart();
   while( cop ) {
      os << (CoarsenOperator<DIM>*) cop() << std::endl;
      cop++;
   }

   os << "Refine operator list: " << std::endl;
   typename tbox::List< tbox::Pointer< RefineOperator<DIM> > >::Iterator
      rop = d_refine_operators.listStart();
   while( rop ) {
      os << (RefineOperator<DIM>*) rop() << std::endl;
      rop++;
   }

   os << "Time interpolate operator list: " << std::endl;
   typename tbox::List< tbox::Pointer< TimeInterpolateOperator<DIM> > >::Iterator
      top = d_time_operators.listStart();
   while( top ) {
      os << (TimeInterpolateOperator<DIM>*) top() << std::endl;
      top++;
   }

   hier::GridGeometry<DIM>::printClassData(os);
}

}
}
#endif
