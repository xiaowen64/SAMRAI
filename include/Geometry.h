// 
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/transfer/operators/Geometry.h $
// Package:	SAMRAI transfer 
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description: Base class for interface between transfer ops and geometry.
//

#ifndef included_xfer_Geometry
#define included_xfer_Geometry

#include "SAMRAI_config.h"
#include "GridGeometry.h"
#include "Variable.h"
#include "tbox/List.h"
#include "tbox/Pointer.h"
#ifndef included_String
#include <string>
#define included_String
#endif
#include "CoarsenOperator.h"
#include "RefineOperator.h"
#include "TimeInterpolateOperator.h"

namespace SAMRAI { 
    namespace xfer {

/**
 * Class Geometry is the base class for SAMRAI geometry classes;
 * it is derived from hier::GridGeometry and is intended to serve as
 * a base class for geometry classes that define specific coordinate
 * system operations for an AMR grid hierarchy.  This transfer geometry class
 * provides a lookup mechanism to search for time interpolation and spatial
 * coarsening/refining operators.  That is, algorithms and applications
 * that manage communication on an AMR hierarchy may query the transfer
 * geometry object for operators that may be applied to specific variables.
 * Typically, the operators are assigned to the transfer geometry object in
 * the constructor of the geometry object that defines the mesh coordinate
 * system (and which is derived from this transfer geometry class).  Additional
 * operators may be added to a tranfer geometry object at any time during
 * program execution. However, each operator must be added BEFORE it is
 * requested or an unrecoverable exception will be thrown and the program will
 * abort.  Also note that operators are added to the heads of the operator
 * lists so that the most recently added operator will be returned if more
 * than one operator satisfies a given request.  See the time interpolation,
 * spatial coarsening, and spatial refinement operator base classes for
 * more information about adding new operators for either new patch data
 * types or new operators for pre-existing patch data types.
 * 
 * @see hier::GridGeometry
 * @see xfer::RefineOperator
 * @see xfer::CoarsenOperator
 * @see xfer::TimeInterpolateOperator
 */

template<int DIM> class Geometry : public hier::GridGeometry<DIM>
{
public:
   /**
    * Constructor for Geometry class just passes the object_name 
    * to the hier::GridGeometry parent class.
    */
   Geometry(const std::string &object_name);

   /**
    * The virtual destructor for the geometry base class does
    * nothing interesting.
    */
   virtual ~Geometry<DIM>();

   /**
    * Add concrete spatial coarsening operator instance to appropriate
    * lookup list.  Note that each concrete operator must implement a
    * lookup function through which it can be identified.
    */
   virtual void addSpatialCoarsenOperator(
      tbox::Pointer< CoarsenOperator<DIM> > coarsen_op);

   /**
    * Add concrete spatial refinement operator instance to appropriate
    * lookup list.  Note that each concrete operator must implement a
    * lookup function through which it can be identified.
    */
   virtual void addSpatialRefineOperator(
      tbox::Pointer< RefineOperator<DIM> > refine_op);

   /**
    * Add concrete time interpolation operator instance to appropriate
    * lookup list.  Note that each concrete operator must implement a
    * lookup function through which it can be identified.
    */
   virtual void addTimeInterpolateOperator(
      tbox::Pointer< TimeInterpolateOperator<DIM> > time_op);

   /**
    * Search list for the spatial coarsening operator matching the
    * request for the given variable.  If the operator is found, a
    * pointer to it will be returned.  Otherwise, an unrecoverable
    * error will result and the program will abort.
    */
   virtual tbox::Pointer< CoarsenOperator<DIM> >
   lookupCoarsenOperator(const tbox::Pointer< hier::Variable<DIM> >& var,
                         const std::string& op_name) const;

   /**
    * Search list for the spatial refinement operator matching the
    * request for the given variable.  If the operator is found, a
    * pointer to it will be returned.  Otherwise, an unrecoverable
    * error will result and the program will abort.
    */
   virtual tbox::Pointer< RefineOperator<DIM> >
   lookupRefineOperator(const tbox::Pointer< hier::Variable<DIM> >& var,
                        const std::string& op_name) const;

   /**
    * Search list for the time interpolation operator matching the
    * request for the given variable.  If the operator is found, a
    * pointer to it will be returned.  Otherwise, an unrecoverable
    * error will result and the program will abort.
    */
   virtual tbox::Pointer< TimeInterpolateOperator<DIM> >
   lookupTimeInterpolateOperator(const tbox::Pointer< hier::Variable<DIM> >& var,
                                 const std::string& op_name
                                       = "STD_LINEAR_TIME_INTERPOLATE") const;

   /**
    * Print class data representation.
    */
   virtual void printClassData(std::ostream& os) const; 

private:
   /*
    * The list of spatial coarsening operators is maintained to lookup
    * operators for specific variables as requested by algorithms and/or
    * applications using this transfer geometry object.  Standard concrete
    * coarsening operators can be found in the patchdata package.
    * Additional operators may be added to this list at any time
    * (see addSpatialCoarsenOperator() function).
    */
   tbox::List< tbox::Pointer< CoarsenOperator<DIM> > > d_coarsen_operators;

   /*
    * The list of spatial refinement operators is maintained to lookup
    * operators for specific variables as requested by algorithms and/or
    * applications using this transfer geometry object.  Standard concrete
    * refinement operators can be found in the patchdata package.
    * Additional operators may be added to this list at any time
    * (see addSpatialRefineOperator() function).
    */
   tbox::List< tbox::Pointer< RefineOperator<DIM> > > d_refine_operators;

   /*
    * The list of time interpolation operators is maintained to lookup
    * operators for specific variables as requested by algorithms and/or
    * applications using this transfer geometry object.  Standard concrete
    * time interpolation operators can be found in the patchdata package.
    * Additional operators may be added to this list at any time
    * (see addTimeInterpolateOperator() function).
    */
   tbox::List< tbox::Pointer< TimeInterpolateOperator<DIM> > > d_time_operators;

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "Geometry.C"
#endif
