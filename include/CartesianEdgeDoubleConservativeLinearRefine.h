//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/geometry/cartesian/operators/edge/CartesianEdgeDoubleConservativeLinearRefine.h $
// Package:	SAMRAI geometry
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description: Conservative linear refine operator for edge-centered 
//              double data on a Cartesian mesh.
//

#ifndef included_geom_CartesianEdgeDoubleConservativeLinearRefine
#define included_geom_CartesianEdgeDoubleConservativeLinearRefine

#include "SAMRAI_config.h"
#include "Box.h"
#include "IntVector.h"
#include "Patch.h"
#include "tbox/Pointer.h"
#ifndef included_String
#include <string>
#define included_String
#endif
#include "RefineOperator.h"

namespace SAMRAI {
    namespace geom {

/**
 * Class CartesianEdgeDoubleConservativeLinearRefine implements
 * conservative linear interpolation for edge-centered double patch data 
 * defined over a Cartesian mesh.  It is derived from the base class 
 * xfer::RefineOperator<DIM>.  The numerical operations for the interpolation 
 * use FORTRAN numerical routines.
 *
 * The findRefineOperator() operator function returns true if the input 
 * variable is edge-centered double, and the string is 
 * "CONSERVATIVE_LINEAR_REFINE".
 * 
 * @see xfer::RefineOperator
 */

template<int DIM> class CartesianEdgeDoubleConservativeLinearRefine 
: public xfer::RefineOperator<DIM>
{
public:
   /**
    * Uninteresting default constructor.
    */
   CartesianEdgeDoubleConservativeLinearRefine();

   /**
    * Uninteresting virtual destructor.
    */
   virtual ~CartesianEdgeDoubleConservativeLinearRefine<DIM>();

   /**
    * Return true if the variable and name string match edge-centered 
    * double conservative linear interpolation; otherwise, return false.
    */
   bool findRefineOperator(const tbox::Pointer< hier::Variable<DIM> >& var,
                           const std::string &op_name) const; 

   /**
    * Return name string identifier of this refinement operation.
    */
   const std::string& getOperatorName() const;

   /**
    * The priority of edge-centered double conservative linear is 0.
    * It will be performed before any user-defined interpolation operations. 
    */
   int getOperatorPriority() const;

   /**
    * The stencil width of the conservative linear interpolation operator is 
    * the vector of ones.
    */
   hier::IntVector<DIM> getStencilWidth() const;

   /**
    * Refine the source component on the coarse patch to the destination
    * component on the fine patch using the edge-centered double conservative
    * linear interpolation operator.  Interpolation is performed on the 
    * intersection of the destination patch and the fine box.   It is assumed 
    * that the coarse patch contains sufficient data for the stencil width of 
    * the refinement operator.
    */
   void refine(hier::Patch<DIM>& fine,
               const hier::Patch<DIM>& coarse,
               const int dst_component,
               const int src_component,
               const hier::Box<DIM>& fine_box,
               const hier::IntVector<DIM>& ratio) const;

private:
   std::string d_name_id;

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "CartesianEdgeDoubleConservativeLinearRefine.C"
#endif
