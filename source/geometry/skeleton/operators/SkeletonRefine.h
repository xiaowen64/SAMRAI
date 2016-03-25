//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/geometry/skeleton/operators/SkeletonRefine.h $
// Package:	SAMRAI geometry
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description: Constant refine operator for cell-centered double data on 
//              a Moving mesh.
//

#ifndef included_geom_SkeletonRefine
#define included_geom_SkeletonRefine

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
 * Class SkeletonRefine implements a dummy refine operator
 * for the skeleton geometry type.  It does nothing but provide   
 * basic implementations of the pure virtual functions in the refine
 * operator interface.  The findCoarsenOperator() operator function 
 * returns true if the the argument string is "SKELETON_REFINE".
 * 
 * @see xfer::RefineOperator
 */

template<int DIM> class SkeletonRefine 
: public xfer::RefineOperator<DIM>
{
public:
   /**
    * Uninteresting default constructor.
    */
   SkeletonRefine();

   /**
    * Uninteresting virtual destructor.
    */
   virtual ~SkeletonRefine<DIM>();

   /**
    * Return true if the variable and name string match cell-centered 
    * double constant interpolation; otherwise, return false.
    */
   bool findRefineOperator(const tbox::Pointer< hier::Variable<DIM> >& var,
                           const std::string &op_name) const; 

   /**
    * Return name string identifier of this refinement operator.
    */
   const std::string& getOperatorName() const;

   /**
    * The priority of cell-centered double constant interpolation is 0.
    * It will be performed before any user-defined interpolation operations. 
    */
   int getOperatorPriority() const;

   /**
    * The stencil width of the constant interpolation operator is the vector 
    * of zeros.  That is, its stencil does not extend outside the fine box.
    */
   hier::IntVector<DIM> getStencilWidth() const;

   /**
    * Refine the source component on the coarse patch to the destination
    * component on the fine patch using the cell-centered double constant
    * interpolation operator.  Interpolation is performed on the intersection 
    * of the destination patch and the fine box.   It is assumed that the
    * coarse patch contains sufficient data for the stencil width of the 
    * refinement operator.
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
#include "SkeletonRefine.C"
#endif
