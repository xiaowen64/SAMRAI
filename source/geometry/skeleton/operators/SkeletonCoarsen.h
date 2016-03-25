//
// File:	SkeletonCoarsen.h
// Package:	SAMRAI geometry
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 601 $
// Modified:	$Date: 2005-09-06 11:23:15 -0700 (Tue, 06 Sep 2005) $
// Description: Weighted averaging operator for cell-centered double data on 
//              a Moving mesh.
//

#ifndef included_geom_SkeletonCoarsen
#define included_geom_SkeletonCoarsen

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_hier_Box
#include "Box.h"
#endif
#ifndef included_hier_IntVector
#include "IntVector.h"
#endif
#ifndef included_hier_Patch
#include "Patch.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif
#ifndef included_String
#include <string>
using namespace std;
#define included_String
#endif
#ifndef included_xfer_CoarsenOperator
#include "CoarsenOperator.h"
#endif

namespace SAMRAI {
    namespace geom {

/**
 * Class SkeletonCoarsen implements a dummy coarsen operator
 * for the skeleton geometry type.  It does nothing but provide  
 * basic implementations of the pure virtual functions in the coarsen
 * operator interface.  The findCoarsenOperator() operator function 
 * returns true if the the argument string is "SKELETON_COARSEN".
 * 
 * @see xfer::CoarsenOperator
 */

template<int DIM> class SkeletonCoarsen 
: public xfer::CoarsenOperator<DIM>
{
public:
   /**
    * Uninteresting default constructor.
    */
   SkeletonCoarsen();

   /**
    * Uninteresting virtual destructor.
    */
   virtual ~SkeletonCoarsen<DIM>();

   /**
    * Return true if the variable and name string match cell-centered 
    * double weighted averaging; otherwise, return false.
    */
   bool findCoarsenOperator(const tbox::Pointer< hier::Variable<DIM> >& var,
                            const string &op_name) const; 

   /**
    * Return name string identifier of this coarsening operator.
    */
   const string& getOperatorName() const;

   /**
    * The priority of cell-centered double weighted averaging is 0.
    * It will be performed before any user-defined coarsen operations. 
    */
   int getOperatorPriority() const;

   /**
    * The stencil width of the weighted averaging operator is the vector of 
    * zeros.  That is, its stencil does not extend outside the fine box.
    */
   hier::IntVector<DIM> getStencilWidth() const;

   /**
    * Coarsen the source component on the fine patch to the destination
    * component on the coarse patch using the cell-centered double weighted 
    * averaging operator.  Coarsening is performed on the intersection of 
    * the destination patch and the coarse box.  It is assumed that the 
    * fine patch contains sufficient data for the stencil width of the 
    * coarsening operator.
    */
   void coarsen(hier::Patch<DIM>& coarse,
                const hier::Patch<DIM>& fine,
                const int dst_component,
                const int src_component,
                const hier::Box<DIM>& coarse_box,
                const hier::IntVector<DIM>& ratio) const;

private:
   string d_name_id;

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "SkeletonCoarsen.C"
#endif
