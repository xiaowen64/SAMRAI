//
// File:	CoarsenOperator.h
// Package:	SAMRAI transfer
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Abstract base class for spatial coarsening operators.
//

#ifndef included_xfer_CoarsenOperator
#define included_xfer_CoarsenOperator

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
#ifndef included_hier_Variable
#include "Variable.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif
#ifndef included_String
#include <string>
using namespace std;
#define included_String
#endif

namespace SAMRAI {
    namespace xfer {

/**
 * Class CoarsenOperator<DIM> is an abstract base class for each
 * spatial coarsening operator used in the SAMRAI framework.  This class
 * defines the interface between numerical coarsening routines and the
 * rest of the framework.  Each concrete coarsening operator subclass
 * must provide four operations:
 * 


 * - \b (1) {an implementation of the coarsening operation
 *            appropriate for its corresponding patch data type.}
 * - \b (2) {a function that determines whether or not the operator
 *             matches an arbitrary request for a coarsening operator.}
 * - \b (3) {a function that returns the stencil width of the operator
 *             (i.e., the number of ghost cells needed by the operator).}
 * - \b (4) {a function that returns an integer stating the priority of the
 *             operator with respect to other coarsening operators.}
 * 


 * To add a new coarsening operator (either for a new patch data type
 * or for a new time coarsening routine on an existing type), define
 * the operator by inheriting from this abstract base class.  The operator
 * subclass must implement the coarsening operation in the coarsen()
 * function, and provide a response to a general operator request in the
 * findCoarsenOperator() function.  The stencil width and operator priority
 * must be returned from the getStencilWidth() and getOperatorPriority()
 * functions, respectively.  Then, the new operator must be added to the
 * operator list for the appropriate transfer geometry object using the
 * Geometry<DIM>::addSpatialCOarsenOperator() function.
 *
 * Since spatial coarsening operators usually depend on patch data centering
 * and data type as well as the mesh coordinate system, they are defined
 * in the {\it geometry} package.
 *
 * @see xfer::Geometry
 */

template<int DIM> class CoarsenOperator : public tbox::DescribedClass
{
public:
   /**
    * The default constructor for the coarsening operator does
    * nothing interesting.
    */
   CoarsenOperator();

   /**
    * The virtual destructor for the coarsening operator does
    * nothing interesting.
    */
   virtual ~CoarsenOperator<DIM>();

   /**
    * Return true if the coarsening operation matches the variable and
    * name string identifier request; false, otherwise.
    */
   virtual bool findCoarsenOperator(
      const tbox::Pointer< hier::Variable<DIM> >& var,
      const string &op_name) const = 0;

   /**
    * Return name string identifier of the coarsening operation.
    */
   virtual const string& getOperatorName() const = 0;

   /**
    * Return the priority of this operator relative to other coarsening
    * operators.  The SAMRAI transfer routines guarantee that coarsening
    * using operators with lower priority will be performed before those
    * with higher priority.
    */
   virtual int getOperatorPriority() const = 0;

   /**
    * Return the stencil width associated with the coarsening operator.
    * The SAMRAI transfer routines guarantee that the source patch will
    * contain sufficient ghost cell data surrounding the interior to
    * satisfy the stencil width requirements for each coarsening operator.
    */
   virtual hier::IntVector<DIM> getStencilWidth() const = 0;

   /**
    * Coarsen the source component on the fine patch to the destination
    * component on the coarse patch. The coarsening operation is performed
    * on the intersection of the destination patch and the coarse box.
    * The fine patch is guaranteed to contain sufficient data for the
    * stencil width of the coarsening operator.
    */
   virtual void coarsen(
      hier::Patch<DIM>& coarse,
      const hier::Patch<DIM>& fine,
      const int dst_component,
      const int src_component,
      const hier::Box<DIM>& coarse_box,
      const hier::IntVector<DIM>& ratio) const = 0;

private:
   CoarsenOperator(const CoarsenOperator<DIM>&);	// not implemented
   void operator=(const CoarsenOperator<DIM>&);	// not implemented

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "CoarsenOperator.C"
#endif
