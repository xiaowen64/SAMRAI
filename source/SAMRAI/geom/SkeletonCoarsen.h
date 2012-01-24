/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Weighted averaging operator for cell-centered double data on
 *                a Moving mesh.
 *
 ************************************************************************/

#ifndef included_geom_SkeletonCoarsen
#define included_geom_SkeletonCoarsen

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/CoarsenOperator.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"

#include <boost/shared_ptr.hpp>
#include <string>

namespace SAMRAI {
namespace geom {

/**
 * Class SkeletonCoarsen implements a dummy coarsen operator
 * for the skeleton geometry type.  It does nothing but provide
 * basic implementations of the pure virtual functions in the coarsen
 * operator interface.  The findCoarsenOperator() operator function
 * returns true if the the argument std::string is "SKELETON_COARSEN".
 *
 * @see hier::CoarsenOperator
 */

class SkeletonCoarsen:
   public hier::CoarsenOperator
{
public:
   /**
    * Uninteresting default constructor.
    */
   explicit SkeletonCoarsen(
      const tbox::Dimension& dim);

   /**
    * Uninteresting virtual destructor.
    */
   virtual ~SkeletonCoarsen();

   /**
    * Return true if the variable and name std::string match cell-centered
    * double weighted averaging; otherwise, return false.
    */
   bool
   findCoarsenOperator(
      const boost::shared_ptr<hier::Variable>& var,
      const std::string& op_name) const;

   /**
    * The priority of cell-centered double weighted averaging is 0.
    * It will be performed before any user-defined coarsen operations.
    */
   int
   getOperatorPriority() const;

   /**
    * The stencil width of the weighted averaging operator is the vector of
    * zeros.  That is, its stencil does not extend outside the fine box.
    */
   hier::IntVector
   getStencilWidth() const;

   /**
    * Coarsen the source component on the fine patch to the destination
    * component on the coarse patch using the cell-centered double weighted
    * averaging operator.  Coarsening is performed on the intersection of
    * the destination patch and the coarse box.  It is assumed that the
    * fine patch contains sufficient data for the stencil width of the
    * coarsening operator.
    */
   void
   coarsen(
      hier::Patch& coarse,
      const hier::Patch& fine,
      const int dst_component,
      const int src_component,
      const hier::Box& coarse_box,
      const hier::IntVector& ratio) const;
};

}
}
#endif
