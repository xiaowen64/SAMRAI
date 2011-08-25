/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Constant refine operator for cell-centered double data on
 *                a Moving mesh.
 *
 ************************************************************************/

#ifndef included_geom_SkeletonRefine
#define included_geom_SkeletonRefine

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/RefineOperator.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/tbox/Pointer.h"

#include <string>

namespace SAMRAI {
namespace geom {

/**
 * Class SkeletonRefine implements a dummy refine operator
 * for the skeleton geometry type.  It does nothing but provide
 * basic implementations of the pure virtual functions in the refine
 * operator interface.  The findCoarsenOperator() operator function
 * returns true if the the argument std::string is "SKELETON_REFINE".
 *
 * @see hier::RefineOperator
 */

class SkeletonRefine:
   public hier::RefineOperator
{
public:
   /**
    * Uninteresting default constructor.
    */
   explicit SkeletonRefine(
      const tbox::Dimension& dim);

   /**
    * Uninteresting virtual destructor.
    */
   virtual ~SkeletonRefine();

   /**
    * Return true if the variable and name std::string match cell-centered
    * double constant interpolation; otherwise, return false.
    */
   bool
   findRefineOperator(
      const tbox::Pointer<hier::Variable>& var,
      const std::string& op_name) const;

   /**
    * The priority of cell-centered double constant interpolation is 0.
    * It will be performed before any user-defined interpolation operations.
    */
   int
   getOperatorPriority() const;

   /**
    * The stencil width of the constant interpolation operator is the vector
    * of zeros.  That is, its stencil does not extend outside the fine box.
    */
   hier::IntVector
   getStencilWidth() const;

   /**
    * No-op implementation of refine()
    */
   void
   refine(
      hier::Patch& fine,
      const hier::Patch& coarse,
      const int dst_component,
      const int src_component,
      const hier::BoxOverlap& fine_overlap,
      const hier::IntVector& ratio) const;
};

}
}
#endif
