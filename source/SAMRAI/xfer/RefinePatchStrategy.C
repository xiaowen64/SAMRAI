/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Strategy interface to user routines for refining AMR data.
 *
 ************************************************************************/

#ifndef included_xfer_RefinePatchStrategy_C
#define included_xfer_RefinePatchStrategy_C

#include "SAMRAI/xfer/RefinePatchStrategy.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/BoxContainerConstIterator.h"

namespace SAMRAI {
namespace xfer {

/*
 *************************************************************************
 *
 * The default constructor and virtual destructor do nothing
 * particularly interesting.
 *
 *************************************************************************
 */

RefinePatchStrategy::RefinePatchStrategy(
   const tbox::Dimension& dim):
   d_dim(dim)
{
   registerObject();
}

RefinePatchStrategy::~RefinePatchStrategy()
{
   unregisterObject();
}

/*
 *************************************************************************
 *
 * Loop over all fill boxes and call the user-defined preprocesses.
 *
 *************************************************************************
 */

void RefinePatchStrategy::preprocessRefineBoxes(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const hier::BoxContainer& fine_boxes,
   const hier::IntVector& ratio)
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(fine, coarse, ratio);

   for (hier::BoxContainer::ConstIterator b(fine_boxes); b != fine_boxes.end(); ++b) {
      this->preprocessRefine(fine, coarse, b(), ratio);
   }
}

/*
 *************************************************************************
 *
 * Loop over all fill boxes and call the user-defined postprocesses.
 *
 *************************************************************************
 */

void RefinePatchStrategy::postprocessRefineBoxes(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const hier::BoxContainer& fine_boxes,
   const hier::IntVector& ratio)
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS3(d_dim, fine, coarse, ratio);

   for (hier::BoxContainer::ConstIterator b(fine_boxes); b != fine_boxes.end(); ++b) {
      this->postprocessRefine(fine, coarse, b(), ratio);
   }
}

/*
 *************************************************************************
 * Register this in the static registry.
 *************************************************************************
 */
void RefinePatchStrategy::registerObject()
{
   std::set<RefinePatchStrategy *>& current_objects =
      RefinePatchStrategy::getCurrentObjects();
   TBOX_DIM_ASSERT_CHECK_DIM(d_dim);
   current_objects.insert(this);
}

/*
 *************************************************************************
 * Unregister this from the static registry.
 *************************************************************************
 */
void RefinePatchStrategy::unregisterObject()
{
   std::set<RefinePatchStrategy *>& current_objects =
      RefinePatchStrategy::getCurrentObjects();
   TBOX_DIM_ASSERT_CHECK_DIM(d_dim);
   current_objects.erase(this);
}

/*
 *************************************************************************
 * Return the static registry.
 *************************************************************************
 */
std::set<RefinePatchStrategy *>& RefinePatchStrategy::getCurrentObjects()
{
   static std::set<RefinePatchStrategy *> current_objects;
   return current_objects;
}

/*
 *************************************************************************
 *************************************************************************
 */
const tbox::Dimension& RefinePatchStrategy::getDim() const
{
   return d_dim;
}

/*
 *************************************************************************
 * Compute the max refine stencil width from all constructed
 * refine patch strategies.
 *************************************************************************
 */
hier::IntVector
RefinePatchStrategy::getMaxRefineOpStencilWidth(
   const tbox::Dimension& dim)
{
   hier::IntVector max_width(dim, 0);

   std::set<RefinePatchStrategy *>& current_objects =
      RefinePatchStrategy::getCurrentObjects();
   for (std::set<RefinePatchStrategy *>::const_iterator
        si = current_objects.begin(); si != current_objects.end(); ++si) {
      const RefinePatchStrategy* op = *si;
      if (op->getDim() == dim) {
         max_width.max(op->getRefineOpStencilWidth());
      }
   }

   return max_width;
}

}
}
#endif
