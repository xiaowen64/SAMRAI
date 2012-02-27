/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Strategy interface to user routines for refining AMR data.
 *
 ************************************************************************/

namespace SAMRAI {
namespace xfer {

/*
 *************************************************************************
 *
 * Loop over all fill boxes and call the user-defined preprocesses.
 *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
RefinePatchStrategy::preprocessRefineBoxes(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const hier::BoxContainer& fine_boxes,
   const hier::IntVector& ratio)
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(fine, coarse, ratio);

   for (hier::BoxContainer::ConstIterator b(fine_boxes); b != fine_boxes.end(); ++b) {
      preprocessRefine(fine, coarse, b(), ratio);
   }
}

/*
 *************************************************************************
 *
 * Loop over all fill boxes and call the user-defined postprocesses.
 *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
RefinePatchStrategy::postprocessRefineBoxes(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const hier::BoxContainer& fine_boxes,
   const hier::IntVector& ratio)
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS3(d_dim, fine, coarse, ratio);

   for (hier::BoxContainer::ConstIterator b(fine_boxes); b != fine_boxes.end(); ++b) {
      postprocessRefine(fine, coarse, b(), ratio);
   }
}

/*
 *************************************************************************
 * Register this in the static registry.
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
RefinePatchStrategy::registerObject()
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
SAMRAI_INLINE_KEYWORD
void
RefinePatchStrategy::unregisterObject()
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
SAMRAI_INLINE_KEYWORD
std::set<RefinePatchStrategy *>&
RefinePatchStrategy::getCurrentObjects()
{
   static std::set<RefinePatchStrategy *> current_objects;
   return current_objects;
}

/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
const tbox::Dimension&
RefinePatchStrategy::getDim() const
{
   return d_dim;
}

}
}
