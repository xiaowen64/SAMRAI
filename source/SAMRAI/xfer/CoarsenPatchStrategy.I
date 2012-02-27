/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Strategy interface to user routines for coarsening AMR data.
 *
 ************************************************************************/

namespace SAMRAI {
namespace xfer {

/*
 *************************************************************************
 * Register this in the static registry.
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
CoarsenPatchStrategy::registerObject()
{
   std::set<CoarsenPatchStrategy *>& current_objects =
      CoarsenPatchStrategy::getCurrentObjects();
   current_objects.insert(this);
}

/*
 *************************************************************************
 * Return the static registry.
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
std::set<CoarsenPatchStrategy *>&
CoarsenPatchStrategy::getCurrentObjects()
{
   static std::set<CoarsenPatchStrategy *> current_objects;
   return current_objects;
}

SAMRAI_INLINE_KEYWORD
const tbox::Dimension&
CoarsenPatchStrategy::getDim() const
{
   return d_dim;
}

}
}
