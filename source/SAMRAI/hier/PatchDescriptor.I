/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Factory class for patch data objects that live on a patch
 *
 ************************************************************************/

namespace SAMRAI {
namespace hier {

SAMRAI_INLINE_KEYWORD
boost::shared_ptr<PatchDataFactory>
PatchDescriptor::getPatchDataFactory(
   const int id) const
{
   TBOX_ASSERT((id >= 0) && (id < d_max_number_registered_components));

   return d_factories[id];
}

SAMRAI_INLINE_KEYWORD
int
PatchDescriptor::getMaxNumberRegisteredComponents() const
{
   return d_max_number_registered_components;
}

SAMRAI_INLINE_KEYWORD
const std::string&
PatchDescriptor::mapIndexToName(
   const int id) const
{
   TBOX_ASSERT((id >= 0) && (id < d_max_number_registered_components));

   return d_names[id];
}

/*
 *************************************************************************
 * Set the mininum value to be returned by getMaxGhostWidth().
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
PatchDescriptor::setMinGhostWidth(
   const IntVector& min_value)
{
   d_min_gcw[min_value.getDim().getValue() - 1] = min_value;
}

}
}
