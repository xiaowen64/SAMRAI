/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Abstract fill pattern class to provide interface for stencils
 *
 ************************************************************************/

namespace SAMRAI {
namespace xfer {

SAMRAI_INLINE_KEYWORD
bool
PatchLevelInteriorFillPattern::needsToCommunicateDestinationFillBoxes() const
{
   return false;
}

SAMRAI_INLINE_KEYWORD
bool
PatchLevelInteriorFillPattern::doesSourceLevelCommunicateToDestination() const
{
   return true;
}

SAMRAI_INLINE_KEYWORD
bool
PatchLevelInteriorFillPattern::fillingCoarseFineGhosts() const
{
   return false;
}

SAMRAI_INLINE_KEYWORD
bool
PatchLevelInteriorFillPattern::fillingEnhancedConnectivityOnly() const
{
   return false;
}

SAMRAI_INLINE_KEYWORD
int
PatchLevelInteriorFillPattern::getMaxFillBoxes() const
{
   return d_max_fill_boxes;
}

}
}
