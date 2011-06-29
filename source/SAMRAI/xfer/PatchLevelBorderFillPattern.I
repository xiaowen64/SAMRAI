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
bool PatchLevelBorderFillPattern::needsToCommunicateDestinationFillBoxes()
const
{
   return true;
}

SAMRAI_INLINE_KEYWORD
bool PatchLevelBorderFillPattern::doesSourceLevelCommunicateToDestination()
const
{
   return false;
}

SAMRAI_INLINE_KEYWORD
bool PatchLevelBorderFillPattern::fillingCoarseFineGhosts() const
{
   return true;
}

SAMRAI_INLINE_KEYWORD
bool PatchLevelBorderFillPattern::fillingEnhancedConnectivityOnly() const
{
   return false;
}

SAMRAI_INLINE_KEYWORD
int PatchLevelBorderFillPattern::getMaxFillBoxes() const
{
   return d_max_fill_boxes;
}

}
}
