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
bool PatchLevelEnhancedFillPattern::needsToCommunicateDestinationFillBoxes()
const
{
   return true;
}

SAMRAI_INLINE_KEYWORD
bool PatchLevelEnhancedFillPattern::doesSourceLevelCommunicateToDestination()
const
{
   return false;
}

SAMRAI_INLINE_KEYWORD
bool PatchLevelEnhancedFillPattern::fillingCoarseFineGhosts() const
{
   return true;
}

SAMRAI_INLINE_KEYWORD
bool PatchLevelEnhancedFillPattern::fillingEnhancedConnectivityOnly() const
{
   return true;
}

SAMRAI_INLINE_KEYWORD
int PatchLevelEnhancedFillPattern::getMaxFillBoxes() const
{
   return d_max_fill_boxes;
}

}
}
