/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Common MappedBox operations for MappedBox containers. 
 *
 ************************************************************************/
#include "SAMRAI/hier/MappedBox.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace hier {

SAMRAI_INLINE_KEYWORD
void MappedBoxContainerUtils::growMappedBoxVectorBoxes(
   std::vector<MappedBox>& mapped_box_vector,
   const IntVector& growth)
{
   for (std::vector<MappedBox>::iterator ni = mapped_box_vector.begin();
      ni != mapped_box_vector.end(); ++ni) {
      ni->getBox().grow(growth);
   }
}

SAMRAI_INLINE_KEYWORD
void MappedBoxContainerUtils::refineMappedBoxVectorBoxes(
   std::vector<MappedBox>& mapped_box_vector,
   const IntVector& ratio)
{
   for (std::vector<MappedBox>::iterator ni = mapped_box_vector.begin();
      ni != mapped_box_vector.end(); ++ni) {
      ni->getBox().refine(ratio);
   }
}

SAMRAI_INLINE_KEYWORD
void MappedBoxContainerUtils::coarsenMappedBoxVectorBoxes(
   std::vector<MappedBox>& mapped_box_vector,
   const IntVector& ratio)
{
   for (std::vector<MappedBox>::iterator ni = mapped_box_vector.begin();
      ni != mapped_box_vector.end(); ++ni) {
      ni->getBox().coarsen(ratio);
   }
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif
