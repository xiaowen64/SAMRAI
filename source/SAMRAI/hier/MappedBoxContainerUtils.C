/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Common MappedBox operations for MappedBox containers. 
 *
 ************************************************************************/
#ifndef included_hier_MappedBoxContainerUtils_C
#define included_hier_MappedBoxContainerUtils_C

#include "SAMRAI/hier/MappedBoxContainerUtils.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"
#include "SAMRAI/tbox/TimerManager.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/hier/MappedBoxContainerUtils.I"
#endif

namespace SAMRAI {
namespace hier {

/*
 * Constructor does nothing because the objects are stateless.
 */

MappedBoxContainerUtils::MappedBoxContainerUtils() {
}

/*
 ***************************************************************************
 ***************************************************************************
 */

void MappedBoxContainerUtils::recursivePrintMappedBoxVector(
   const std::vector<MappedBox>& mapped_boxes,
   std::ostream& os,
   const std::string& border,
   int detail_depth)
{
   (void)detail_depth;
   os << border;
   for (std::vector<MappedBox>::const_iterator ni = mapped_boxes.begin();
        ni != mapped_boxes.end();
        ++ni) {
      os << "  " << *ni;
   }
}

void MappedBoxContainerUtils::convertBoxListToMappedBoxVector(
   const BoxList& box_list,
   std::vector<MappedBox>& mapped_box_vector,
   const BlockId& block_id)
{
   mapped_box_vector.reserve(mapped_box_vector.size() + box_list.size());
   LocalId last_used_id(-1);
   for (BoxList::Iterator bi(box_list); bi; bi++) {
      mapped_box_vector.push_back(MappedBox(*bi, ++last_used_id, 0, block_id));
   }
}

void MappedBoxContainerUtils::convertMappedBoxVectorToBoxList(
   const std::vector<MappedBox>& mapped_box_vector,
   BoxList& box_list)
{
   for (std::vector<MappedBox>::const_iterator ni = mapped_box_vector.begin();
        ni != mapped_box_vector.end(); ++ni) {
      box_list.appendItem(ni->getBox());
   }
}

}
}
#endif
