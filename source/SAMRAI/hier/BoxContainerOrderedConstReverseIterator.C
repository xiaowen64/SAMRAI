/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   An iterator over containers of boxes.
 *
 ************************************************************************/

#ifndef included_hier_BoxContainerOrderedConstReverseIterator_C
#define included_hier_BoxContainerOrderedConstReverseIterator_C

#include "SAMRAI/hier/BoxContainerOrderedConstReverseIterator.h"
#include "SAMRAI/hier/BoxContainer.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/hier/BoxContainerOrderedConstReverseIterator.I"
#endif

namespace SAMRAI {
namespace hier {

BoxContainerOrderedConstReverseIterator::BoxContainerOrderedConstReverseIterator(
   const BoxContainer& container,
   bool from_back):
   d_box_container(&container),
   d_set_iter(from_back ? container.d_set.rbegin() :
              container.d_set.rend())
{
}

BoxContainerOrderedConstReverseIterator& BoxContainerOrderedConstReverseIterator::operator = (
   const BoxContainerOrderedConstReverseIterator& rhs)
{
   if (this != &rhs) {
      d_box_container = rhs.d_box_container;
      d_set_iter = rhs.d_set_iter;
   }
   return *this;
}

}
}

#endif
