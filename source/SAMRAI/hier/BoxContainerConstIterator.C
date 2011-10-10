/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   A const iterator over containers of boxes.
 *
 ************************************************************************/

#ifndef included_hier_BoxContainerConstIterator_C
#define included_hier_BoxContainerConstIterator_C

#include "SAMRAI/hier/BoxContainerConstIterator.h"
#include "SAMRAI/hier/BoxContainer.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/hier/BoxContainerConstIterator.I"
#endif

namespace SAMRAI {
namespace hier {

BoxContainerConstIterator::BoxContainerConstIterator(
   const BoxContainer& container,
   bool from_start):
   d_box_container(&container),
   d_list_iter(from_start ? container.d_list.begin() :
               container.d_list.end()),
   d_set_iter(from_start ? container.d_set.begin() :
              container.d_set.end()),
   d_ordered(container.d_ordered)
{
}

BoxContainerConstIterator&
BoxContainerConstIterator::operator = (
   const BoxContainerConstIterator& rhs)
{
   if (this != &rhs) {
      d_box_container = rhs.d_box_container;
      d_list_iter = rhs.d_list_iter;
      d_set_iter = rhs.d_set_iter;
      d_ordered = rhs.d_ordered;
   }
   return *this;
}

}
}

#endif
