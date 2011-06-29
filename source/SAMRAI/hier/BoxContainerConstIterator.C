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
   d_list_iter(from_start ? container.begin().d_list_iter :
                            container.end().d_list_iter)
{
}

BoxContainerConstIterator&
BoxContainerConstIterator::operator = (
   const BoxContainerConstIterator& rhs)
{
   if (this != &rhs) {
      d_box_container = rhs.d_box_container;
      d_list_iter = rhs.d_list_iter;
   }
   return *this;
}

}
}

#endif
