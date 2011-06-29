/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Base class that describes intersections between AMR boxes 
 *
 ************************************************************************/

#ifndef included_hier_BoxOverlap_C
#define included_hier_BoxOverlap_C

#include "SAMRAI/hier/BoxOverlap.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/hier/BoxOverlap.I"
#endif

namespace SAMRAI {
namespace hier {

BoxOverlap::~BoxOverlap()
{
}

void BoxOverlap::print(
   std::ostream& os) const
{
   os << "print() method not implemented for this overlap type" << std::endl;
}

}
}
#endif
