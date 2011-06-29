/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Globally unique identifier that can be locally determined.
 *
 ************************************************************************/
#ifndef included_hier_GlobalId_C
#define included_hier_GlobalId_C

#include "SAMRAI/hier/GlobalId.h"

#include <iostream>

#ifndef SAMRAI_INLINE
#include "SAMRAI/hier/GlobalId.I"
#endif

namespace SAMRAI {
namespace hier {

std::ostream& operator << (
   std::ostream& co,
   const GlobalId& r)
{
   co << r.d_owner_rank << '#' << r.d_local_id;
   return co;
}

}
}
#endif
