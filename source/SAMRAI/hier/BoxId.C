/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Identifier for a Box.
 *
 ************************************************************************/

#ifndef included_hier_BoxId_C
#define included_hier_BoxId_C

#include "SAMRAI/hier/BoxId.h"

#include <iostream>

#ifndef SAMRAI_INLINE
#include "SAMRAI/hier/BoxId.I"
#endif

namespace SAMRAI {
namespace hier {

/*
 *******************************************************************************
 * Stream-insert operator.
 *******************************************************************************
 */
std::ostream& operator << (
   std::ostream& co,
   const BoxId& r)
{
   co << r.d_global_id.getOwnerRank()
   << ':' << r.d_block_id
   << '#' << r.d_global_id.getLocalId()
   << '/' << r.d_periodic_id;
   return co;
}

}
}
#endif
