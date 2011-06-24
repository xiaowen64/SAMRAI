/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Periodic shift identifier in periodic domain.
 *
 ************************************************************************/

#ifndef included_hier_PeriodicId_C
#define included_hier_PeriodicId_C

#include "SAMRAI/hier/PeriodicId.h"

#include <iostream>

#ifndef SAMRAI_INLINE
#include "SAMRAI/hier/PeriodicId.I"
#endif

namespace SAMRAI {
namespace hier {

const PeriodicId PeriodicId::s_invalid_id(-1);
const PeriodicId PeriodicId::s_zero_id(0);



/*
********************************************************************************
********************************************************************************
*/
std::ostream &operator << (
   std::ostream &co,
   const PeriodicId &r)
{
   co << r.d_value;
   return co;
}



}
}
#endif
