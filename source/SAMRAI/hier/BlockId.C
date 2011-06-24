/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Block identifier in multiblock domain.
 *
 ************************************************************************/

#ifndef included_hier_BlockId_C
#define included_hier_BlockId_C

#include "SAMRAI/hier/BlockId.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include <iostream>

#ifndef SAMRAI_INLINE
#include "SAMRAI/hier/BlockId.I"
#endif

namespace SAMRAI {
namespace hier {

const BlockId BlockId::s_invalid_id(tbox::MathUtilities<int>::getMax());
const BlockId BlockId::s_zero_id(0);



/*
********************************************************************************
********************************************************************************
*/
std::ostream &operator << (
   std::ostream &co,
   const BlockId &r)
{
   co << r.d_value;
   return co;
}



}
}
#endif
