/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Factory abstract base class for creating patch data objects 
 *
 ************************************************************************/

#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace hier {

SAMRAI_INLINE_KEYWORD
PatchDataFactory::PatchDataFactory(
   const IntVector& ghosts):
   d_ghosts(ghosts)
{
   TBOX_ASSERT(ghosts.min() >= 0);
}

SAMRAI_INLINE_KEYWORD
const tbox::Dimension& PatchDataFactory::getDim() const
{
   return d_ghosts.getDim();
}

}
}
