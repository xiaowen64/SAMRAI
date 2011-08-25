/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Factory abstract base class for creating patch data objects
 *
 ************************************************************************/

#ifndef included_hier_PatchDataFactory_C
#define included_hier_PatchDataFactory_C

#include "SAMRAI/hier/PatchDataFactory.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/hier/PatchDataFactory.I"
#endif
namespace SAMRAI {
namespace hier {

PatchDataFactory::~PatchDataFactory()
{
}

const hier::IntVector&
PatchDataFactory::getGhostCellWidth() const
{
   return d_ghosts;
}

/**********************************************************************
* Default implementation                                             *
**********************************************************************/

MultiblockDataTranslator *
PatchDataFactory::getMultiblockDataTranslator()
{
   return (MultiblockDataTranslator *)NULL;
}

}
}
#endif
