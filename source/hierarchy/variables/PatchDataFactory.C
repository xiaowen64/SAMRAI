//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/hierarchy/variables/PatchDataFactory.C $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1776 $
// Modified:	$LastChangedDate: 2007-12-13 16:40:01 -0800 (Thu, 13 Dec 2007) $
// Description:	Factory abstract base class for creating patch data objects
//

#ifndef included_hier_PatchDataFactory_C
#define included_hier_PatchDataFactory_C

#include "PatchDataFactory.h"

#ifdef DEBUG_NO_INLINE
#include "PatchDataFactory.I"
#endif
namespace SAMRAI {
    namespace hier {

template<int DIM>  PatchDataFactory<DIM>::~PatchDataFactory()
{
}

template<int DIM>
const hier::IntVector<DIM>&
PatchDataFactory<DIM>::getGhostCellWidth() const
{
   return(d_ghosts);
}

/**********************************************************************
 * Default implementation                                             *
 **********************************************************************/
template<int DIM>
MultiblockDataTranslator<DIM>*
PatchDataFactory<DIM>::getMultiblockDataTranslator()
{
   return (MultiblockDataTranslator<DIM>*) NULL;
}


}
}
#endif
