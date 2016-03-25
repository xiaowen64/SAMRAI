//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/trunk/source/multiblock/MultiblockGridGeometry.C $
// Package:	SAMRAI multiblock package
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 878 $
// Modified:	$LastChangedDate: 2006-01-09 16:55:30 -0800 (Mon, 09 Jan 2006) $
// Description:	data translator for Multiblock.
//

#ifndef included_hier_MultiblockDataTranslator_C
#define included_hier_MultiblockDataTranslator_C 

#include "MultiblockDataTranslator.h"

namespace SAMRAI {
    namespace hier {

/*
*************************************************************************
*									*
* The default constructor and virtual destructor do nothing             *
* particularly interesting.		                                *
*									*
*************************************************************************
*/

template<int DIM>
MultiblockDataTranslator<DIM>::MultiblockDataTranslator()
{
}

template<int DIM>
MultiblockDataTranslator<DIM>::~MultiblockDataTranslator()
{
}


}
}
#endif
