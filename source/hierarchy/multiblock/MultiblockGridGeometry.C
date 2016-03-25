//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/trunk/source/multiblock/MultiblockGridGeometry.C $
// Package:	SAMRAI multiblock package
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 878 $
// Modified:	$LastChangedDate: 2006-01-09 16:55:30 -0800 (Mon, 09 Jan 2006) $
// Description:	Strategy interface to user routines for refining AMR data.
//

#ifndef included_hier_MultiblockGridGeometry_C
#define included_hier_MultiblockGridGeometry_C 

#include "MultiblockGridGeometry.h"

#include "tbox/Utilities.h"

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
MultiblockGridGeometry<DIM>::MultiblockGridGeometry(
   tbox::Array< tbox::Pointer< hier::GridGeometry<DIM> > >& block_geoms)
{
   d_block_geometry = block_geoms;
}

template<int DIM>
MultiblockGridGeometry<DIM>::~MultiblockGridGeometry()
{
}


}
}
#endif
