/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   AMR hierarchy generation and regridding routines. 
 *
 ************************************************************************/

#ifndef included_mesh_BaseGriddingAlgorithm_C
#define included_mesh_BaseGriddingAlgorithm_C

#include "SAMRAI/mesh/BaseGriddingAlgorithm.h"

namespace SAMRAI {
namespace mesh {

/*
 *************************************************************************
 *                                                                       *
 * Constructor and destructor for BaseGriddingAlgorithm.            *
 *                                                                       *
 *************************************************************************
 */
BaseGriddingAlgorithm::BaseGriddingAlgorithm(
   const tbox::Pointer<hier::PatchHierarchy> &hierarchy )
   : d_hierarchy(hierarchy)
{
}

/*
 *************************************************************************
 *                                                                       *
 * Destructor tells the tbox_RestartManager to remove this object from   *
 * the list of restart items.                                            *
 *                                                                       *
 *************************************************************************
 */
BaseGriddingAlgorithm::~BaseGriddingAlgorithm()
{
}

/*
 *************************************************************************
 *
 * Provide default for override_mapping.
 *
 *************************************************************************
 */
void BaseGriddingAlgorithm::makeCoarsestLevel(
   const double level_time)

{
   makeCoarsestLevel(
      level_time,
      hier::MappedBoxLevel(d_hierarchy->getDim()));
}

}
}
#endif
