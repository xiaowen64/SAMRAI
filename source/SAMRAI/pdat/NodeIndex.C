/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_NodeIndex_C
#define included_pdat_NodeIndex_C

#include "SAMRAI/pdat/NodeIndex.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/pdat/NodeIndex.I"
#endif
namespace SAMRAI {
namespace pdat {

std::vector<hier::IntVector> NodeIndex::s_offsets[tbox::Dimension::
                                                  MAXIMUM_DIMENSION_VALUE];
bool NodeIndex::s_offsets_are_set[tbox::Dimension::MAXIMUM_DIMENSION_VALUE] = { false };

}
}
#endif
