/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Interface for the AMR Index object
 *
 ************************************************************************/

#ifndef included_hier_Index_C
#define included_hier_Index_C

#include "SAMRAI/hier/Index.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/hier/Index.I"
#endif

namespace SAMRAI {
namespace hier {

Index * Index::s_zeros[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
Index * Index::s_ones[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

Index * Index::s_mins[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
Index * Index::s_maxs[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

tbox::StartupShutdownManager::Handler
Index::s_initialize_finalize_handler(
   Index::initializeCallback,
   0,
   0,
   Index::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);

void Index::initializeCallback()
{
   for (unsigned short d = 0; d < tbox::Dimension::MAXIMUM_DIMENSION_VALUE; ++d) {
      s_zeros[d] = new Index(tbox::Dimension(static_cast<unsigned short>(d + 1)), 0);
      s_ones[d] = new Index(tbox::Dimension(static_cast<unsigned short>(d + 1)), 1);

      s_mins[d] = new Index(tbox::Dimension(static_cast<unsigned short>(d + 1)),
            tbox::MathUtilities<int>::getMin());
      s_maxs[d] = new Index(tbox::Dimension(static_cast<unsigned short>(d + 1)),
            tbox::MathUtilities<int>::getMax());
   }
}

void Index::finalizeCallback()
{
   for (int d = 0; d < tbox::Dimension::MAXIMUM_DIMENSION_VALUE; ++d) {
      delete s_zeros[d];
      delete s_ones[d];

      delete s_mins[d];
      delete s_maxs[d];
   }
}

}
}

#endif
