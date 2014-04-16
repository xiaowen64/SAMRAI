/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2014 Lawrence Livermore National Security, LLC
 * Description:   Interface to application-specific patch functions to support
 *                MethodOfLines integration algorithm
 *
 ************************************************************************/

#ifndef included_algs_MethodOfLinesPatchStrategy_C
#define included_algs_MethodOfLinesPatchStrategy_C

#include "SAMRAI/algs/MethodOfLinesPatchStrategy.h"
#include "SAMRAI/hier/VariableDatabase.h"

namespace SAMRAI {
namespace algs {

/*
 *************************************************************************
 *
 * Note: hier::Variable contexts should be consistent with those in
 *       MethodOfLinesIntegrator class.
 *
 *************************************************************************
 */

MethodOfLinesPatchStrategy::MethodOfLinesPatchStrategy():
   xfer::RefinePatchStrategy(),
   xfer::CoarsenPatchStrategy()
{
   d_interior_with_ghosts.reset();
   d_interior.reset();
}

MethodOfLinesPatchStrategy::~MethodOfLinesPatchStrategy()
{
}

}
}
#endif
