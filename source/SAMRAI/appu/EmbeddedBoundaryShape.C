/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Base class for analytic embedded Boundaries 
 *
 ************************************************************************/

#ifndef included_appu_EmbeddedBoundaryShape_C
#define included_appu_EmbeddedBoundaryShape_C

#include "SAMRAI/appu/EmbeddedBoundaryShape.h"

#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace appu {

/*
 *******************************************************************
 *
 *  Empty constructor and destructor
 *
 *******************************************************************
 */
EmbeddedBoundaryShape::EmbeddedBoundaryShape()
{
}

EmbeddedBoundaryShape::~EmbeddedBoundaryShape()
{
}

/*
 *************************************************************************
 *                                                                       *
 * Default virtual function implementations.                             *
 *                                                                       *
 *************************************************************************
 */
bool
EmbeddedBoundaryShape::isInside(
   const double* xyz) const
{
   NULL_USE(xyz);
   TBOX_ERROR("EmbeddedBoundaryShape::isInside(): "
      << "\nNo implementation provided for this shape."
      << std::endl);
   return false;
}

void
EmbeddedBoundaryShape::isInside(
   const int* nx,
   const double* dx,
   const double* origin,
   int* inout) const
{
   NULL_USE(nx);
   NULL_USE(dx);
   NULL_USE(origin);
   NULL_USE(inout);

   TBOX_ERROR("EmbeddedBoundaryShape::isInside(): "
      << "\nNo implementation provided for this shape. "
      << std::endl);
}

}
}
#endif
