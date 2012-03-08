//
// File:        $URL$
// Package:     SAMRAI algorithms
// Copyright:   (c) 1997-2000 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision$
// Modified:    $LastChangedDate$
// Description: Routines for summing node data at patch boundaries.
//

namespace SAMRAI {
namespace algs {

/*
 *************************************************************************
 *
 * Static functions to determine number of patch data slots needed
 * for PatchBoundaryNodeSum objects.
 *
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
int
PatchBoundaryNodeSum::getNumSharedPatchDataSlots(
   int max_variables_to_register)
{
   // node boundary sum requires two internal outernode variables
   // (source and destination) for each registered variable.

   return 2 * max_variables_to_register;
}

SAMRAI_INLINE_KEYWORD
int
PatchBoundaryNodeSum::getNumUniquePatchDataSlots(
   int max_variables_to_register)
{
   NULL_USE(max_variables_to_register);
   // all patch data slots used by node boundary sum are static
   // and shared among all objects.

   return 0;
}

/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
const std::string&
PatchBoundaryNodeSum::getObjectName() const
{
   return d_object_name;
}

}
}
