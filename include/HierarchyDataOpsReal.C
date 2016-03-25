//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/mathops/hierarchy/HierarchyDataOpsReal.C $
// Package:     SAMRAI mathops
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Interface to templated operations for real data on hierarchy.
//

#ifndef included_math_HierarchyDataOpsReal_C
#define included_math_HierarchyDataOpsReal_C

#include "HierarchyDataOpsReal.h"

namespace SAMRAI {
    namespace math {

template<int DIM, class TYPE>
HierarchyDataOpsReal<DIM,TYPE>::HierarchyDataOpsReal()
{
}

template<int DIM, class TYPE> 
HierarchyDataOpsReal<DIM,TYPE>::~HierarchyDataOpsReal()
{
}

/*
*************************************************************************
*                                                                       *
* The following are private and cannot be used, but they are defined    *
* here for compilers that require that every template declaration have  *
* a definition (a stupid requirement, if you ask me).                   *
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
HierarchyDataOpsReal<DIM,TYPE>::HierarchyDataOpsReal(
   const HierarchyDataOpsReal<DIM,TYPE>& foo)
{
   (void) foo;  // private and not used (but included for stupid compilers)
}

template<int DIM, class TYPE>
void HierarchyDataOpsReal<DIM,TYPE>::operator=(
   const HierarchyDataOpsReal<DIM,TYPE>& foo)
{
   (void) foo;  // private and not used (but included for stupid compilers)
}

}
}
#endif
