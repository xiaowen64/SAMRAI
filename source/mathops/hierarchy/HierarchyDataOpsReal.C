//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-1/source/mathops/hierarchy/HierarchyDataOpsReal.C $
// Package:     SAMRAI mathops
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1846 $
// Modified:    $LastChangedDate: 2008-01-11 09:51:05 -0800 (Fri, 11 Jan 2008) $
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
   NULL_USE(foo);
}

template<int DIM, class TYPE>
void HierarchyDataOpsReal<DIM,TYPE>::operator=(
   const HierarchyDataOpsReal<DIM,TYPE>& foo)
{
   NULL_USE(foo);
}

}
}
#endif
