//
// File:	HierarchyDataOpsReal.C
// Package:     SAMRAI mathops
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
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
