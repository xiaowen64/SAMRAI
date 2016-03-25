//
// File:        VisDerivedDataStrategy.C
// Package:     SAMRAI application utilities
// Copyright:   (c) 1997-2003 The Regents of the University of California
// Revision:    $Revision: 47 $
// Modified:    $Date: 2004-12-09 16:08:57 -0800 (Thu, 09 Dec 2004) $
// Description: Interface for writing user-defined data to either VisIt or
//              Vizamrai dump file
//

#ifndef included_appu_VisDerivedDataStrategy_C
#define included_appu_VisDerivedDataStrategy_C

#include "VisDerivedDataStrategy.h"

namespace SAMRAI {
    namespace appu {

template<int DIM>  VisDerivedDataStrategy<DIM>::VisDerivedDataStrategy()
{
}

template<int DIM>  VisDerivedDataStrategy<DIM>::~VisDerivedDataStrategy()
{
}

}
}

#endif
