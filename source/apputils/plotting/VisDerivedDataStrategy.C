//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/apputils/plotting/VisDerivedDataStrategy.C $
// Package:     SAMRAI application utilities
// Copyright:   (c) 1997-2003 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
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
