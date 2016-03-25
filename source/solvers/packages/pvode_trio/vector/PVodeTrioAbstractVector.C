//
// File:        solv::PVodeTrioAbstractVector.C
// Package:     SAMRAI solvers package
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Interface to C++ vector implementation for PVodeTrio package.
//

#include "PVodeTrioAbstractVector.h"

#if HAVE_KINSOL || HAVE_PVODE

namespace SAMRAI {
   namespace solv {

PVodeTrioAbstractVector::PVodeTrioAbstractVector()
{
}

PVodeTrioAbstractVector::~PVodeTrioAbstractVector()
{
}


}
}

#endif
