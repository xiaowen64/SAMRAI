//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-1/source/solvers/packages/sundials/vector/solv_NVector.C $
// Package:     SAMRAI solvers
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1892 $
// Modified:    $LastChangedDate: 2008-01-23 08:29:36 -0800 (Wed, 23 Jan 2008) $
// Description: Interface to C++ vector implementation for Sundials package.
//

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_tbox_Utilities
#include "tbox/Utilities.h"
#endif

#ifdef HAVE_SUNDIALS

#ifndef included_NVector_SAMRAI
#include "solv_NVector.h"
#endif

#ifndef included_solv_SundialsAbstractVector
#include "SundialsAbstractVector.h"
#endif

#define SABSVEC_CAST(v) (static_cast<SAMRAI::solv::SundialsAbstractVector *>(v->content))

extern "C" {

void N_VPrint_SAMRAI(N_Vector v) {
   SABSVEC_CAST(v) -> printVector();
}

}

#endif
