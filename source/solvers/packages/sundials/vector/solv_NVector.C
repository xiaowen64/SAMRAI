//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/solvers/packages/sundials/vector/solv_NVector.C $
// Package:     SAMRAI solvers
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2132 $
// Modified:    $LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description: Interface to C++ vector implementation for Sundials package.
//

#include "SAMRAI_config.h"

#include "tbox/Utilities.h"

#ifdef HAVE_SUNDIALS

#include "solv_NVector.h"

#include "SundialsAbstractVector.h"

#define SABSVEC_CAST(v) (static_cast<SAMRAI::solv::SundialsAbstractVector *>(v->content))

extern "C" {

void N_VPrint_SAMRAI(N_Vector v) {
   SABSVEC_CAST(v) -> printVector();
}

}

#endif
