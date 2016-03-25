/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/solvers/packages/sundials/vector/solv_NVector.h $
 * Package:     SAMRAI solvers
 * Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 2132 $
 * Modified:    $LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
 * Description: C interface to C++ vector implementation for Sundials package.
 */

#ifndef included_NVector_SAMRAI
#define included_NVector_SAMRAI
#endif

#include "SAMRAI_config.h"

#ifdef HAVE_SUNDIALS

#include "sundials/sundials_nvector.h"


extern "C" {

   /**
    * @brief Helper funtion for printing SAMRAI N_Vector.
    *
    */
   void N_VPrint_SAMRAI(N_Vector v);

}

#endif
