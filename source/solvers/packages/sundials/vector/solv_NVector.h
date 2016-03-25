/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/solvers/packages/sundials/vector/solv_NVector.h $
 * Package:     SAMRAI solvers
 * Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 1765 $
 * Modified:    $LastChangedDate: 2007-12-11 15:15:21 -0800 (Tue, 11 Dec 2007) $
 * Description: C interface to C++ vector implementation for Sundials package.
 */

#ifndef included_NVector_SAMRAI
#define included_NVector_SAMRAI
#endif

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

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
