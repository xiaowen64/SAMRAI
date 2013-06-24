/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2013 Lawrence Livermore National Security, LLC
 * Description:   Utility functions for using OpenMP.
 *
 ************************************************************************/

#ifndef included_tbox_OpenMPUtilities
#define included_tbox_OpenMPUtilities

#ifdef _OPENMP

#include "SAMRAI/SAMRAI_config.h"

#include <omp.h>

#define TBOX_omp_lock_t omp_lock_t

#define TBOX_omp_init_lock(LOCK_PTR) omp_init_lock(LOCK_PTR)
#define TBOX_omp_destroy_lock(LOCK_PTR) omp_destroy_lock(LOCK_PTR)

#define TBOX_omp_set_lock(LOCK_PTR) omp_set_lock(LOCK_PTR)
#define TBOX_omp_unset_lock(LOCK_PTR) omp_unset_lock(LOCK_PTR)

#define TBOX_omp_get_num_threads() omp_get_num_threads()
#define TBOX_omp_get_max_threads() omp_get_max_threads()

#else

#define TBOX_omp_lock_t int

#define TBOX_omp_init_lock(LOCK_PTR)
#define TBOX_omp_destroy_lock(LOCK_PTR)

#define TBOX_omp_set_lock(LOCK_PTR)
#define TBOX_omp_unset_lock(LOCK_PTR)

#define TBOX_omp_get_num_threads() (1)
#define TBOX_omp_get_max_threads() (1)

#endif

#endif
