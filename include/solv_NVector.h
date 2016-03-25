/*
 * File:        NVector.h
 * Package:     SAMRAI solvers
 * Copyright:   (c) 1997-2005 The Regents of the University of California
 * Revision:    $Revision: 173 $
 * Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
 * Description: Interface to C++ vector implementation for PVodeTrio package.
 */

#ifndef included_NVector
#define included_NVector

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

/**
 * The current implementation of PVodeTrio does not allow one to incorporate
 * a new vector kernal implementation at link time.  The PVodeTrio package
 * must be compiled so that it sees the declarations of the vector structure 
 * and operations.  This header file contains the external definitions 
 * that must be seen by PVodeTrio include files for a new vector kernel to 
 * be used with that package.  Specifically, this file must be included by 
 * the header file "nvector.h" in the PVodeTrio package.
 *
 * NOTE: These function definitions assume that the vector data will be
 *       of type double.  Also, the return type int is used here whereas
 *       in PVodeTrio "boole" is used (boole is #define'd to be int).
 * 
 * @see solv::PVodeTrioAbstractVector
 */

#ifdef __cplusplus

#ifndef included_solv_PVodeTrioAbstractVector
#include "PVodeTrioAbstractVector.h"
#endif

typedef struct SAMRAI::solv::PVodeTrioAbstractVector* N_Vector;

#else
struct PVodeTrioAbstractVector;
typedef struct PVodeTrioAbstractVector* N_Vector;
#endif

extern N_Vector N_VNew(int N, void* machEnv);
extern void     N_VFree(N_Vector v);
extern void     N_VPrint(const N_Vector v);
extern void     N_VConst(double c, N_Vector v );
extern void     N_VScale(const double c, const N_Vector x,
                         N_Vector z);
extern void     N_VLinearSum(const double a, const N_Vector x,
                             const double b, const N_Vector y,
                             N_Vector z);
extern void     N_VProd(const N_Vector x, const N_Vector y,
                        N_Vector z);
extern void     N_VDiv(const N_Vector x, const N_Vector y,
                       N_Vector z);
extern void     N_VAbs(const N_Vector x,
                       N_Vector z);
extern void     N_VInv(const N_Vector x,
                       N_Vector z);
extern void     N_VAddConst(const N_Vector x, const double b,
                            N_Vector z);
extern double   N_VDotProd(const N_Vector x, const N_Vector y);
extern double   N_VMaxNorm(const N_Vector x);
extern double   N_VL1Norm(const N_Vector x);
extern double   N_VWL2Norm(const N_Vector x, const N_Vector w);
extern double   N_VWrmsNorm(const N_Vector x, const N_Vector w);
extern double   N_VMin(const N_Vector x);
extern int      N_VConstrProdPos(const N_Vector c, const N_Vector x);
extern void     N_VCompare(const double c, const N_Vector x,
                           N_Vector z);
extern int      N_VInvTest(const N_Vector x,
                           N_Vector z);

#endif
