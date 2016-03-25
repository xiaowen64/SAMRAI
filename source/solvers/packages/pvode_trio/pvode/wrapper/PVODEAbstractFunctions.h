//
// File:        PVODEAbstractFunctions.h
// Package:     SAMRAI solvers package
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Interface to user-specified functions for PVODE package
//

#ifndef included_solv_PVODEAbstractFunctions
#define included_solv_PVODEAbstractFunctions

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_solv_PVodeTrioAbstractVector
#include "PVodeTrioAbstractVector.h"
#endif

#ifdef HAVE_PVODE

namespace SAMRAI {
   namespace solv {

/**
 * Class PVODEAbstractFunctions is an abstract base class that defines
 * an interface for the user-supplied RHSFunction and preconditioner 
 * routines to be used with PVODE and CVSpgmr via the C++ wrapper 
 * class PVODESolver.  To use PVODE with the C++ wrapper one must 
 * derive a subclass of this base class and pass it into the PVODESolver 
 * constructor.  The pure virtual member functions in this interface are 
 * used by PVODE and CVSpgmr during the ODE integration process.  The 
 * complete argument lists in the function signatures defined by PVODE 
 * for the user-supplied routines have been preserved for the most part.  
 * In a few cases, some arguments do not appear in the function signatures 
 * below since they are superfluous via this interface.
 *
 * @see solv::PVODESolver
 * @see solv::PVodeTrioAbstractVector
 */

class PVODEAbstractFunctions
{
public:
   /**
    * The constructor and destructor for PVODEAbstractFunctions
    * is empty.
    */
   PVODEAbstractFunctions();
   virtual ~PVODEAbstractFunctions();

   /**
    * User-supplied right-hand side function evaluation.
    *
    * The function arguments are:
    * 


    * - \b t        (INPUT) {current value of the independent variable}
    * - \b y        (INPUT) {current value of dependent variable vector}
    * - \b y_dot   (OUTPUT){current value of the derivative of y}
    * 


    *
    * IMPORTANT: This function must not modify the vector y. 
    */
   virtual void evaluateRHSFunction(double t,
                                    PVodeTrioAbstractVector* y,
                                    PVodeTrioAbstractVector* y_dot) = 0;

   /**
    * User-supplied function for setting up the preconditioner 
    * to be used in the solution of the linear system that arises
    * during Newton iteration.
    */
   virtual int CVSpgmrPrecondSet(int neq,
                                 double t,
                                 PVodeTrioAbstractVector* y,
                                 PVodeTrioAbstractVector* fy,
                                 int jok,
                                 int *jcurPtr,
                                 double gamma,
                                 PVodeTrioAbstractVector* ewt,
                                 double h,
                                 double mach_roundoff,
                                 long int *nfePtr,
                                 void *P_data,
                                 PVodeTrioAbstractVector* vtemp1,
                                 PVodeTrioAbstractVector* vtemp2,
                                 PVodeTrioAbstractVector* vtemp3) = 0; 

   /**
    * User-supplied function for setting up the preconditioner 
    * to be used in the solution of the linear system that arises
    * during Newton iteration.
    */
   virtual int CVSpgmrPrecondSolve(int neq,
                                   double t,
                                   PVodeTrioAbstractVector* y,
                                   PVodeTrioAbstractVector* fy,
                                   PVodeTrioAbstractVector* vtemp,
                                   double gamma,
                                   PVodeTrioAbstractVector* ewt,
                                   double delta,
                                   long int *nfePtr,
                                   PVodeTrioAbstractVector* r,
                                   int lr,
                                   void *P_data,
                                   PVodeTrioAbstractVector* z) = 0;

};

}
}

#endif

#endif
