//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/solvers/packages/sundials/cvode/CVODEAbstractFunctions.h $
// Package:     SAMRAI solvers package
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1750 $
// Modified:    $LastChangedDate: 2007-12-10 15:19:17 -0800 (Mon, 10 Dec 2007) $
// Description: Interface to user-specified functions for CVODE package
//

#ifndef included_solv_CVODEAbstractFunctions
#define included_solv_CVODEAbstractFunctions

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_solv_SundialsAbstractVector
#include "SundialsAbstractVector.h"
#endif

#ifdef HAVE_SUNDIALS

namespace SAMRAI {
   namespace solv {

/**
 * Class CVODEAbstractFunctions is an abstract base class that defines
 * an interface for the user-supplied RHSFunction and preconditioner 
 * routines to be used with CVODE and CVSpgmr via the C++ wrapper 
 * class CVODESolver.  To use CVODE with the C++ wrapper one must 
 * derive a subclass of this base class and pass it into the CVODESolver 
 * constructor.  The pure virtual member functions in this interface are 
 * used by CVODE and CVSpgmr during the ODE integration process.  The 
 * complete argument lists in the function signatures defined by CVODE 
 * for the user-supplied routines have been preserved for the most part.  
 * In a few cases, some arguments do not appear in the function signatures 
 * below since they are superfluous via this interface.
 *
 * @see solv::CVODESolver
 * @see solv::SundialsAbstractVector
 */

class CVODEAbstractFunctions
{
public:
   /**
    * The constructor and destructor for CVODEAbstractFunctions
    * is empty.
    */
   CVODEAbstractFunctions();
   virtual ~CVODEAbstractFunctions();

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
   virtual int evaluateRHSFunction(double t,
				   SundialsAbstractVector* y,
				   SundialsAbstractVector* y_dot) = 0;

   /**
    * User-supplied function for setting up the preconditioner 
    * to be used in the solution of the linear system that arises
    * during Newton iteration.
    */
   virtual int CVSpgmrPrecondSet(double t,
                                 SundialsAbstractVector* y,
                                 SundialsAbstractVector* fy,
                                 int jok,
                                 int *jcurPtr,
                                 double gamma,
                                 SundialsAbstractVector* vtemp1,
                                 SundialsAbstractVector* vtemp2,
                                 SundialsAbstractVector* vtemp3) = 0; 

   /**
    * User-supplied function for setting up the preconditioner 
    * to be used in the solution of the linear system that arises
    * during Newton iteration.
    */
   virtual int CVSpgmrPrecondSolve(double t,
                                   SundialsAbstractVector* y,
                                   SundialsAbstractVector* fy,
                                   SundialsAbstractVector* r,
                                   SundialsAbstractVector* z,
                                   double gamma,
                                   double delta,
                                   int lr,
				   SundialsAbstractVector* vtemp) = 0;

};

}
}

#endif

#endif
