//
// File:        SNESAbstractFunctions.h
// Package:     SAMRAI solvers
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 601 $
// Modified:    $Date: 2005-09-06 11:23:15 -0700 (Tue, 06 Sep 2005) $
// Description: Interface to user functions for SAMRAI-based PETSc SNES context
//

#ifndef included_solv_SNESAbstractFunctions
#define included_solv_SNESAbstractFunctions

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifdef REQUIRES_CMATH
#include <cmath>
#endif

/*
************************************************************************
*  THIS CLASS WILL BE UNDEFINED IF THE LIBRARY IS BUILT WITHOUT PETSC
************************************************************************
*/
#ifdef HAVE_PETSC
#ifndef included_petsc_snes
#define included_petsc_snes
extern "C" {
#ifdef PETSC2028
#include "snes.h"
#else
#include "petscsnes.h"
#endif
}
#endif

namespace SAMRAI {
    namespace solv {

/*!
 * @brief Abstract base class that declares
 * the functions to be used with the PETSc SNES nonlinear solver package.  
 *
 * Class SNESAbstractFunctions is an abstract base class that declares
 * the functions to be used with the PETSc SNES nonlinear solver package.  
 * This class works in cooperation with the SNES_SAMRAIContext class.  
 * To provide these functions to the PETSc SNES solver package, a subclass
 * of this base class must be instantiated and be supplied to the 
 * SNES_SAMRAIContext constructor.  Pointers to these functions will 
 * be stored in a SNES context and invoked from within the nonlinear solver.
 * Note that the virtual members of this class are all protected.  They should
 * not be used outside of a subclass of this class.
 *
 * @see solv::SNES_SAMRAIContext
 */

class SNESAbstractFunctions
{
public:
   /*!
    * Uninteresting constructor for SNESAbstractFunctions.
    */
   SNESAbstractFunctions();

   /*!
    * Uninteresting destructor for SNESAbstractFunctions.
    */
   virtual ~SNESAbstractFunctions();

   /*!
    * User-supplied nonlinear function evaluation.  Returns 0 if successful.
    * Arguments:
    * 
    * @param xcur (IN) the current iterate for the nonlinear system
    * @param fcur (OUT) current function value
    *
    * IMPORTANT: This function must not modify xcur.
    */
   virtual int evaluateNonlinearFunction(Vec xcur, Vec fcur) = 0;

   /*!
    * Optional user-supplied routine to evaluate the Jacobian of the
    * system.  This function can be empty if the matrix-free option
    * has been selected.  Returns 0 if successful.  Arguments:
    * 
    * @param x (IN) current Newton iterate.
    *
    * IMPORTANT: This function must not modify x.
    */
   virtual int evaluateJacobian(Vec x) = 0;

   /*!
    * Optional user-supplied Jacobian-vector product routine.  This 
    * function can be empty if the matrix-free option has been selected.
    * Returns 0 if successful.  Arguments:
    * 
    *  @param x (IN) vector to be multiplied by the Jacobian.
    *  @param y (OUT) the product of the Jacobian and vector
    *
    * IMPORTANT: This function must not modify x.
    */
   virtual int jacobianTimesVector(Vec x, Vec y) = 0;

   /*!
    * User-supplied preconditioner setup function.  The setup   
    * function is called to provide matrix data for the subsequent
    * call(s) to applyPreconditioner().  The integer return value
    * is a flag indicating success if 0 is returned, and failure
    * otherwise.  Together setupPreconditioner() and applyPreconditioner() 
    * form a right preconditoner for the PETSc Krylov solver.
    * Returns 0 if successful.  
    *
    */
   virtual int setupPreconditioner(Vec xcur) = 0;

   /*!
    * User-supplied preconditioner solve function.  This function must
    * solve \f$M z = r\f$, where \f$M\f$ is the right preconditioner
    * matrix formed by setupPreconditioner(). The integer return value
    * is a flag indicating success if 0 is returned, and failure otherwise.  
    * Arguments:
    * 
    * @param r (IN) right-hand side of preconditioning system
    * @param z (OUT) result of applying preconditioner to right-hand side
    *
    * IMPORTANT: This function must not modify r.
    */
   virtual int applyPreconditioner(Vec r, Vec z) = 0;
};

}
}

#endif
#endif
