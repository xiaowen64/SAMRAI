//
// File:        NonlinearSolverStrategy.h
// Package:     SAMRAI solvers
// Copyright:   (c) 1997-2000 The Regents of the University of California
// Revision:    $Revision: 47 $
// Modified:    $Date: 2004-12-09 16:08:57 -0800 (Thu, 09 Dec 2004) $
// Description: Interface between implicit integrator and nonlinear solver.
//

#ifndef included_solv_NonlinearSolverStrategy
#define included_solv_NonlinearSolverStrategy


#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_solv_SAMRAIVectorReal
#include "SAMRAIVectorReal.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif

namespace SAMRAI {
    namespace solv {

/**
 * @brief Abstract base class defining interface between an
 * algs::ImplicitIntegrator object and 
 * a nonlinear solver used to advance the solution in time.
 *
 * The interface follows the Strategy design pattern.
 * The methods declared in the interface are provided in a
 * concrete solver derived from this base class.
 *
 * @see algs::ImplicitIntegrator
 */

template<int DIM> class NonlinearSolverStrategy
{
public:
   /**
    * Empty constructor for algs::NonlinearSolverStrategy.
    */
   NonlinearSolverStrategy();

   /**
    * Empty constructor for algs::NonlinearSolverStrategy.
    */
   virtual ~NonlinearSolverStrategy<DIM>();

   /**
    * Initialize the solver state.  The vector argument represents the
    * solution of the nonlinear system.  In general, this routine must
    * be called before the solve() routine is invoked.
    */
   virtual void initialize(
      const tbox::Pointer< SAMRAIVectorReal<DIM,double> > solution) = 0;

   /**
    * Solve the nonlinear problem and return the integer code defined by the
    * particular nonlinear solver package in use (e.g., indicating success 
    * or failure of solution process).  In general, the initialize() routine 
    * must be called before this solve function.
    */
   virtual int solve() = 0;

};

}
}

#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "NonlinearSolverStrategy.C"
#endif
