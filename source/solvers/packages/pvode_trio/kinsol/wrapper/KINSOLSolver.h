//
// File:        KINSOLSolver.h
// Package:     SAMRAI solvers
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Wrapper class for KINSOL solver function calls and data
//

#ifndef included_solv_KINSOLSolver
#define included_solv_KINSOLSolver

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

/*
************************************************************************
*  THIS CLASS WILL BE UNDEFINED IF THE LIBRARY IS BUILT WITHOUT KINSOL
************************************************************************
*/
#ifdef HAVE_KINSOL

#ifndef included_String
#include <string>
using namespace std;
#define included_String
#endif
#ifndef included_tbox_IOStream
#include "tbox/IOStream.h"
#endif
#ifndef included_solv_PVodeTrioAbstractVector
#include "PVodeTrioAbstractVector.h"
#endif
#ifndef included_solv_KINSOLAbstractFunctions
#include "KINSOLAbstractFunctions.h"
#endif

// KINSOL includes
#ifndef included_kinsol_h
#define included_kinsol_h
extern "C" {
#include "kinsol.h"
}
#endif
#ifndef included_kinspgmr_h
#define included_kinspgmr_h
extern "C" {
#include "kinspgmr.h"
}
#endif

namespace SAMRAI {
    namespace solv {

/**
 * Class KINSOLSolver serves as a C++ wrapper for the KINSOL nonlinear
 * algebraic equation solver package and its data structures.  It is intended 
 * to be sufficiently generic to be used independently of the SAMRAI framework.
 * This class declares four private static member functions to link 
 * user-defined routines for nonlinear residual calculation, preconditioner
 * setup and solve, and Jacobian-vector product.  The implementation of these
 * functions is defined by the user in a subclass of the abstract base class
 * KINSOLAbstractFunctions.  The vector objects used within the solver 
 * are given in a subclass of the abstract class PVodeTrioAbstractVector. 
 * The PVodeTrioAbstractVector class defines the vector kernel operations 
 * required by the KINSOL package so that they may be easily supplied 
 * by a user who opts not to use the vector kernel supplied by the KINSOL
 * package.
 * 
 * Note that this class provides no input or restart capabilities and 
 * relies on KINSOL for output reporting.  When using KINSOL in an
 * application using SAMRAI, it is straightforward to include this
 * functionality in the entity using this solver class.
 *
 * KINSOL was developed in the Center for Applied Scientific Computing (CASC) 
 * at Lawrence Livermore National Laboratory (LLNL).  For more information 
 * about KINSOL and a complete description of the operations and data 
 * structures used by this class, see A.G. Taylor and A.C. Hindmarsh, 
 * "User documentation for KINSOL, a nonlinear solver for sequential and 
 * parallel computers", UCRL-ID-131185, Lawrence Livermore National 
 * Laboratory, 1998. 
 *
 * @see solv::KINSOLAbstractFunctions
 * @see solv::PVodeTrioAbstractVector 
 */

class KINSOLSolver
{
public:
   /**
    * Constructor for KINSOLSolver sets default KINSOL parameters 
    * and initializes the solver package with user-supplied functions.  Solver 
    * parameters may be changed later using member functions described
    * below.  The integer flags indicate whether user-supplied preconditioner 
    * and Jacobian-vector product function should be used.  Zero indicates
    * no user function; otherwise, user function will be used by the nonlinear
    * solver.
    *
    * Important note:  The solution vector is not passed into the constructor.
    * Before the solver can be used, the initialize() function must be called.
    *
    * When assertion checking is active, an unrecoverable assertion will
    * result if pointer to functions is null or string is empty.
    */
   KINSOLSolver(const string& object_name,
                        KINSOLAbstractFunctions* my_functions,
                        const int uses_preconditioner,
                        const int uses_jac_times_vector);

   /**
    * Virtual destructor for KINSOLSolver.
    */
   ~KINSOLSolver();

   /**
    * Initialize solver with solution vector.  The solution vector is
    * required to initialize the memory record used internally within
    * KINSOL.  This routine must be called before the solver can be used.
    * 
    * When assertion checking is active, an unrecoverable assertion will
    * result if vector pointer is null or solution vector has already been
    * set.
    */
   void initialize(PVodeTrioAbstractVector* solution); 

   /**
    * Solve nonlinear problem and return integer termination code defined
    * by KINSOL.  The default return value is KINSOL_SUCCESS (= 1) 
    * indicating success.  Return values which indicate non-recoverable
    * nonlinear solver behavior are KINSOL_NO_MEM (= -1), 
    * KINSOL_INPUT_ERROR (= -2), and KINSOL_LSOLV_NO_MEM (= -3).  
    * Return values PRECONDSET_FAILURE (= 9), and PRECONDSOLVE_FAILURE (= 10)
    * generally indicate non-recoverable behavior in the preconditioner.
    * See kinsol.h header file for more information about return values.  
    *
    * If KINSOL requires re-initialization, it is automatically done before 
    * the solve.  This may be required if any of the KINSOL data parameters 
    * have changed since the last call to the solver. 
    */
   int solve();

   /**
    * Accessory function for setting KINSOL output log file name and output
    * printing options.  Output file name and options may be changed
    * throughout run as desired.
    *
    * KINSOL printing options are:
    * 


    * - \b 0 {no statistics printed}
    * - \b 1 {output iteration count, residual norm, number function calls}
    * - \b 2 {same as 1, but with statistics on globalization process}
    * - \b 3 {same as 2, but with more Krylov iteration statistics}
    * 


    * The default is no output (i.e., 0).  If the file name string is empty
    * the default file name "kinsol.log" is used.
    *
    * See KINSOL documentation for more information.
    */
   void setLogFileData(const string& log_fname,
                       const int flag);

   /**
    * Accessory functions for passing user-defined function information
    * to KINSOL.   
    *
    * my_functions is a pointer to the abstract function subclass object
    * that defines the residual calculation and preconditioner functions.
    *
    * uses_preconditioner turns user preconditioner on or off.   
    *
    * uses_jac_times_vector turns user Jacobian-vector product on or off.
    *
    * Flags use "TRUE"/"FALSE" values defined in KINSOL.  See KINSOL 
    * documentation for more information.
    */
   void setKINSOLFunctions(KINSOLAbstractFunctions* my_functions,
                           const int uses_preconditioner,
                           const int uses_jac_times_vector);

   ///
   void setPreconditioner(const int uses_preconditioner);

   ///
   void setJacobianTimesVector(const int uses_jac_times_vector); 

   /**
    * Return pointer to object that provides user-defined functions for KINSOL.
    */
   KINSOLAbstractFunctions* getKINSOLFunctions() const;

   /**
    * Set vector used by KINSOL to scale either nonlinear solution vector
    * or nonlinear residual vector.  The elements of the scaling vectors 
    * must be positive.  In either case, the scaling vector should be 
    * defined so that the vector formed by taking the element-wise product 
    * of the solution/residual vector and scaling vector has all elements 
    * roughly the same magnitude when the solution vector IS/IS NOT NEAR 
    * a root of the nonlinear function.  
    *
    * See KINSOL documentation for more information.
    */
   void setSolutionScaleVector(PVodeTrioAbstractVector* uscale);

   ///
   void setResidualScaleVector(PVodeTrioAbstractVector* fscale);

   /**
    * Set constraints on nonlinear solution.  By default the constraint
    * vector is null.
    *
    * The constraints are applied in KINSOL as follows:
    * 


    * - \b {if constraints[i] > 0.0, then the constraint is solution[i]>0.0}
    * - \b {if constraints[i] < 0.0, then the constraint is solution[i]<0.0}
    * - \b {if constraints[i] = 0.0, then no constraint on solution[i]}
    * 


    * 
    * See KINSOL documentation for more information.
    */
   void setConstraintVector(PVodeTrioAbstractVector* constraints);

   /**
    * Accessory functions for setting nonlinear solver parameters.
    * Parameters and default values are:
    *
    * Residual stopping tolerance is tolerarnce on max_norm(fscale * residual),
    * where product of vectors is another vector each element of which is
    * the product of the corresponding entries in the original vectors. 
    * The default is \f$machine_epsilon^(1/3)\f$.
    *
    * Default maximum nonlinear iterations is 200.
    * 
    * Default maximum Krylov dimension is 1.
    *
    * Options for global Newton method are: INEXACT_NEWTON = 0, LINESEARCH = 1.
    * The default is INEXACT_NEWTON.
    *
    * Default maximum Newton step is 1000*max(norm(uscale*u_0), norm(uscale)),
    * where u_0 is the initial guess at the solution.
    *
    * Default scaled step tolerarnce between successive nonlinear iterates is
    * \f$machine_epsilon^(2/3)\f$.
    * 
    * Default relative error for nonlinear function is set to machine_epsilon.
    *
    * Scalar update constraint value restricts update of solution to 
    * del(u)/u < constraint_value.  Here, vector ratio is another vector
    * each element of which is the ratio of the corresponding entries in 
    * the original vectors.  The default is no constraint.
    * 
    * See KINSOL documentation for more information.
    */
   void setResidualStoppingTolerance(const double tol);

   ///
   void setMaxIterations(const int maxits);

   ///
   void setMaxKrylovDimension(const int kdim);

   ///
   void setGlobalStrategy(const int global);

   ///
   void setMaxNewtonStep(const double maxstep);

   ///
   void setNonlinearStepTolerance(const double tol);

   ///
   void setRelativeFunctionError(const double reserr);

   ///
   void setSolutionUpdateConstraint(const double constraint);

   /**
    * Accessory functions for setting convergence tests for inner linear
    * solvers within an inexact Newton method.  In general, the linear
    * solver attempts to produce a step p, satisfying:
    * norm(F(u) + J(u)*p) <= (eta + u_round)*norm(F(u)), where the norm
    * is a scaled L2-norm.
    *
    * The convergence test indicates the value for eta; options are:
    * 


    * - \b 0 == ETACHOICE1{Choice 1 of Eisenstat and Walker}
    * - \b 1 == ETACHOICE2{Choice 2 of Eisenstat and Walker}
    * - \b 2 == ETACONSTANT{use constant value for eta}.
    * 

  
    * The default option is ETACONSTANT.
    *
    * The default constant value for eta is 0.1.
    *
    * For choice ETACHOICE2, alpha = 2.0 and gamma = 0.9 are defaults.
    *
    * See KINSOL documentation for more information.
    */
   void setLinearSolverConvergenceTest(const int conv);

   ///
   void setLinearSolverConstantTolerance(const double tol);

   ///
   void setEisenstatWalkerParameters(const double alpha, const double gamma);

   /**
    * Accessory functions for setting preconditioner parameters.
    *
    * Set preconditioner setup flag options are:
    * 


    * - \b 0 {force call to precond setup routine on each call to solver}
    * - \b 1 {prevent initial call to precond setup routine on call to solver}
    * 


    * Typically, one chooses 1 only after beginning the first of a series 
    * of calls with 0 value.  The default is 0 (i.e., precond setup routine 
    * always called).
    *
    * Default maximum number of steps calling the preconditioner without
    * calling the preconditioner setup routine is 1. 
    *
    * Default maximum number of linear solver restarts is 0 (no restarts).
    *
    * See KINSOL documentation for more information.
    */
   void setPrecondSetupFlag(const int flag);
   
   ///
   void setMaxStepsWithNoPrecondSetup(const int maxsolv);

   ///
   void setMaxLinearSolveRestarts(const int restarts);

   /**
    * Accessory functions to retrieve information fom KINSOL.
    *
    * See KINSOL documentation for more information.
    */
   int getTotalNumberOfNonlinearIterations() const;

   /// 
   int getTotalNumberOfFunctionCalls() const;

   ///
   int getTotalNumberOfBetaConditionFailures() const;

   ///
   int getTotalNumberOfBacktracks() const;

   ///
   double getScaledResidualNorm() const;

   ///
   double getNewtonStepLength() const;

   /**
    * Print out all data members for this object.
    */
   void printClassData(ostream& os) const;

private:
   /*
    * Static member functions for linkage with KINSOL routines.
    * See header file for KINSOLAbstractFunctions for more information.
    */ 
   static void KINSOLFuncEval(int neq,
                              PVodeTrioAbstractVector* soln,
                              PVodeTrioAbstractVector* fval,
                              void* my_solver);

   static int KINSOLPrecondSet(int neq,
                               PVodeTrioAbstractVector* soln, 
                               PVodeTrioAbstractVector* soln_scale, 
                               PVodeTrioAbstractVector* fval, 
                               PVodeTrioAbstractVector* fval_scale, 
                               PVodeTrioAbstractVector* vtemp1, 
                               PVodeTrioAbstractVector* vtemp2, 
                               SysFn sys_func, 
                               double mach_roundoff, 
                               long int* nfePtr, 
                               void* my_solver);

   static int KINSOLPrecondSolve(int neq,
                                 PVodeTrioAbstractVector* soln,
                                 PVodeTrioAbstractVector* soln_scale,
                                 PVodeTrioAbstractVector* fval,
                                 PVodeTrioAbstractVector* fval_scale,
                                 PVodeTrioAbstractVector* rhs,
                                 PVodeTrioAbstractVector* vtemp,
                                 SysFn sys_func,
                                 double mach_roundoff,
                                 long int* nfePtr,
                                 void* my_solver);

   static int KINSOLJacobianTimesVector(void* my_solver,
                                        PVodeTrioAbstractVector* vector,
                                        PVodeTrioAbstractVector* prod,
                                        int* flags,
                                        PVodeTrioAbstractVector* soln);

   /*
    * Open KINSOL log file, allocate main memory for KINSOL and initialize
    * KINSOL memory record.  KINSOL is initialized based on current state
    * of solver parameter data members.  If any solver parameters have 
    * changed since last initialization, this function will be automatically
    * invoked at next call to solver. 
    */
   void initializeKINSOL();

   string d_object_name;

   /*
    * The following data members are input or set to default values in
    * the KINSOLSolver constructor.  Many of these can be altered at
    * any time through class member functions.  When this occurs,
    * KINSOL may need to be re-initialized (e.g., if Krylov dimension 
    * changes, KINSOL must change its memeory record).  Then the
    * initializeKINSOL() member function will be invoked when 
    * nonlinear solve function is called next.
    */

   /*
    * Nonlinear solution vector.
    */
   PVodeTrioAbstractVector* d_solution_vector;

   /*
    * Pointer to object which provides user-supplied functions to KINSOL.
    */
   KINSOLAbstractFunctions* d_KINSOL_functions;

   /*
    * Boolean flags used during KINSOL initialization to provide correct 
    * static function linkage with KINSOL package.
    */
   bool d_uses_preconditioner;
   bool d_uses_jac_times_vector;

   /*
    * KINSOL input and initialization parameters.
    */
   KINMem    d_kin_mem;                       // KINSOL memory structure
   FILE*     d_kinsol_log_file;               // KINSOL message log file
   string    d_kinsol_log_file_name;          // KINSOL log file name

   /*
    * Nonlinear solution and residual scaling vectors, and integer flags
    * to determine ownership of scaling vectors.
    */
   PVodeTrioAbstractVector* d_soln_scale;
   bool d_my_soln_scale_vector;
   PVodeTrioAbstractVector* d_fval_scale;
   bool d_my_fval_scale_vector;

   /*
    * Constraints on nonlinear solution vector.
    */
   PVodeTrioAbstractVector* d_constraints; 

   /*
    * Integer flag indicating whether KINSOL needs initialization 
    * when solver is called.
    */
   int d_KINSOL_needs_initialization;

   /*
    * KINSOL nonlinear and linear solver parameters
    */
// Krylov method parameters
   int d_krylov_dimension;  // maximum krylov dimension
   int d_max_restarts;      // max. num. of linear solver restarts allowed
   int d_max_solves_no_set; // max. num. of steps calling preconditioner
                            // without resetting preconditioner
   int d_neq;               // number of eqns. in algebraic system
                            // I really don't know what this is for??????

// Nonlinear solver inputs
   int    d_global_strategy;     // globalization method for Newton steps.
   double d_residual_tol;        // stop tol. on scaled nonlinear residual
   double d_step_tol;            // stop tol. on consecutive step difference

   /*
    * Optional KINSOL nonlinear solver parameters.
    */
   long int  d_int_optional_input[OPT_SIZE];
   double    d_real_optional_input[OPT_SIZE];

};

}
}
#endif
#endif
