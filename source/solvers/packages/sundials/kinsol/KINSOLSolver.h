//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/solvers/packages/sundials/kinsol/KINSOLSolver.h $
// Package:     SAMRAI solvers
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1765 $
// Modified:    $LastChangedDate: 2007-12-11 15:15:21 -0800 (Tue, 11 Dec 2007) $
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
#ifdef HAVE_SUNDIALS

#ifndef included_String
#include <string>
#define included_String
#endif
#ifndef included_tbox_IOStream
#include "tbox/IOStream.h"
#endif
#ifndef included_solv_SundialsAbstractVector
#include "SundialsAbstractVector.h"
#endif
#ifndef included_solv_KINSOLAbstractFunctions
#include "KINSOLAbstractFunctions.h"
#endif

// KINSOL includes
#ifndef included_kinsol_h
#define included_kinsol_h
#include "kinsol/kinsol.h"
#endif
#ifndef included_kinspgmr_h
#define included_kinspgmr_h
#include "kinsol/kinsol_spgmr.h"
#endif

#ifndef LACKS_SSTREAM
#define KINSOL_SAMRAI_ERROR(ierr) do {						\
      if (ierr != KIN_SUCCESS) {                                   				\
         std::ostringstream tboxos;							\
         SAMRAI::tbox::Utilities::abort(tboxos.str().c_str(), __FILE__, __LINE__);	\
      } 									\
} while (0)
#else
#define KINSOL_SAMRAI_ERROR(ierr) do {						\
      if (ierr != KIN_SUCCESS) {                                   				\
         std::ostrstream tboxos;							\
         CHKERRCONTINUE(ierr); 							\
         SAMRAI::tbox::Utilities::abort(tboxos.str(), __FILE__, __LINE__);	        \
      } 									\
} while (0)
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
 * are given in a subclass of the abstract class SundialsAbstractVector. 
 * The SundialsAbstractVector class defines the vector kernel operations 
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
 * @see solv::SundialsAbstractVector 
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
   KINSOLSolver(const std::string& object_name,
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
   void initialize(SundialsAbstractVector* solution); 

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
   void setLogFileData(const std::string& log_fname,
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
   void setSolutionScaleVector(SundialsAbstractVector* uscale);

   ///
   void setResidualScaleVector(SundialsAbstractVector* fscale);

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
   void setConstraintVector(SundialsAbstractVector* constraints);

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

   ///
   void setMaxStepsWithNoPrecondSetup(const int maxsolv);

   ///
   void setMaxLinearSolveRestarts(const int restarts);

   /** 
    * The number of nonlinear iterations between checks by the
    * nonlinear residual monitoring algorithm (specifies lenght of
    * subinterval) NOTE: should be a multiple of
    * MaxStepsWithNoPrecondSetup
    */
   void setMaxSubSetupCalls(const int maxsub);
   
   /**
    * Set values of omega_min and omega_max scalars used by nonlinear
    * residual monitoring algorithm.
    *
    *  Defaults is [0.00001 and 0.9]
    */
   void setResidualMonitoringParams(const double omega_min, const double omega_max);

   /**
    * Set constant value used by residual monitoring algorithm. If
    * omega=0, then it is estimated using omega_min and
    * omega_max.
    *
    * Default is 0.0.
    */
   void setResidualMonitoringConstant(const double omega);

   /** 
    * Set flag controlling whether or not the value * of eps is
    * bounded below by 0.01*fnormtol.
    * 
    *       FALSE constrain value of eps by setting to the following:
    *                eps = MAX{0.01*fnormtol, eps}
    *
    *       TRUE do notconstrain value of eps 
    *
    * Default is FALSE
   */
   void setNoMinEps(const bool flag);

   
   /**
    * Set maximum number of beta condition failures in the line search algorithm.
    * 
    * Default is [MXNBCF_DEFAULT] (defined in kinsol_impl.h)
    */
   void setMaxBetaFails(const int max_beta_fails);

   /**
    * Flag controlling whether or not the KINSol routine makes an
    * initial call to the linearl solver setup routine.
    * Default is false.
    */
   void setNoInitialSetup(const bool flag);


   /**
    * Flag controlling whether or not the nonlinear residual
    * monitoring schemes is used to control Jacobian updating Default
    * is FALSE.
    */
   void setNoResidualMonitoring(const bool flag);

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
   virtual void printClassData(std::ostream& os) const;

   /*
    * Open KINSOL log file, allocate main memory for KINSOL and initialize
    * KINSOL memory record.  KINSOL is initialized based on current state
    * of solver parameter data members.  If any solver parameters have 
    * changed since last initialization, this function will be automatically
    * invoked at next call to solver. 
    */
   void initializeKINSOL();


private:
   /*
    * Static member functions for linkage with KINSOL routines.
    * See header file for KINSOLAbstractFunctions for more information.
    */ 
   static int KINSOLFuncEval(N_Vector soln,
                              N_Vector fval,
                              void* my_solver);


   static int KINSOLPrecondSet(N_Vector uu, 
			       N_Vector uscale,
			       N_Vector fval, 
			       N_Vector fscale,
			       void *my_solver, 
			       N_Vector vtemp1,
			       N_Vector vtemp2);


   static int KINSOLPrecondSolve(N_Vector uu, 
				 N_Vector uscale, 
				 N_Vector fval, 
				 N_Vector fscale, 
                                 N_Vector vv, 
				 void *my_solver,
				 N_Vector vtemp);

   static int KINSOLJacobianTimesVector(N_Vector v, 
					N_Vector Jv,
					N_Vector uu, 
					int *new_uu, 
					void *my_solver);


   std::string d_object_name;

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
   SundialsAbstractVector *d_solution_vector;

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
   void*       d_kin_mem;                       // KINSOL memory structure
   FILE*       d_kinsol_log_file;               // KINSOL message log file
   std::string d_kinsol_log_file_name;          // KINSOL log file name

   /*
    * Nonlinear solution and residual scaling vectors, and integer flags
    * to determine ownership of scaling vectors.
    */
   SundialsAbstractVector *d_soln_scale;
   bool d_my_soln_scale_vector;
   SundialsAbstractVector *d_fval_scale;
   bool d_my_fval_scale_vector;

   /*
    * Constraints on nonlinear solution vector.
    */
   SundialsAbstractVector *d_constraints; 

   /*
    * Integer flag indicating whether KINSOL needs initialization 
    * when solver is called.
    */
   int d_KINSOL_needs_initialization;

   /*
    * KINSOL nonlinear and linear solver parameters
    */
   int d_krylov_dimension;  // maximum krylov dimension
   int d_max_restarts;      // max. num. of linear solver restarts allowed
   int d_max_solves_no_set; // max. num. of steps calling preconditioner
                            // without resetting preconditioner

   int d_max_iter;          // maximum number of nonlinear iterations
   int d_max_newton_step;   // maximum scaled length of Newton step

   int    d_global_strategy;     // globalization method for Newton steps.
   double d_residual_tol;        // stop tol. on scaled nonlinear residual
   double d_step_tol;            // stop tol. on consecutive step difference

   double d_maxsub;              // number of nonlinear iterations
				 // between checks by the nonlinear
				 // residual monitoring alg

   
   int d_no_initial_setup;       // intitial setup
   int d_no_residual_monitoring; // residual monitoring to control Jacobian update


   double d_omega_min;           // residual monitoring alg params
   double d_omega_max;
   double d_omega;

   int d_no_min_eps;             //  eps is bounded
   
   int d_max_beta_fails;         // maximum number of beta condition failures 

   // flag indicating which method to use to compute the value of the
   // eta coefficient used in the calculation of the linear solver
   // convergence tolerance:
   int d_eta_choice;

   // KINSetEtaConstValue constant value of eta - use with
   // KIN_ETACONSTANT option
   double d_eta_constant;

   // values of eta_gamma (egamma) and eta_alpha (ealpha) coefficients
   // use with KIN_ETACHOICE2
   double d_eta_gamma;
   double d_eta_alpha;

   // real scalar equal to realative error in computing F(u)
   double d_relative_function_error;

   // level of verbosity of output
   int d_print_level; 

};

}
}
#endif
#endif
