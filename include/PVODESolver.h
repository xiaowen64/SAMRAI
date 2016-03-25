//
// File:        PVODESolver.h
// Package:     SAMRAI solvers
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Wrapper class for PVODE solver function calls and data
//

#ifndef included_solv_PVODESolver
#define included_solv_PVODESolver

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifdef HAVE_PVODE

#ifndef included_String
#include <string>
using namespace std;
#define included_String
#endif

#ifndef included_tbox_IOStream
#define included_tbox_IOStream
#include "tbox/IOStream.h"
#endif
#ifndef included_solv_PVODEAbstractFunctions
#include "PVODEAbstractFunctions.h"
#endif
#ifndef included_solv_PVodeTrioAbstractVector
#include "PVodeTrioAbstractVector.h"
#endif

// PVODE includes
#ifndef included_cvode_h
#define included_cvode_h
extern "C" {
#include "cvode.h"
}
#endif


namespace SAMRAI {
   namespace solv {


/*!
 * @brief Class PVODESolver serves as a C++ wrapper for the PVODE 
 * ordinary differential equation solver package.  
 *
 * It is intended to be 
 * sufficiently generic to be used independently of the SAMRAI framework.
 * This class declares one private static member function to link the
 * user-defined routine for right-hand side function evaluation and
 * two private statice member functions to link the user-defined 
 * preconditioner setup and solve routines.  The implementation of these
 * functions is defined by the user in a subclass of the abstract base 
 * class PVODEAbstractFunctions.  The vector objects used within the 
 * solver are given in a subclass of the abstract class 
 * PVodeTrioAbstractVector. The PVodeTrioAbstractVector 
 * class defines the vector kernel operations required by the PVODE
 * package so that they may be easily supplied by a user who opts not 
 * to use the vector kernel supplied by the PVODE package.  (It should be
 * noted that the vector kernel used by PVODE is the same as the one
 * used by the other packages in the PVodeTrio of solvers).
 * 
 * Note that this class provides no input or restart capabilities and 
 * relies on PVODE for output reporting.  
 *
 * PVODESolver Usage:
 * 
 *
 *    -  In order to use the PVODESolver, the user must provide a
 *           concrete subclass of PVODEAbstractFunctions abstract
 *           base class which defines the evaluateRHSFunction(),
 *           CVSpgmrPrecondSet(), and CVSpgmrPrecondSolve() methods.
 *
 *    -  Solving a system of ODEs using this PVODE C++ interface 
 *           requires four main stages.  First, a PVODESolver
 *           object is created with a user-specified name and
 *           PVODEAbstractFunctions object.  Second, the
 *           user must specify the integration parameters that s/he
 *           wishes to use.  Next, the user must call the PVODESolver
 *           method initialize(solution_vector) with the 
 *           PVodeTrioAbstractVector that s/he wants to put the solution 
 *           in.  Finally, the solve() method is invoked to solve the 
 *           system of ODEs to the specified value of the independent 
 *           variable.
 *           
 *    -  The following is a list of integration parameters that
 *           must be specified by the user before calling the solve()
 *           method:
 *        
 * 
 *            - Number of equations - setNumberOfEquations(neq)
 * 
 *            - Either relative or absolute tolerance must
 *                  be set - setRelativeTolerance(relative_tolerance),
 *                  setAbsoluteTolerance(absolute_tolerance)
 * 
 *            - Initial value of independent variable -
 *                  setInitialValueOfIndependentVariable(init_time)
 *            - Final value of independent variable -
 *                  setFinalValueOfIndependentVariable(final_time
 *                      pvode_needs_initialization)
 *            - Initial condition vector - 
 *                  setInitialConditionVector(ic_vector)
 *
 *        
 *           
 *    -  The following is a list of default values for integration
 *           parameters:
 * 
 *        
 * 
 *           - @b Linear Multistep Method               
 *                BDF
 * 
 *           - @b Iteration Type                        
 *                FUNCTIONAL
 * 
 *           - @b Tolerance Type                        
 *                SS (scalar relative and scalar absolute tolerances)
 * 
 *           - @b Relative Tolerance                    
 *                0.0
 * 
 *           - @b Scalar Absolute Tolerance             
 *                0.0
 * 
 *           - @b Vector Absolute Tolerance             
 *                NULL
 * 
 *           - @b Stepping Method                       
 *                NORMAL
 * 
 *           - @b Maximum Order for Multistep Method    
 *                12 for ADAMS, 5 for BDF
 * 
 *           - @b Maximum Number of Internal Steps      
 *                500
 * 
 *           - @b Maximum Number of NIL Step Warnings   
 *                10
 * 
 *           - @b Initial Step Size                     
 *                determined by PVODE
 * 
 *           - @b Maximum Absolute Value of Step Size   
 *                infinity
 * 
 *           - @b Minimum Absolute Value of Step Size   
 *                0.0
 * 
 *           - @b CVSpgmr Preconditioning Type          
 *                NONE
 * 
 *           - @b CVSpgmr Gram Schmidt Algorithm        
 *                MODIFIED_GS
 * 
 *           - @b CVSpgmr Maximum Krylov Dimension      
 *                MIN(num_equations, CVSPGMR_MAXL=5)
 * 
 *           - @b CVSpgmr Tolerance Scale Factor        
 *                CVSPGMR_DELT = 0.05.
 * 
 *        
 *
 * 
 *
 * PVODE was developed in the Center for Applied Scientific Computing (CASC) 
 * at Lawrence Livermore National Laboratory (LLNL).  Many of the comments
 * in this class were taken verbatim from PVODE header files.  For more 
 * information about PVODE and a complete description of the operations 
 * and data structures used by this class, see S.D. Cohen and A.C. Hindmarsh, 
 * "PVODE User Guide", UCRL-MA-118618, Lawrence Livermore National 
 * Laboratory, 1994. 
 *
 * @see solv::PVODEAbstractFunctions
 * @see solv::PVodeTrioAbstractVector 
 */

class PVODESolver
{
public:
   /**
    * Constructor for PVODESolver sets default PVODE parameters 
    * and initializes the solver package with user-supplied functions
    * PVODESolver parameters may be changed later using member 
    * functions described below.  
    *
    * Notes:
    * 


    *
    *    -
    *        The solution vector is not passed into the constructor.
    *        Before the solver can be used, the initialize() function must 
    *        be called.
    *
    * 


    *
    * Assertion checks: 
    * 


    *
    *    -
    *        my_functions must not be null
    * 
    *    -
    *        object_name must not be empty.
    * 
    * 


    * 
    */
   PVODESolver(const string& object_name,
                    PVODEAbstractFunctions* my_functions,
                    const bool uses_preconditioner);

   /**
    * Virtual destructor for PVODESolver closes the 
    * PVODE log file and frees the memory allocated for the 
    * PVODE memory record.
    */
   ~PVODESolver();

   /**
    * Initialize solver with solution vector.  The solution vector is
    * required to initialize the memory record used internally within
    * PVODE.  This routine must be called before the solver can be used.
    * 
    * Assertion checks:
    * 


    *
    *    -
    *        the solution vector must not be null
    * 
    *    -
    *        the solution vector must not have already been set
    * 
    * 


    */
   void initialize(PVodeTrioAbstractVector* solution); 

   /**
    * Integrate ODE system specified t_f.  The integer return value is  
    * a termination code defined by PVODE.  The following is a table
    * of termination codes and a brief description of their meanings.
    * 
    * PVODE Termination Codes:
    * 


    *
    *    - @b SUCCESS (=0)            
    *        CVode succeeded.
    *
    *    - @b PVODE_NO_MEM (=-1)    
    *        The cvode_mem argument was NULL.
    *
    *    - @b ILL_INPUT (=-2)        
    *        One of the inputs to CVode is illegal. This    
    *        includes the situation when a component of the 
    *        error weight vectors becomes < 0 during
    *        internal time-stepping. The ILL_INPUT flag 
    *        will also be returned if the linear solver 
    *        routine CV--- (called by the user after
    *        calling CVodeMalloc) failed to set one of the
    *        linear solver-related fields in cvode_mem or 
    *        if the linear solver's init routine failed. In
    *        any case, the user should see the printed   
    *        error message for more details.
    *
    *    - @b TOO_MUCH_WORK (=-3)   
    *        The solver took maxstep internal steps but 
    *        could not reach t_f. The default value for  
    *        mxstep is MXSTEP_DEFAULT = 500.
    *
    *    - @b TOO_MUCH_ACC (=-4)    
    *        The solver could not satisfy the accuracy 
    *        demanded by the user for some internal step. 
    *
    *    - @b ERR_FAILURE (=-5)      
    *        Error test failures occurred too many times 
    *        (= MXNEF = 7) during one internal time step or 
    *        occurred with |h| = hmin.
    *
    *    - @b CONV_FAILURE (=-6)     
    *        Convergence test failures occurred too many  
    *        times (= MXNCF = 10) during one internal time
    *        step or occurred with |h| = hmin.
    *
    *    - @b SETUP_FAILURE (=-7)    
    *        The linear solver's setup routine failed in an
    *                 unrecoverable manner.
    *
    *    - @b SOLVE_FAILURE (=-8)    
    *        The linear solver's solve routine failed in an 
    *                 unrecoverable manner.
    * 
    * 


    * 
    * See cvode.h header file for more information about return values.  
    *
    * If PVODE or CVSpgmr requires re-initialization, it is 
    * automatically done before the solve.  This may be required if any 
    * of the PVODE or CVSpgmr data parameters have changed since the 
    * last call to the solver.  
    *
    * Assertion checks:
    * 


    *
    *     - 
    *        The user specified final value for the independent variable t
    *        must be greater than the specified initial value.
    *
    * 


    */
   int solve();

   /**
    * Accessor function for setting PVODE output log file name and output
    * printing options.  Output file name and options may be changed
    * throughout run as desired.
    *
    * If the file name string is empty the default file name "cvode.log" 
    * is used.
    */
   void setLogFileData(const string& log_fname = string());

   /**
    * Set PVODESolver to use my_functions as the concrete subclass 
    * of the PVODEAbstractFunctions class that defines the 
    * right-hand side evaluation and preconditioner functions.  The
    * uses_preconditioner argument indicates whether or not the
    * the user has defined preconditioner routines in their concrete
    * subclass of the PVODEAbstractFunctions class.
    *
    * Assertion checks:
    * 


    *
    *    -
    *        my_function must not be a null pointer
    *
    * 


    */
   void setPVODEFunctions(PVODEAbstractFunctions* my_functions,
                          const bool uses_preconditioner);

   /**
    * Return pointer to object that provides user-defined functions for 
    * PVODE and CVSpgmr.
    */
   PVODEAbstractFunctions* getPVODEFunctions() const;

   // Methods for setting PVODE parameters.

   /**
    * Set number of equations in system of ODEs to be solved.
    * 
    * Assertion checks:
    * 


    *
    *    -
    *        neq must be positive
    *
    * 


    */
   void setNumberOfEquations(int neq);

   /**
    * Set linear multistep method.  The user can specify either
    * ADAMS or BDF (backward differentiation formula) methods 
    * The BDF method is recommended  for stiff problems, and 
    * the ADAMS method is recommended for nonstiff problems.
    *
    * Assertion checks:
    * 


    *
    *    -
    *        linear_multistep_method must be one of ADAMS or BDF.
    *
    * 


    *
    * Note: the enumeration constants ADAMS and BDF are defined in cvode.h.
    */
   void setLinearMultistepMethod(int linear_multistep_method);

   /**
    * Set iteration type.  The user can specify either FUNCTIONAL
    * iteration, which does not require linear algebra, or a 
    * NEWTON iteration, which requires the solution of linear 
    * systems. In the NEWTON case, the user must also specify a 
    * PVODE linear solver. NEWTON is recommended in case of 
    * stiff problems.  
    *
    * Assertion checks:
    * 


    *
    *    -
    *        iteration_type must be one of FUNCTIONAL or NEWTON
    *
    * 


    *
    * Note: the enumeration constants FUNCTIONAL and NEWTON are defined 
    * in cvode.h.
    */
   void setIterationType(int iteration_type);

   /**
    * Set tolerance type.  This parameter specifies the relative 
    * and absolute tolerance types to be used. The SS tolerance type 
    * means a scalar relative and absolute tolerance, while the SV 
    * tolerance type means a scalar relative tolerance and a 
    * vector absolute tolerance (a potentially different 
    * absolute tolerance for each vector component).    
    *
    * Assertion checks:
    * 


    *
    *    -
    *        tolerance_type must be one of SS or SV
    *
    * 


    *
    * Note: the enumeration constants SS and SV are defined in cvode.h.
    */
   void setToleranceType(int tolerance_type);

   /**
    * Set the relative tolerance level.  
    *
    * Assertion checks:
    * 


    *
    *    -
    *        relative_tolerance must be greater than or equal to 0.0
    *
    * 


    *
    * Note that pure absolute tolerance can be used by
    * setting the relative tolerance to 0.  However, 
    * it is an error to simultaneously set relative and 
    * absolute tolerances to 0.
    */
   void setRelativeTolerance(double relative_tolerance);

   /**
    * Set the scalar absolute tolerance level.
    *
    * Assertion checks:
    * 


    *
    *    -
    *        absolute_tolerance must be greater than or equal to 0.0
    *
    * 


    *
    * Note that pure relative tolerance can be used by
    * setting the absolute tolerance to 0.  However, 
    * it is an error to simultaneously set relative and 
    * absolute tolerances to 0.
    */
   void setAbsoluteTolerance(double absolute_tolerance);

   /**
    * Set the vector absolute tolerance level.
    *
    * Assertion checks:
    * 


    *
    *    -
    *        absolute_tolerance must not be a null pointer
    *
    *    -
    *        each component of absolute_tolerance must be 
    *        greater than or equal to 0.0
    *
    * 


    *
    * Note that pure relative tolerance can be used by
    * setting the absolute tolerance to 0.  However, 
    * it is an error to simultaneously set relative and 
    * absolute tolerances to 0.
    */
   void setAbsoluteTolerance(PVodeTrioAbstractVector* absolute_tolerance);

   /**
    * Set stepping method to use for integration.  There are 
    * stepping methods: NORMAL and ONE_STEP.  The NORMAL
    * method has the solver take internal steps until 
    * it has reached or just passed the user specified t_f
    * parameter. The solver then interpolates in order to 
    * return an approximate value of y(t_f). The ONE_STEP 
    * option tells the solver to just take one internal step 
    * and return the solution at the point reached by that 
    * step.                       
    *
    * Assertion checks:
    * 


    *
    *    -
    *        stepping_method must be one of NORMAL or ONE_STEP
    *
    * 


    *
    * Note: the enumeration constants NORMAL and ONE_STEP are 
    * defined in cvode.h.
    */
   void setSteppingMethod(int stepping_method);

   /**
    * Set initial value for independent variable.
    */
   void setInitialValueOfIndependentVariable(double t_0);

   /**
    * Set final value for independent variable (i.e. the value of
    * independent variable to integrate the system to).  The boolean
    * argument specifies whether PVODE should be re-initialized (i.e.
    * on first step) or if we are taking subsequent steps in a 
    * sequence, in which case it is not initialized.
    */
   void setFinalValueOfIndependentVariable(double t_f,
      bool pvode_needs_initialization);

   /**
    * Set initial condition vector.
    *
    * Assertion checks:
    * 


    *
    *    -
    *        ic_vector must not be null
    *
    * 


    */
   void setInitialConditionVector(PVodeTrioAbstractVector* ic_vector);

   /**
    * Set maximum order for the linear multistep method.
    * By default, this is set to 12 for ADAMS methods and 5 for BDF 
    * methods.
    *
    * Assertion checks:
    * 


    *
    *    -
    *        max_order must be greater than or equal to 0
    *
    * 


    */
   void setMaximumLinearMultistepMethodOrder(int max_order);

   /**
    * Set maximum number of internal steps to be taken by 
    * the solver in its attempt to reach t_f.
    * By default, this is set to 500.
    *
    * Assertion checks:
    * 


    *
    *    -
    *        max_num_internal_steps must be greater than or equal to 0
    *
    * 


    */
   void setMaximumNumberOfInternalSteps(int max_num_internal_steps);

   /**
    * Set maximum number of warning messages issued by the solver
    * that (t + h == t) on the next internal step.  By default, 
    * this is set to 10.
    *
    * Assertion checks:
    * 


    *
    *    -
    *        max_num_warnings must be greater than or equal to 0
    *
    * 


    */
   void setMaximumNumberOfNilStepWarnings(int max_num_warnings);

   /**
    * Set initial step size.
    *
    * Assertion checks:
    * 


    *
    *    -
    *        init_step_size must be greater than or equal to 0.0
    *
    * 


    */
   void setInitialStepSize(double init_step_size);

   /**
    * Set maximum absolute value of step size allowed.
    * By default, there is no upper bound on the absolute value
    * of step size.
    *
    * Assertion checks:
    * 


    *
    *    -
    *        max_step_size must be greater than or equal to 0.0
    *
    * 


    */
   void setMaximumAbsoluteStepSize(double max_step_size);

   /**
    * Set minimum absolute value of step size allowed.
    * By default, this is set to 0.0.
    *
    * Assertion checks:
    * 


    *
    *    -
    *        min_step_size must be greater than or equal to 0.0
    *
    * 


    */
   void setMinimumAbsoluteStepSize(double min_step_size);

   // Methods for setting CVSpgmr parameters.

   /**
    * Set the preconditioning type to be used by CVSpgmr.
    * This must be one of the four enumeration constants
    * NONE, LEFT, RIGHT, or BOTH defined in iterativ.h.
    * These correspond to no preconditioning, left preconditioning only,
    * right preconditioning only, and both left and right
    * preconditioning, respectively.
    *
    * Assertion Checks:
    * 


    *
    *	 -
    *	    precondition_type must be one of NONE, LEFT, RIGHT, or BOTH.
    *
    * 


    */
   void setPreconditioningType(int precondition_type);

   /**
    * Set the Gram-Schmidt orthogonalization type to be used by CVSpgmr.
    * This must be one of the two enumeration constants MODIFIED_GS
    * or CLASSICAL_GS defined in iterativ.h. These correspond to
    * using modified Gram-Schmidt and classical Gram-Schmidt, respectively.
    *
    * Assertion Checks:
    * 


    *
    *	 -
    *	    gs_type must be one of CLASSICAL_GS or MODIFIED_GS.
    *
    * 


    */
   void setGramSchmidtType(int gs_type);

   /**
    * Set the maximum Krylov dimension to be used by CVSpgmr.
    * This is an optional input to the CVSPGMR solver. Pass 0 to
    * use the default value MIN(num_equations, CVSPGMR_MAXL=5).
    *
    * Assertion Checks:
    * 


    *
    *	 -
    *	    max_krylov_dim must be nonnegative
    *
    * 


    */
   void setMaxKrylovDimension(int max_krylov_dim);

   /**
    * Set the factor by which the tolerance on the nonlinear
    * iteration is multiplied to get a tolerance on the linear iteration.
    * This is an optional input to the CVSPGMR solver. Pass 0 to
    * use the default value CVSPGMR_DELT = 0.05.
    *
    * Assertion Checks:
    * 


    *
    *	 -
    *	    tol_scale_factor must be nonnegative
    *
    * 


    */
   void setCVSpgmrToleranceScaleFactor(double tol_scale_factor);

   /**
    * Get solution vector.
    */
   PVodeTrioAbstractVector* getSolutionVector() const;

   /**
    * Get k-th derivative vector at the specified value of the
    * independent variable, t.  The integer return value is 
    * return code the PVODE CVodeDky() function.  The following is a table
    * of termination codes and a brief description of their meanings.
    * 
    * CVodeDky Return Codes:
    * 


    *
    *    - @b OKAY (=0)           
    *        CVodeDky succeeded.
    *
    *    - @b BAD_K (=-1)        
    *        
    *    - @b BAD_T (=-2)        
    *        
    *    - @b BAD_DKY (=-3)      
    *        
    *    - @b DKY_NO_MEM (=-4)  
    *        
    * 


    * 
    * Important Notes:
    * 


    *
    *    -
    *       t must lie in the interval [t_cur - h, t_cur] 
    *       where t_cur is the current internal time reached
    *       and h is the last internal step size successfully 
    *       used by the solver.
    *
    *    -
    *       k may take on value 0, 1, . . . q where q is the order 
    *       of the current linear multistep method being used. 
    * 
    *    -
    *       the dky vector must be allocated by the user.       
    *
    *    -
    *       it is only leagal to call this method after a 
    *       successful return from the solve() method.
    * 
    * 


    *
    */
   int getDkyVector(double t, int k, PVodeTrioAbstractVector* dky) const;

   /**
    * Get actual value of the independent variable that PVODE integrated
    * to (i.e. the value of t that actually corresponds to the solution
    * vector y).
    */
   double getActualFinalValueOfIndependentVariable() const;

   /**
    * Set PVODE to collect statistics.  This must be called before
    * solve() is invoked.
    */
   void turnOnPVODEStatisticsCollection();

   /**
    * Print PVODE and CVSpgmr statistics.
    */
   void printStatistics(ostream& os) const;

   /**
    * Print PVODE statistics to the stream.
    *
    * The abbreviations printed out refer to the following
    * quantities:
    * 


    * 
    *    - @b lenrw            
    *       size (in double words) of memory used for doubles 
    * 
    *    - @b leniw            
    *       size (in integer words) of memory used for integers
    * 
    *    - @b nst              
    *       cumulative number of internal steps taken by solver
    * 
    *    - @b nfe              
    *       number of right-hand side function evaluations
    * 
    *    - @b nni              
    *       number of NEWTON iterations performed 
    * 
    *    - @b nsetups          
    *       number of calls made to linear solver's setup routine
    * 
    *    - @b netf             
    *       number of local error test failures
    * 
    *    - @b ncfn             
    *       number of nonlinear convergence failures
    * 
    *    - @b qu               
    *       order used during the last internal step
    * 
    *    - @b qcur             
    *       order to be used on the next internal step
    * 
    *    - @b hu               
    *       step size for the last internal step
    * 
    *    - @b hcur             
    *       step size to be attempted on the next internal step
    * 
    *    - @b tcur             
    *       current internal value of t reached by the solver
    * 
    *    - @b tolsf            
    *       suggested tolerance scaling factor
    *
    * 


    */
   void printPVODEStatistics(ostream& os) const;

   // PVODE optional return values.

   /**
    * Return the cumulative number of internal steps taken by
    * the solver.  
    *
    * Note: if the solver was not set to collect statistics, 
    * a value of -1 is returned.
    */
   int getNumberOfInternalStepsTaken() const;

   /**
    * Return the number of calls to the right-hand side function.
    *
    * Note: if the solver was not set to collect statistics, 
    * a value of -1 is returned.
    */
   int getNumberOfRHSFunctionCalls() const;

   /**
    * Return the number of calls made to linear solver setup
    * routines.
    *
    * Note: if the solver was not set to collect statistics, 
    * a value of -1 is returned.
    */
   int getNumberOfLinearSolverSetupCalls() const;

   /**
    * Return the number of NEWTON iterations performed.
    *
    * Note: if the solver was not set to collect statistics, 
    * a value of -1 is returned.
    */
   int getNumberOfNewtonIterations() const;

   /**
    * Return the number of nonlinear convergence failures that have 
    * occurred.
    *
    * Note: if the solver was not set to collect statistics, 
    * a value of -1 is returned.
    */
   int getNumberOfNonlinearConvergenceFailures() const;

   /**
    * Return the number of local error test failures.
    *
    * Note: if the solver was not set to collect statistics, 
    * a value of -1 is returned.
    */
   int getNumberOfLocalErrorTestFailures() const;

   /**
    * Return the order of the linear multistep method used during
    * the last internal step.
    *
    * Note: if the solver was not set to collect statistics, 
    * a value of -1 is returned.
    */
   int getOrderUsedDuringLastInternalStep() const;

   /**
    * Return the order of the linear multistep method to be used during
    * the next internal step.
    *
    * Note: if the solver was not set to collect statistics, 
    * a value of -1 is returned.
    */
   int getOrderToBeUsedDuringNextInternalStep() const;

   /**
    * Return the size (in LLNL_REAL words) of memory used
    * for LLNL_REALS.
    *
    * Note: if the solver was not set to collect statistics, 
    * a value of -1 is returned.
    */
   int getPVODEMemoryUsageForDoubles() const;

   /**
    * Return the size (in integer words) of memory used
    * for integers.
    *
    * Note: if the solver was not set to collect statistics, 
    * a value of -1 is returned.
    */
   int getPVODEMemoryUsageForIntegers() const;

   /**
    * Return the step size for the last internal step.
    *
    * Note: if the solver was not set to collect statistics, 
    * a value of -1 is returned.
    */
   double getStepSizeForLastInternalStep() const;

   /**
    * Return the step size to be used in the next internal step.
    *
    * Note: if the solver was not set to collect statistics, 
    * a value of -1 is returned.
    */
   double getStepSizeForNextInternalStep() const;

   /**
    * Return the current internal value of the independent
    * variable reached by the solver.
    *
    * Note: if the solver was not set to collect statistics, 
    * the minimum double value (as defined in float.h) is 
    * returned.
    */
   double getCurrentInternalValueOfIndependentVariable() const;

   /**
    * Return the suggested tolerance scaling factor.
    *
    * Note: if the solver was not set to collect statistics, 
    * a value of -1 is returned.
    */
   double getPVODESuggestedToleranceScalingFactor() const;

   // CVSpgmr optional return values.

   /**
    * Print CVSpgmr statistics to the stream.
    *
    * The abbreviations printed out refer to the following
    * quantities:
    * 


    *
    *    - @b spgmr_lrw        
    *      size (in double words) of memory used for doubles
    *
    *    - @b spgmr_liw        
    *       size (in integer words) of memory used for integers
    *
    *    - @b nli              
    *       number of linear iterations
    *
    *    - @b ncfl             
    *       number of linear convergence failures
    *
    *    - @b npe              
    *       number of preconditioner evaluations
    *
    *    - @b nps              
    *       number of calls to CVSpgmrPrecondSolve()
    *
    * 


    */
   void printCVSpgmrStatistics(ostream& os) const;

   /**
    * Return the number of preconditioner evaluations.
    */
   int getNumberOfPreconditionerEvaluations() const;

   /**
    * Return the number of linear iterations.
    */
   int getNumberOfLinearIterations() const;

   /**
    * Return the number of CVSpgmrPrecondSolve() calls.
    */
   int getNumberOfPrecondSolveCalls() const;

   /**
    * Return the number of linear convergence failures.
    */
   int getNumberOfLinearConvergenceFailures() const;

   /**
    * Return the size (in double words) of memory used for doubles.
    */
   int getCVSpgmrMemoryUsageForDoubles() const;

   /**
    * Return the size (in integer words) of memory used for integers.
    */
   int getCVSpgmrMemoryUsageForIntegers() const;

   /**
    * Print out all data members for this object.
    */
   void printClassData(ostream& os) const;

private:
   /*
    * Static member function for linkage with PVODE routines.
    */ 
   static void PVODERHSFuncEval(int neq,
                                double t,
                                PVodeTrioAbstractVector* y,
                                PVodeTrioAbstractVector* y_dot,
                                void* my_solver);

   /*
    * Static member functions for linkage with CVSpgmr routines.
    */ 
   static int CVSpgmrPrecondSet(int neq,
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
				void *my_solver,
				PVodeTrioAbstractVector* vtemp1,
				PVodeTrioAbstractVector* vtemp2,
				PVodeTrioAbstractVector* vtemp3);

   static int CVSpgmrPrecondSolve(int neq,
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
                                  void *my_solver,
                                  PVodeTrioAbstractVector* z);


   /*
    * Open PVODE log file, allocate main memory for PVODE and initialize
    * PVODE memory record.  PVODE is initialized based on current state
    * of solver parameter data members.  If any solver parameters have 
    * changed since last initialization, this function will be automatically
    * invoked at next call to the solve() method.  Also, if NEWTON iteration 
    * is specified, this method also initializes the CVSpgmr linear solver. 
    *
    * Assertion checks:
    * 


    *
    *    -
    *       the solution vector must have already been set.
    *
    * 


    *
    */
   void initializePVODE();

   string d_object_name;

   /*
    * The following data members are input or set to default values in
    * the PVODESolver constructor.  Many of these can be altered at
    * any time through class member functions.  When this occurs,
    * PVODE may need to be re-initialized (e.g., if the linear solver
    * changes, PVODE must change its memory record).  In this case,
    * the initializePVODE() member function is invoked in the next 
    * call to solve().
    */

   /*
    * Solution vector.
    */
   PVodeTrioAbstractVector* d_solution_vector;

   /*
    * Pointer to object which provides user-supplied functions to PVODE
    * and CVSpgmr.
    */
   PVODEAbstractFunctions* d_pvode_functions;

   /*
    * PVODE memory record.
    */
   CVodeMem  d_pvode_mem;                    // PVODE memory structure

   /*
    * PVODE log file information.
    */
   FILE*     d_pvode_log_file;               // PVODE message log file
   string    d_pvode_log_file_name;          // PVODE log file name

   /*
    * ODE parameters.
    */
   int d_neq;           // number of dependent variables in ODE system 
   double d_t_0;        // initial value for independent variable
   double d_user_t_f;   // user-specified final value for independent variable
   double d_actual_t_f; // actual final value of indep. variable after a step
   PVodeTrioAbstractVector* d_ic_vector;

   /*
    * ODE integration parameters.
    */
   int d_linear_multistep_method;
   int d_iteration_type;
   int d_tolerance_type;
   double d_relative_tolerance;
   bool d_use_scalar_absolute_tolerance;
   double d_absolute_tolerance_scalar;
   PVodeTrioAbstractVector* d_absolute_tolerance_vector;
   int d_stepping_method;

   /*
    * Optional PVODE parameters.
    */
   bool      d_use_optional_data;   
   long int  d_int_optional_data[OPT_SIZE];
   double    d_real_optional_data[OPT_SIZE];

   /*
    * CVSpgmr parameters
    */
   int d_precondition_type;
   int d_gram_schmidt_type;
   int d_max_krylov_dim;
   double d_tol_scale_factor;

   /*
    * Boolean flag indicating whether PVODE needs initialization 
    * when solver is called.
    */
   bool d_PVODE_needs_initialization;

   /*
    * Boolean flag indicating whether user-supplied preconditioner
    * routines are provided in the concrete subclass of 
    * PVODEAbstractFunctions.
    */
   bool d_uses_preconditioner;
};


}
}

#endif
#endif
