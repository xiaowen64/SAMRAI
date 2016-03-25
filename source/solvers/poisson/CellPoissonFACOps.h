/*
 * File:        CellPoissonFACOps.h
 * Copyright:   (c) 1997-2005 The Regents of the University of California
 * Revision:    $Revision: 453 $
 * Modified:    $Date: 2005-06-16 10:19:28 -0700 (Thu, 16 Jun 2005) $
 * Description: Operator class for cell-centered scalar Poisson using FAC
 */

#ifndef included_solv_CellPoissonFACOps
#define included_solv_CellPoissonFACOps

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_hier_CoarseFineBoundary
#include "CoarseFineBoundary.h"
#endif

#ifndef included_hier_Patch
#include "Patch.h"
#endif

#ifndef included_hier_PatchHierarchy
#include "PatchHierarchy.h"
#endif

#ifndef included_hier_PatchLevel
#include "PatchLevel.h"
#endif

#ifndef included_hier_IntVector
#include "IntVector.h"
#endif

#ifndef included_hier_Box
#include "Box.h"
#endif

#ifndef included_math_HierarchyCellDataOpsReal
#include "HierarchyCellDataOpsReal.h"
#endif

#ifndef included_math_HierarchySideDataOpsReal
#include "HierarchySideDataOpsReal.h"
#endif

#ifndef included_pdat_CellData
#include "CellData.h"
#endif

#ifndef included_pdat_CellVariable
#include "CellVariable.h"
#endif

#ifndef included_pdat_CellDoubleConstantRefine
#include "CellDoubleConstantRefine.h"
#endif

#ifndef included_pdat_OutersideData
#include "OutersideData.h"
#endif

#ifndef include_pdat_OutersideVariable
#include "OutersideVariable.h"
#endif

#ifndef included_pdat_SideData
#include "SideData.h"
#endif

#ifndef included_pdat_SideVariable
#include "SideVariable.h"
#endif

#ifndef included_solv_CartesianRobinBcHelper
#include "CartesianRobinBcHelper.h"
#endif

#ifndef included_solv_FACPreconditioner
#include "FACPreconditioner.h"
#endif

#ifndef included_solv_FACOperatorStrategy
#include "FACOperatorStrategy.h"
#endif

#ifndef included_solv_RobinBcCoefStrategy
#include "RobinBcCoefStrategy.h"
#endif

#ifndef included_solv_CellPoissonHypreSolver
#include "CellPoissonHypreSolver.h"
#endif

#ifndef included_solv_PoissonSpecifications
#include "PoissonSpecifications.h"
#endif

#ifndef included_xfer_CoarsenAlgorithm
#include "CoarsenAlgorithm.h"
#endif

#ifndef included_xfer_CoarsenOperator
#include "CoarsenOperator.h"
#endif

#ifndef included_xfer_RefineAlgorithm
#include "RefineAlgorithm.h"
#endif

#ifndef included_xfer_RefineOperator
#include "RefineOperator.h"
#endif

#ifndef included_solv_SAMRAIVectorReal
#include "SAMRAIVectorReal.h"
#endif

#ifndef included_hier_VariableContext
#include "VariableContext.h"
#endif

#ifndef included_tbox_Database
#include "tbox/Database.h"
#endif

#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif

#ifndef included_tbox_Timer
#include "tbox/Timer.h"
#endif

#ifndef included_String
#include <string>
using namespace std;
#define included_String
#endif

#ifndef included_xfer_CoarsenSchedule
#include "CoarsenSchedule.h"
#endif

#ifndef included_xfer_RefineSchedule
#include "RefineSchedule.h"
#endif

namespace SAMRAI {
    namespace solv {


/*!
 * @brief FAC operator class to solve Poisson's equation on a SAMR grid,
 * using cell-centered, second-order finite-volume method, with Robin
 * boundary conditions.
 *
 * This class provides operators that are used by the FAC
 * preconditioner FACPreconditioner<DIM>.
 * It is used to solve the scalar Poisson's equation using a cell-centered
 * second-order finite-volume discretization.
 * It is designed to provide all operations specific to
 * the scalar Poisson's equation,
 * @f[ \nabla \cdot D \nabla u + C u = f @f]
 * (see PoissonSpecifications) where
 * - C, D and f are indpendent of u
 * - C is a cell-centered scalar field
 * - D is the @em diffusion @em coefficients, stored on faces
 * - f is a cell-centered scalar function
 *
 * You are left to provide the source function, initial guess, etc.,
 * by specifying them in specific forms.
 *
 * This class provides:
 * -# 5-point (second order), cell-centered stencil operations
 *    for the discrete Laplacian.
 * -# Red-black Gauss-Seidel smoothing.
 * -# Provisions for working Robin boundary conditions
 *    (see RobinBcCoefStrategy).
 *
 * This class is meant to provide the Poisson-specific operator
 * used by the FAC preconditioner, FACPreconditioner<DIM>.
 * To use the preconditioner with this class, you will have to provide:
 * -# The solution vector SAMRAIVectorReal<DIM>,
 *    with appropriate norm weighting for the cell-centered AMR mesh.
 *    This class provides the function computeVectorWeights()
 *    to help with computing the appropriate weights.
 *    Since this is for a scalar equation, only the first depth
 *    of the first component of the vectors are used.
 *    All other parts are ignored.
 * -# The source vector SAMRAIVectorReal<DIM> for f.
 * -# A PoissonSpecifications objects to specify
 *    the cell-centered scalar field C and the side-centered
 *    diffusion coefficients D
 * -# The boundary condition specifications in terms of the coefficients
 *    @f$ \alpha @f$, @f$ \beta @f$ and @f$ \gamma @f$ in the
 *    Robin formula @f$  \alpha u + \beta u_n = \gamma @f$ applied on the
 *    boundary faces.  See RobinBcCoefStrategy<DIM>.
 *
 * This class allocates and deallocates only its own scratch data.
 * Other data that it manipuates are passed in as function arguments.
 * Hence, it owns none of the solution vectors, error vectors,
 * diffusion coefficient data, or any such things.
 *
 * Input Examples
 * @verbatim
 * coarse_solver_choice = "hypre"    // see setCoarsestLevelSolverChoice()
 * coarse_solver_tolerance = 1e-14   // see setCoarsestLevelSolverTolerance()
 * coarse_solver_max_iterations = 10 // see setCoarsestLevelSolverMaxIterations()
 * smoothing_choice = "redblack"     // see setSmoothingChoice()
 * cf_discretization = "Ewing"       // see setCoarseFineDiscretization()
 * prolongation_method = "LINEAR_REFINE" // see setProlongationMethod()
 * hypre_solver = { ... }            // tbox::Database for initializing Hypre solver
 * @endverbatim
 */
template<int DIM> class CellPoissonFACOps
  : public FACOperatorStrategy<DIM>
{

public:

   /*!
    * @brief Constructor.
    *
    * If you want standard output and logging,
    * pass in valid pointers for those streams.
    * @param object_name Ojbect name
    * @param database Input database
    */
   CellPoissonFACOps(
      const string &object_name=string() ,
      tbox::Pointer<tbox::Database> database=NULL );

   /*!
    * @brief Destructor.
    *
    * Deallocate internal data.
    */
   ~CellPoissonFACOps<DIM>(void);



   /*!
    * @brief Set the scalar Poisson equation specifications.
    */
   void setPoissonSpecifications(
      const PoissonSpecifications &spec );


   /*!
    * @brief Enable logging.
    *
    * By default, logging is disabled.  The logging flag is
    * propagated to the major components used by this class.
    */
   void enableLogging( bool enable_logging );



   //@{
   /*!
    * @name Functions for setting solver mathematic algorithm controls
    */

   /*!
    * @brief Set the choice of smoothing algorithms.
    *
    * Current smoothing choices are:
    * - "redblack": Red-black Gauss-Seidel smoothing.
    */
   void setSmoothingChoice( const string &smoothing_choice );

   /*!
    * @brief Set coarse level solver.
    *
    * Select from these:
    * - @c "redblack" (red-black smoothing until convergence--very slow!)
    * - @c "hypre" (only if the HYPRE library is available).
    */
   void setCoarsestLevelSolverChoice( const string &choice );

   /*!
    * @brief Set tolerance for coarse level solve.
    *
    * If the coarse level solver requires a tolerance (currently, they all do),
    * the specified value is used.
    */
   void setCoarsestLevelSolverTolerance( double tol );

   /*!
    * @brief Set max iterations for coarse level solve.
    *
    * If the coarse level solver requires a max iteration limit
    * (currently, they all do), the specified value is used.
    */
   void setCoarsestLevelSolverMaxIterations( int max_iterations );

   /*!
    * @brief Set the coarse-fine boundary discretization method.
    *
    * Specify the @c op_name string which will be passed to
    * xfer::Geometry<DIM>::lookupRefineOperator() to get the operator
    * for setting fine grid ghost cells from the coarse grid.
    * Note that chosing this operator implicitly choses the
    * discretization method at the coarse-fine boundary.
    *
    * There is one important instance where this string is
    * @em not passed to xfer::Geometry<DIM>::lookupRefineOperator.
    * If this variable is set to "Ewing", Ewing's coarse-fine
    * discretization is used (a constant refinement is performed,
    * and the flux is later corrected to result in Ewing's scheme).
    * For a reference to Ewing's discretization method, see
    * "Local Refinement Techniques for Elliptic Problems on Cell-Centered
    * Grids, I. Error Analysis", Mathematics of Computation, Vol. 56, No. 194,
    * April 1991, pp. 437-461.
    *
    * @param method String selecting the coarse-fine discretization method.
    */
   void setCoarseFineDiscretization( const string &coarsefine_method );

   /*!
    * @brief Set the name of the prolongation method.
    *
    * Specify the @c op_name string which will be passed to
    * xfer::Geometry<DIM>::lookupRefineOperator() to get the operator
    * for prolonging the coarse-grid correction.
    *
    * By default, "CONSTANT_REFINE" is used.  "LINEAR_REFINE" seems to 
    * to lead to faster convergence, but it does NOT satisfy the Galerkin
    * condition.
    *
    * Prolonging using linear refinement requires a Robin bc
    * coefficient implementation that is capable of delivering
    * coefficients for non-hierarchy data, because linear refinement
    * requires boundary conditions to be set on temporary levels.
    *
    * @param prolongation_method String selecting the coarse-fine
    *        discretization method.
    */
   void setProlongationMethod( const string &prolongation_method );

#ifdef HAVE_HYPRE
   /*!
    * @brief Set whether to use Hypre's PFMG algorithm instead of the
    * SMG algorithm.
    *
    * This flag affects the Hypre solver (used to solve the coarsest level).
    * The flag is used to select which of HYPRE's linear solver algorithms
    * to use if true, the semicoarsening multigrid algorithm is used, and if
    * false, the ``PF'' multigrid algorithm is used.
    * By default, the SMG algorithm is used.
    *
    * This setting has effect only when Hypre is chosen for the coarsest
    * level solver.  See setCoarsestLevelSolverChoice().
    *
    * Changing the algorithm must be done before initializing the solver
    * state and must NOT be done while the state is initialized
    * (the program will exit), as that would corrupt the state.
    */
   void setUseSMG( bool use_smg );
#endif

   //@}


   //@{
   /*!
    * @name Functions for setting patch data indices and coefficients
    */

   /*!
    * @brief Set the scratch patch data index for the flux.
    *
    * The use of this function is optional.
    * The patch data index should be a pdat::SideData<DIM> type of variable.
    * If the flux id is -1 (the default initial value), scratch space
    * for the flux is allocated as needed and immediately deallocated
    * afterward, level by level.  If you have space preallocated for
    * flux and you would like that to be used, set flux id to the
    * patch data index of that space.
    */
   void setFluxId( int flux_id );

   //@}


   /*!
    * @brief Provide an implementation for getting the
    * physical bc coefficients
    *
    * If your solution is fixed at the physical boundary
    * ghost cell centers AND those cells have the correct
    * values before entering solveSystem(), you may use a
    * GhostCellRobinBcCoefs<DIM> object.
    *
    * If your solution is @b not fixed at the ghost cell centers,
    * the ghost cell values will change as the interior
    * cell values change.  In those cases, the flexible
    * Robin boundary conditions are applied.  You must
    * call this function to provide the implementation for
    * determining the boundary condition coefficients.
    *
    * @param physical_bc_coef tbox::Pointer to an object that can
    *        set the Robin bc coefficients.
    */
   void setPhysicalBcCoefObject(
      const RobinBcCoefStrategy<DIM> *physical_bc_coef );


   //@{

   /*!
    * @name Functions for checking validity and correctness of state.
    */
  
   /*!
    * @brief Check validity and correctness of input patch data indices.
    *
    * Descriptors checked:
    * -# Diffusion coefficient (see setDiffcoefId())
    * -# Flux (see setFluxId())
    * -# Source (see setScalarFieldId())
    */
   void checkInputPatchDataIndices() const;

   //@}


   /*!
    * @brief Set weight appropriate for computing vector norms.
    *
    * If you this function to set the weights used when you
    * SAMRAIVectorReal<DIM>::addComponent, you can use the
    * vector norm functions of SAMRAIVectorReal<DIM>, and
    * the weights will be used to blank out coarse grid
    * regions under fine grids.
    *
    * The weights computed are specific to the cell-centered
    * discretization used by this class.  The weight is equal
    * to the cell volume if the cell has not been refined,
    * and zero if it has.
    *
    * This function is state-independent.  All inputs are in
    * the argument list.
    *
    * @param hierarchy Hierarchy configuration to compute weights for
    * @param weight_id hier::Patch data index of the weight
    * @param coarsest_ln Coarsest level number.  Must be included
    *        in hierarchy.  Must not be greater than @c finest_ln.
    *        Default to 0.
    * @param finest_ln Finest level number.  Must be included
    *        in hierarchy.  Must not be less than @c coarsest_ln.
    *        Default to finest level in @c hierarchy.
    */
   void computeVectorWeights(
      tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy ,
      int weight_id ,
      int coarsest_ln=-1 ,
      int finest_ln=-1
   ) const;


   /*!
    * @brief Set the FAC preconditioner that will be using this object.
    *
    * The FAC preconditioner is accessed to get convergence data during
    * the cycle postprocessing step.  It is optional.
    */
   void setPreconditioner( const FACPreconditioner<DIM> *preconditioner );


   /*!
    * @brief function to compute flux, using general diffusion
    * coefficient data.
    *
    * Recall that this solver class discretizes the PDE 
    * @f[ \nabla \cdot D \nabla u + C u = f @f] on an AMR grid.  This member
    * function allows users of this solver class to compute gradient 
    * terms, @f[ D \nabla w @f], in their code in a manner consistent with the
    * solver discretization.   In particular, when solving PDE systems, it may 
    * be necessary to discretize the gradient operator appearing in equations 
    * not treated by the solver class in the same way as those treated by this 
    * class.  These funtions allow users to do this easily.  The divergence
    * operator used in this solver is the standard sum of centered differences
    * involving flux terms on the cell sides computed by these routines.
    *
    * Note that the patch must exist on a level in an AMR hierarchy so that
    * the discretization can be computed properly at the coarse-fine interface. 
    * Poisson coefficients C and D must exist on the patch, if they are variable.
    * Also, calling this function does not affect the internal solver state in any 
    * way.  However, the solver must be fully initialized before it is called and care 
    * should be exercised to pass arguments so that the solver solution quantity and
    * other internal solver quantities are not adversely affected.
    *
    * @param patch patch on which computation will take place
    * @param ratio_to_coarser_level refinement ratio from coarser level to level
    *                               on which patch lives; if current patch level
    *                               is level zero, this is ignored
    * @param w_data cell-centered data 
    * @param Dgradw_data side-centered flux data (i.e., D (grad w))
    */
   void computeFluxOnPatch (
      const hier::Patch<DIM> &patch ,
      const hier::IntVector<DIM> &ratio_to_coarser_level,
      const pdat::CellData<DIM,double> &w_data ,
      pdat::SideData<DIM,double> &Dgradw_data ) const;



   //@{ @name FACOperatorStrategy<DIM> virtuals

   virtual void restrictSolution(
      const SAMRAIVectorReal<DIM,double> &s ,
      SAMRAIVectorReal<DIM,double> &d ,
      int dest_ln );
   virtual void restrictResidual(
      const SAMRAIVectorReal<DIM,double> &s ,
      SAMRAIVectorReal<DIM,double> &d ,
      int dest_ln );

   virtual void prolongErrorAndCorrect(
      const SAMRAIVectorReal<DIM,double> &s ,
      SAMRAIVectorReal<DIM,double> &d ,
      int dest_ln );

   virtual void smoothError(
      SAMRAIVectorReal<DIM,double> &error ,
      const SAMRAIVectorReal<DIM,double> &residual ,
      int ln ,
      int num_sweeps );

   virtual void solveCoarsestLevel(
      SAMRAIVectorReal<DIM,double> &error ,
      const SAMRAIVectorReal<DIM,double> &residual ,
      int ln );

   virtual void computeCompositeResidualOnLevel(
      SAMRAIVectorReal<DIM,double> &residual ,
      const SAMRAIVectorReal<DIM,double> &solution ,
      const SAMRAIVectorReal<DIM,double> &rhs ,
      int ln ,
      bool error_equation_indicator );

   virtual double computeResidualNorm(
      const SAMRAIVectorReal<DIM,double> &residual ,
      int fine_ln ,
      int coarse_ln );

   virtual void initializeOperatorState(
      const SAMRAIVectorReal<DIM,double> &solution ,
      const SAMRAIVectorReal<DIM,double> &rhs );

   virtual void deallocateOperatorState();

   virtual void postprocessOneCycle(
      int iteration_num ,
      const SAMRAIVectorReal<DIM,double> &current_soln ,
      const SAMRAIVectorReal<DIM,double> &residual );

   //@}

private:

   //@{
   /*!
    * @name Private workhorse functions.
    */

   /*!
    * @brief Red-black Gauss-Seidel error smoothing on a level.
    *
    * Smoothes on the residual equation @f$ Ae=r @f$ on a level.
    *
    * @param error error vector
    * @param residual residual vector
    * @param ln level number
    * @param num_sweeps number of sweeps
    * @param residual_tolerance the maximum residual considered to be
    *        converged
    */
   void smoothErrorByRedBlack(
      SAMRAIVectorReal<DIM,double> &error ,
      const SAMRAIVectorReal<DIM,double> &residual ,
      int ln ,
      int num_sweeps ,
      double residual_tolerance=-1.0 );

   /*!
    * @brief Solve the coarsest level using HYPRE
    */
   void solveCoarsestLevel_HYPRE(
      SAMRAIVectorReal<DIM,double> &error ,
      const SAMRAIVectorReal<DIM,double> &residual ,
      int ln );

   /*!
    * @brief Fix flux per Ewing's coarse-fine boundary treatment.
    *
    * Ewing's coarse-fine boundary treatment can be implemented
    * using a constant refinement into the fine-grid ghost boundary,
    * naively computing the flux using the constant-refined data then
    * fixing up the flux to correct the error.
    *
    * To use this function
    * -# you must use constant refinement to fill the fine level ghost cells
    * -# the flux must first be computed and stored
    *
    * @param patch patch
    * @param soln_data cell-centered solution data
    * @param flux_data side-centered flux data
    * @param diffcoef_data side-centered diffusion coefficient data
    * @param cfb coarse-fine boundary object for the level
    *        in which patch resides
    * @param ratio_to_coarser Refinement ratio to the next coarser level.
    */
   void ewingFixFlux(
      const hier::Patch<DIM> &patch ,
      const pdat::CellData<DIM,double> &soln_data ,
      pdat::SideData<DIM,double> &flux_data ,
      const hier::IntVector<DIM> &ratio_to_coarser ) const;


   /*!
    * @brief AMR-unaware function to compute residual on a single patch,
    * with variable scalar field.
    *
    * @param patch patch
    * @param flux_data side-centered flux data
    * @param soln_data cell-centered solution data
    * @param rhs_data cell-centered rhs data
    * @param residual_data cell-centered residual data
    */
   void computeResidualOnPatch(
      const hier::Patch<DIM> &patch ,
      const pdat::SideData<DIM,double> &flux_data ,
      const pdat::CellData<DIM,double> &soln_data ,
      const pdat::CellData<DIM,double> &rhs_data ,
      pdat::CellData<DIM,double> &residual_data ) const;


   /*!
    * @brief AMR-unaware function to red or black smoothing on a single patch,
    * for variable diffusion coefficient and variable scalar field.
    *
    * @param patch patch
    * @param flux_data side-centered flux data
    * @param rhs_data cell-centered rhs data
    * @param scalar_field_data
    *        cell-centered scalar field data
    * @param soln_data cell-centered solution data
    * @param red_or_black red-black switch.  Set to 'r' or 'b'.
    * @param p_maxres max residual output.  Set to NULL to avoid computing.
    */
   void redOrBlackSmoothingOnPatch(
      const hier::Patch<DIM> &patch ,
      const pdat::SideData<DIM,double> &flux_data ,
      const pdat::CellData<DIM,double> &rhs_data ,
      pdat::CellData<DIM,double> &soln_data ,
      char red_or_black ,
      double *p_maxres=NULL ) const;

   //@}


   //@{ @name For executing, caching and resetting communication schedules.

   /*!
    * @brief Execute a refinement schedule
    * for prolonging cell data.
    *
    * General notes regarding internal objects for communication:
    * We maintain objects to support caching schedules to improve
    * efficiency.  Communication is needed in 5 distinct tasks.
    *   -# Prolongation
    *   -# Restriction
    *   -# Flux coarsening.  Changing the coarse grid flux to the
    *      composite grid flux by coarsening the fine grid flux
    *      at the coarse-fine boundaries.
    *   -# Fill boundary data from other patches in the same level
    *      and physical boundary condition.
    *   -# Fill boundary data from same level, coarser levels
    *      and physical boundary condition.
    *
    * For each task, we maintain a refine or coarsen operator,
    * a cache (array) of communication schedules and a list of keys
    * assocating the schedules to the exact data being transfered.
    * A key is a single integer uniquely identifying the data ids
    * and destination level of the transfer, allowing a simple 1D
    * lookup of the cached schedules.
    *
    * The 5 member functions named @c xeqSchedule... return
    * communication schedules appropriate for five specific tasks.
    * They use a cached schedule if possible or create and cache
    * a new schedule if needed.  These functions and the data
    * they manipulate are as follows:
    * <ol>
    *   <li> xeqScheduleProlongation():
    *        d_prolongation_refine_operator
    *        d_prolongation_refine_schedules
    *   <li> xeqScheduleURestriction():
    *        d_restriction_coarsen_operator,
    *        d_urestriction_coarsen_schedules.
    *   <li> xeqScheduleRRestriction():
    *        d_restriction_coarsen_operator,
    *        d_rrestriction_coarsen_schedules.
    *   <li> xeqScheduleFluxCoarsen():
    *        d_flux_coarsen_operator,
    *        d_flux_coarsen_schedules.
    *   <li> xeqScheduleGhostFill():
    *        d_ghostfill_refine_operator,
    *        d_ghostfill_refine_schedules.
    *   <li> xeqScheduleGhostFillNoCoarse():
    *        d_ghostfill_nocoarse_refine_operator,
    *        d_ghostfill_nocoarse_refine_schedules.
    * </ol>
    *
    * @return refinement schedule for prolongation
    */
   void xeqScheduleProlongation(
      int dst_id,
      int src_id,
      int scr_id,
      int dest_ln);


   /*!
    * @brief Execute schedule for restricting solution to the specified
    * level or reregister an existing one.
    *
    * See general notes for xeqScheduleProlongation().
    *
    * @return coarsening schedule for restriction
    */
   void xeqScheduleURestriction(
      int dst_id,
      int src_id,
      int dest_ln);


   /*!
    * @brief Execute schedule for restricting residual to the specified
    * level or reregister an existing one.
    *
    * See general notes for xeqScheduleProlongation().
    *
    * @return coarsening schedule for restriction
    */
   void xeqScheduleRRestriction(
      int dst_id,
      int src_id,
      int dest_ln);


   /*!
    * @brief Execute schedule for coarsening flux to the specified
    * level or reregister an existing one.
    *
    * See general notes for xeqScheduleProlongation().
    *
    * @return coarsening schedule for setting composite grid flux at
    * coarse-fine boundaries.
    */
   void xeqScheduleFluxCoarsen(
      int dst_id,
      int src_id,
      int dest_ln);


   /*!
    * @brief Execute schedule for filling ghosts on the specified
    * level or reregister an existing one.
    *
    * See general notes for xeqScheduleProlongation().
    *
    * @return refine schedule for filling ghost data from coarser level
    * and physical bc.
    */
   void xeqScheduleGhostFill(
      int dst_id,
      int dest_ln);


   /*!
    * @brief Execute schedule for filling ghosts on the specified
    * level or reregister an existing one.
    * This version does not get data from coarser levels.
    *
    * See general notes for xeqScheduleProlongation().
    *
    * This function is used for the bottom solve level, since it does
    * not access data from any coarser level.  (Ghost data obtained
    * from coarser level must have been placed there before solve begins!)
    *
    * @return refine schedule for filling ghost data from same level
    * and physical bc.
    */
   void xeqScheduleGhostFillNoCoarse(
      int dst_id,
      int dest_ln);

   //@}

   //! @brief Return the patch data index for cell scratch data.
   int registerCellScratch() const;
   //! @brief Return the patch data index for flux scratch data.
   int registerFluxScratch() const;
   //! @brief Return the patch data index for outerflux scratch data.
   int registerOfluxScratch() const;

   //! @brief Free static variables at shutdown time.
   static void freeVariables();


   /*!
    * @brief Object name.
    */
   string d_object_name;



   //@{ @name Hierarchy-dependent objects.

   /*!
    * @brief Reference hierarchy
    *
    * This variable is non-null between the initializeOperatorState()
    * and deallocateOperatorState() calls.  It is not truly needed,
    * because the hierarchy is obtainable through variables in most
    * function argument lists.  We use it to enforce working on one
    * hierarchy at a time.
    */
   tbox::Pointer< hier::PatchHierarchy<DIM> > d_hierarchy;

   /*!
    * @brief Coarsest level for solve.
    */
   int d_ln_min;

   /*!
    * @brief Finest level for solve.
    */
   int d_ln_max;

   /*!
    * @brief Description of coarse-fine boundaries.
    *
    * There is one coarse-fine boundary object for each level.
    * d_coarse_fine_boundary[i] is the description of
    * the coarse-fine boundary between level i and level i-1.
    * The coarse-fine boundary does not exist at the coarsest level,
    * although the hier::CoarseFineBoundary<DIM> object still exists (it
    * should not contain any boxes).
    *
    * This array is initialized in initializeOperatorState() and
    * deallocated in deallocateOperatorState().  When allocated,
    * it is allocated for the index range [0,d_ln_max], though
    * the range [0,d_ln_min-1] is not used.  This is okay because
    *hier::CoarseFineBoundary<DIM> is a light object before
    * it is set for a level.
    */
   tbox::Array< hier::CoarseFineBoundary<DIM> > d_cf_boundary;


   //@}


   //@{
   /*!
    * @name Private state variables for solution process.
    */

   /*!
    * @brief Scalar Poisson equations specifications.
    * @see setPoissonSpecifications().
    */
   PoissonSpecifications d_poisson_spec;

   /*!
    * @brief Smoothing choice.
    * @see setSmoothingChoice.
    */
   string d_smoothing_choice;

   /*!
    * @brief Coarse level solver.
    * @see setCoarsestLevelSolverChoice
    */
   string d_coarse_solver_choice;

   /*!
    * @brief Coarse-fine discretization method.
    * @see setCoarseFineDiscretization().
    */
   string d_cf_discretization;

   /*!
    * @brief Coarse-fine discretization method.
    *
    * The name of the refinement operator used to prolong the
    * coarse grid correction.
    *
    * @see setProlongationMethod()
    */
   string d_prolongation_method;

   /*!
    * @brief Tolerance specified to coarse solver
    * @see setCoarsestLevelSolverTolerance()
    */
   double d_coarse_solver_tolerance;

   /*!
    * @brief Coarse level solver iteration limit.
    * @see setCoarsestLevelSolverMaxIterations()
    */
   int d_coarse_solver_max_iterations;

   /*!
    * @brief Residual tolerance to govern smoothing.
    *
    * When we use one of the internal error smoothing functions
    * and want to terminate the smoothing sweeps at a certain
    * level of residual, this will be set to > 0.  If it is
    * < 0, the smoothing function effectively ignores it.
    *
    * This variable is needed because some coarse-level solver
    * simply runs the smoothing function until convergence.
    * It sets this variable to > 0, calls the smoothing function,
    * then resets it to < 0.
    */
   double d_residual_tolerance_during_smoothing;

   /*!
    * @brief Id of the flux.
    *
    * If set to -1, create and delete storage space on the fly.
    * Else, user has provided space for flux.
    *
    * @see setFluxId
    */
   int d_flux_id;

#ifdef HAVE_HYPRE
   /*!
    * @brief HYPRE coarse-level solver object.
    */
   CellPoissonHypreSolver<DIM> d_hypre_solver;
#endif

   /*!
    * @brief Externally provided physical boundary condition object.
    *
    * see setPhysicalBcCoefObject()
    */
   const RobinBcCoefStrategy<DIM> *d_physical_bc_coef;

   //@}



   //@{ @name Internal context and scratch data

   static tbox::Pointer<pdat::CellVariable<DIM,double> > s_cell_scratch_var;

   static tbox::Pointer<pdat::SideVariable<DIM,double> > s_flux_scratch_var;

   static tbox::Pointer<pdat::OutersideVariable<DIM,double> > s_oflux_scratch_var;

   /*!
    * @brief Default context of internally maintained hierarchy data.
    */
   tbox::Pointer<hier::VariableContext> d_context;

   /*!
    * @brief ID of the solution-like scratch data.
    *
    * Set in constructor and never changed.
    * Corresponds to a pdat::CellVariable<DIM,double> named
    * @c d_object_name+"::cell_scratch".
    * Scratch data is allocated and removed as needed
    * to reduce memory usage.
    */
   int d_cell_scratch_id;

   /*!
    * @brief ID of the side-centered scratch data.
    *
    * Set in constructor and never changed.
    * Corresponds to a pdat::SideVariable<DIM,double> named
    * @c d_object_name+"::flux_scratch".
    *
    * This data is allocated only as needed and deallocated
    * immediately after use.
    */
   int d_flux_scratch_id;

   /*!
    * @brief ID of the outerside-centered scratch data.
    *
    * Set in constructor and never changed.
    * Corresponds to a pdat::OutersideVariable<DIM,double> named
    * @c d_object_name+"::oflux_scratch".
    */
   int d_oflux_scratch_id;

   //@}


   //@{
   /*!
    * @name Various refine and coarsen objects used internally.
    */

   //! @brief Error prolongation (refinement) operator.
   tbox::Pointer<xfer::RefineOperator<DIM> > d_prolongation_refine_operator;
   tbox::Pointer<xfer::RefineAlgorithm<DIM> > d_prolongation_refine_algorithm;
   tbox::Array<tbox::Pointer<xfer::RefineSchedule<DIM> > >
      d_prolongation_refine_schedules;

   //! @brief Solution restriction (coarsening) operator.
   tbox::Pointer<xfer::CoarsenOperator<DIM> > d_urestriction_coarsen_operator;
   tbox::Pointer<xfer::CoarsenAlgorithm<DIM> > d_urestriction_coarsen_algorithm;
   tbox::Array<tbox::Pointer<xfer::CoarsenSchedule<DIM> > >
      d_urestriction_coarsen_schedules;

   //! @brief Residual restriction (coarsening) operator.
   tbox::Pointer<xfer::CoarsenOperator<DIM> > d_rrestriction_coarsen_operator;
   tbox::Pointer<xfer::CoarsenAlgorithm<DIM> > d_rrestriction_coarsen_algorithm;
   tbox::Array<tbox::Pointer<xfer::CoarsenSchedule<DIM> > >
      d_rrestriction_coarsen_schedules;

   //! @brief Coarsen operator for outerflux-to-flux
   tbox::Pointer<xfer::CoarsenOperator<DIM> > d_flux_coarsen_operator;
   tbox::Pointer<xfer::CoarsenAlgorithm<DIM> > d_flux_coarsen_algorithm;
   tbox::Array<tbox::Pointer<xfer::CoarsenSchedule<DIM> > >
      d_flux_coarsen_schedules;

   //! @brief Refine operator for cell-like data from coarser level.
   tbox::Pointer<xfer::RefineOperator<DIM> > d_ghostfill_refine_operator;
   tbox::Pointer<xfer::RefineAlgorithm<DIM> > d_ghostfill_refine_algorithm;
   tbox::Array<tbox::Pointer<xfer::RefineSchedule<DIM> > >
      d_ghostfill_refine_schedules;

   //! @brief Refine operator for cell-like data from same level.
   tbox::Pointer<xfer::RefineOperator<DIM> > d_ghostfill_nocoarse_refine_operator;
   tbox::Pointer<xfer::RefineAlgorithm<DIM> > d_ghostfill_nocoarse_refine_algorithm;
   tbox::Array<tbox::Pointer<xfer::RefineSchedule<DIM> > >
      d_ghostfill_nocoarse_refine_schedules;

   //@}


   /*!
    * @brief Utility object employed in setting ghost cells and providing
    * xfer::RefinePatchStrategy<DIM> implementation.
    *
    * Since this class deals only in scalar variables having
    * Robin boundary conditions, we take advantage of the corresponding
    * implementation in CartesianRobinBcHelper<DIM>.  Whenever
    * we need an implementation of xfer::RefinePatchStrategy<DIM>,
    * this object is used.  Note that in the code, before we
    * use this object to set ghost cell values, directly or
    * indirectly by calling xfer::RefineSchedule::fillData(),
    * we must tell d_bc_helper the patch data index we want
    * to set and whether we are setting data with homogeneous
    * boundary condition.
    */
   CartesianRobinBcHelper<DIM> d_bc_helper;




   //@{
   /*!
    * @name Non-essential objects used in outputs and debugging.
    */

   /*!
    * @brief Logging flag.
    */
   bool d_enable_logging;

   /*!
    * @brief Preconditioner using this object.
    *
    * See setPreconditioner().
    */
   const FACPreconditioner<DIM> *d_preconditioner;

   /*!
    * @brief Hierarchy cell operator used in debugging.
    */
   tbox::Pointer<math::HierarchyCellDataOpsReal<DIM,double> > d_hopscell;

   /*!
    * @brief Hierarchy side operator used in debugging.
    */
   tbox::Pointer<math::HierarchySideDataOpsReal<DIM,double> > d_hopsside;

   /*!
    * @brief Timers for performance measurement.
    */
   tbox::Pointer<tbox::Timer> t_restrict_solution;
   tbox::Pointer<tbox::Timer> t_restrict_residual;
   tbox::Pointer<tbox::Timer> t_prolong;
   tbox::Pointer<tbox::Timer> t_smooth_error;
   tbox::Pointer<tbox::Timer> t_solve_coarsest;
   tbox::Pointer<tbox::Timer> t_compute_composite_residual;
   tbox::Pointer<tbox::Timer> t_compute_residual_norm;

};

}
}

#ifndef DEBUG_NO_INLINE
#include "CellPoissonFACOps.I"
#endif

#endif // included_solv_CellPoissonFACOps

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "CellPoissonFACOps.C"
#endif
