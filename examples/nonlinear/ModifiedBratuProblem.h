//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/examples/nonlinear/ModifiedBratuProblem.h $
// Package:     SAMRAI application
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2132 $
// Modified:    $LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description: Class containing numerical routines for modified Bratu problem
//

#ifndef included_ModifiedBratuProblem
#define included_ModifiedBratuProblem


#include "SAMRAI_config.h"

#if !defined(HAVE_PETSC) || !defined(HAVE_SUNDIALS) || !defined(HAVE_HYPRE)

/*
*************************************************************************
* If the library is not compiled with PETSC -and- KINSOL, print an error.
* If we're running autotests, skip the error and compile an empty
* class.
*************************************************************************
*/
#if (TESTING != 1)
#error "This example requires SAMRAI be compiled with KINSOL, PETSC, and HYPRE."
#endif

#else

#include "tbox/AbstractStream.h"
#include "tbox/Array.h"
#include "Box.h"
#include "BoxArray.h"
#include "BoxList.h"
#include "CartesianGridGeometry.h"
#include "CellData.h"
#include "CellVariable.h"
#include "CoarsenAlgorithm.h"
#include "CoarsenOperator.h"
#include "CoarsenPatchStrategy.h"
#include "CoarsenSchedule.h"
#include "ComponentSelector.h"
#include "tbox/Database.h"
#include "FaceVariable.h"
#include "ImplicitEquationStrategy.h"
#include "IntVector.h"
#include "OutersideVariable.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "tbox/Pointer.h"
#include "RefineAlgorithm.h"
#include "RefineOperator.h"
#include "RefinePatchStrategy.h"
#include "RefineSchedule.h"
#include "SAMRAIVectorReal.h"
#include "tbox/Serializable.h"
#include "SideVariable.h"
#include "StandardTagAndInitStrategy.h"
#include <string>
using namespace std;
#define included_String
#include "VariableContext.h"
#include "CartesianVizamraiDataWriter.h"
#include "VisItDataWriter.h"

#include "CellPoissonFACSolver.h"

#include "KINSOLAbstractFunctions.h"
#include "SundialsAbstractVector.h"
#include "PETScAbstractVectorReal.h"
#include "SNESAbstractFunctions.h"

// Dimension macros.

using namespace SAMRAI;
using namespace xfer;

/**
 * Class ModifiedBratuProblem class provides operations needed to solve
 *
 *     du/dt = div( D(x,t)*grad(u) ) + lambda * exp(u) + f(u,x,t)
 *
 * using implicit time integration and either KINSOL or PETSc to solve the 
 * nonlinear system at each step.  Specifically, it provides operations 
 * needed by the algs::ImplicitIntegrator<NDIM> class as well as those defined by the 
 * interfaces to KINSOL and PETSc; i.e., KINSOLAbstractFunctions and 
 * SNESAbstractFunctions respectively.
 *
 * This example is implemented only for 2D, 2:1 refinement ratios only.
 */

class ModifiedBratuProblem : 
   public algs::ImplicitEquationStrategy<NDIM>, 
   public mesh::StandardTagAndInitStrategy<NDIM>,
   public solv::SNESAbstractFunctions,
   public solv::KINSOLAbstractFunctions,
   public RefinePatchStrategy<NDIM>, 
   public CoarsenPatchStrategy<NDIM>,
   public tbox::Serializable
{
public: 

   /**
    * Constructor for ModifiedBratuProblem class creates problem variables
    * to represent the solution and other quantities on the patch hierarchy.
    * It initializes data members to default values and sets others based
    * on input and/or restart values.  The constructor also sets up algorithms
    * for communicating data between patches on the hierarchy.
    */
   ModifiedBratuProblem(
      const string& object_name,
      tbox::Pointer<tbox::Database> input_db,
      tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geometry,
      tbox::Pointer<appu::CartesianVizamraiDataWriter<NDIM> > viz_data_writer = NULL,
      tbox::Pointer<appu::VisItDataWriter<NDIM> > visit_data_writer = NULL);

   /**
    * Destructor for ModifiedBratuProblem class does nothing.
    */
   ~ModifiedBratuProblem();


   //@{
   /*!
     @name Implicit integrator interfaces
   */

   /**
    * Set vector weights on the hierarchy.
    */
   void setVectorWeights( tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy );

   /**
    * Set the nonlinear solution vector so that the new solution data is
    * solved for when the nonlinear solver advances the solution.
    *
    * Function overloaded from algs::ImplicitEquationStrategy<NDIM>.
    */
   void setupSolutionVector(tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> > solution);

   /**
    * Return time increment for advancing the solution at the first timestep.
    *
    * Function overloaded from algs::ImplicitEquationStrategy<NDIM>.
    */
   double getInitialDt();

   /**
    * Return the next time increment through which to advance the solution.
    * The good_solution is the value returned by a call to checkNewSolution(),
    * which determines whether the computed solution is acceptable or not.
    * The integer solver_retcode is the return code generated by the
    * nonlinear solver.   This value must be interpreted in a manner
    * consistant with the solver in use.
    *
    * Function overloaded from algs::ImplicitEquationStrategy<NDIM>.
    */
   double getNextDt(const bool good_solution,
                    const int solver_retcode);

   /**
    * Set the initial guess for the time advanced solution at the start
    * of the nonlinear iteration.  The boolean argument first_step
    * indicates whether we are at the first step on the current hierarchy
    * configuration.  This is true when the hierarchy is constructed
    * initially and after regridding.  In these cases, setting the initial
    * iterate using extrapolation, for example, may not be possible.
    *
    * Function overloaded from algs::ImplicitEquationStrategy<NDIM>.
    */
   void setInitialGuess(const bool first_step,
                        const double current_time,
                        const double current_dt,
                        const double old_dt);

   /**
    * Check the computed solution and return true if it is acceptable;
    * otherwise return false.  The integer solver_retcode is the return
    * code generated by the nonlinear solver.  This value must be
    * interpreted in a manner consistent with the solver in use.
    *
    * Function overloaded from algs::ImplicitEquationStrategy<NDIM>.
    */
   bool checkNewSolution(const int solver_retcode);

   /**
    * Update solution storage and dependent quantities after computing an
    * acceptable time advanced solution.   The new_time value is the new
    * solution time.
    *
    * Function overloaded from algs::ImplicitEquationStrategy<NDIM>.
    */
   void updateSolution(const double new_time);

   //@}


   //@{
   /*!
    * @name Functions overloaded from mesh::StandardTagAndInitStrategy<NDIM>.
    */
   
   void initializeLevelData(const tbox::Pointer<hier::BasePatchHierarchy<NDIM> > hierarchy,
                            const int level_number,
                            const double time,
                            const bool can_be_refined,
                            const bool initial_time,
                            const tbox::Pointer<hier::BasePatchLevel<NDIM> > old_level = NULL,
                            const bool allocate_data = true );

   void resetHierarchyConfiguration(const tbox::Pointer<hier::BasePatchHierarchy<NDIM> > hierarchy,
                                    const int coarsest_level,
                                    const int finest_level);
   //@}

   //@{
   /*!
    * @name Interface functions overloaded from KINSOLAbstractFunctions.
    */

   void evaluateNonlinearFunction(solv::SundialsAbstractVector* soln,
                                  solv::SundialsAbstractVector* fval);

#if 0
   // SGS need to fix; make note for users
   int precondSetup(solv::SundialsAbstractVector* soln,
                    solv::SundialsAbstractVector* soln_scale,
                    solv::SundialsAbstractVector* fval,
                    solv::SundialsAbstractVector* fval_scale,
                    solv::SundialsAbstractVector* vtemp1,
                    solv::SundialsAbstractVector* vtemp2,
                    double mach_roundoff,
                    int& num_feval);
   int precondSolve(solv::SundialsAbstractVector* soln,
                    solv::SundialsAbstractVector* soln_scale,
                    solv::SundialsAbstractVector* fval,
                    solv::SundialsAbstractVector* fval_scale,
                    solv::SundialsAbstractVector* rhs,
                    solv::SundialsAbstractVector* vtemp,
                    double mach_roundoff,
                    int& num_feval);
#endif
   
   int precondSetup(solv::SundialsAbstractVector* soln,
		    solv::SundialsAbstractVector* soln_scale,
		    solv::SundialsAbstractVector* fval,
		    solv::SundialsAbstractVector* fval_scale,
		    solv::SundialsAbstractVector* vtemp1,
		    solv::SundialsAbstractVector* vtemp2,
		    int& num_feval);

   int precondSolve(solv::SundialsAbstractVector* soln,
		    solv::SundialsAbstractVector* soln_scale,
		    solv::SundialsAbstractVector* fval,
		    solv::SundialsAbstractVector* fval_scale,
		    solv::SundialsAbstractVector* rhs,
		    solv::SundialsAbstractVector* vtemp,
		    int& num_feval);

   int jacobianTimesVector(solv::SundialsAbstractVector* vector,
                           solv::SundialsAbstractVector* product,
                           const bool soln_changed,
                           solv::SundialsAbstractVector* soln);


   //@}

   //@{
   /*!
    * @name Interface functions overloaded from SNESAbstractFunctions.
    */
   
   int evaluateNonlinearFunction(Vec xcur,
                                 Vec fcur);

   int evaluateJacobian(Vec x);
 
   int jacobianTimesVector(Vec xin,
                           Vec xout);

   int setupPreconditioner(Vec x);

   int applyPreconditioner(Vec r,
                           Vec z);
   //@}

   /*!
    * @brief Set solution ghost cell values along physical boundaries.
    *
    * Function is overloaded from RefinePatchStrategy.
    */

   void setPhysicalBoundaryConditions(hier::Patch<NDIM>& patch,
                                  const double time,
                                  const hier::IntVector<NDIM>& ghost_width_to_fill);

   //@{
   /*!
    * @name Empty functions for applying user-defined data refine operations
    */

   /*
    * These are overloaded from RefinePatchStrategy.
    * There are no such user-defined operations here.
    */

   void preprocessRefine(hier::Patch<NDIM>& fine,
                         const hier::Patch<NDIM>& coarse,
                         const hier::Box<NDIM>& fine_box,
                         const hier::IntVector<NDIM>& ratio)
   {
      (void) fine;
      (void) coarse;
      (void) fine_box;
      (void) ratio;
   }

   void postprocessRefine(hier::Patch<NDIM>& fine,
                          const hier::Patch<NDIM>& coarse,
                          const hier::Box<NDIM>& fine_box,
                          const hier::IntVector<NDIM>& ratio)
   {
      (void) fine;
      (void) coarse;
      (void) fine_box;
      (void) ratio;
   }

   hier::IntVector<NDIM> getRefineOpStencilWidth() const 
   {
       return(hier::IntVector<NDIM>(0));
   }

   //@}

   //@{
   /*!
    * @name Empty functions for applying user-defined data coarsen operations
    */
   /*
    * These are overloaded from CoarsenPatchStrategy.
    * There are no such user-defined operations here.
    */

   void preprocessCoarsen(hier::Patch<NDIM>& coarse,
                          const hier::Patch<NDIM>& fine,
                          const hier::Box<NDIM>& coarse_box,
                          const hier::IntVector<NDIM>& ratio)
   {
      (void) coarse;
      (void) fine;
      (void) coarse_box;
      (void) ratio;
   }

   void postprocessCoarsen(hier::Patch<NDIM>& coarse,
                           const hier::Patch<NDIM>& fine,
                           const hier::Box<NDIM>& coarse_box,
                           const hier::IntVector<NDIM>& ratio)
   {
      (void) coarse;
      (void) fine;
      (void) coarse_box;
      (void) ratio;
   }

   hier::IntVector<NDIM> getCoarsenOpStencilWidth() const 
   {
       return(hier::IntVector<NDIM>(0));
   }

   //@}

   /**
    * Write data members to given data base for restart.
    *
    * Overloaded from tbox::Serializable.
    */
   void putToDatabase(tbox::Pointer<tbox::Database> db);

   /**
    * Write class data to given output stream.
    */
   void printClassData(ostream& os) const;

private:
   /*
    * Functions to read data from input and restart databases. If the boolean
    * flag is true, all data members are read from restart.  They can
    * later be overwritten from values in the input file.  When the flag
    * is false, all data values are set from thos given in input.
    *
    * An assertion results if the database pointer is null.
    */
   void getFromInput(tbox::Pointer<tbox::Database> db,
                     bool is_from_restart);

   /* 
    * Functions for fixing up flux computations along coarse/fine interfaces
    * when ghost cells are filled with CONSTANT_REFINE refinement operators.
    */

   void getLevelEdges(hier::BoxList<NDIM>& boxes,
                      tbox::Pointer<hier::Patch<NDIM> > patch,
                      tbox::Pointer<hier::PatchLevel<NDIM> > level,
                      const int dim,
                      const int face);

   void correctLevelFlux(tbox::Pointer<hier::PatchLevel<NDIM> > level);

   void correctPatchFlux(tbox::Pointer<hier::PatchLevel<NDIM> > level,
                         tbox::Pointer<hier::Patch<NDIM> > patch,
                         tbox::Pointer< pdat::CellData<NDIM,double> > u);

   //@{
   /*!
    * @name Numerical routines specific to modified Bratu problem
    */
   /*
    * These are needed by the nonlinear solvers.
    * They are called by the interface
    * routines after the vectors and other data has been appropriately
    * unwrapped so that these routines are solver-independent.
    */

   void evaluateBratuFunction(tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> > x,
                              tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> > f);

   /*!
     @brief Compute A(x)*x.

     The A(x) used is the one computed in evaluateBratuJacobian()
     and stored at d_jacobian_a_id and d_jacobian_b_id.
   */
   int jacobianTimesVector(tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> > vector,
                           tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> > product);

   void setupBratuPreconditioner(tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> > x);

   int applyBratuPreconditioner(tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> > r,
                                tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> > z);

   /*!
     @brief Recompute the jacobian A(x).

     The diagonal of A(x) is placed at d_jacobian_a_id.
     The off-diagonals are not computed here, because they are
     independent of x, and it is easier to not explicitly compute them.
   */
   void evaluateBratuJacobian(tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> > x);

   //@}

   /*
    * The object name is used as a handle to databases stored in
    * restart files and for error reporting purposes.
    */
   string d_object_name;

   /*
    * We cache a pointer to the grid geometry object to set up initial
    * data and set physical boundary conditions.
    */
   tbox::Pointer<geom::CartesianGridGeometry<NDIM> > d_grid_geometry;

   /*
    * Parameters read from input.
    */

   double d_lambda;      // factor multiplying exponential term
   double d_input_dt;    // time increment

   /*
    *hier::Variable<NDIM> data management.  
    *
    * Contexts are labels to describe the way variables are used.
    */
   tbox::Pointer<hier::VariableContext> d_current;
   tbox::Pointer<hier::VariableContext> d_new;
   tbox::Pointer<hier::VariableContext> d_scratch;

   /*
    * Variables for the discrete problem; see comments above class constructor.
    */
   tbox::Pointer< pdat::CellVariable<NDIM,double> > d_solution;
   tbox::Pointer< pdat::CellVariable<NDIM,double> > d_source_term;
   tbox::Pointer< pdat::CellVariable<NDIM,double> > d_exponential_term;
   tbox::Pointer< pdat::SideVariable<NDIM,double> > d_diffusion_coef;
   tbox::Pointer< pdat::SideVariable<NDIM,double> > d_flux;
   tbox::Pointer< pdat::OutersideVariable<NDIM,double> > d_coarse_fine_flux;

   /*
    * For storing Jacobian A(x) stuff and computing Jacobian-vector
    * multiply A(x)*v.
    */
   tbox::Pointer< pdat::CellVariable<NDIM,double> > d_jacobian_a;
   tbox::Pointer< pdat::FaceVariable<NDIM,double> > d_jacobian_b;
   int d_jacobian_a_id;
   int d_jacobian_b_id;
   hier::ComponentSelector d_jacobian_data;

   /*
    * For storing Jacobian A(x) stuff in setting up and applying
    * the preconditioner A(x)*z=r.
    */
   tbox::Pointer< pdat::CellVariable<NDIM,double> > d_precond_a;
   tbox::Pointer< pdat::FaceVariable<NDIM,double> > d_precond_b;
   int d_precond_a_id;
   int d_precond_b_id;
   hier::ComponentSelector d_precond_data;

   int d_soln_scratch_id;
   int d_flux_id;
   int d_coarse_fine_flux_id;
   int d_function_id;

   hier::ComponentSelector d_problem_data;
   hier::ComponentSelector d_new_patch_problem_data;


   hier::IntVector<NDIM> d_nghosts;

   /*
    * The nonlinear solution process requires a solution vector; we cache
    * a pointer to it here.  A variable is used to define weights for the
    * solution vector entries on a composite grid.
    */
   tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> > d_solution_vector;
   tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> > d_current_soln_vector;

   tbox::Pointer< pdat::CellVariable<NDIM,double> > d_weight;

   int d_weight_id;

   /*
    * Communication algorithms and schedules used for filling ghost cells
    * and moving data between levels.
    * Schedules stored in arrays are indexed by the destination
    * level number in the transfer.  They are cached to save the cost
    * of generating them multiple times for the same hierarchy
    * configuration.
    */
   RefineAlgorithm<NDIM>  d_fill_new_level;
   RefineAlgorithm<NDIM> d_soln_fill;
   tbox::Array<tbox::Pointer<RefineSchedule<NDIM> > > d_soln_fill_schedule;
   CoarsenAlgorithm<NDIM> d_flux_coarsen;
   tbox::Array<tbox::Pointer<CoarsenSchedule<NDIM> > > d_flux_coarsen_schedule;
   CoarsenAlgorithm<NDIM> d_soln_coarsen;
   tbox::Array<tbox::Pointer<CoarsenSchedule<NDIM> > > d_soln_coarsen_schedule;
   CoarsenAlgorithm<NDIM> d_scratch_soln_coarsen;
   tbox::Array<tbox::Pointer<CoarsenSchedule<NDIM> > > d_scratch_soln_coarsen_schedule;

   tbox::Pointer<RefineOperator<NDIM> >  d_soln_refine_op;
   tbox::Pointer<CoarsenOperator<NDIM> > d_soln_coarsen_op;

   /*
    * Current solution time and time increment used in the solution process.
    * New time is current time + current dt.
    */
   double d_current_time;
   double d_new_time;
   double d_current_dt;

   /*
    * Preconditioner and parameters used for Jacobian system.
    *
    * The FAC solver manages the composite grid solution procedure.
    * The Poisson level strategy solves the problem on each level
    * in the hierarchy.  
    */
   bool d_use_old_solver;
   tbox::Pointer<solv::CellPoissonFACSolver<NDIM> > d_FAC_solver;

   int    d_max_precond_its;
   double d_precond_tol;



   static tbox::Pointer<tbox::Timer> s_copy_timer;
   static tbox::Pointer<tbox::Timer> s_pc_timer;

};

#endif
#endif
