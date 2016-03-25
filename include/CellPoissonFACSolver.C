#ifndef included_solv_CellPoissonFACSolver_C
#define included_solv_CellPoissonFACSolver_C

/*
 * File:         $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/solvers/poisson/CellPoissonFACSolver.C $
 * Package:      SAMRAI solvers
 * Copyright:    (c) 1997-2007 Lawrence Livermore National Security, LLC
 * Revision:     $LastChangedRevision: 1704 $
 * Modified:     $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
 * Description: High-level solver (wrapper) for scalar poisson equation.
 */

#include "CellVariable.h"
#include "CellPoissonFACSolver.h"
#include "tbox/PIO.h"
#include "tbox/Utilities.h"

#include IOMANIP_HEADER_FILE


#ifdef DEBUG_NO_INLINE
#include "CellPoissonFACSolver.I"
#endif

namespace SAMRAI {
    namespace solv {

/*
*************************************************************************
*                                                                       *
* Constructor sets uninitialized solver state.                          *
* Set default iteration and convergence parameters.                     *
*                                                                       *
* By default settings:                                                  *
*   - Poisson equation specified has D=1, C=0.                          *
*   - State is uninitialized                                            *
*   - Logging is disabled                                               *
*   - Context for internal data is set based on object name.            *
*                                                                       *
*************************************************************************
*/

template<int DIM>  CellPoissonFACSolver<DIM>::CellPoissonFACSolver (
   const std::string &object_name,
   tbox::Pointer<tbox::Database> database )
   :
   d_object_name(object_name),
   d_poisson_spec(object_name+"::poisson_spec"),
   d_fac_ops(object_name+"::fac_ops"),
   d_fac_precond(object_name+"::fac_precond",d_fac_ops),
   d_bc_object(NULL),
   d_simple_bc(object_name+"::bc"),
   d_hierarchy(NULL),
   d_ln_min(-1),
   d_ln_max(-1),
   d_context(hier::VariableDatabase<DIM>::getDatabase()
             ->getContext(object_name+"::CONTEXT")) ,
   d_uv(NULL),
   d_fv(NULL),
   d_weight_id(hier::VariableDatabase<DIM>::getDatabase()->
      registerVariableAndContext(
         new pdat::CellVariable<DIM,double>(object_name+"::weight"),
         d_context )),
   d_solver_is_initialized(false),
   d_enable_logging(false)
{

   setMaxCycles(10);
   setResidualTolerance(1e-6);
   setPresmoothingSweeps(1);
   setPostsmoothingSweeps(1);
   setCoarseFineDiscretization("Ewing");
#ifdef HAVE_HYPRE
   setCoarsestLevelSolverChoice("hypre");
   setCoarsestLevelSolverTolerance(1e-10);
   setCoarsestLevelSolverMaxIterations(20);
   setUseSMG(true);
#else
   setCoarsestLevelSolverChoice("redblack");
   setCoarsestLevelSolverTolerance(1e-8);
   setCoarsestLevelSolverMaxIterations(500);
#endif

   /*
    * The default RobinBcCoefStrategy<DIM> used,
    * SimpleCellRobinBcCoefs<DIM> only works with constant refine
    * for prolongation.  So we use constant refinement
    * for prolongation by default.
    */
   setProlongationMethod("CONSTANT_REFINE");

   /*
    * The FAC operator optionally uses the preconditioner
    * to get data for logging.
    */
   d_fac_ops.setPreconditioner((const FACPreconditioner<DIM>*)(&d_fac_precond));


   if ( database ) {
      getFromInput(database);
   }

   return;
}

/*
*************************************************************************
*                                                                       *
* Destructor for CellPoissonFACSolver<DIM>.                            *
* Deallocate internal data.                                             *
*                                                                       *
*************************************************************************
*/

template<int DIM>  CellPoissonFACSolver<DIM>::~CellPoissonFACSolver()
{
   deallocateSolverState();
   return;
}


/*
********************************************************************
* Set state from database                                          *
*                                                                  *
* Do not allow FAC preconditioner and Poisson FAC operators to be  *
* set from database, as that may cause them to be inconsistent     *
* with this object if user does not coordinate the inputs          *
* correctly.  This is also why we don't allow direct access to     *
* those objects.  The responsibility for maintaining consistency   *
* lies in the public functions to set parameters, so use them      *
* instead of setting the parameters directly in this function.     *
********************************************************************
*/

template<int DIM> void CellPoissonFACSolver<DIM>::getFromInput(
   tbox::Pointer<tbox::Database> database )
{
   if ( database ) {
      if ( database->isBool("enable_logging") ) {
         bool logging = database->getBool("enable_logging");
         enableLogging(logging);
      }
      if ( database->isInteger("max_cycles") ) {
         int max_cycles = database->getInteger("max_cycles");
         setMaxCycles(max_cycles);
      }
      if ( database->isDouble("residual_tol") ) {
         double residual_tol = database->getDouble("residual_tol");
         setResidualTolerance(residual_tol);
      }
      if ( database->isInteger("num_pre_sweeps") ) {
         int num_pre_sweeps = database->getInteger("num_pre_sweeps");
         setPresmoothingSweeps(num_pre_sweeps);
      }
      if ( database->isInteger("num_post_sweeps") ) {
         int num_post_sweeps = database->getInteger("num_post_sweeps");
         setPostsmoothingSweeps(num_post_sweeps);
      }
      if ( database->isString("coarse_fine_discretization") ) {
         std::string s = database->getString("coarse_fine_discretization");
         setCoarseFineDiscretization(s);
      }
      if ( database->isString("prolongation_method") ) {
         std::string s = database->getString("prolongation_method");
         setProlongationMethod(s);
      }
      if ( database->isString("coarse_solver_choice") ) {
         std::string s = database->getString("coarse_solver_choice");
         setCoarsestLevelSolverChoice(s);
      }
      if ( database->isDouble("coarse_solver_tolerance") ) {
         double tol = database->getDouble("coarse_solver_tolerance");
         setCoarsestLevelSolverTolerance(tol);
      }
      if ( database->isInteger("coarse_solver_max_iterations") ) {
         int itr = database->getInteger("coarse_solver_max_iterations");
         setCoarsestLevelSolverMaxIterations(itr);
      }
#ifdef HAVE_HYPRE
      if ( database->isBool("use_smg") ) {
         bool smg = database->getBool("use_smg");
         setUseSMG(smg);
      }
#endif
   }
   return;
}





/*
*************************************************************************
*                                                                       *
* Prepare internal data for solve.                                      *
* Allocate scratch data.  Create vectors for u and f                    *
* required by the FACPreconditioner<DIM> interface.                    *
* Set up internal boundary condition object.                            *
* Share data to coordinate with FAC preconditioner and                  *
* Poisson FAC operator.                                                 *
*                                                                       *
*************************************************************************
*/

template<int DIM> void CellPoissonFACSolver<DIM>::initializeSolverState(
   const int solution ,
   const int rhs ,
   tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
   const int coarse_level,
   const int fine_level )
{
   if ( d_bc_object == NULL ) {
      TBOX_ERROR(d_object_name << ": No BC coefficient strategy object!\n"
                 << "Use either setBoundaries or setPhysicalBcCoefObject\n"
                 << "to specify the boundary conidition.\n");
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   if ( solution < 0 || rhs < 0 ) {
      TBOX_ERROR(d_object_name << ": Bad patch data id.\n");
   }
#endif

#ifdef DEBUG_CHECK_ASSERTIONS
   if ( !hierarchy ) {
      TBOX_ERROR(d_object_name << ": NULL hierarchy pointer not allowed\n"
                 << "in inititialization.");
   }
#endif
   d_hierarchy = hierarchy;

   d_ln_min = coarse_level;
   d_ln_max = fine_level;
   if ( d_ln_min == -1 ) {
      d_ln_min = 0;
   }
   if ( d_ln_max == -1 ) {
      d_ln_max = d_hierarchy->getFinestLevelNumber();
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   if ( d_ln_min < 0 || d_ln_max < 0 || d_ln_min > d_ln_max ) {
      TBOX_ERROR(d_object_name << ": Bad range of levels in\n"
                 << "inititialization.\n");
   }
#endif

   int ln;
   for ( ln=d_ln_min; ln<=d_ln_max; ++ln ) {
      d_hierarchy->getPatchLevel(ln)->allocatePatchData(d_weight_id);
   }

   d_fac_ops.computeVectorWeights( d_hierarchy,
                                   d_weight_id,
                                   d_ln_min,
                                   d_ln_max );

   if ( d_bc_object == &d_simple_bc ) {
      d_simple_bc.setHierarchy(d_hierarchy,
                               d_ln_min,
                               d_ln_max);
      if ( d_poisson_spec.dIsConstant() ) {
         d_simple_bc.setDiffusionCoefConstant( d_poisson_spec.getDConstant() );
      }
      else {
         d_simple_bc.setDiffusionCoefId( d_poisson_spec.getDPatchDataId() );
      }
   }

   d_fac_ops.setPoissonSpecifications(d_poisson_spec);

   createVectorWrappers(solution,rhs);

   d_fac_precond.initializeSolverState(*d_uv, *d_fv);

   d_solver_is_initialized = true;

   return;
}




template<int DIM> void CellPoissonFACSolver<DIM>::deallocateSolverState()
{
   if ( d_hierarchy ) {

      d_fac_precond.deallocateSolverState();

      /*
       * Delete internally managed data.
       */
      int ln;
      for ( ln=d_ln_min; ln<=d_ln_max; ++ln ) {
         d_hierarchy->getPatchLevel(ln)->deallocatePatchData(d_weight_id);
      }

      d_hierarchy.setNull();
      d_ln_min = -1;
      d_ln_max = -1;
      d_solver_is_initialized = false;

      destroyVectorWrappers();

   }
   return;
}




/*
*************************************************************************
* Enable logging and propagate logging flag to major components.        *
*************************************************************************
*/

template<int DIM> void CellPoissonFACSolver<DIM>::enableLogging( bool logging )
{
   d_enable_logging = logging;
   d_fac_precond.enableLogging(d_enable_logging);
   d_fac_ops.enableLogging(d_enable_logging);
   return;
}





template<int DIM> void CellPoissonFACSolver<DIM>::setBoundaries(const std::string& boundary_type,
                                               const int fluxes,
                                               const int flags,
                                               int* bdry_types)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if ( d_bc_object != NULL && d_bc_object != &d_simple_bc ) {
      TBOX_ERROR(d_object_name << ": Bad attempt to set boundary condition\n"
                 << "by using default bc object after it has been overriden.\n");
   }
#endif
   d_simple_bc.setBoundaries(boundary_type,
                             fluxes,
                             flags,
                             bdry_types);
   d_bc_object = &d_simple_bc;
   d_fac_ops.setPhysicalBcCoefObject(d_bc_object);
   return;
}





template<int DIM> void CellPoissonFACSolver<DIM>::setBcObject(
   const RobinBcCoefStrategy<DIM> *bc_object )
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if ( !bc_object ) {
      TBOX_ERROR(d_object_name << ": NULL pointer for boundary condition\n"
                 << "object.\n");
   }
#endif
   d_bc_object = bc_object;
   d_fac_ops.setPhysicalBcCoefObject(d_bc_object);
   return;
}



/*
*************************************************************************
*                                                                       *
* Solve the linear system and report whether iteration converged.       *
*                                                                       *
* This version is for an initialized solver state.                      *
* Before solving, set the final piece of the boundary condition,        *
* which is not known until now, and initialize some internal            *
* solver quantities.                                                    *
*                                                                       *
*************************************************************************
*/

template<int DIM> bool CellPoissonFACSolver<DIM>::solveSystem(
   const int u,
   const int f)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if ( !d_solver_is_initialized ) {
      TBOX_ERROR(d_object_name << ".solveSystem(int,int): uninitialized\n"
                 << "solver state.  You must call initializeSolverState()\n"
                 << "before using this function.  Or you can use\n"
                 << "solveSystem(int,int,...) to initialize the solver,\n"
                 << "solve and deallocate the solver.\n");
   }
   if ( u < 0 || f < 0 ) {
      TBOX_ERROR(d_object_name << ": Bad patch data id.\n");
   }
#endif
   if ( d_bc_object == &d_simple_bc ) {
      /*
       * Knowing that we are using the SimpelCellRobinBcCoefsX
       * implementation of RobinBcCoefStrategy, we must save
       * the ghost data in u before solving.
       * The solver overwrites it, but SimpleCellRobinBcCoefs<DIM>
       * needs to get to access it repeatedly.
       */
      d_simple_bc.cacheDirichletData(u);
   }

   createVectorWrappers(u,f);
   bool solver_rval;
   solver_rval = d_fac_precond.solveSystem( *d_uv, *d_fv );

   if ( d_bc_object == &d_simple_bc ) {
      /*
       * Restore the Dirichlet cell data that were overwritten by the
       * solve process.  We do this to be backward compatible with the
       * user code.
       */
      d_simple_bc.restoreDirichletData(u);
   }

   return solver_rval;
}



/*
*************************************************************************
*                                                                       *
* Solve the linear system and report whether iteration converged.       *
*                                                                       *
* This version is for an uninitialized solver state.                    *
* 1. Initialize the (currently uninitialized) solver state.             *
* 2. Solve.                                                             *
* 3. Deallocate the solver state.                                       *
*                                                                       *
*************************************************************************
*/

template<int DIM> bool CellPoissonFACSolver<DIM>::solveSystem(
   const int u,
   const int f,
   tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
   int coarse_ln,
   int fine_ln)
{
   if ( d_enable_logging ) {
      tbox::plog << "CellPoissonFACSolver<DIM>::solveSystem (" << d_object_name
           << ")\n";
      d_poisson_spec.printClassData(tbox::plog);
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   if ( d_solver_is_initialized ) {
      TBOX_ERROR(d_object_name << ".solveSystem(int,int,...): initialized\n"
                 << "solver state.  This function can only used when the\n"
                 << "solver state is uninitialized.  You should deallocate\n"
                 << "the solver state or use solveSystem(int,int).\n");
   }
   if ( !hierarchy ) {
      TBOX_ERROR(d_object_name << ".solveSystem(): Null hierarchy\n"
                 << "specified.\n");
   }
#endif
   initializeSolverState( u, f, hierarchy, coarse_ln, fine_ln );

   bool solver_rval;
   solver_rval = solveSystem( u, f );

   deallocateSolverState();

   return solver_rval;
}




template<int DIM> void CellPoissonFACSolver<DIM>::createVectorWrappers( int u, int f ) {

   hier::VariableDatabase<DIM> &vdb(*hier::VariableDatabase<DIM>::getDatabase());
   tbox::Pointer< hier::Variable<DIM> > variable;

   if ( !d_uv || d_uv->getComponentDescriptorIndex(0) != u ) {
      d_uv.setNull();
      d_uv = new SAMRAIVectorReal<DIM,double>(d_object_name+"::uv",
                                                d_hierarchy,
                                                d_ln_min,
                                                d_ln_max);
      vdb.mapIndexToVariable(u, variable);
#ifdef DEBUG_CHECK_ASSERTIONS
      if ( !variable ) {
         TBOX_ERROR(d_object_name << ": No variable for patch data index "
                    << u << "\n");
      }
      tbox::Pointer<pdat::CellVariable<DIM,double> > cell_variable = variable;
      if ( !cell_variable ) {
         TBOX_ERROR(d_object_name << ": hier::Patch data index " << u
                    << " is not a cell-double variable.\n");
      }
#endif
      d_uv->addComponent( variable, u, d_weight_id );
   }

   if ( !d_fv || d_fv->getComponentDescriptorIndex(0) != f ) {
      d_fv.setNull();
      d_fv = new SAMRAIVectorReal<DIM,double>(d_object_name+"::fv",
                                                d_hierarchy,
                                                d_ln_min,
                                                d_ln_max);
      vdb.mapIndexToVariable(f, variable);
#ifdef DEBUG_CHECK_ASSERTIONS
      if ( !variable ) {
         TBOX_ERROR(d_object_name << ": No variable for patch data index "
                    << f << "\n");
      }
      tbox::Pointer<pdat::CellVariable<DIM,double> > cell_variable = variable;
      if ( !cell_variable ) {
         TBOX_ERROR(d_object_name << ": hier::Patch data index " << f
                    << " is not a cell-double variable.\n");
      }
#endif
      d_fv->addComponent( variable, f, d_weight_id );
   }

   return;
}



/*
***********************************************************************
* Delete the vector wrappers.  Do not freeVectorComponents because    *
* we do not control their data allocation.  The user does that.       *
***********************************************************************
*/
template<int DIM> void CellPoissonFACSolver<DIM>::destroyVectorWrappers() {
   d_uv.setNull();
   d_fv.setNull();
   return;
}


}
}
#endif
