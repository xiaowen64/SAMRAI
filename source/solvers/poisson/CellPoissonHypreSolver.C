#ifndef included_solv_CellPoissonHypreSolver_C
#define included_solv_CellPoissonHypreSolver_C

/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-0/source/solvers/poisson/CellPoissonHypreSolver.C $
 * Package:     SAMRAI solvers
 * Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 2294 $
 * Modified:    $LastChangedDate: 2008-07-14 10:34:59 -0700 (Mon, 14 Jul 2008) $
 * Description: Hypre solver interface for diffusion-like elliptic problems.
 */

#include "CellPoissonHypreSolver.h"

#ifdef HAVE_HYPRE

#include <stdlib.h>
#include "BoundaryBoxUtils.h"
#include "CartesianPatchGeometry.h"
#include "CartesianGridGeometry.h"
#include "PatchGeometry.h"
#include "ArrayDataBasicOps.h"
#include "PatchSideDataBasicOps.h"
#include "ArrayData.h"
#include "CellIndex.h"
#include "CellIterator.h"
#include "FaceIndex.h"
#include "SideData.h"
#include "SideIndex.h"
#include "SideVariable.h"
#include "OuterfaceData.h"
#include "OutersideData.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/PIO.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/ShutdownRegistry.h"
#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"


#ifdef DEBUG_NO_INLINE
#include "CellPoissonHypreSolver.I"
#endif

extern "C" {

  void compdiagvariablec2d_(
     double *diag ,
     const double *c ,
     const double *offdiagi ,
     const double *offdiagj ,
     const int *ifirst ,
     const int *ilast ,
     const int *jfirst ,
     const int *jlast ,
     const double *cscale ,
     const double *dscale );
  void compdiagscalarc2d_(
     double *diag ,
     const double *c ,
     const double *offdiagi ,
     const double *offdiagj ,
     const int *ifirst ,
     const int *ilast ,
     const int *jfirst ,
     const int *jlast ,
     const double *cscale ,
     const double *dscale );
  void compdiagzeroc2d_(
     double *diag ,
     const double *offdiagi ,
     const double *offdiagj ,
     const int *ifirst ,
     const int *ilast ,
     const int *jfirst ,
     const int *jlast ,
     const double *cscale ,
     const double *dscale );
  void adjbdry2d_(
     double *diag ,
     const double *offdiagi ,
     const double *offdiagj ,
     const int *pifirst , const int *pilast ,
     const int *pjfirst , const int *pjlast ,
     const double *acoef ,
     const double *bcoef ,
     const int *aifirst , const int *ailast ,
     const int *ajfirst , const int *ajlast ,
     const double *Ak0 ,
     const int *kifirst , const int *kilast ,
     const int *kjfirst , const int *kjlast ,
     const int *lower, const int *upper ,
     const int *location ,
     const double *h );
  void adjbdryconstoffdiags2d_(
     double *diag ,
     const double *offdiag ,
     const int *pifirst ,
     const int *pilast ,
     const int *pjfirst ,
     const int *pjlast ,
     const double *acoef ,
     const int *aifirst ,
     const int *ailast ,
     const int *ajfirst ,
     const int *ajlast ,
     const double *Ak0 ,
     const int *kifirst ,
     const int *kilast ,
     const int *kjfirst ,
     const int *kjlast ,
     const int *lower, const int *upper ,
     const int *location ,
     const double *h );
  void adjustrhs2d_( double *rhs ,
                     const int *rifirst ,
                     const int *rilast ,
                     const int *rjfirst ,
                     const int *rjlast ,
                     const double *Ak0 ,
                     const int *kifirst ,
                     const int *kilast ,
                     const int *kjfirst ,
                     const int *kjlast ,
                     const double *gcoef ,
                     const int *aifirst ,
                     const int *ailast ,
                     const int *ajfirst ,
                     const int *ajlast ,
                     const int *lower, const int *upper ,
                     const int *location );

  void compdiagvariablec3d_(
     double *diag ,
     const double *c ,
     const double *offdiagi ,
     const double *offdiagj ,
     const double *offdiagk ,
     const int *ifirst ,
     const int *ilast ,
     const int *jfirst ,
     const int *jlast ,
     const int *kfirst ,
     const int *klast ,
     const double *cscale ,
     const double *dscale );
  void compdiagscalarc3d_(
     double *diag ,
     const double *c ,
     const double *offdiagi ,
     const double *offdiagj ,
     const double *offdiagk ,
     const int *ifirst ,
     const int *ilast ,
     const int *jfirst ,
     const int *jlast ,
     const int *kfirst ,
     const int *klast ,
     const double *cscale ,
     const double *dscale );
  void compdiagzeroc3d_(
     double *diag ,
     const double *offdiagi ,
     const double *offdiagj ,
     const double *offdiagk ,
     const int *ifirst ,
     const int *ilast ,
     const int *jfirst ,
     const int *jlast ,
     const int *kfirst ,
     const int *klast ,
     const double *cscale ,
     const double *dscale );
  void adjbdry3d_(
     double *diag ,
     const double *offdiagi ,
     const double *offdiagj ,
     const double *offdiagk ,
     const int *pifirst ,
     const int *pilast ,
     const int *pjfirst ,
     const int *pjlast ,
     const int *pkfirst ,
     const int *pklast ,
     const double *acoef ,
     const double *bcoef ,
     const int *aifirst ,
     const int *ailast ,
     const int *ajfirst ,
     const int *ajlast ,
     const int *akfirst ,
     const int *aklast ,
     const double *Ak0 ,
     const int *kifirst ,
     const int *kilast ,
     const int *kjfirst ,
     const int *kjlast ,
     const int *kkfirst ,
     const int *kklast ,
     const int *lower, const int *upper ,
     const int *location ,
     const double *h );
  void adjbdryconstoffdiags3d_(
     double *diag ,
     const double *offdiag ,
     const int *pifirst ,
     const int *pilast ,
     const int *pjfirst ,
     const int *pjlast ,
     const int *pkfirst ,
     const int *pklast ,
     const double *acoef ,
     const int *aifirst ,
     const int *ailast ,
     const int *ajfirst ,
     const int *ajlast ,
     const int *akfirst ,
     const int *aklast ,
     const double *Ak0 ,
     const int *kifirst ,
     const int *kilast ,
     const int *kjfirst ,
     const int *kjlast ,
     const int *kkfirst ,
     const int *kklast ,
     const int *lower, const int *upper ,
     const int *location ,
     const double *h );
  void adjustrhs3d_( double *rhs ,
                     const int *rifirst ,
                     const int *rilast ,
                     const int *rjfirst ,
                     const int *rjlast ,
                     const int *rkfirst ,
                     const int *rklast ,
                     const double *Ak0 ,
                     const int *kifirst ,
                     const int *kilast ,
                     const int *kjfirst ,
                     const int *kjlast ,
                     const int *kkfirst ,
                     const int *kklast ,
                     const double *gcoef ,
                     const int *aifirst ,
                     const int *ailast ,
                     const int *ajfirst ,
                     const int *ajlast ,
                     const int *akfirst ,
                     const int *aklast ,
                     const int *lower, const int *upper ,
                     const int *location );

}

namespace SAMRAI {
    namespace solv {

template<int DIM>
tbox::Pointer<pdat::OutersideVariable<DIM,double> >
CellPoissonHypreSolver<DIM>::s_Ak0_var;

#ifndef NULL
#define NULL (0)
#endif

/*
*************************************************************************
* Constructor                                                           *
*************************************************************************
*/

template<int DIM>  CellPoissonHypreSolver<DIM>::CellPoissonHypreSolver(
   const std::string& object_name,
   tbox::Pointer<tbox::Database> database )
: d_object_name(object_name) ,
  d_hierarchy(NULL) ,
  d_ln(-1) ,
  d_context( hier::VariableDatabase<DIM>::getDatabase()->
             getContext(object_name+"::context") ) ,
  d_physical_bc_coef_strategy(&d_physical_bc_simple_case) ,
  d_physical_bc_variable(NULL),
  d_physical_bc_simple_case(d_object_name+"::simple bc"),
  d_cf_bc_coef(object_name+"::coarse-fine bc coefs") ,
  d_coarsefine_bc_variable(NULL),
  d_Ak0_id(-1) ,
  d_soln_depth(0) ,
  d_rhs_depth(0) ,
  d_max_iterations(10) ,
  d_relative_residual_tol(1e-10) ,
  d_number_iterations(-1) ,
  d_num_pre_relax_steps(1) ,
  d_num_post_relax_steps(1) ,
  d_relative_residual_norm(-1.0) ,
  d_use_smg(false) ,
  d_grid(NULL) ,
  d_stencil(NULL) ,
  d_matrix(NULL) ,
  d_linear_rhs(NULL) ,
  d_linear_sol(NULL) ,
  d_mg_data(NULL) ,
  d_print_solver_info(false)
{

   if (DIM == 1 || DIM > 3) 
   {
      TBOX_ERROR(" CellPoissonHypreSolver : DIM == 1 or > 3 not implemented");
   }
   
   t_solve_system = tbox::TimerManager::getManager()->
      getTimer("solv::CellPoissonHypreSolver::solveSystem()");
   t_set_matrix_coefficients = tbox::TimerManager::getManager()->
      getTimer("solv::CellPoissonHypreSolver::setMatrixCoefficients()");

   hier::VariableDatabase<DIM> *vdb = hier::VariableDatabase<DIM>::getDatabase();
   if ( s_Ak0_var.isNull() ) {
      s_Ak0_var = new
         pdat::OutersideVariable<DIM,double>(d_object_name + "::Ak0", 1);
      tbox::ShutdownRegistry::registerShutdownRoutine(freeVariables, 254);
   }
   d_Ak0_id =
      vdb->registerVariableAndContext( s_Ak0_var,
                                       d_context,
                                       hier::IntVector<DIM>(0) );
   if ( database ) {
      getFromInput(database);
   }
   return;
}


/*
********************************************************************
* Set state from database                                          *
********************************************************************
*/

template<int DIM> void CellPoissonHypreSolver<DIM>::getFromInput(
   tbox::Pointer<tbox::Database> database )
{
   if ( database ) {
      d_print_solver_info = database->getBoolWithDefault("print_solver_info",
                                                         d_print_solver_info);
      d_max_iterations = database->getIntegerWithDefault("max_iterations",
                                                         d_max_iterations);
      d_relative_residual_tol = database->getDoubleWithDefault("relative_residual_tol",
                                                               d_relative_residual_tol);
      if ( database->isDouble("residual_tol") ) {
         TBOX_ERROR("CellPoissonHypreSolver input error.\n"
                    <<"The parameter 'residual_tol' has been replaced\n"
                    <<"by 'relative_residual_tol' to be more descriptive.\n"
                    <<"Please change the parameter name in the input database.");
      }
      d_num_pre_relax_steps =
        database->getIntegerWithDefault("num_pre_relax_steps",
                                        d_num_pre_relax_steps);
      if ( d_num_pre_relax_steps < 0 ) {
	 TBOX_ERROR(d_object_name<< ": Number of relaxation steps must be\n"
		    <<"non-negative.\n");
      }
      d_num_post_relax_steps =
        database->getIntegerWithDefault("num_post_relax_steps",
                                        d_num_post_relax_steps);
      if ( d_num_post_relax_steps < 0 ) {
	 TBOX_ERROR(d_object_name<< ": Number of relaxation steps must be\n"
		    <<"non-negative.\n");
      }
      if ( database->isBool("use_smg") ) {
         bool use_smg = database->getBool("use_smg");
         if ( use_smg != d_use_smg ) {
            setUseSMG(use_smg);
         }
      }
   }
   return;
}

/*
********************************************************************
* Initialize internal data for a given hierarchy level             *
* After setting internal data, propagate the information           *
* to the major algorithm objects.  Allocate data for               *
* storing boundary condition-dependent quantities for              *
* adding to souce term before solving.                             *
********************************************************************
*/

template<int DIM> void CellPoissonHypreSolver<DIM>::initializeSolverState(
   tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy ,
   int ln )
{
   deallocateSolverState();

   d_hierarchy = hierarchy;
   d_ln = ln;

   hier::IntVector<DIM> max_gcw(1);
   d_cf_boundary.computeFromHierarchy(*d_hierarchy, d_ln, max_gcw);

   d_physical_bc_simple_case.setHierarchy(d_hierarchy, d_ln, d_ln);

   d_number_iterations = -1;
   d_relative_residual_norm = -1.0;

   tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(d_ln);
   level->allocatePatchData(d_Ak0_id);
   allocateHypreData();

   return;
}


/*
********************************************************************
* Deallocate data initialized by initializeSolverState             *
********************************************************************
*/

template<int DIM> void CellPoissonHypreSolver<DIM>::deallocateSolverState()
{
   if ( d_hierarchy.isNull() ) return;
   d_cf_boundary.clear();
   tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(d_ln);
   level->deallocatePatchData(d_Ak0_id);
   deallocateHypreData();
   d_hierarchy.setNull();
   d_ln = -1;
   return;
}

/*
*************************************************************************
*                                                                       *
* Allocate the HYPRE data structures that depend only on the level      *
* and will not change (grid, stencil, matrix, and vectors).             *
*                                                                       *
*************************************************************************
*/
template<int DIM> void CellPoissonHypreSolver<DIM>::allocateHypreData()
{
#ifdef HAVE_MPI
   MPI_Comm communicator = tbox::SAMRAI_MPI::getCommunicator();
#else
   MPI_Comm communicator;
#endif

   /*
    * Set up the grid data - only set grid data for local boxes
    */

   tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(d_ln);
   tbox::Pointer< geom::CartesianGridGeometry<DIM> > grid_geometry
     = d_hierarchy->getGridGeometry();
   const hier::IntVector<DIM> ratio = level->getRatio();
   hier::IntVector<DIM> periodic_shift =
      grid_geometry->getPeriodicShift(ratio);

   int periodic_flag[DIM];
   int d;
   bool is_periodic = false;
   for ( d=0; d<DIM; ++d ) {
      periodic_flag[d] = periodic_shift[d] != 0;
      is_periodic = is_periodic || periodic_flag[d];
   }

   HYPRE_StructGridCreate(communicator, DIM, &d_grid);
   for (typename hier::PatchLevel<DIM>::Iterator p(level); p; p++) {
      const hier::Box<DIM>& box = level->getPatch(p())->getBox();
      hier::Index<DIM> lower = box.lower();
      hier::Index<DIM> upper = box.upper();
      HYPRE_StructGridSetExtents(d_grid, lower, upper);
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   if ( is_periodic ) {
      const hier::BoxArray<DIM> &level_domain = level->getPhysicalDomain();
      hier::Box<DIM> domain_bound(level_domain[0]);
      for ( int i=1; i<level_domain.size(); ++i ) {
         domain_bound.lower().min( level_domain[i].lower() );
         domain_bound.upper().max( level_domain[i].upper() );
      }
      for ( d=0; d<DIM; ++d ) {
         if ( periodic_flag[d] == true ) {
            int tmpi=1;
            unsigned int p_of_two;
            for ( p_of_two=0; p_of_two<8*sizeof(p_of_two)-1; ++p_of_two ) {
               if ( tmpi == domain_bound.numberCells(d) ) {
                  break;
               }
               if ( tmpi > domain_bound.numberCells(d) ) {
                  TBOX_ERROR(d_object_name << ": Hypre currently requires\n"
                             <<"that grid size in periodic directions be\n"
                             <<"powers of two.  (This requirement may go\n"
                             <<"away in future versions of hypre.)\n"
                             <<"Size problem in direction " << d << "\n"
                             <<"Domain bound is " << domain_bound << ",\n"
                             <<"Size of " << domain_bound.numberCells() << "\n");
               }
               tmpi = tmpi ? tmpi << 1 : 1;
            }
         }
      }
   }
#endif

   HYPRE_StructGridSetPeriodic( d_grid, periodic_shift );
   HYPRE_StructGridAssemble(d_grid);

   {
      /*
       * Allocate stencil data and set stencil offsets
       */
      
      if (DIM == 1) {
	 const int stencil_size = 2; 
	 int stencil_offsets[2][1] = {
	    { -1 }, { 0 }
	 };
	 HYPRE_StructStencilCreate(DIM, stencil_size, &d_stencil);
	 for (int s = 0; s < stencil_size; s++) {
	    HYPRE_StructStencilSetElement(d_stencil, s,
					  stencil_offsets[s]);
	 }
      } else if (DIM == 2) {
	 const int stencil_size = 3; 
	 int stencil_offsets[3][2] = {
	    { -1, 0 }, { 0, -1}, { 0, 0 }
	 };
	 HYPRE_StructStencilCreate(DIM, stencil_size, &d_stencil);
	 for (int s = 0; s < stencil_size; s++) {
	    HYPRE_StructStencilSetElement(d_stencil, s,
                                       stencil_offsets[s]);
	 }
      } else if (DIM == 3) {
	 const int stencil_size = 4;  
	 int stencil_offsets[4][3] = {
	    { -1,  0,  0}, { 0,  -1,  0}, { 0,  0,  -1}, { 0,  0,  0}
	 };
	 HYPRE_StructStencilCreate(DIM, stencil_size, &d_stencil);
	 for (int s = 0; s < stencil_size; s++) {
	    HYPRE_StructStencilSetElement(d_stencil, s,
					  stencil_offsets[s]);
	 }
      }
   }
   
   
   {
      int full_ghosts1[2*3] = { 1, 1, 0, 0, 0, 0 };
      int no_ghosts1  [2*3] = { 0, 0, 0, 0, 0, 0 };

      int full_ghosts2[2*3] = { 1, 1, 1, 1, 0, 0 };
      int no_ghosts2  [2*3] = { 0, 0, 0, 0, 0, 0 };

      int full_ghosts3[2*3] = { 1, 1, 1, 1, 1, 1 };
      int no_ghosts3  [2*3] = { 0, 0, 0, 0, 0, 0 };

      /*
       * Allocate the structured matrix
       */

      int *full_ghosts;
      int *no_ghosts;

      if (DIM == 1) {
	 full_ghosts = full_ghosts1;
	 no_ghosts = no_ghosts1;
      } else if (DIM == 2) {
	 full_ghosts = full_ghosts2;
	 no_ghosts = no_ghosts2;
      } else if (DIM == 3) {
	 full_ghosts = full_ghosts3;
	 no_ghosts = no_ghosts3;
      }

      HYPRE_StructMatrixCreate(communicator,
                               d_grid,
                               d_stencil,
                               &d_matrix);
      HYPRE_StructMatrixSetNumGhost(d_matrix, full_ghosts);
      HYPRE_StructMatrixSetSymmetric(d_matrix, 1);
      HYPRE_StructMatrixInitialize(d_matrix);

      HYPRE_StructVectorCreate(communicator,
                               d_grid,
                               &d_linear_rhs);
      HYPRE_StructVectorSetNumGhost(d_linear_rhs, no_ghosts);
      HYPRE_StructVectorInitialize(d_linear_rhs);

      HYPRE_StructVectorCreate(communicator,
                               d_grid,
                               &d_linear_sol);
      HYPRE_StructVectorSetNumGhost(d_linear_sol, full_ghosts);
      HYPRE_StructVectorInitialize(d_linear_sol);
   }

   return;
}

/*
*************************************************************************
*                                                                       *
* The destructor deallocates solver data.                               *
*                                                                       *
*************************************************************************
*/

template<int DIM>  CellPoissonHypreSolver<DIM>::~CellPoissonHypreSolver()
{
   deallocateHypreData();

   if ( ! d_hierarchy.isNull() ) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(0);
      level->deallocatePatchData(d_Ak0_id);
   }
   hier::VariableDatabase<DIM> *vdb =
      hier::VariableDatabase<DIM>::getDatabase();
   vdb->removePatchDataIndex(d_Ak0_id);
   return;
}

/*
*************************************************************************
*                                                                       *
* Deallocate HYPRE data and solver.  HYPRE requires that we             *
* check whether HYPRE has already deallocated this data.                *
* Note that the HYPRE solver, d_mg_data, was created at                 *
* the end of setMatrixCoefficients.                                     *
*                                                                       *
*************************************************************************
*/

template<int DIM> void CellPoissonHypreSolver<DIM>::deallocateHypreData()
{
   if (d_stencil) {
      HYPRE_StructStencilDestroy(d_stencil);
   }
   if (d_grid) {
      HYPRE_StructGridDestroy(d_grid);
   }
   if (d_matrix) {
      HYPRE_StructMatrixDestroy(d_matrix);
   }
   if (d_linear_rhs) {
      HYPRE_StructVectorDestroy(d_linear_rhs);
   }
   if (d_linear_sol) {
      HYPRE_StructVectorDestroy(d_linear_sol);
   }
   destroyHypreSolver();
   d_grid = NULL;
   d_stencil = NULL;
   d_matrix = NULL;
   d_linear_rhs = NULL;
   d_linear_sol = NULL;

   return;
}

/*
*************************************************************************
*                                                                       *
* Copy data into the HYPRE vector structures.                           *
*                                                                       *
*************************************************************************
*/

template<int DIM> void CellPoissonHypreSolver<DIM>::copyToHypre(
   HYPRE_StructVector vector,
   pdat::CellData<DIM,double> &src,
   int depth,
   const hier::Box<DIM> &box)
{
// Tracer t("CellPoissonHypreSolver<DIM>::copyToHypre");
   for (pdat::CellIterator<DIM> c(box); c; c++) {
      hier::IntVector<DIM> ic = c();
      HYPRE_StructVectorSetValues(vector, ic, src(c(),depth));
   }
   return;
}

/*
*************************************************************************
*                                                                       *
* Copy data out of the HYPRE vector structures.                         *
*                                                                       *
*************************************************************************
*/

template<int DIM> void CellPoissonHypreSolver<DIM>::copyFromHypre(
   pdat::CellData<DIM,double> &dst,
   int depth,
   HYPRE_StructVector vector,
   const hier::Box<DIM> box)
{
// Tracer t("CellPoissonHypreSolver<DIM>::copyFromHypre");

   for (pdat::CellIterator<DIM> c(box); c; c++) {
      double value;
      hier::IntVector<DIM> ic = c();
      HYPRE_StructVectorGetValues(vector, ic, &value);
      dst(c(),depth) = value;
   }
   return;
}





/*
*************************************************************************
*                                                                       *
* Set the matrix coefficients for the linear system.                    *
* The matrix coefficients are dependent on the problem                  *
* specification described by the PoissonSpecificiations            *
* object and by the boundary condition.                                 *
*                                                                       *
*************************************************************************
*/

template<int DIM> void CellPoissonHypreSolver<DIM>::setMatrixCoefficients(
   const PoissonSpecifications &spec )
{
   if ( d_physical_bc_coef_strategy == NULL ) {
      TBOX_ERROR(d_object_name << ": No BC coefficient strategy object!\n"
		 << "Use either setBoundaries or setPhysicalBcCoefObject\n"
		 << "to specify the boundary conidition.  Do it before\n"
		 << "calling setMatrixCoefficients.");
   }

   t_set_matrix_coefficients->start();

   int i=0;
 
   tbox::Pointer< pdat::CellData<DIM,double> > C_data;
   tbox::Pointer< pdat::SideData<DIM,double> > D_data;

   /*
    * Some computations can be done using high-level math objects.
    * Define the math objects.
    */
   math::ArrayDataBasicOps<DIM,double> array_math;
   math::PatchSideDataBasicOps<DIM,double> patch_side_math;

   /*
    * The value of the ghost cell based on the Robin boundary condition
    * can be written as the sum of a constant, k0, plus a multiple of the
    * internal cell value, k1*ui.  k1*ui depends on the value of u so it
    * contributes to the product Au,
    * while the constant k0 contributes the right hand side f.
    * We save Ak0 = A*k0(a) to add to f when solving.
    * We assume unit g here because we will multiply it in just before
    * solving, thus allowing everything that does not affect A to change
    * from solve to solve.
    */
   tbox::Pointer< pdat::OutersideData<DIM,double> > Ak0;

   /*
    * Loop over patches and set matrix entries for each patch.
    */
   tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(d_ln);
   const hier::IntVector<DIM> no_ghosts(0);
   for (typename hier::PatchLevel<DIM>::Iterator pi(*level); pi; pi++ ) {

      hier::Patch<DIM> &patch = *level->getPatch(*pi);

      tbox::Pointer< geom::CartesianPatchGeometry<DIM> > pg =
         patch.getPatchGeometry();

      const double *h = pg->getDx();
 
      const hier::Box<DIM> patch_box = patch.getBox();
      const hier::Index<DIM> patch_lo = patch_box.lower();
      const hier::Index<DIM> patch_up = patch_box.upper();

      if ( !spec.cIsZero() && !spec.cIsConstant() ) {
         C_data = patch.getPatchData( spec.getCPatchDataId() );
         if ( C_data.isNull() ) {
            TBOX_ERROR(d_object_name << ": Invalid cell variable index "
                       <<  spec.getCPatchDataId()
                       <<" for the C parameter.  It is not\n"
                       <<"cell-centered double data." );
         }
      }

      if ( ! spec.dIsConstant() ) {
         D_data = patch.getPatchData( spec.getDPatchDataId() );
         if ( D_data.isNull() ) {
            TBOX_ERROR(d_object_name << ": Invalid cell variable index "
                       <<  spec.getDPatchDataId()
                       << " for diffusion coefficient.  It is not\n"
                       << "side-centered double data." );
         }
      }

      Ak0 = patch.getPatchData(d_Ak0_id);

      Ak0->fillAll(0.0); 
       
      pdat::CellData<DIM,double> diagonal(patch_box, 1, no_ghosts);

      /*
       * Set diagonals to zero so we can accumulate to it.
       * Accumulation is used at boundaries to shift weights
       * for ghost cells onto the diagonal.
       */
      diagonal.fillAll(0.0);

      const tbox::Pointer< geom::CartesianPatchGeometry<DIM> > 
         geometry = patch.getPatchGeometry();

      const hier::Index<DIM> ifirst = patch_box.lower();
      const hier::Index<DIM> ilast = patch_box.upper();

      /*
       * Storage for off-diagonal entries,
	which can be variable or constant.
       */
      pdat::SideData<DIM,double> off_diagonal(patch_box, 1, no_ghosts);

      /*
       * Compute all off-diagonal entries with no regard to BCs.
       * These off-diagonal entries are simply D/(h*h), according
	to our central difference formula.
       */
      if ( spec.dIsConstant() ) {
        for ( i=0; i<DIM; ++i ) {
          double dhh = spec.getDConstant() / (h[i]*h[i]);
          off_diagonal.fill( dhh );
        }
      }
      else {
        for ( i=0; i<DIM; ++i ) {
          hier::Box<DIM> sbox(patch_box);
          sbox.growUpper( i, 1 );
          array_math.scale( off_diagonal.getArrayData(i) ,
                            1.0/(h[i]*h[i]) ,
                            D_data->getArrayData(i) ,
                            sbox );
        }
      }

      /*
       * Compute diagonal entries using off-diagonal contributions.
       */
      if ( spec.cIsZero() ) {
         computeDiagonalEntries( diagonal,
                                 off_diagonal,
                                 patch_box );
      }
      else if ( spec.cIsConstant() ) {
         computeDiagonalEntries( diagonal,
                                 spec.getCConstant(),
                                 off_diagonal,
                                 patch_box );
      }
      else {
         computeDiagonalEntries( diagonal,
                                 *C_data,
                                 off_diagonal,
                                 patch_box );
      }

      /*
       * Walk physical domain boundaries and adjust off-diagonals
       * before computation of diagonal entries.  
       * The exterior cell's value is
       * uo = ( h*gamma + ui*(beta-h*alpha/2) )/( beta+h*alpha/2 )
       *   = k0 + k1*ui
       * where k0 = h*gamma/( beta+h*alpha/2 )
       * k1 = ( beta-h*alpha/2 )/( beta+h*alpha/2 )
       * Split coupling between interior-exterior cells
       * into two parts: interior-interior coupling (k1)
       * and rhs contribution (k0).
       */
      {
         const tbox::Array< hier::BoundaryBox<DIM> > &surface_boxes =
            pg->getCodimensionBoundaries(1);
         const int n_bdry_boxes = surface_boxes.getSize();
         for ( int n=0; n<n_bdry_boxes; ++n ) {

            const hier::BoundaryBox<DIM> &boundary_box = surface_boxes[n];
            if ( boundary_box.getBoundaryType() != 1 ) {
               TBOX_ERROR(d_object_name << ": Illegal boundary type in "
                          << "CellPoissonHypreSolver<DIM>::setMatrixCoefficients\n");
            }
            const hier::BoundaryBoxUtils<DIM> bbu(boundary_box);
            const int location_index = boundary_box.getLocationIndex();
            const hier::BoundaryBox<DIM> trimmed_boundary_box =
               bbu.trimBoundaryBox( patch.getBox() );
            const hier::Box<DIM> bccoef_box =
               bbu.getSurfaceBoxFromBoundaryBox();
            tbox::Pointer<pdat::ArrayData<DIM,double> >
               acoef_data = new pdat::ArrayData<DIM,double>( bccoef_box, 1 );
            tbox::Pointer<pdat::ArrayData<DIM,double> >
               bcoef_data = new pdat::ArrayData<DIM,double>( bccoef_box, 1 );
            tbox::Pointer<pdat::ArrayData<DIM,double> >
               gcoef_data = NULL;
            static const double fill_time = 0.0;
            d_physical_bc_coef_strategy->setBcCoefs( acoef_data,
                                                     bcoef_data,
                                                     gcoef_data,
                                                     d_physical_bc_variable,
                                                     patch,
                                                     boundary_box,
                                                     fill_time );
            pdat::ArrayData<DIM,double> &Ak0_data =
               Ak0->getArrayData(location_index/2,
                                 location_index%2);
            adjustBoundaryEntries( diagonal,
                                   off_diagonal,
                                   patch_box,
                                   *acoef_data,
                                   *bcoef_data,
                                   bccoef_box,
                                   Ak0_data,
                                   trimmed_boundary_box,
                                   h );
         }
      }

      /*
       * Walk coarse-fine boundaries and adjust off-diagonals
       * according data in ghost cells.
       */
      if ( d_ln > 0 ) {
         /*
          * There are potentially coarse-fine boundaries to deal with.
          */

         tbox::Array< hier::BoundaryBox<DIM> > surface_boxes;

	 if (DIM == 2) {
            surface_boxes = d_cf_boundary.getEdgeBoundaries(*pi);
	 } else if (DIM == 3) {
            surface_boxes = d_cf_boundary.getFaceBoundaries(*pi);
	 } 

         const int n_bdry_boxes = surface_boxes.getSize();
         for ( int n=0; n<n_bdry_boxes; ++n ) {

            const hier::BoundaryBox<DIM> &boundary_box = surface_boxes[n];
            if ( boundary_box.getBoundaryType() != 1 ) {
               TBOX_ERROR(d_object_name << ": Illegal boundary type in "
                          << "CellPoissonHypreSolver<DIM>::setMatrixCoefficients\n");
            }
            const int location_index = boundary_box.getLocationIndex();
            const hier::BoundaryBoxUtils<DIM> bbu(boundary_box);
            const hier::BoundaryBox<DIM> trimmed_boundary_box =
               bbu.trimBoundaryBox( patch.getBox() );
            const hier::Box<DIM> bccoef_box =
               bbu.getSurfaceBoxFromBoundaryBox();
            tbox::Pointer<pdat::ArrayData<DIM,double> >
               acoef_data = new pdat::ArrayData<DIM,double>( bccoef_box, 1 );
            tbox::Pointer<pdat::ArrayData<DIM,double> >
               bcoef_data = new pdat::ArrayData<DIM,double>( bccoef_box, 1 );
            tbox::Pointer<pdat::ArrayData<DIM,double> >
               gcoef_data = NULL;
            static const double fill_time = 0.0;
            /*
	     * Reset invalid ghost data id to help detect use in setBcCoefs.
             */
            d_cf_bc_coef.setGhostDataId( -1 );
            d_cf_bc_coef.setBcCoefs( acoef_data,
                                     bcoef_data,
                                     gcoef_data,
                                     d_coarsefine_bc_variable,
                                     patch,
                                     boundary_box,
                                     fill_time );
            pdat::ArrayData<DIM,double> &Ak0_data =
               Ak0->getArrayData(location_index/2,
                                 location_index%2);
            adjustBoundaryEntries( diagonal,
                                   off_diagonal,
                                   patch_box,
                                   *acoef_data,
                                   *bcoef_data,
                                   bccoef_box,
                                   Ak0_data,
                                   trimmed_boundary_box,
                                   h );
         }
      }

      /*
       * Copy matrix entries to HYPRE matrix structure.  Note that
       * we translate our temporary diagonal/off-diagonal storage into the
       * HYPRE symmetric storage scheme for the stencil specified earlier.
       */
      const int stencil_size = DIM+1;
      int stencil_indices[stencil_size];
      double mat_entries[stencil_size];

      for ( i=0; i<stencil_size; i++ ) stencil_indices[i] = i;

      pdat::CellIterator<DIM> ic(patch_box);

      /*
	To do: This loop uses inefficient high-level syntax.
	See if it can be replaced by a Fortran loop or if we
	can set matrix entries for an entire box at once.
       */
      for ( ; ic; ic++) {

        hier::IntVector<DIM> icell = ic();
        pdat::SideIndex<DIM>  ixlower(ic(),
                                 pdat::SideIndex<DIM>::X,
                                 pdat::SideIndex<DIM>::Lower);
        mat_entries[0] = (off_diagonal)(ixlower);
   
	if (DIM > 1) {
	   pdat::SideIndex<DIM>  iylower(ic(),
					 pdat::SideIndex<DIM>::Y,
					 pdat::SideIndex<DIM>::Lower);
	   mat_entries[1] = (off_diagonal)(iylower);
	}

	if (DIM > 2) {
	   pdat::SideIndex<DIM>  izlower(ic(),
					 pdat::SideIndex<DIM>::Z,
					 pdat::SideIndex<DIM>::Lower);
	   // The "funny" indexing prevents a warning when compiling for 
	   // DIM < 2.  This code is only reached if DIM > 2 when 
	   // executing.
	   mat_entries[DIM > 2 ? 2 : 0] = (off_diagonal)(izlower);
	}

        mat_entries[DIM] = (diagonal)(ic());
        HYPRE_StructMatrixSetValues(d_matrix, icell,
                                        stencil_size, stencil_indices,
                                        mat_entries);
      } // end cell loop

   } // end patch loop

   if (d_print_solver_info) {
      HYPRE_StructMatrixPrint("mat_bA.out",d_matrix,1);
   }
   
   HYPRE_StructMatrixAssemble(d_matrix);

   if (d_print_solver_info) {
      HYPRE_StructMatrixPrint("mat_aA.out",d_matrix,1);
   }

   t_set_matrix_coefficients->stop();

   setupHypreSolver();

   return;
}




/*
**********************************************************************
* Add g*A*k0(a) from physical boundaries to rhs.                     *
* This operation is done for physical as well as cf boundaries,      *
* so it is placed in a function.                                     *
**********************************************************************
*/

template<int DIM> void CellPoissonHypreSolver<DIM>::add_gAk0_toRhs(
   const hier::Patch<DIM> &patch,
   const tbox::Array< hier::BoundaryBox<DIM> > &bdry_boxes,
   const RobinBcCoefStrategy<DIM> *robin_bc_coef,
   pdat::CellData<DIM,double> &rhs )
{
   /*
    * g*A*k0(a) is the storage for adjustments to be made to the rhs
    * when we solve. This is the value of the weight of the ghost cell
    * value for the interior cell, times k0.  It is independent of u,
    * and so is moved to the rhs.  Before solving, g*A*k0(a) is added
    * to rhs.
    */
   tbox::Pointer< pdat::OutersideData<DIM,double> > Ak0;

   tbox::Pointer< geom::CartesianPatchGeometry<DIM> > pg =
      patch.getPatchGeometry();

   Ak0 = patch.getPatchData(d_Ak0_id);

   const int n_bdry_boxes = bdry_boxes.getSize();
   for ( int n=0; n<n_bdry_boxes; ++n ) {

      const hier::BoundaryBox<DIM> &boundary_box = bdry_boxes[n];
#ifdef DEBUG_CHECK_ASSERTIONS
      if ( boundary_box.getBoundaryType() != 1 ) {
         TBOX_ERROR(d_object_name << ": Illegal boundary type in "
                    << "CellPoissonHypreSolver<DIM>::add_gAk0_toRhs\n");
      }
#endif
      const int location_index = boundary_box.getLocationIndex();
      const hier::BoundaryBoxUtils<DIM> bbu(boundary_box);
      const hier::BoundaryBox<DIM> trimmed_boundary_box =
         bbu.trimBoundaryBox( patch.getBox() );
      const hier::Index<DIM> &lower = trimmed_boundary_box.getBox().lower();
      const hier::Index<DIM> &upper = trimmed_boundary_box.getBox().upper();
      const hier::Box<DIM> &rhsbox = rhs.getArrayData().getBox();
      const hier::Box<DIM> &Ak0box = Ak0->getArrayData(location_index/2,
                                                  location_index%2).getBox();
      const hier::Box<DIM> bccoef_box = bbu.getSurfaceBoxFromBoundaryBox();
      tbox::Pointer<pdat::ArrayData<DIM,double> >
         acoef_data = NULL;
      tbox::Pointer<pdat::ArrayData<DIM,double> >
         bcoef_data = NULL;
      tbox::Pointer<pdat::ArrayData<DIM,double> >
         gcoef_data = new pdat::ArrayData<DIM,double>( bccoef_box, 1 );
      static const double fill_time = 0.0;
      robin_bc_coef->setBcCoefs( acoef_data,
                                 bcoef_data,
                                 gcoef_data,
                                 d_physical_bc_variable,
                                 patch,
                                 boundary_box,
                                 fill_time );
      /*
       * Nomenclature for indices: cel=first-cell, gho=ghost,
       * beg=beginning, end=ending.
       */
      if (DIM == 2) {
	 adjustrhs2d_( rhs.getPointer(d_rhs_depth),
		       &rhsbox.lower()[0],
		       &rhsbox.upper()[0],
		       &rhsbox.lower()[1],
		       &rhsbox.upper()[1],
		       Ak0->getPointer(location_index/2,location_index%2),
		       &Ak0box.lower()[0],
		       &Ak0box.upper()[0],
		       &Ak0box.lower()[1],
		       &Ak0box.upper()[1],
		       gcoef_data->getPointer() ,
		       &bccoef_box.lower()[0] ,
		       &bccoef_box.upper()[0] ,
		       &bccoef_box.lower()[1] ,
		       &bccoef_box.upper()[1] ,
		       lower, upper,
		       &location_index );
      } else if (DIM == 3) { 
	 adjustrhs3d_( rhs.getPointer(d_rhs_depth),
		       &rhsbox.lower()[0],
		       &rhsbox.upper()[0],
		       &rhsbox.lower()[1],
		       &rhsbox.upper()[1],
		       &rhsbox.lower()[2],
		       &rhsbox.upper()[2],
		       Ak0->getPointer(location_index/2,location_index%2),
		       &Ak0box.lower()[0],
		       &Ak0box.upper()[0],
		       &Ak0box.lower()[1],
		       &Ak0box.upper()[1],
		       &Ak0box.lower()[2],
		       &Ak0box.upper()[2],
		       gcoef_data->getPointer() ,
		       &bccoef_box.lower()[0] ,
		       &bccoef_box.upper()[0] ,
		       &bccoef_box.lower()[1] ,
		       &bccoef_box.upper()[1] ,
		       &bccoef_box.lower()[2] ,
		       &bccoef_box.upper()[2] ,
		       lower, upper,
		       &location_index );
      }
   }
   return;
}




/*
*************************************************************************
* Create the hypre solver and set it according to the current state.    *
*************************************************************************
*/
template<int DIM> void CellPoissonHypreSolver<DIM>::setupHypreSolver()
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( d_mg_data == NULL );
#endif

#ifdef HAVE_MPI
   MPI_Comm communicator = tbox::SAMRAI_MPI::getCommunicator();
#else
   MPI_Comm communicator;
#endif

   if ( d_use_smg ) {
      HYPRE_StructSMGCreate(communicator, &d_mg_data);
      HYPRE_StructSMGSetMemoryUse(d_mg_data, 0);
      HYPRE_StructSMGSetMaxIter(d_mg_data, d_max_iterations);
      HYPRE_StructSMGSetTol(d_mg_data, d_relative_residual_tol);
      HYPRE_StructSMGSetLogging(d_mg_data, 1);
      HYPRE_StructSMGSetNumPreRelax(d_mg_data,
                                    d_num_pre_relax_steps);
      HYPRE_StructSMGSetNumPostRelax(d_mg_data,
                                     d_num_post_relax_steps);
      HYPRE_StructSMGSetup(d_mg_data,
                           d_matrix,
                           d_linear_rhs,
                           d_linear_sol);
   } 
   else {
      HYPRE_StructPFMGCreate(communicator, &d_mg_data);
      HYPRE_StructPFMGSetMaxIter(d_mg_data, d_max_iterations);
      HYPRE_StructPFMGSetTol(d_mg_data, d_relative_residual_tol);
      HYPRE_StructPFMGSetLogging(d_mg_data, 1);
      HYPRE_StructPFMGSetNumPreRelax(d_mg_data,
                                     d_num_pre_relax_steps);
      HYPRE_StructPFMGSetNumPostRelax(d_mg_data,
                                      d_num_post_relax_steps);
      HYPRE_StructPFMGSetup(d_mg_data,
                            d_matrix,
                            d_linear_rhs,
                            d_linear_sol);
   }
   return;
}




template<int DIM> void CellPoissonHypreSolver<DIM>::destroyHypreSolver()
{
   if ( d_mg_data != NULL ) {
      if (d_use_smg) {
         HYPRE_StructSMGDestroy(d_mg_data);
      }
      else {
         HYPRE_StructPFMGDestroy(d_mg_data);
      }
      d_mg_data = NULL;
   }
   return;
}





/*
*************************************************************************
*                                                                       *
* Solve the linear system.  This routine assumes that the boundary      *
* conditions and the matrix coefficients have been specified.           *
*                                                                       *
*************************************************************************
*/

template<int DIM> int CellPoissonHypreSolver<DIM>::solveSystem( const int u ,
                                               const int f ,
                                               bool homogeneous_bc )
{
   if ( d_physical_bc_coef_strategy == NULL ) {
      TBOX_ERROR(d_object_name << ": No BC coefficient strategy object!\n"
		 << "Use either setBoundaries or setPhysicalBcCoefObject\n"
		 << "to specify the boundary conidition.  Do it before\n"
		 << "calling solveSystem.");
   }
   // Tracer t("CellPoissonHypreSolver<DIM>::solveSystem");

   t_solve_system->start();

   tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(d_ln);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(u >= 0);
   TBOX_ASSERT(u < level->getPatchDescriptor()->getMaxNumberRegisteredComponents());
   TBOX_ASSERT(f >= 0);
   TBOX_ASSERT(f < level->getPatchDescriptor()->getMaxNumberRegisteredComponents());
#endif


   if ( d_physical_bc_coef_strategy == &d_physical_bc_simple_case ) {
      /*
       * If we are using the simple bc implementation, the final piece
       * of information it requires is the Dirichlet boundary value
       * set in the ghost cells.  Now that we have the ghost cell data,
       * we can complete the boundary condition setup.
       */
      d_physical_bc_simple_case.cacheDirichletData(u);
   }
   

   /*
    * Modify right-hand-side to account for boundary conditions and
    * copy solution and right-hand-side to HYPRE structures.  
    */

   const hier::IntVector<DIM> no_ghosts(0);
   const hier::IntVector<DIM> ghosts(1);

   /*
    * At coarse-fine boundaries, we expect ghost cells to have correct
    * values to be used in our bc, so u provides the ghost cell data.
    * Assume that the user only provided data for the immediate first
    * ghost cell, so pass zero for the number of extensions fillable.
    */
   d_cf_bc_coef.setGhostDataId( u, hier::IntVector<DIM>(0) );

   for (typename hier::PatchLevel<DIM>::Iterator p(level); p; p++) {
      tbox::Pointer< hier::Patch<DIM> > patch = level->getPatch(p());

      const hier::Box<DIM> box = patch->getBox();

      /*
       * Set up variable data needed to prepare linear system solver.
       */
      tbox::Pointer<pdat::CellData<DIM,double> > u_data_ = patch->getPatchData(u);
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT( !u_data_.isNull() );
#endif
      pdat::CellData<DIM,double> &u_data = *u_data_;
      pdat::CellData<DIM,double> rhs_data(box,1,no_ghosts);

      /*
       * Copy rhs and solution from the hierarchy into HYPRE structures.
       * For rhs, add in the contribution from boundary conditions, if
       * needed.  If boundary condition is homogenous, this only adds
	zero, so we skip it.
       */
      copyToHypre(d_linear_sol, u_data, d_soln_depth, box );
      rhs_data.copy( *(patch->getPatchData(f)) );
      if ( ! homogeneous_bc ) {
         /*
	   Add g*A*k0(a) from physical and coarse-fine boundaries to rhs.
          */
         add_gAk0_toRhs(*patch,
                        patch->getPatchGeometry()->getCodimensionBoundaries(1),
                        d_physical_bc_coef_strategy,
                        rhs_data);
         add_gAk0_toRhs(*patch,
                        d_cf_boundary.getBoundaries(patch->getPatchNumber(),1),
                        &d_cf_bc_coef,
                        rhs_data);
      }
      copyToHypre(d_linear_rhs, rhs_data, d_rhs_depth, box);

   } // end patch loop

   /*
    * Reset invalid ghost data id to help detect erroneous further use.
    */
   d_cf_bc_coef.setGhostDataId( -1 );

   /*
    * Finish assembly of the vectors
    */
   HYPRE_StructVectorAssemble(d_linear_sol);

   HYPRE_StructVectorAssemble(d_linear_rhs);

   /*
    * Solve the system - zero means convergence
    * Solve takes the same arguments as Setup
    */

   if (d_print_solver_info) {
      HYPRE_StructVectorPrint("sol0.out", d_linear_sol, 1);
      HYPRE_StructMatrixPrint("mat0.out",d_matrix,1);
      HYPRE_StructVectorPrint("rhs.out",d_linear_rhs,1);
   }

   if ( d_use_smg ) {
      // HYPRE_StructSMGSetMaxIter(d_mg_data, d_max_iterations);
      HYPRE_StructSMGSetTol(d_mg_data, d_relative_residual_tol);
      /* converge = */ HYPRE_StructSMGSolve(d_mg_data,
                                      d_matrix,
                                      d_linear_rhs,
                                      d_linear_sol);
   } 
   else {
      // HYPRE_StructPFMGSetMaxIter(d_mg_data, d_max_iterations);
      HYPRE_StructPFMGSetTol(d_mg_data, d_relative_residual_tol);
      /* converge = */ HYPRE_StructPFMGSolve(d_mg_data,
                                       d_matrix,
                                       d_linear_rhs,
                                       d_linear_sol);
   }

   if (d_print_solver_info) {
      HYPRE_StructMatrixPrint("mat.out",d_matrix,1);
      HYPRE_StructVectorPrint("sol.out",d_linear_sol,1);
   }

   if ( d_use_smg ) {
      HYPRE_StructSMGGetNumIterations(d_mg_data,
                                      &d_number_iterations);
      HYPRE_StructSMGGetFinalRelativeResidualNorm(d_mg_data,
                                                  &d_relative_residual_norm);
   } 
   else {
      HYPRE_StructPFMGGetNumIterations(d_mg_data,
                                       &d_number_iterations);
      HYPRE_StructPFMGGetFinalRelativeResidualNorm(d_mg_data, 
                                                   &d_relative_residual_norm);
   }

   /*
    * Pull the solution vector out of the HYPRE structures
    */
   for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++) {
      tbox::Pointer< hier::Patch<DIM> > patch = level->getPatch(ip());
      tbox::Pointer<pdat::CellData<DIM,double> > u_data_ = patch->getPatchData(u);
      pdat::CellData<DIM,double> &u_data = *u_data_;
      copyFromHypre(u_data,
                    d_soln_depth,
                    d_linear_sol,
                    patch->getBox());
   }

   t_solve_system->stop();

   return( d_relative_residual_norm <= d_relative_residual_tol );
}



template<int DIM> void CellPoissonHypreSolver<DIM>::computeDiagonalEntries(
   pdat::CellData<DIM,double> &diagonal,
   const pdat::CellData<DIM,double> &C_data,
   const pdat::SideData<DIM,double> &off_diagonal,
   const hier::Box<DIM> &patch_box )
{
   const hier::Index<DIM> patch_lo = patch_box.lower();
   const hier::Index<DIM> patch_up = patch_box.upper();
   const double c=1.0, d=1.0;
   if (DIM == 2) {
      compdiagvariablec2d_( diagonal.getPointer(),
			    C_data.getPointer(),
			    off_diagonal.getPointer(0),
			    off_diagonal.getPointer(1),
			    &patch_lo[0], &patch_up[0],
			    &patch_lo[1], &patch_up[1],
			    &c, &d );
   } else if (DIM == 3) {
      compdiagvariablec3d_( diagonal.getPointer(),
			    C_data.getPointer(),
			    off_diagonal.getPointer(0),
			    off_diagonal.getPointer(1),
                         off_diagonal.getPointer(2),
			    &patch_lo[0], &patch_up[0],
			    &patch_lo[1], &patch_up[1],
			    &patch_lo[2], &patch_up[2],
			    &c, &d );
   }

   return;
}



template<int DIM> void CellPoissonHypreSolver<DIM>::computeDiagonalEntries(
   pdat::CellData<DIM,double> &diagonal,
   const double C,
   const pdat::SideData<DIM,double> &off_diagonal,
   const hier::Box<DIM> &patch_box )
{
   const hier::Index<DIM> patch_lo = patch_box.lower();
   const hier::Index<DIM> patch_up = patch_box.upper();
   const double c=1.0, d=1.0;
   if (DIM == 2) {
      compdiagscalarc2d_( diagonal.getPointer(),
			  &C,
			  off_diagonal.getPointer(0),
			  off_diagonal.getPointer(1),
			  &patch_lo[0], &patch_up[0],
			  &patch_lo[1], &patch_up[1],
			  &c, &d );
   } else if (DIM == 3) {
      compdiagscalarc3d_( diagonal.getPointer(),
			  &C,
			  off_diagonal.getPointer(0),
                       off_diagonal.getPointer(1),
			  off_diagonal.getPointer(2),
			  &patch_lo[0], &patch_up[0],
			  &patch_lo[1], &patch_up[1],
			  &patch_lo[2], &patch_up[2],
			  &c, &d );
   } else {
      TBOX_ERROR("CellPoissonHypreSolver error...\n"
		       << "DIM > 3 not supported." << std::endl);
   }

   return;
}



template<int DIM> void CellPoissonHypreSolver<DIM>::computeDiagonalEntries(
   pdat::CellData<DIM,double> &diagonal,
   const pdat::SideData<DIM,double> &off_diagonal,
   const hier::Box<DIM> &patch_box )
{
   const hier::Index<DIM> patch_lo = patch_box.lower();
   const hier::Index<DIM> patch_up = patch_box.upper();
   const double c=1.0, d=1.0;
   if (DIM == 2) {
      compdiagzeroc2d_( diagonal.getPointer(),
			off_diagonal.getPointer(0),
			off_diagonal.getPointer(1),
			&patch_lo[0], &patch_up[0],
			&patch_lo[1], &patch_up[1],
			&c, &d );
   } else if (DIM == 3) {
      compdiagzeroc3d_( diagonal.getPointer(),
			off_diagonal.getPointer(0),
			off_diagonal.getPointer(1),
			off_diagonal.getPointer(2),
			&patch_lo[0], &patch_up[0],
			&patch_lo[1], &patch_up[1],
			&patch_lo[2], &patch_up[2],
			&c, &d );
   } else {
      TBOX_ERROR("CellPoissonHypreSolver error...\n"
		 << "DIM > 3 not supported." << std::endl);
   }

   return;
}




template<int DIM> void CellPoissonHypreSolver<DIM>::adjustBoundaryEntries(
   pdat::CellData<DIM,double> &diagonal,
   const pdat::SideData<DIM,double> &off_diagonal,
   const hier::Box<DIM> &patch_box,
   const pdat::ArrayData<DIM,double> &acoef_data,
   const pdat::ArrayData<DIM,double> &bcoef_data,
   const hier::Box<DIM> bccoef_box,
   pdat::ArrayData<DIM,double> &Ak0_data,
   const hier::BoundaryBox<DIM> &trimmed_boundary_box,
   const double h[DIM] )
{
   const hier::Index<DIM> patch_lo = patch_box.lower();
   const hier::Index<DIM> patch_up = patch_box.upper();
   const int location_index = trimmed_boundary_box.getLocationIndex();
   const hier::Index<DIM> &lower = trimmed_boundary_box.getBox().lower();
   const hier::Index<DIM> &upper = trimmed_boundary_box.getBox().upper();
   const hier::Box<DIM> &Ak0_box = Ak0_data.getBox();
   if (DIM == 2) {
      adjbdry2d_( diagonal.getPointer() ,
		  off_diagonal.getPointer(0) ,
		  off_diagonal.getPointer(1) ,
		  &patch_lo[0], &patch_up[0] ,
		  &patch_lo[1], &patch_up[1] ,
		  acoef_data.getPointer() ,
		  bcoef_data.getPointer() ,
		  &bccoef_box.lower()[0] ,
		  &bccoef_box.upper()[0] ,
		  &bccoef_box.lower()[1] ,
		  &bccoef_box.upper()[1] ,
		  Ak0_data.getPointer() ,
		  &Ak0_box.lower()[0] ,
		  &Ak0_box.upper()[0] ,
		  &Ak0_box.lower()[1] ,
		  &Ak0_box.upper()[1] ,
		  lower, upper ,
		  &location_index, h );
   } else if (DIM == 3) {
      adjbdry3d_( diagonal.getPointer() ,
		  off_diagonal.getPointer(0) ,
		  off_diagonal.getPointer(1) ,
		  off_diagonal.getPointer(2) ,
		  &patch_lo[0], &patch_up[0] ,
		  &patch_lo[1], &patch_up[1] ,
		  &patch_lo[2], &patch_up[2] ,
		  acoef_data.getPointer() ,
		  bcoef_data.getPointer() ,
		  &bccoef_box.lower()[0] ,
		  &bccoef_box.upper()[0] ,
		  &bccoef_box.lower()[1] ,
		  &bccoef_box.upper()[1] ,
		  &bccoef_box.lower()[2] ,
		  &bccoef_box.upper()[2] ,
		  Ak0_data.getPointer() ,
		  &Ak0_box.lower()[0] ,
		  &Ak0_box.upper()[0] ,
		  &Ak0_box.lower()[1] ,
		  &Ak0_box.upper()[1] ,
		  &Ak0_box.lower()[2] ,
		  &Ak0_box.upper()[2] ,
		  lower, upper ,
		  &location_index, h );
   } else {
      TBOX_ERROR("CellPoissonHypreSolver error...\n"
		 << "DIM > 3 not supported." << std::endl);
   }
}




template<int DIM> void
CellPoissonHypreSolver<DIM>::freeVariables()
{
   s_Ak0_var.setNull();
   return;
}

}
}

#endif
#endif
