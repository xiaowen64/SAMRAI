#ifndef included_solv_CellPoissonFACOps_C
#define included_solv_CellPoissonFACOps_C

/*
 * File:        CellPoissonFACOps.C
 * Copyright:   (c) 1997-2005 The Regents of the University of California
 * Revision:    $Revision: 719 $
 * Modified:    $Date: 2005-11-10 11:45:11 -0800 (Thu, 10 Nov 2005) $
 * Description: Operator class for cell-centered scalar Poisson using FAC
 */

#include "CellPoissonFACOps.h"

#include IOMANIP_HEADER_FILE

#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "Index.h"
#include "Variable.h"
#include "VariableDatabase.h"
#include "CellDoubleConstantRefine.h"
#include "CellVariable.h"
#include "OutersideData.h"
#include "OutersideVariable.h"
#include "PatchData.h"
#include "SideVariable.h"
#include "FACPreconditioner.h"
#include "CellPoissonHypreSolver.h"
#include "tbox/Array.h"
#include "tbox/ShutdownRegistry.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"
#include "CoarsenAlgorithm.h"
#include "CoarsenOperator.h"
#include "CoarsenSchedule.h"
#include "RefineAlgorithm.h"
#include "RefineOperator.h"
#include "RefineSchedule.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

#ifdef DEBUG_NO_INLINE
#include "CellPoissonFACOps.I"
#endif

namespace SAMRAI {

    namespace solv {

template<int DIM>
tbox::Pointer<pdat::CellVariable<DIM,double> >
CellPoissonFACOps<DIM>::s_cell_scratch_var;

template<int DIM>
tbox::Pointer<pdat::SideVariable<DIM,double> >
CellPoissonFACOps<DIM>::s_flux_scratch_var;

template<int DIM>
tbox::Pointer<pdat::OutersideVariable<DIM,double> >
CellPoissonFACOps<DIM>::s_oflux_scratch_var;


extern "C" {

   void compfluxvardc2d_(
      double *xflux ,
      double *yflux ,
      const int *fluxgi ,
      const int *fluxgj ,
      const double *xdiff_coef ,
      const double *ydiff_coef ,
      const int *dcgi ,
      const int *dcgj ,
      const double *soln ,
      const int *solngi ,
      const int *solngj ,
      const int *ifirst ,
      const int *ilast ,
      const int *jfirst ,
      const int *jlast ,
      const double *dx );
   void compfluxcondc2d_(
      double *xflux ,
      double *yflux ,
      const int *fluxgi ,
      const int *fluxgj ,
      const double &diff_coef ,
      const double *soln ,
      const int *solngi ,
      const int *solngj ,
      const int *ifirst ,
      const int *ilast ,
      const int *jfirst ,
      const int *jlast ,
      const double *dx );
   void rbgswithfluxmaxvardcvarsf2d_(
      const double *xflux ,
      const double *yflux ,
      const int *fluxgi ,
      const int *fluxgj ,
      const double *xdiff_coef ,
      const double *ydiff_coef ,
      const int *dcgi ,
      const int *dcgj ,
      const double *rhs ,
      const int *rhsgi ,
      const int *rhsgj ,
      const double *scalar_field ,
      const int *scalar_field_gi ,
      const int *scalar_field_gj ,
      double *soln ,
      const int *solngi ,
      const int *solngj ,
      const int *ifirst ,
      const int *ilast ,
      const int *jfirst ,
      const int *jlast ,
      const double *dx ,
      const int *offset ,
      const double *maxres );
   void rbgswithfluxmaxcondcvarsf2d_(
      const double *xflux ,
      const double *yflux ,
      const int *fluxgi ,
      const int *fluxgj ,
      const double &dc ,
      const double *rhs ,
      const int *rhsgi ,
      const int *rhsgj ,
      const double *scalar_field ,
      const int *scalar_field_gi ,
      const int *scalar_field_gj ,
      double *soln ,
      const int *solngi ,
      const int *solngj ,
      const int *ifirst ,
      const int *ilast ,
      const int *jfirst ,
      const int *jlast ,
      const double *dx ,
      const int *offset ,
      const double *maxres );
   void rbgswithfluxmaxvardcconsf2d_(
      const double *xflux ,
      const double *yflux ,
      const int *fluxgi ,
      const int *fluxgj ,
      const double *xdiff_coef ,
      const double *ydiff_coef ,
      const int *dcgi ,
      const int *dcgj ,
      const double *rhs ,
      const int *rhsgi ,
      const int *rhsgj ,
      const double &scalar_field ,
      double *soln ,
      const int *solngi ,
      const int *solngj ,
      const int *ifirst ,
      const int *ilast ,
      const int *jfirst ,
      const int *jlast ,
      const double *dx ,
      const int *offset ,
      const double *maxres );
   void rbgswithfluxmaxcondcconsf2d_(
      const double *xflux ,
      const double *yflux ,
      const int *fluxgi ,
      const int *fluxgj ,
      const double &dc ,
      const double *rhs ,
      const int *rhsgi ,
      const int *rhsgj ,
      const double &scalar_field ,
      double *soln ,
      const int *solngi ,
      const int *solngj ,
      const int *ifirst ,
      const int *ilast ,
      const int *jfirst ,
      const int *jlast ,
      const double *dx ,
      const int *offset ,
      const double *maxres );
   void compresvarsca2d_(
      const double *xflux ,
      const double *yflux ,
      const int *fluxgi ,
      const int *fluxgj ,
      const double *rhs ,
      const int *rhsgi ,
      const int *rhsgj ,
      double *residual ,
      const int *residualgi ,
      const int *residualgj ,
      const double *scalar_field ,
      const int *scalar_field_gi ,
      const int *scalar_field_gj ,
      const double *soln ,
      const int *solngi ,
      const int *solngj ,
      const int *ifirst ,
      const int *ilast ,
      const int *jfirst ,
      const int *jlast ,
      const double *dx );
   void compresconsca2d_(
      const double *xflux ,
      const double *yflux ,
      const int *fluxgi ,
      const int *fluxgj ,
      const double *rhs ,
      const int *rhsgi ,
      const int *rhsgj ,
      double *residual ,
      const int *residualgi ,
      const int *residualgj ,
      const double &scalar_field ,
      const double *soln ,
      const int *solngi ,
      const int *solngj ,
      const int *ifirst ,
      const int *ilast ,
      const int *jfirst ,
      const int *jlast ,
      const double *dx );
   void ewingfixfluxvardc2d_(
      const double *xflux ,
      const double *yflux ,
      const int *fluxgi ,
      const int *fluxgj ,
      const double *xdiff_coef ,
      const double *ydiff_coef ,
      const int *dcgi ,
      const int *dcgj ,
      const double *soln ,
      const int *solngi ,
      const int *solngj ,
      const int *ifirst ,
      const int *ilast ,
      const int *jfirst ,
      const int *jlast ,
      const int *location_index,
      const int *ratio_to_coarser,
      const int *blower,
      const int *bupper,
      const double *dx );
   void ewingfixfluxcondc2d_(
      const double *xflux ,
      const double *yflux ,
      const int *fluxgi ,
      const int *fluxgj ,
      const double &diff_coef ,
      const double *soln ,
      const int *solngi ,
      const int *solngj ,
      const int *ifirst ,
      const int *ilast ,
      const int *jfirst ,
      const int *jlast ,
      const int *location_index,
      const int *ratio_to_coarser,
      const int *blower,
      const int *bupper,
      const double *dx );

   void compfluxvardc3d_(
      double *xflux ,
      double *yflux ,
      double *zflux ,
      const int *fluxgi ,
      const int *fluxgj ,
      const int *fluxgk ,
      const double *xdiff_coef ,
      const double *ydiff_coef ,
      const double *zdiff_coef ,
      const int *dcgi ,
      const int *dcgj ,
      const int *dcgk ,
      const double *soln ,
      const int *solngi ,
      const int *solngj ,
      const int *solngk ,
      const int *ifirst ,
      const int *ilast ,
      const int *jfirst ,
      const int *jlast ,
      const int *kfirst ,
      const int *klast ,
      const double *dx );
   void compfluxcondc3d_(
      double *xflux ,
      double *yflux ,
      double *zflux ,
      const int *fluxgi ,
      const int *fluxgj ,
      const int *fluxgk ,
      const double &diff_coef ,
      const double *soln ,
      const int *solngi ,
      const int *solngj ,
      const int *solngk ,
      const int *ifirst ,
      const int *ilast ,
      const int *jfirst ,
      const int *jlast ,
      const int *kfirst ,
      const int *klast ,
      const double *dx );
   void rbgswithfluxmaxvardcvarsf3d_(
      const double *xflux ,
      const double *yflux ,
      const double *zflux ,
      const int *fluxgi ,
      const int *fluxgj ,
      const int *fluxgk ,
      const double *xdiff_coef ,
      const double *ydiff_coef ,
      const double *zdiff_coef ,
      const int *dcgi ,
      const int *dcgj ,
      const int *dcgk ,
      const double *rhs ,
      const int *rhsgi ,
      const int *rhsgj ,
      const int *rhsgk ,
      const double *scalar_field ,
      const int *scalar_field_gi ,
      const int *scalar_field_gj ,
      const int *scalar_field_gk ,
      double *soln ,
      const int *solngi ,
      const int *solngj ,
      const int *solngk ,
      const int *ifirst ,
      const int *ilast ,
      const int *jfirst ,
      const int *jlast ,
      const int *kfirst ,
      const int *klast ,
      const double *dx ,
      const int *offset ,
      const double *maxres );
   void rbgswithfluxmaxcondcvarsf3d_(
      const double *xflux ,
      const double *yflux ,
      const double *zflux ,
      const int *fluxgi ,
      const int *fluxgj ,
      const int *fluxgk ,
      const double &dc ,
      const double *rhs ,
      const int *rhsgi ,
      const int *rhsgj ,
      const int *rhsgk ,
      const double *scalar_field ,
      const int *scalar_field_gi ,
      const int *scalar_field_gj ,
      const int *scalar_field_gk ,
      double *soln ,
      const int *solngi ,
      const int *solngj ,
      const int *solngk ,
      const int *ifirst ,
      const int *ilast ,
      const int *jfirst ,
      const int *jlast ,
      const int *kfirst ,
      const int *klast ,
      const double *dx ,
      const int *offset ,
      const double *maxres );
   void rbgswithfluxmaxvardcconsf3d_(
      const double *xflux ,
      const double *yflux ,
      const double *zflux ,
      const int *fluxgi ,
      const int *fluxgj ,
      const int *fluxgk ,
      const double *xdiff_coef ,
      const double *ydiff_coef ,
      const double *zdiff_coef ,
      const int *dcgi ,
      const int *dcgj ,
      const int *dcgk ,
      const double *rhs ,
      const int *rhsgi ,
      const int *rhsgj ,
      const int *rhsgk ,
      const double &scalar_field ,
      double *soln ,
      const int *solngi ,
      const int *solngj ,
      const int *solngk ,
      const int *ifirst ,
      const int *ilast ,
      const int *jfirst ,
      const int *jlast ,
      const int *kfirst ,
      const int *klast ,
      const double *dx ,
      const int *offset ,
      const double *maxres );
   void rbgswithfluxmaxcondcconsf3d_(
      const double *xflux ,
      const double *yflux ,
      const double *zflux ,
      const int *fluxgi ,
      const int *fluxgj ,
      const int *fluxgk ,
      const double &dc ,
      const double *rhs ,
      const int *rhsgi ,
      const int *rhsgj ,
      const int *rhsgk ,
      const double &scalar_field ,
      double *soln ,
      const int *solngi ,
      const int *solngj ,
      const int *solngk ,
      const int *ifirst ,
      const int *ilast ,
      const int *jfirst ,
      const int *jlast ,
      const int *kfirst ,
      const int *klast ,
      const double *dx ,
      const int *offset ,
      const double *maxres );
   void compresvarsca3d_(
      const double *xflux ,
      const double *yflux ,
      const double *zflux ,
      const int *fluxgi ,
      const int *fluxgj ,
      const int *fluxgk ,
      const double *rhs ,
      const int *rhsgi ,
      const int *rhsgj ,
      const int *rhsgk ,
      double *residual ,
      const int *residualgi ,
      const int *residualgj ,
      const int *residualgk ,
      const double *scalar_field ,
      const int *scalar_field_gi ,
      const int *scalar_field_gj ,
      const int *scalar_field_gk ,
      const double *soln ,
      const int *solngi ,
      const int *solngj ,
      const int *solngk ,
      const int *ifirst ,
      const int *ilast ,
      const int *jfirst ,
      const int *jlast ,
      const int *kfirst ,
      const int *klast ,
      const double *dx );
   void compresconsca3d_(
      const double *xflux ,
      const double *yflux ,
      const double *zflux ,
      const int *fluxgi ,
      const int *fluxgj ,
      const int *fluxgk ,
      const double *rhs ,
      const int *rhsgi ,
      const int *rhsgj ,
      const int *rhsgk ,
      double *residual ,
      const int *residualgi ,
      const int *residualgj ,
      const int *residualgk ,
      const double &scalar_field ,
      const double *soln ,
      const int *solngi ,
      const int *solngj ,
      const int *solngk ,
      const int *ifirst ,
      const int *ilast ,
      const int *jfirst ,
      const int *jlast ,
      const int *kfirst ,
      const int *klast ,
      const double *dx );
   void ewingfixfluxvardc3d_(
      const double *xflux ,
      const double *yflux ,
      const double *zflux ,
      const int *fluxgi ,
      const int *fluxgj ,
      const int *fluxgk ,
      const double *xdiff_coef ,
      const double *ydiff_coef ,
      const double *zdiff_coef ,
      const int *dcgi ,
      const int *dcgj ,
      const int *dcgk ,
      const double *soln ,
      const int *solngi ,
      const int *solngj ,
      const int *solngk ,
      const int *ifirst ,
      const int *ilast ,
      const int *jfirst ,
      const int *jlast ,
      const int *kfirst ,
      const int *klast ,
      const int *location_index,
      const int *ratio_to_coarser,
      const int *blower,
      const int *bupper,
      const double *dx );
   void ewingfixfluxcondc3d_(
      const double *xflux ,
      const double *yflux ,
      const double *zflux ,
      const int *fluxgi ,
      const int *fluxgj ,
      const int *fluxgk ,
      const double &diff_coef ,
      const double *soln ,
      const int *solngi ,
      const int *solngj ,
      const int *solngk ,
      const int *ifirst ,
      const int *ilast ,
      const int *jfirst ,
      const int *jlast ,
      const int *kfirst ,
      const int *klast ,
      const int *location_index,
      const int *ratio_to_coarser,
      const int *blower,
      const int *bupper,
      const double *dx );

}


/*
********************************************************************
* Constructor.                                                     *
********************************************************************
*/
template<int DIM>  CellPoissonFACOps<DIM>::CellPoissonFACOps(
   const string &object_name ,
   tbox::Pointer<tbox::Database> database
) :
   d_object_name(object_name) ,
   d_hierarchy() ,
   d_ln_min(-1) ,
   d_ln_max(-1) ,
   d_cf_boundary() ,
   d_poisson_spec(object_name+"::Poisson specs") ,
   d_smoothing_choice("redblack") ,
   d_coarse_solver_choice(
#ifdef HAVE_HYPRE
                          "hypre"
#else
                          "redblack"
#endif


                          ) ,
   d_cf_discretization("Ewing") ,
   d_prolongation_method("CONSTANT_REFINE") ,
   d_coarse_solver_tolerance(1.e-8) ,
   d_coarse_solver_max_iterations(10) ,
   d_residual_tolerance_during_smoothing(-1.0) ,
   d_flux_id(-1) ,
#ifdef HAVE_HYPRE
   d_hypre_solver(object_name+"::hypre_solver",
                  database && database->isDatabase("hypre_solver") ?
                  database->getDatabase("hypre_solver") :
                  tbox::Pointer<tbox::Database>(NULL)),
#endif
   d_physical_bc_coef(NULL) ,
   d_context(hier::VariableDatabase<DIM>::getDatabase()
             ->getContext(object_name+"::PRIVATE_CONTEXT")) ,
   d_cell_scratch_id(-1),
   d_flux_scratch_id(-1),
   d_oflux_scratch_id(-1),
   d_prolongation_refine_operator() ,
   d_prolongation_refine_algorithm() ,
   d_prolongation_refine_schedules() ,
   d_urestriction_coarsen_operator() ,
   d_urestriction_coarsen_algorithm() ,
   d_urestriction_coarsen_schedules() ,
   d_rrestriction_coarsen_operator() ,
   d_rrestriction_coarsen_algorithm() ,
   d_rrestriction_coarsen_schedules() ,
   d_flux_coarsen_operator() ,
   d_flux_coarsen_algorithm() ,
   d_flux_coarsen_schedules() ,
   d_ghostfill_refine_operator() ,
   d_ghostfill_refine_algorithm() ,
   d_ghostfill_refine_schedules() ,
   d_ghostfill_nocoarse_refine_operator() ,
   d_ghostfill_nocoarse_refine_algorithm() ,
   d_ghostfill_nocoarse_refine_schedules() ,
   d_bc_helper( d_object_name+"::bc helper" ),
   d_enable_logging(false) ,
   d_preconditioner(NULL) ,
   d_hopscell() ,
   d_hopsside()
{

   t_restrict_solution = tbox::TimerManager::getManager()->
      getTimer("solv::CellPoissonFACOps::restrictSolution()");
   t_restrict_residual = tbox::TimerManager::getManager()->
      getTimer("solv::CellPoissonFACOps::restrictResidual()");
   t_prolong = tbox::TimerManager::getManager()->
      getTimer("solv::CellPoissonFACOps::prolongErrorAndCorrect()");
   t_smooth_error = tbox::TimerManager::getManager()->
      getTimer("solv::CellPoissonFACOps::smoothError()");
   t_solve_coarsest = tbox::TimerManager::getManager()->
      getTimer("solv::CellPoissonFACOps::solveCoarsestLevel()");
   t_compute_composite_residual = tbox::TimerManager::getManager()->
      getTimer("solv::CellPoissonFACOps::computeCompositeResidualOnLevel()");
   t_compute_residual_norm = tbox::TimerManager::getManager()->
      getTimer("solv::CellPoissonFACOps::computeResidualNorm()");

   if (DIM == 1) {
      TBOX_ERROR(d_object_name << ": 1D not implemented yet.\n");
   }

   if ( s_cell_scratch_var.isNull() ) {
#ifdef DEBUG_CHECK_ASSERTIONS
      assert( s_cell_scratch_var.isNull() );
      assert( s_cell_scratch_var.isNull() );
#endif
      s_cell_scratch_var = new pdat::CellVariable<DIM,double>
         ("CellPoissonFACOps::private_cell_scratch");
      s_flux_scratch_var = new pdat::SideVariable<DIM,double>
         ("CellPoissonFACOps::private_flux_scratch");
      s_oflux_scratch_var = new pdat::OutersideVariable<DIM,double>
         ("CellPoissonFACOps::private_oflux_scratch");
      tbox::ShutdownRegistry::registerShutdownRoutine(freeVariables, 254);
   }

   hier::VariableDatabase<DIM> *vdb = hier::VariableDatabase<DIM>::getDatabase();
   d_cell_scratch_id = vdb->
      registerVariableAndContext( s_cell_scratch_var,
                                  d_context ,
                                  hier::IntVector<DIM>(1) );
   d_flux_scratch_id = vdb->
      registerVariableAndContext( s_flux_scratch_var,
                                  d_context );
   d_oflux_scratch_id = vdb->
      registerVariableAndContext( s_oflux_scratch_var,
                                  d_context );


   /*
    * Some variables initialized by default are overriden by input.
    */
   if ( database ) {

      d_coarse_solver_choice =
         database->getStringWithDefault("coarse_solver_choice",
                                        d_coarse_solver_choice);
      d_coarse_solver_tolerance =
         database->getDoubleWithDefault("coarse_solver_tolerance",
                                        d_coarse_solver_tolerance);
      d_coarse_solver_max_iterations =
         database->getIntegerWithDefault("coarse_solver_max_iterations",
                                         d_coarse_solver_max_iterations);
      d_smoothing_choice =
         database->getStringWithDefault("smoothing_choice",
                                        d_smoothing_choice);

      d_cf_discretization =
         database->getStringWithDefault( "cf_discretization" ,
                                         d_cf_discretization );

      d_prolongation_method =
         database->getStringWithDefault( "prolongation_method" ,
                                         d_prolongation_method );


      d_enable_logging =
         database->getBoolWithDefault( "enable_logging" ,
                                       d_enable_logging );

   }


   /*
    * Check input validity and correctness.
    */
   checkInputPatchDataIndices();


   return;
}



template<int DIM>  CellPoissonFACOps<DIM>::~CellPoissonFACOps(void)
{
   return;
}




/*
************************************************************************
* FACOperatorStrategy<DIM> virtual initializeOperatorState function.  *
*                                                                      *
* Set internal variables to correspond to the solution passed in.      *
* Look up transfer operators.                                          *
************************************************************************
*/

template<int DIM> void CellPoissonFACOps<DIM>::initializeOperatorState (
   const SAMRAIVectorReal<DIM,double> &solution ,
   const SAMRAIVectorReal<DIM,double> &rhs )
{
   deallocateOperatorState();
   int ln;
   hier::VariableDatabase<DIM> *vdb = hier::VariableDatabase<DIM>::getDatabase();

   d_hierarchy = solution.getPatchHierarchy();
   d_ln_min = solution.getCoarsestLevelNumber();
   d_ln_max = solution.getFinestLevelNumber();
   d_hopscell = new math::HierarchyCellDataOpsReal<DIM,double>( d_hierarchy,
                                                            d_ln_min,
                                                            d_ln_max );
   d_hopsside = new math::HierarchySideDataOpsReal<DIM,double>( d_hierarchy,
                                                            d_ln_min,
                                                            d_ln_max );

#ifdef DEBUG_CHECK_ASSERTIONS

   int pn;

   if ( d_physical_bc_coef == NULL ) {
      /*
       * It's an error not to have bc object set.
       * Note that the bc object cannot be passed in through
       * the argument because the interface is inherited.
       */
      TBOX_ERROR(d_object_name << ": No physical bc object in\n"
                 << "CellPoissonFACOps<DIM>::initializeOperatorState\n"
                 << "You must use "
                 << "CellPoissonFACOps<DIM>::setPhysicalBcCoefObject\n"
                 << "to set one before calling initializeOperatorState\n");
   }


   if (  solution.getNumberOfComponents() != 1 ) {
      TBOX_WARNING(d_object_name
                   << ": Solution vector has multiple components.\n"
                   << "Solver is for component 0 only.\n");
   }
   if ( rhs.getNumberOfComponents() != 1 ) {
      TBOX_WARNING(d_object_name
                   << ": RHS vector has multiple components.\n"
                   << "Solver is for component 0 only.\n");
   }


   /*
    * Make sure that solution and rhs data
    *   are of correct type
    *   are allocated
    *   has sufficient ghost width
    */
   tbox::Pointer< hier::Variable<DIM> > var;
   {
      vdb->mapIndexToVariable( rhs.getComponentDescriptorIndex(0),
                               var );
      if ( !var ) {
         TBOX_ERROR(d_object_name << ": RHS component does not\n"
                    << "correspond to a variable.\n");
      }
      tbox::Pointer<pdat::CellVariable<DIM,double> > cell_var = var;
      if ( !cell_var ) {
         TBOX_ERROR(d_object_name
                    << ": RHS variable is not cell-centered double\n");
      }
   }
   {
      vdb->mapIndexToVariable( solution.getComponentDescriptorIndex(0),
                               var );
      if ( !var ) {
         TBOX_ERROR(d_object_name << ": Solution component does not\n"
                    << "correspond to a variable.\n");
      }
      tbox::Pointer<pdat::CellVariable<DIM,double> > cell_var = var;
      if ( !cell_var ) {
         TBOX_ERROR(d_object_name
                    << ": Solution variable is not cell-centered double\n");
      }
   }
   for ( ln=d_ln_min; ln<=d_ln_max; ++ln ) {
      tbox::Pointer<hier::PatchLevel<DIM> > level_ptr =
         d_hierarchy->getPatchLevel(ln);
      hier::PatchLevel<DIM> &level = *level_ptr;
      typename hier::PatchLevel<DIM>::Iterator pi;
      for ( pi.initialize(level); pi; pi++ ) {
         pn = pi();
         hier::Patch<DIM> &patch = *level.getPatch(pn);
         tbox::Pointer< hier::PatchData<DIM> > fd =
            patch.getPatchData( rhs.getComponentDescriptorIndex(0) );
         if ( fd ) {
            /*
              Some data checks can only be done if the data already exists.
            */
            tbox::Pointer<pdat::CellData<DIM,double> > cd = fd;
            if ( !cd ) {
               TBOX_ERROR(d_object_name
                          << ": RHS data is not cell-centered double\n");
            }
            if ( cd->getDepth() > 1 ) {
               TBOX_WARNING(d_object_name
                           << ": RHS data has multiple depths.\n"
                           << "Solver is for depth 0 only.\n");
            }
         }
         tbox::Pointer< hier::PatchData<DIM> > ud =
            patch.getPatchData( solution.getComponentDescriptorIndex(0) );
         if ( ud ) {
            /*
              Some data checks can only be done if the data already exists.
            */
            tbox::Pointer<pdat::CellData<DIM,double> > cd = ud;
            if ( !cd ) {
               TBOX_ERROR(d_object_name
                          << ": Solution data is not cell-centered double\n");
            }
            if ( cd->getDepth() > 1 ) {
               TBOX_WARNING(d_object_name
                           << ": Solution data has multiple depths.\n"
                           << "Solver is for depth 0 only.\n");
            }
            if ( cd->getGhostCellWidth() < hier::IntVector<DIM>(1) ) {
               TBOX_ERROR(d_object_name
                          << ": Solution data has insufficient ghost width\n");
            }
         }
      }
   }

   /*
    * Solution and rhs must have some similar properties.
    */
   if( rhs.getPatchHierarchy() != d_hierarchy
       || rhs.getCoarsestLevelNumber() != d_ln_min
       || rhs.getFinestLevelNumber() != d_ln_max ) {
      TBOX_ERROR(d_object_name << ": solution and rhs do not have\n"
                 << "the same set of patch levels.\n" );
   }

#endif

   /*
    * Initialize the coarse-fine boundary description for the
    * hierarchy.
    */
   d_cf_boundary.resizeArray( d_hierarchy->getNumberLevels() );

   hier::IntVector<DIM> max_gcw(1);
   for ( ln=d_ln_min; ln<=d_ln_max; ++ln ) {
      d_cf_boundary[ln].computeFromHierarchy(*d_hierarchy,
                                              ln,
                                              max_gcw);
   }
#if HAVE_HYPRE
   if ( d_coarse_solver_choice == "hypre" ) {
      d_hypre_solver.initializeSolverState( d_hierarchy, d_ln_min );
      /*
       * Share the boundary condition object with the hypre solver
       * to make sure that boundary condition settings are consistent
       * between the two objects.
       */
      d_hypre_solver.setPhysicalBcCoefObject( d_physical_bc_coef );
      d_hypre_solver.setMatrixCoefficients( d_poisson_spec );
   }
#endif

   /*
    * Get the transfer operators.
    * Flux coarsening is conservative.
    * Cell (solution, error, etc) coarsening is conservative.
    * Cell refinement from same level is constant refinement.
    * Cell refinement from coarser level is chosen by the
    *   choice of coarse-fine discretization, d_cf_discretization,
    *   which should be set to either "Ewing" or one of the
    *   acceptable strings for looking up the refine operator.
    */
   tbox::Pointer< geom::CartesianGridGeometry<DIM> > geometry
      = d_hierarchy->getGridGeometry();
   tbox::Pointer< hier::Variable<DIM> > variable;

   vdb->mapIndexToVariable( d_cell_scratch_id, variable );
   d_prolongation_refine_operator =
      geometry->lookupRefineOperator( variable,
                                      d_prolongation_method );

   vdb->mapIndexToVariable( d_cell_scratch_id, variable );
   d_urestriction_coarsen_operator =
   d_rrestriction_coarsen_operator =
      geometry->lookupCoarsenOperator(variable,
                                      "CONSERVATIVE_COARSEN");

   vdb->mapIndexToVariable( d_oflux_scratch_id, variable );
   d_flux_coarsen_operator =
      geometry->lookupCoarsenOperator(variable,
                                      "CONSERVATIVE_COARSEN");

   vdb->mapIndexToVariable( d_cell_scratch_id, variable );
   d_ghostfill_refine_operator =
      geometry->lookupRefineOperator(variable,
                                     d_cf_discretization == "Ewing" ?
                                     "CONSTANT_REFINE" : d_cf_discretization);

   vdb->mapIndexToVariable( d_cell_scratch_id, variable );
   d_ghostfill_nocoarse_refine_operator =
      geometry->lookupRefineOperator(variable,
                                     "CONSTANT_REFINE");

#ifdef DEBUG_CHECK_ASSERTIONS
   if ( !d_prolongation_refine_operator ) {
      TBOX_ERROR(d_object_name
                 << ": Cannot find prolongation refine operator");
   }
   if ( !d_urestriction_coarsen_operator ) {
      TBOX_ERROR(d_object_name
                 << ": Cannot find restriction coarsening operator");
   }
   if ( !d_rrestriction_coarsen_operator ) {
      TBOX_ERROR(d_object_name
                 << ": Cannot find restriction coarsening operator");
   }
   if ( !d_flux_coarsen_operator ) {
      TBOX_ERROR(d_object_name
                 << ": Cannot find flux coarsening operator");
   }
   if ( !d_ghostfill_refine_operator ) {
      TBOX_ERROR(d_object_name
                 << ": Cannot find ghost filling refinement operator");
   }
   if ( !d_ghostfill_nocoarse_refine_operator ) {
      TBOX_ERROR(d_object_name
                 << ": Cannot find ghost filling refinement operator");
   }
#endif

   for ( ln=d_ln_min+1; ln<=d_ln_max; ++ln ) {
      d_hierarchy->getPatchLevel(ln)->
         allocatePatchData(d_oflux_scratch_id);
   }

   /*
     Make space for saving communication schedules.
     There is no need to delete the old schedules first
     because we have deallocated the solver state above.
   */
   d_prolongation_refine_schedules.resizeArray(d_ln_max+1);
   d_ghostfill_refine_schedules.resizeArray(d_ln_max+1);
   d_ghostfill_nocoarse_refine_schedules.resizeArray(d_ln_max+1);
   d_urestriction_coarsen_schedules.resizeArray(d_ln_max+1);
   d_rrestriction_coarsen_schedules.resizeArray(d_ln_max+1);
   d_flux_coarsen_schedules.resizeArray(d_ln_max+1);

   d_prolongation_refine_algorithm = new xfer::RefineAlgorithm<DIM>;
   d_urestriction_coarsen_algorithm = new xfer::CoarsenAlgorithm<DIM>;
   d_rrestriction_coarsen_algorithm = new xfer::CoarsenAlgorithm<DIM>;
   d_flux_coarsen_algorithm = new xfer::CoarsenAlgorithm<DIM>;
   d_ghostfill_refine_algorithm = new xfer::RefineAlgorithm<DIM>;
   d_ghostfill_nocoarse_refine_algorithm = new xfer::RefineAlgorithm<DIM>;

   d_prolongation_refine_algorithm->
      registerRefine( d_cell_scratch_id ,
                      solution.getComponentDescriptorIndex(0) ,
                      d_cell_scratch_id ,
                      d_prolongation_refine_operator );
   d_urestriction_coarsen_algorithm->
      registerCoarsen( solution.getComponentDescriptorIndex(0) ,
                       solution.getComponentDescriptorIndex(0) ,
                       d_urestriction_coarsen_operator );
   d_rrestriction_coarsen_algorithm->
      registerCoarsen( rhs.getComponentDescriptorIndex(0) ,
                       rhs.getComponentDescriptorIndex(0) ,
                       d_rrestriction_coarsen_operator );
   d_ghostfill_refine_algorithm->
      registerRefine( solution.getComponentDescriptorIndex(0) ,
                      solution.getComponentDescriptorIndex(0) ,
                      solution.getComponentDescriptorIndex(0) ,
                      d_ghostfill_refine_operator );
   d_flux_coarsen_algorithm->
      registerCoarsen( ( ( d_flux_id != -1 ) ? d_flux_id : d_flux_scratch_id ) ,
                       d_oflux_scratch_id ,
                       d_flux_coarsen_operator );
   d_ghostfill_nocoarse_refine_algorithm->
      registerRefine( solution.getComponentDescriptorIndex(0) ,
                      solution.getComponentDescriptorIndex(0) ,
                      solution.getComponentDescriptorIndex(0) ,
                      d_ghostfill_nocoarse_refine_operator );

   for ( int dest_ln=d_ln_min+1; dest_ln<=d_ln_max; ++dest_ln ) {
      d_prolongation_refine_schedules[dest_ln] =
         d_prolongation_refine_algorithm->
         createSchedule( d_hierarchy->getPatchLevel(dest_ln) ,
                         tbox::Pointer<hier::PatchLevel<DIM> >() ,
                         dest_ln-1 ,
                         d_hierarchy ,
                         &d_bc_helper );
      if ( ! d_prolongation_refine_schedules[dest_ln] ) {
         TBOX_ERROR(d_object_name
                    << ": Cannot create a refine schedule for prolongation!\n");
      }
      d_ghostfill_refine_schedules[dest_ln] =
         d_ghostfill_refine_algorithm->
         createSchedule( d_hierarchy->getPatchLevel(dest_ln) ,
                         dest_ln-1 ,
                         d_hierarchy ,
                         &d_bc_helper );
      if ( ! d_ghostfill_refine_schedules[dest_ln] ) {
         TBOX_ERROR(d_object_name
                    << ": Cannot create a refine schedule for ghost filling!\n");
      }
      d_ghostfill_nocoarse_refine_schedules[dest_ln] =
         d_ghostfill_nocoarse_refine_algorithm->
         createSchedule( d_hierarchy->getPatchLevel(dest_ln) ,
                         &d_bc_helper );
      if ( ! d_ghostfill_nocoarse_refine_schedules[dest_ln] ) {
         TBOX_ERROR(d_object_name
                    << ": Cannot create a refine schedule for ghost filling on bottom level!\n");
      }
   }
   for ( int dest_ln=d_ln_min; dest_ln<d_ln_max; ++dest_ln ) {
      d_urestriction_coarsen_schedules[dest_ln] =
         d_urestriction_coarsen_algorithm->
         createSchedule(d_hierarchy->getPatchLevel(dest_ln),
                        d_hierarchy->getPatchLevel(dest_ln+1));
      if ( ! d_urestriction_coarsen_schedules[dest_ln] ) {
         TBOX_ERROR(d_object_name
                    << ": Cannot create a coarsen schedule for U restriction!\n");
      }
      d_rrestriction_coarsen_schedules[dest_ln] =
         d_rrestriction_coarsen_algorithm->
         createSchedule(d_hierarchy->getPatchLevel(dest_ln),
                        d_hierarchy->getPatchLevel(dest_ln+1));
      if ( ! d_rrestriction_coarsen_schedules[dest_ln] ) {
         TBOX_ERROR(d_object_name
                    << ": Cannot create a coarsen schedule for R restriction!\n");
      }
      d_flux_coarsen_schedules[dest_ln] =
         d_flux_coarsen_algorithm->
         createSchedule(d_hierarchy->getPatchLevel(dest_ln),
                        d_hierarchy->getPatchLevel(dest_ln+1));
      if ( ! d_flux_coarsen_schedules[dest_ln] ) {
         TBOX_ERROR(d_object_name
                    << ": Cannot create a coarsen schedule for flux transfer!\n");
      }
   }
   d_ghostfill_nocoarse_refine_schedules[d_ln_min] =
      d_ghostfill_nocoarse_refine_algorithm->
      createSchedule( d_hierarchy->getPatchLevel(d_ln_min) ,
                      &d_bc_helper );
   if ( ! d_ghostfill_nocoarse_refine_schedules[d_ln_min] ) {
      TBOX_ERROR(d_object_name
                 << ": Cannot create a refine schedule for ghost filling on bottom level!\n");
   }

   return;
}



/*
********************************************************************
* FACOperatorStrategy<DIM> virtual deallocateOperatorState        *
* function.  Deallocate internal hierarchy-dependent data.         *
* State is allocated iff hierarchy is set.                         *
********************************************************************
*/

template<int DIM> void CellPoissonFACOps<DIM>::deallocateOperatorState()
{
   if ( d_hierarchy ) {
      int ln;
      for ( ln=d_ln_min+1; ln<=d_ln_max; ++ln ) {
         d_hierarchy->getPatchLevel(ln)->
            deallocatePatchData(d_oflux_scratch_id);
      }
      d_cf_boundary.resizeArray(0);
#if HAVE_HYPRE
      d_hypre_solver.deallocateSolverState();
#endif
      d_hierarchy.setNull();
      d_ln_min = -1;
      d_ln_max = -1;

      d_prolongation_refine_algorithm.setNull();
      d_prolongation_refine_schedules.setNull();

      d_urestriction_coarsen_algorithm.setNull();
      d_urestriction_coarsen_schedules.setNull();

      d_rrestriction_coarsen_algorithm.setNull();
      d_rrestriction_coarsen_schedules.setNull();

      d_flux_coarsen_algorithm.setNull();
      d_flux_coarsen_schedules.setNull();

      d_ghostfill_refine_algorithm.setNull();
      d_ghostfill_refine_schedules.setNull();

      d_ghostfill_nocoarse_refine_algorithm.setNull();
      d_ghostfill_nocoarse_refine_schedules.setNull();

   }
   return;
}




/*
********************************************************************
* FACOperatorStrategy<DIM> virtual postprocessOneCycle function.  *
********************************************************************
*/

template<int DIM> void CellPoissonFACOps<DIM>::postprocessOneCycle(
   int fac_cycle_num ,
   const SAMRAIVectorReal<DIM,double> &current_soln ,
   const SAMRAIVectorReal<DIM,double> &residual )
{
   if ( d_enable_logging ) {
      if ( d_preconditioner ) {
         /*
          * Output convergence progress.  This is probably only appropriate
          * if the solver is NOT being used as a preconditioner.
          */
         double avg_factor, final_factor;
         d_preconditioner->getConvergenceFactors( avg_factor, final_factor );
         tbox::plog
            << "iter=" << setw(4) << fac_cycle_num
            << " resid=" << d_preconditioner->getResidualNorm()
            << " net conv=" << d_preconditioner->getNetConvergenceFactor()
            << " final conv=" << d_preconditioner->getNetConvergenceFactor()
            << " avg conv=" << d_preconditioner->getAvgConvergenceFactor()
            << endl;
      }
   }
   return;
}



/*
********************************************************************
* FACOperatorStrategy<DIM> virtual restrictSolution function.     *
* After restricting solution, update ghost cells of the affected   *
* level.                                                           *
********************************************************************
*/

template<int DIM> void CellPoissonFACOps<DIM>::restrictSolution(
   const SAMRAIVectorReal<DIM,double> &s ,
   SAMRAIVectorReal<DIM,double> &d ,
   int dest_ln ) {

   t_restrict_solution->start();

   xeqScheduleURestriction(d.getComponentDescriptorIndex(0),
                           s.getComponentDescriptorIndex(0),
                           dest_ln);

   d_bc_helper.setHomogeneousBc(false);
   d_bc_helper.setTargetDataId(d.getComponentDescriptorIndex(0));

   if ( dest_ln == d_ln_min ) {
      xeqScheduleGhostFillNoCoarse(d.getComponentDescriptorIndex(0),
                                   dest_ln);
   }
   else {
      xeqScheduleGhostFill(d.getComponentDescriptorIndex(0),
                           dest_ln);
   }

   t_restrict_solution->stop();
   return;
}




/*
********************************************************************
* FACOperatorStrategy<DIM> virtual restrictresidual function.     *
********************************************************************
*/

template<int DIM> void CellPoissonFACOps<DIM>::restrictResidual(
   const SAMRAIVectorReal<DIM,double> &s ,
   SAMRAIVectorReal<DIM,double> &d ,
   int dest_ln ) {

   t_restrict_residual->start();

   xeqScheduleRRestriction(d.getComponentDescriptorIndex(0),
                           s.getComponentDescriptorIndex(0),
                           dest_ln);

   t_restrict_residual->stop();
   return;
}




/*
***********************************************************************
* FACOperatorStrategy<DIM> virtual prolongErrorAndCorrect function.  *
* After the prolongation, we set the physical boundary condition      *
* for the correction, which is zero.  Other ghost cell values,        *
* which are preset to zero, need not be set.                          *
***********************************************************************
*/

template<int DIM> void CellPoissonFACOps<DIM>::prolongErrorAndCorrect(
   const SAMRAIVectorReal<DIM,double> &s ,
   SAMRAIVectorReal<DIM,double> &d ,
   int dest_ln )  {

   t_prolong->start();

#ifdef DEBUG_CHECK_ASSERTIONS
   if( s.getPatchHierarchy() != d_hierarchy
       || d.getPatchHierarchy() != d_hierarchy ) {
      TBOX_ERROR(d_object_name << ": Vector hierarchy does not match\n"
                 "internal state hierarchy.");
   }
#endif

   tbox::Pointer< hier::PatchLevel<DIM> > coarse_level
      = d_hierarchy->getPatchLevel(dest_ln-1);
   tbox::Pointer< hier::PatchLevel<DIM> > fine_level
      = d_hierarchy->getPatchLevel(dest_ln);

   /*
    * Data is prolonged into the scratch space corresponding
    * to index d_cell_scratch_id and allocated here.
    */
   fine_level->allocatePatchData(d_cell_scratch_id);

   /*
    * Refine solution into scratch space to fill the fine level
    * interior in the scratch space, then use that refined data
    * to correct the fine level error.
    */
   d_bc_helper.setTargetDataId(d_cell_scratch_id);
   d_bc_helper.setHomogeneousBc(true);
   const int src_index = s.getComponentDescriptorIndex(0);
   xeqScheduleProlongation(d_cell_scratch_id,
                           src_index,
                           d_cell_scratch_id,
                           dest_ln);

   /*
    * Add the refined error in the scratch space
    * to the error currently residing in the destination level.
    */
   math::HierarchyCellDataOpsReal<DIM,double>
      hierarchy_math_ops( d_hierarchy, dest_ln, dest_ln );
   const int dst_index = d.getComponentDescriptorIndex(0);
   hierarchy_math_ops.add( dst_index, dst_index, d_cell_scratch_id );

   fine_level->deallocatePatchData(d_cell_scratch_id);

   t_prolong->stop();

   return;
}



/*
********************************************************************
********************************************************************
*/

template<int DIM> void CellPoissonFACOps<DIM>::smoothError(
   SAMRAIVectorReal<DIM,double> &data ,
   const SAMRAIVectorReal<DIM,double> &residual ,
   int ln ,
   int num_sweeps )
{

   t_smooth_error->start();

   checkInputPatchDataIndices();
   if ( d_smoothing_choice == "redblack" ) {
      smoothErrorByRedBlack ( data ,
                              residual ,
                              ln ,
                              num_sweeps ,
                              d_residual_tolerance_during_smoothing );
   }
   else {
      TBOX_ERROR(d_object_name << ": Bad smoothing choice '"
                 << d_smoothing_choice
                 << "' in CellPoissonFACOps<DIM>.");
   }

   t_smooth_error->stop();
   return;
}



/*
********************************************************************
* Workhorse function to smooth error using red-black               *
* Gauss-Seidel iterations.                                         *
********************************************************************
*/

template<int DIM> void CellPoissonFACOps<DIM>::smoothErrorByRedBlack(
   SAMRAIVectorReal<DIM,double> &data ,
   const SAMRAIVectorReal<DIM,double> &residual ,
   int ln ,
   int num_sweeps ,
   double residual_tolerance )
{

   checkInputPatchDataIndices();

#ifdef DEBUG_CHECK_ASSERTIONS
   if( data.getPatchHierarchy() != d_hierarchy
       || residual.getPatchHierarchy() != d_hierarchy ) {
      TBOX_ERROR(d_object_name << ": Vector hierarchy does not match\n"
                 "internal hierarchy.");
   }
#endif
   tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);

   const int data_id = data.getComponentDescriptorIndex(0);

   const int flux_id = ( d_flux_id != -1 ) ? d_flux_id : d_flux_scratch_id;

   d_bc_helper.setTargetDataId(data_id);
   d_bc_helper.setHomogeneousBc(true);
   xeqScheduleGhostFillNoCoarse(data_id, ln);

   if ( ln > d_ln_min ) {
      /*
       * Perform a one-time transfer of data from coarser level,
       * to fill ghost boundaries that will not change through
       * the smoothing loop.
       */
      xeqScheduleGhostFill(data_id, ln);
   }

   /*
    * Smooth the number of sweeps specified or until
    * the convergence is satisfactory.
    */
   int isweep;
   double red_maxres, blk_maxres, maxres=0;
   red_maxres = blk_maxres = residual_tolerance+1;
   /*
    * Instead of checking residual convergence globally,
    * we check the not_converged flag.  This avoids possible
    * round-off errors affecting different processes differently,
    * leading to disagreement on whether to continue smoothing.
    */
   int not_converged = 1;
   for ( isweep=0; isweep<num_sweeps && not_converged; ++isweep ) {
      red_maxres = blk_maxres = 0;

      // Red sweep.
      xeqScheduleGhostFillNoCoarse(data_id, ln);
      int pn;
      typename hier::PatchLevel<DIM>::Iterator pi;
      for ( pi.initialize(*level); pi; pi++ ) {
         pn = *pi;
         tbox::Pointer< hier::Patch<DIM> > patch = level->getPatch(pn);

         bool deallocate_flux_data_when_done=false;
         if ( flux_id == d_flux_scratch_id ) {
            /*
             * Using internal temporary storage for flux.
             * For each patch, make sure the internal
             * side-centered data is allocated and note
             * whether that data should be deallocated when done.
             */
            if ( !patch->checkAllocated(flux_id) ) {
               patch->allocatePatchData(flux_id);
               deallocate_flux_data_when_done = true;
            }
         }

         tbox::Pointer<pdat::CellData<DIM,double> >
            scalar_field_data = d_poisson_spec.cIsVariable() ?
            patch->getPatchData( d_poisson_spec.getCPatchDataId() ) :
            tbox::Pointer< hier::PatchData<DIM> >(NULL);
         tbox::Pointer<pdat::CellData<DIM,double> >
            err_data = data.getComponentPatchData ( 0 , *patch );
         tbox::Pointer<pdat::CellData<DIM,double> >
            residual_data = residual.getComponentPatchData ( 0 , *patch );
         tbox::Pointer<pdat::SideData<DIM,double> >
            flux_data = patch->getPatchData( flux_id );

         computeFluxOnPatch(
                            *patch ,
                            level->getRatioToCoarserLevel() ,
                            *err_data ,
                            *flux_data );

         redOrBlackSmoothingOnPatch( *patch ,
                                     *flux_data ,
                                     *residual_data ,
                                     *err_data ,
                                     'r' ,
                                     &red_maxres );

         if ( deallocate_flux_data_when_done ) {
            patch->deallocatePatchData(flux_id);
         }
      }        // End patch number pn
      xeqScheduleGhostFillNoCoarse(data_id, ln);

      // Black sweep.
      for ( pi.initialize(*level); pi; pi++ ) {
         pn = *pi;
         tbox::Pointer< hier::Patch<DIM> > patch = level->getPatch(pn);

         bool deallocate_flux_data_when_done=false;
         if ( flux_id == d_flux_scratch_id ) {
            /*
             * Using internal temporary storage for flux.
             * For each patch, make sure the internal
             * side-centered data is allocated and note
             * whether that data should be deallocated when done.
             */
            if ( !patch->checkAllocated(flux_id) ) {
               patch->allocatePatchData(flux_id);
               deallocate_flux_data_when_done = true;
            }
         }

         tbox::Pointer<pdat::CellData<DIM,double> >
            scalar_field_data = d_poisson_spec.cIsVariable() ?
            patch->getPatchData( d_poisson_spec.getCPatchDataId() ) :
            tbox::Pointer< hier::PatchData<DIM> >(NULL);
         tbox::Pointer<pdat::CellData<DIM,double> >
            err_data = data.getComponentPatchData ( 0 , *patch );
         tbox::Pointer<pdat::CellData<DIM,double> >
            residual_data = residual.getComponentPatchData ( 0 , *patch );
         tbox::Pointer<pdat::SideData<DIM,double> >
            flux_data = patch->getPatchData( flux_id );

         computeFluxOnPatch(
                            *patch ,
                            level->getRatioToCoarserLevel() ,
                            *err_data ,
                            *flux_data );

         redOrBlackSmoothingOnPatch( *patch ,
                                     *flux_data ,
                                     *residual_data ,
                                     *err_data ,
                                     'b' ,
                                     &red_maxres );

         if ( deallocate_flux_data_when_done ) {
            patch->deallocatePatchData(flux_id);
         }
      }        // End patch number pn
      xeqScheduleGhostFillNoCoarse(data_id, ln);
      if ( residual_tolerance >= 0.0 ) {
         /*
           Check for early end of sweeps due to convergence
           only if it is numerically possible (user gave a
           non negative value for residual tolerance).
          */
         maxres = tbox::Utilities::dmax(red_maxres,blk_maxres);
         not_converged = maxres>residual_tolerance;
         not_converged = tbox::MPI::maxReduction( not_converged );
      }
   }        // End sweep number isweep
   if ( d_enable_logging ) tbox::plog
      << d_object_name << " RBGS smoothing maxres = " << maxres << "\n"
      << "  after " << isweep << " sweeps.\n";

   return;
}




/*
********************************************************************
* Fix flux on coarse-fine boundaries computed from a               *
* constant-refine interpolation of coarse level data.              *
********************************************************************
*/

template<int DIM> void CellPoissonFACOps<DIM>::ewingFixFlux (
   const hier::Patch<DIM> &patch ,
   const pdat::CellData<DIM,double> &soln_data ,
   pdat::SideData<DIM,double> &flux_data ,
   const hier::IntVector<DIM> &ratio_to_coarser ) const
{
   const int patch_ln = patch.getPatchLevelNumber();
   const int pn = patch.getPatchNumber();
   tbox::Pointer< geom::CartesianPatchGeometry<DIM> > patch_geom
      = patch.getPatchGeometry();
   const double *dx = patch_geom->getDx();
   const hier::Box<DIM> &patch_box( patch.getBox() );
   const hier::Index<DIM> &plower = patch_box.lower();
   const hier::Index<DIM> &pupper = patch_box.upper();

   const tbox::Array< hier::BoundaryBox<DIM> > &bboxes =
      d_cf_boundary[patch_ln].getBoundaries(pn,1);
   int bn, nboxes=bboxes.getSize();

   if ( d_poisson_spec.dIsVariable() ) {

      tbox::Pointer<pdat::SideData<DIM,double> > diffcoef_data;
      diffcoef_data = patch.getPatchData( d_poisson_spec.getDPatchDataId() );

      for ( bn=0; bn<nboxes; ++bn ) {
         const hier::BoundaryBox<DIM> &boundary_box=bboxes[bn];
#ifdef DEBUG_CHECK_ASSERTIONS
         assert( boundary_box.getBoundaryType() == 1 );
#endif
         const hier::Box<DIM> &bdry_box = boundary_box.getBox();
         const hier::Index<DIM> &blower = bdry_box.lower();
         const hier::Index<DIM> &bupper = bdry_box.upper();
         const int location_index = boundary_box.getLocationIndex();
	 if (DIM == 2) {
	    ewingfixfluxvardc2d_(
	       flux_data.getPointer(0) , flux_data.getPointer(1) ,
	       &flux_data.getGhostCellWidth()[0],
	       &flux_data.getGhostCellWidth()[1] ,
	       diffcoef_data->getPointer(0) , diffcoef_data->getPointer(1) ,
	       &diffcoef_data->getGhostCellWidth()[0],
	       &diffcoef_data->getGhostCellWidth()[1] ,
	       soln_data.getPointer() ,
	       &soln_data.getGhostCellWidth()[0],
	       &soln_data.getGhostCellWidth()[1] ,
	       &plower[0], &pupper[0], &plower[1], &pupper[1] ,
	       &location_index,
	       ratio_to_coarser,
	       blower, bupper,
	       dx );
	 } else if (DIM == 3) {
	    ewingfixfluxvardc3d_(
	       flux_data.getPointer(0) ,
	       flux_data.getPointer(1) ,
	       flux_data.getPointer(2) ,
	       &flux_data.getGhostCellWidth()[0],
	       &flux_data.getGhostCellWidth()[1] ,
	       &flux_data.getGhostCellWidth()[2] ,
	       diffcoef_data->getPointer(0) ,
	       diffcoef_data->getPointer(1) ,
	       diffcoef_data->getPointer(2) ,
	       &diffcoef_data->getGhostCellWidth()[0],
	       &diffcoef_data->getGhostCellWidth()[1] ,
	       &diffcoef_data->getGhostCellWidth()[2] ,
	       soln_data.getPointer() ,
	       &soln_data.getGhostCellWidth()[0],
	       &soln_data.getGhostCellWidth()[1] ,
	       &soln_data.getGhostCellWidth()[2] ,
	       &plower[0], &pupper[0],
	       &plower[1], &pupper[1] ,
	       &plower[2], &pupper[2] ,
	       &location_index,
	       ratio_to_coarser,
	       blower, bupper,
	       dx );
	 }
      }
   }
   else {

      const double diffcoef_constant = d_poisson_spec.getDConstant();

      for ( bn=0; bn<nboxes; ++bn ) {
         const hier::BoundaryBox<DIM> &boundary_box=bboxes[bn];
#ifdef DEBUG_CHECK_ASSERTIONS
         assert( boundary_box.getBoundaryType() == 1 );
#endif
         const hier::Box<DIM> &bdry_box = boundary_box.getBox();
         const hier::Index<DIM> &blower = bdry_box.lower();
         const hier::Index<DIM> &bupper = bdry_box.upper();
         const int location_index = boundary_box.getLocationIndex();
	 if (DIM == 2) {
	    ewingfixfluxcondc2d_(
	       flux_data.getPointer(0) , flux_data.getPointer(1) ,
	       &flux_data.getGhostCellWidth()[0],
	       &flux_data.getGhostCellWidth()[1] ,
	       diffcoef_constant ,
	       soln_data.getPointer() ,
	       &soln_data.getGhostCellWidth()[0],
	       &soln_data.getGhostCellWidth()[1] ,
	       &plower[0], &pupper[0],
	       &plower[1], &pupper[1] ,
	       &location_index,
	       ratio_to_coarser,
	       blower, bupper,
	       dx );
	 } else if (DIM == 3) {
	    ewingfixfluxcondc3d_(
	       flux_data.getPointer(0) ,
	       flux_data.getPointer(1) ,
	       flux_data.getPointer(2) ,
	       &flux_data.getGhostCellWidth()[0],
	       &flux_data.getGhostCellWidth()[1] ,
	       &flux_data.getGhostCellWidth()[2] ,
	       diffcoef_constant ,
	       soln_data.getPointer() ,
	       &soln_data.getGhostCellWidth()[0],
	       &soln_data.getGhostCellWidth()[1] ,
	       &soln_data.getGhostCellWidth()[2] ,
	       &plower[0], &pupper[0],
	       &plower[1], &pupper[1] ,
	       &plower[2], &pupper[2] ,
	       &location_index,
	       ratio_to_coarser,
	       blower, bupper,
	       dx );
	 }
      }
   }

   return;
}




/*
********************************************************************
* FACOperatorStrategy<DIM> virtual solveCoarsestLevel             *
* function                                                         *
********************************************************************
*/

template<int DIM> void CellPoissonFACOps<DIM>::solveCoarsestLevel(
   SAMRAIVectorReal<DIM,double> &data ,
   const SAMRAIVectorReal<DIM,double> &residual ,
   int coarsest_ln )  {

   t_solve_coarsest->start();

   checkInputPatchDataIndices();

   if ( d_coarse_solver_choice == "jacobi" ) {
      d_residual_tolerance_during_smoothing = d_coarse_solver_tolerance;
      smoothError( data ,
                   residual ,
                   coarsest_ln ,
                   d_coarse_solver_max_iterations );
      d_residual_tolerance_during_smoothing = -1.0;
   }
   else if ( d_coarse_solver_choice == "redblack" ) {
      d_residual_tolerance_during_smoothing = d_coarse_solver_tolerance;
      smoothError( data ,
                   residual ,
                   coarsest_ln ,
                   d_coarse_solver_max_iterations );
      d_residual_tolerance_during_smoothing = -1.0;
   }
   else if ( d_coarse_solver_choice == "hypre" ) {
#ifndef HAVE_HYPRE
      TBOX_ERROR(d_object_name << ": Coarse level solver choice '"
                 << d_coarse_solver_choice
                 << "' unavailable in "
                 << "scapCellPoissonOps::solveCoarsestLevel.");
#else
      solveCoarsestLevel_HYPRE( data, residual, coarsest_ln );
#endif
   }
   else {
      TBOX_ERROR(d_object_name << ": Bad coarse level solver choice '"
                 << d_coarse_solver_choice
                 << "' in scapCellPoissonOps::solveCoarsestLevel.");
   }

   xeqScheduleGhostFillNoCoarse(data.getComponentDescriptorIndex(0),
                                coarsest_ln);

   t_solve_coarsest->stop();

   return;
}




/*
********************************************************************
* Solve coarsest level using Hypre                                 *
* We only solve for the error, so we always use homogeneous bc.    *
********************************************************************
*/

template<int DIM> void CellPoissonFACOps<DIM>::solveCoarsestLevel_HYPRE(
   SAMRAIVectorReal<DIM,double> &data ,
   const SAMRAIVectorReal<DIM,double> &residual ,
   int coarsest_ln )  {

#ifndef HAVE_HYPRE
      TBOX_ERROR(d_object_name << ": Coarse level solver choice '"
                 << d_coarse_solver_choice
                 << "' unavailable in "
                 << "CellPoissonFACOps::solveCoarsestLevel.");
#else

   checkInputPatchDataIndices();
   d_hypre_solver.setStoppingCriteria( d_coarse_solver_max_iterations,
                                       d_coarse_solver_tolerance );
   const int solver_ret =
      d_hypre_solver.solveSystem(
         data.getComponentDescriptorIndex(0),
         residual.getComponentDescriptorIndex(0) ,
         true);
   /*
    * Present data on the solve.
    * The Hypre solver returns 0 if converged.
    */
   if( d_enable_logging ) tbox::plog
      << d_object_name << " Hypre solve " << (solver_ret?"":"NOT ")
      << "converged\n"
      << "\titerations: " << d_hypre_solver.getNumberIterations() << "\n"
      << "\tresidual: " << d_hypre_solver.getRelativeResidualNorm() << "\n";

#endif

   return;

}



/*
********************************************************************
* FACOperatorStrategy<DIM> virtual                                *
* computeCompositeResidualOnLevel function                         *
********************************************************************
*/

template<int DIM> void CellPoissonFACOps<DIM>::computeCompositeResidualOnLevel(
   SAMRAIVectorReal<DIM,double> &residual ,
   const SAMRAIVectorReal<DIM,double> &solution ,
   const SAMRAIVectorReal<DIM,double> &rhs ,
   int ln ,
   bool error_equation_indicator )  {

   t_compute_composite_residual->start();

   checkInputPatchDataIndices();
#ifdef DEBUG_CHECK_ASSERTIONS
   if( residual.getPatchHierarchy() != d_hierarchy
       || solution.getPatchHierarchy() != d_hierarchy
       || rhs.getPatchHierarchy() != d_hierarchy ) {
      TBOX_ERROR(d_object_name << ": Vector hierarchy does not match\n"
                 "internal hierarchy.");
   }
#endif
   tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);

   /*
    * Set up the bc helper so that when we use a refine schedule
    * to fill ghosts, the correct data is operated on.
    */
   const int soln_id  = solution.getComponentDescriptorIndex(0);
   d_bc_helper.setTargetDataId(soln_id);
   d_bc_helper.setHomogeneousBc(error_equation_indicator);

   const int flux_id = ( d_flux_id != -1 ) ? d_flux_id : d_flux_scratch_id;

   /*
    * Assumptions:
    * 1. Data does not yet exist in ghost boundaries.
    * 2. Residual data on next finer grid (if any)
    *    has been computed already.
    * 3. Flux data from next finer grid (if any) has
    *    been computed but has not been coarsened to
    *    this level.
    *
    * Steps:
    * S1. Fill solution ghost data by refinement
    *     or setting physical boundary conditions.
    *     This also brings in information from coarser
    *     to form the composite grid flux.
    * S2. Compute flux on ln.
    * S3. If next finer is available,
    *     Coarsen flux data on next finer level,
    *     overwriting flux computed from coarse data.
    * S4. Compute residual data from flux.
    */

   /* S1. Fill solution ghost data. */
   {
      tbox::Pointer< xfer::RefineSchedule<DIM> > ln_refine_schedule;
      if ( ln > d_ln_min ) {
         /* Fill from current, next coarser level and physical boundary */
         xeqScheduleGhostFill(soln_id, ln);
      }
      else {
         /* Fill from current and physical boundary */
         xeqScheduleGhostFillNoCoarse(soln_id, ln);
      }
   }

   /*
    * For the whole level, make sure the internal
    * side-centered data is allocated and note
    * whether that data should be deallocated when done.
    * We do this for the whole level because the data
    * undergoes transfer operations which require the
    * whole level data.
    */
   bool deallocate_flux_data_when_done=false;
   if ( flux_id == d_flux_scratch_id ) {
      if ( !level->checkAllocated(flux_id) ) {
         level->allocatePatchData(flux_id);
         deallocate_flux_data_when_done = true;
      }
   }

   /*
    * S2. Compute flux on patches in level.
    */
   int pn;
   typename hier::PatchLevel<DIM>::Iterator pi;
   for ( pi.initialize(*level); pi; pi++ ) {
      pn = *pi;
      tbox::Pointer< hier::Patch<DIM> > patch = level->getPatch(pn);

      tbox::Pointer<pdat::CellData<DIM,double> >
         soln_data = solution.getComponentPatchData ( 0 , *patch );
      tbox::Pointer<pdat::CellData<DIM,double> >
         rhs_data = rhs.getComponentPatchData ( 0 , *patch );
      tbox::Pointer<pdat::CellData<DIM,double> >
         residual_data = residual.getComponentPatchData( 0 , *patch );
      tbox::Pointer<pdat::SideData<DIM,double> >
         flux_data = patch->getPatchData( flux_id );
      computeFluxOnPatch(
                         *patch ,
                         level->getRatioToCoarserLevel() ,
                         *soln_data ,
                         *flux_data );

   }

   /*
    * S3. Coarsen oflux data from next finer level so that
    * the computed flux becomes the composite grid flux.
    */
   if ( ln < d_ln_max ) {
      xeqScheduleFluxCoarsen(flux_id, d_oflux_scratch_id, ln);
   }

   /*
    * S4. Compute residual on patches in level.
    */
   for ( pi.initialize(*level); pi; pi++ ) {
      pn = *pi;
      tbox::Pointer< hier::Patch<DIM> > patch = level->getPatch(pn);
      tbox::Pointer<pdat::CellData<DIM,double> >
         soln_data = solution.getComponentPatchData ( 0 , *patch );
      tbox::Pointer<pdat::CellData<DIM,double> >
         rhs_data = rhs.getComponentPatchData ( 0 , *patch );
      tbox::Pointer<pdat::CellData<DIM,double> >
         residual_data = residual.getComponentPatchData( 0 , *patch );
      tbox::Pointer<pdat::SideData<DIM,double> >
         flux_data = patch->getPatchData( flux_id );
      computeResidualOnPatch( *patch ,
                              *flux_data ,
                              *soln_data ,
                              *rhs_data ,
                              *residual_data );

      if ( ln > d_ln_min ) {
         /*
          * Save outerflux data so that next coarser level
          *  can compute its coarse-fine composite flux.
          *  This is not strictly needed in this "compute residual"
          *  loop through the patches, but we put it here to
          *  avoid writing another loop for it.
          */
         tbox::Pointer<pdat::OutersideData<DIM,double> >
            oflux_data = patch->getPatchData( d_oflux_scratch_id );
#ifdef DEBUG_CHECK_ASSERTIONS
         assert( oflux_data );
#endif
         oflux_data->copy(*flux_data);
      }
   }


   if ( deallocate_flux_data_when_done ) {
      level->deallocatePatchData(flux_id);
   }

   t_compute_composite_residual->stop();
   return;
}




/*
********************************************************************
* FACOperatorStrategy<DIM> virtual computeResidualNorm             *
* function                                                         *
********************************************************************
*/

template<int DIM> double CellPoissonFACOps<DIM>::computeResidualNorm(
   const SAMRAIVectorReal<DIM,double> &residual ,
   int fine_ln ,
   int coarse_ln )
{

   if ( coarse_ln != residual.getCoarsestLevelNumber() ||
        fine_ln != residual.getFinestLevelNumber() ) {
      TBOX_ERROR("CellPoissonFACOps::computeResidualNorm() is not\n"
                 <<"set up to compute residual except on the range of\n"
                 <<"levels defining the vector.\n");
   }
   t_compute_residual_norm->start();
   /*
    * The residual vector was cloned from vectors that has
    * the proper weights associated with them, so we do not
    * have to explicitly weight the residuals.
    *
    * maxNorm: not good to use because Hypre's norm does not
    *   correspond to it.  Also maybe too sensitive to spikes.
    * L2Norm: maybe good.  But does not correspond to the
    *   scale of the quantity.
    * L1Norm: maybe good.  Correspond to scale of quantity,
    *   but may be too insensitive to spikes.
    * RMSNorm: maybe good.
    */
   double norm = residual.RMSNorm();
   t_compute_residual_norm->stop();
   return norm;
}



/*
********************************************************************
* Compute the vector weight and put it at a specified patch data   *
* index.                                                           *
********************************************************************
*/

template<int DIM> void CellPoissonFACOps<DIM>::computeVectorWeights(
   tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy ,
   int weight_id ,
   int coarsest_ln ,
   int finest_ln ) const {

   if ( coarsest_ln == -1 ) coarsest_ln = 0;
   if ( finest_ln == -1 ) finest_ln = hierarchy->getFinestLevelNumber();
   if ( finest_ln < coarsest_ln ) {
      TBOX_ERROR(d_object_name
                 << ": Illegal level number range.  finest_ln < coarsest_ln.");
   }

   int ln;
   for ( ln=finest_ln; ln >= coarsest_ln; --ln ) {

      /*
       * On every level, first assign cell volume to vector weight.
       */

      tbox::Pointer< hier::PatchLevel<DIM> > level =
         hierarchy->getPatchLevel(ln);
      for (typename hier::PatchLevel<DIM>::Iterator p(level); p; p++) {
         tbox::Pointer< hier::Patch<DIM> > patch = level->getPatch(p());
         tbox::Pointer< geom::CartesianPatchGeometry<DIM> > patch_geometry =
            patch->getPatchGeometry();
         const double* dx = patch_geometry->getDx();
         double cell_vol = dx[0];
	 if (DIM > 1) {
	    cell_vol *= dx[1];
	 }

	 if (DIM > 2) {
	    cell_vol *= dx[2];
	 }

         tbox::Pointer< pdat::CellData<DIM,double> > w =
            patch->getPatchData(weight_id);
         if ( !w ) {
            TBOX_ERROR(d_object_name
                       << ": weight id must refer to a pdat::CellVariable");
         }
         w->fillAll(cell_vol); 
      }

      /*
       * On all but the finest level, assign 0 to vector
       * weight to cells covered by finer cells.
       */

      if (ln < finest_ln) {

         /*
          * First get the boxes that describe index space of the next finer 
          * level and coarsen them to describe corresponding index space 
          * at this level.
          */

         tbox::Pointer< hier::PatchLevel<DIM> > next_finer_level
            = hierarchy->getPatchLevel(ln+1);
         hier::BoxArray<DIM> coarsened_boxes = next_finer_level->getBoxes(); 
         hier::IntVector<DIM> coarsen_ratio = next_finer_level->getRatio();
         coarsen_ratio /= level->getRatio();
         coarsened_boxes.coarsen(coarsen_ratio);

         /*
          * Then set vector weight to 0 wherever there is
          * a nonempty intersection with the next finer level.
          * Note that all assignments are local.
          */

         for (typename hier::PatchLevel<DIM>::Iterator p(level); p; p++) {

            tbox::Pointer< hier::Patch<DIM> > patch = level->getPatch(p());
            for ( int i = 0; i < coarsened_boxes.getNumberOfBoxes(); i++ ) {

               hier::Box<DIM> coarse_box = coarsened_boxes(i);
               hier::Box<DIM> intersection = coarse_box*(patch->getBox());
               if ( !intersection.empty() ) {
                  tbox::Pointer< pdat::CellData<DIM,double> > w =
                     patch->getPatchData(weight_id);
                  w->fillAll(0.0, intersection);

               }  // assignment only in non-empty intersection
            }  // loop over coarsened boxes from finer level
         }  // loop over patches in level
      }  // all levels except finest
   }  // loop over levels
}









/*
********************************************************************
* Check the validity and correctness of input data for this class. *
********************************************************************
*/

template<int DIM> void CellPoissonFACOps<DIM>::checkInputPatchDataIndices() const {
   /*
    * Check input validity and correctness.
    */
   hier::VariableDatabase<DIM> &vdb(*hier::VariableDatabase<DIM>::getDatabase());

   if( !d_poisson_spec.dIsConstant()
       && d_poisson_spec.getDPatchDataId() != -1 ) {
      tbox::Pointer<hier::Variable<DIM> > var;
      vdb.mapIndexToVariable(d_poisson_spec.getDPatchDataId(), var);
      tbox::Pointer<pdat::SideVariable<DIM,double> > diffcoef_var = var;
      if ( !diffcoef_var ) {
         TBOX_ERROR(d_object_name
                    << ": Bad diffusion coefficient patch data index.");
      }
   }

   if( !d_poisson_spec.cIsConstant() && !d_poisson_spec.cIsZero() ) {
      tbox::Pointer<hier::Variable<DIM> > var;
      vdb.mapIndexToVariable(d_poisson_spec.getCPatchDataId(), var);
      tbox::Pointer<pdat::CellVariable<DIM,double> > scalar_field_var = var;
      if ( !scalar_field_var ) {
         TBOX_ERROR(d_object_name << ": Bad linear term patch data index.");
      }
   }

   if( d_flux_id != -1 ) {
      tbox::Pointer<hier::Variable<DIM> > var;
      vdb.mapIndexToVariable(d_flux_id, var);
      tbox::Pointer<pdat::SideVariable<DIM,double> > flux_var = var;
#ifdef DEBUG_CHECK_ASSERTIONS
      assert( flux_var );
#endif
   }


}

/*
*******************************************************************
*                                                                 *
* AMR-unaware patch-centered computational kernels.               *
*                                                                 *
*******************************************************************
*/

template<int DIM> void CellPoissonFACOps<DIM>::computeFluxOnPatch(
   const hier::Patch<DIM> &patch ,
   const hier::IntVector<DIM> &ratio_to_coarser_level,
   const pdat::CellData<DIM,double> &w_data ,
   pdat::SideData<DIM,double> &Dgradw_data ) const 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( patch.inHierarchy() );
   assert( w_data.getGhostCellWidth() >= hier::IntVector<DIM>(1) );
#endif

   tbox::Pointer< geom::CartesianPatchGeometry<DIM> > patch_geom
      = patch.getPatchGeometry();
   const hier::Box<DIM> &box=patch.getBox();
   const int *lower = box.lower();
   const int *upper = box.upper();
   const double *dx = patch_geom->getDx();

   double D_value;
   tbox::Pointer<pdat::SideData<DIM,double> > D_data;
   if ( d_poisson_spec.dIsConstant() ) {
      D_value = d_poisson_spec.getDConstant();
   }
   else {
      D_data = patch.getPatchData( d_poisson_spec.getDPatchDataId() );
   }

   if ( d_poisson_spec.dIsConstant() ) {
      if (DIM == 2) {
	 compfluxcondc2d_(
	    Dgradw_data.getPointer(0) ,
	    Dgradw_data.getPointer(1) ,
	    &Dgradw_data.getGhostCellWidth()[0],
	    &Dgradw_data.getGhostCellWidth()[1] ,
	    D_value ,
	    w_data.getPointer() ,
	    &w_data.getGhostCellWidth()[0],
	    &w_data.getGhostCellWidth()[1] ,
	    &lower[0], &upper[0],
	    &lower[1], &upper[1] ,
	    dx );
      } else if (DIM == 3) {
	 compfluxcondc3d_(
	    Dgradw_data.getPointer(0) ,
	    Dgradw_data.getPointer(1) ,
	    Dgradw_data.getPointer(2) ,
	    &Dgradw_data.getGhostCellWidth()[0] ,
	    &Dgradw_data.getGhostCellWidth()[1] ,
	    &Dgradw_data.getGhostCellWidth()[2] ,
	    D_value ,
	    w_data.getPointer() ,
	    &w_data.getGhostCellWidth()[0] ,
	    &w_data.getGhostCellWidth()[1] ,
	    &w_data.getGhostCellWidth()[2] ,
	    &lower[0], &upper[0],
	    &lower[1], &upper[1] ,
	    &lower[2], &upper[2] ,
	    dx );
      }
   } else {
      if (DIM == 2) {
	 compfluxvardc2d_(
	    Dgradw_data.getPointer(0) ,
	    Dgradw_data.getPointer(1) ,
	    &Dgradw_data.getGhostCellWidth()[0],
	    &Dgradw_data.getGhostCellWidth()[1] ,
	    D_data->getPointer(0) ,
	    D_data->getPointer(1) ,
	    &D_data->getGhostCellWidth()[0],
	    &D_data->getGhostCellWidth()[1] ,
	    w_data.getPointer() ,
	    &w_data.getGhostCellWidth()[0],
	    &w_data.getGhostCellWidth()[1] ,
	    &lower[0], &upper[0],
	    &lower[1], &upper[1] ,
	    dx );
      } if (DIM == 3) {
	 compfluxvardc3d_(
	    Dgradw_data.getPointer(0) ,
	    Dgradw_data.getPointer(1) ,
	    Dgradw_data.getPointer(2) ,
	    &Dgradw_data.getGhostCellWidth()[0] ,
	    &Dgradw_data.getGhostCellWidth()[1] ,
	    &Dgradw_data.getGhostCellWidth()[2] ,
	    D_data->getPointer(0) ,
	    D_data->getPointer(1) ,
	    D_data->getPointer(2) ,
	    &D_data->getGhostCellWidth()[0],
	    &D_data->getGhostCellWidth()[1] ,
	    &D_data->getGhostCellWidth()[2] ,
	    w_data.getPointer() ,
	    &w_data.getGhostCellWidth()[0] ,
	    &w_data.getGhostCellWidth()[1] ,
	    &w_data.getGhostCellWidth()[2] ,
	    &lower[0], &upper[0],
	    &lower[1], &upper[1] ,
	    &lower[2], &upper[2] ,
	    dx );
      }
   }

   const int patch_ln = patch.getPatchLevelNumber();

   if ( d_cf_discretization == "Ewing" && patch_ln > d_ln_min ) {
      ewingFixFlux( patch ,
                    w_data ,
                    Dgradw_data ,
                    ratio_to_coarser_level );
   }

}

template<int DIM> void CellPoissonFACOps<DIM>::computeResidualOnPatch(
   const hier::Patch<DIM> &patch ,
   const pdat::SideData<DIM,double> &flux_data ,
   const pdat::CellData<DIM,double> &soln_data ,
   const pdat::CellData<DIM,double> &rhs_data ,
   pdat::CellData<DIM,double> &residual_data ) const
{

   tbox::Pointer< geom::CartesianPatchGeometry<DIM> > patch_geom
     = patch.getPatchGeometry();
   const hier::Box<DIM> &box=patch.getBox();
   const int *lower = box.lower();
   const int *upper = box.upper();
   const double *dx = patch_geom->getDx();

   tbox::Pointer<pdat::CellData<DIM,double> > scalar_field_data;
   double scalar_field_constant;
   if ( d_poisson_spec.cIsVariable() ) {
      scalar_field_data =
         patch.getPatchData( d_poisson_spec.getCPatchDataId() );
      if (DIM == 2) {
	 compresvarsca2d_(
	    flux_data.getPointer(0) ,
	    flux_data.getPointer(1) ,
	    &flux_data.getGhostCellWidth()[0],
	    &flux_data.getGhostCellWidth()[1] ,
	    rhs_data.getPointer() ,
	    &rhs_data.getGhostCellWidth()[0],
	    &rhs_data.getGhostCellWidth()[1] ,
	    residual_data.getPointer() ,
	    &residual_data.getGhostCellWidth()[0],
	    &residual_data.getGhostCellWidth()[1] ,
	    scalar_field_data->getPointer() ,
	    &scalar_field_data->getGhostCellWidth()[0],
	    &scalar_field_data->getGhostCellWidth()[1] ,
	    soln_data.getPointer() ,
	    &soln_data.getGhostCellWidth()[0],
	    &soln_data.getGhostCellWidth()[1] ,
	    &lower[0], &upper[0], &lower[1], &upper[1] ,
	    dx );
      } else if (DIM == 3) {
	 compresvarsca3d_(
	    flux_data.getPointer(0) ,
	    flux_data.getPointer(1) ,
	    flux_data.getPointer(2) ,
	    &flux_data.getGhostCellWidth()[0],
	    &flux_data.getGhostCellWidth()[1] ,
	    &flux_data.getGhostCellWidth()[2] ,
	    rhs_data.getPointer() ,
	    &rhs_data.getGhostCellWidth()[0],
	    &rhs_data.getGhostCellWidth()[1] ,
	    &rhs_data.getGhostCellWidth()[2] ,
	    residual_data.getPointer() ,
	    &residual_data.getGhostCellWidth()[0],
	    &residual_data.getGhostCellWidth()[1] ,
	    &residual_data.getGhostCellWidth()[2] ,
	    scalar_field_data->getPointer() ,
	    &scalar_field_data->getGhostCellWidth()[0],
	    &scalar_field_data->getGhostCellWidth()[1] ,
	    &scalar_field_data->getGhostCellWidth()[2] ,
	    soln_data.getPointer() ,
	    &soln_data.getGhostCellWidth()[0],
	    &soln_data.getGhostCellWidth()[1] ,
	    &soln_data.getGhostCellWidth()[2] ,
	    &lower[0], &upper[0], &lower[1], &upper[1] , &lower[2], &upper[2] ,
	    dx );
      }
   }
   else if ( d_poisson_spec.cIsConstant() ) {
      scalar_field_constant = d_poisson_spec.getCConstant();
      if (DIM == 2) {
	 compresconsca2d_(
	    flux_data.getPointer(0) ,
	    flux_data.getPointer(1) ,
	    &flux_data.getGhostCellWidth()[0],
	    &flux_data.getGhostCellWidth()[1] ,
	    rhs_data.getPointer() ,
	    &rhs_data.getGhostCellWidth()[0],
	    &rhs_data.getGhostCellWidth()[1] ,
	    residual_data.getPointer() ,
	    &residual_data.getGhostCellWidth()[0],
	    &residual_data.getGhostCellWidth()[1] ,
	    scalar_field_constant ,
	    soln_data.getPointer() ,
	    &soln_data.getGhostCellWidth()[0],
	    &soln_data.getGhostCellWidth()[1] ,
	    &lower[0], &upper[0], &lower[1], &upper[1] ,
	    dx );
      } else if (DIM == 3) { 
	 compresconsca3d_(
	    flux_data.getPointer(0) ,
	    flux_data.getPointer(1) ,
	    flux_data.getPointer(2) ,
	    &flux_data.getGhostCellWidth()[0],
	    &flux_data.getGhostCellWidth()[1] ,
	    &flux_data.getGhostCellWidth()[2] ,
	    rhs_data.getPointer() ,
	    &rhs_data.getGhostCellWidth()[0],
	    &rhs_data.getGhostCellWidth()[1] ,
	    &rhs_data.getGhostCellWidth()[2] ,
	    residual_data.getPointer() ,
	    &residual_data.getGhostCellWidth()[0],
	    &residual_data.getGhostCellWidth()[1] ,
	    &residual_data.getGhostCellWidth()[2] ,
	    scalar_field_constant ,
	    soln_data.getPointer() ,
	    &soln_data.getGhostCellWidth()[0],
	    &soln_data.getGhostCellWidth()[1] ,
	    &soln_data.getGhostCellWidth()[2] ,
	    &lower[0], &upper[0], &lower[1], &upper[1] , &lower[2], &upper[2] ,
	    dx );
      }
   }
   else {
      scalar_field_constant = 0.0;
      if (DIM == 2) {
	 compresconsca2d_(
	    flux_data.getPointer(0) ,
	    flux_data.getPointer(1) ,
	    &flux_data.getGhostCellWidth()[0],
	    &flux_data.getGhostCellWidth()[1] ,
	    rhs_data.getPointer() ,
	    &rhs_data.getGhostCellWidth()[0],
	    &rhs_data.getGhostCellWidth()[1] ,
	    residual_data.getPointer() ,
	    &residual_data.getGhostCellWidth()[0],
	    &residual_data.getGhostCellWidth()[1] ,
	    0.0 ,
	    soln_data.getPointer() ,
	    &soln_data.getGhostCellWidth()[0],
	    &soln_data.getGhostCellWidth()[1] ,
	    &lower[0], &upper[0], &lower[1], &upper[1] ,
	    dx );
      } else if (DIM == 3) {
	 compresconsca3d_(
	    flux_data.getPointer(0) ,
	    flux_data.getPointer(1) ,
	    flux_data.getPointer(2) ,
	    &flux_data.getGhostCellWidth()[0],
	    &flux_data.getGhostCellWidth()[1] ,
	    &flux_data.getGhostCellWidth()[2] ,
	    rhs_data.getPointer() ,
	    &rhs_data.getGhostCellWidth()[0],
	    &rhs_data.getGhostCellWidth()[1] ,
	    &rhs_data.getGhostCellWidth()[2] ,
	    residual_data.getPointer() ,
	    &residual_data.getGhostCellWidth()[0],
	    &residual_data.getGhostCellWidth()[1] ,
	    &residual_data.getGhostCellWidth()[2] ,
	    0.0 ,
	    soln_data.getPointer() ,
	    &soln_data.getGhostCellWidth()[0],
	    &soln_data.getGhostCellWidth()[1] ,
	    &soln_data.getGhostCellWidth()[2] ,
	    &lower[0], &upper[0], &lower[1], &upper[1] , &lower[2], &upper[2] ,
	    dx );
      }
   }

   return;
}



template<int DIM> void CellPoissonFACOps<DIM>::redOrBlackSmoothingOnPatch(
   const hier::Patch<DIM> &patch ,
   const pdat::SideData<DIM,double> &flux_data ,
   const pdat::CellData<DIM,double> &rhs_data ,
   pdat::CellData<DIM,double> &soln_data ,
   char red_or_black ,
   double *p_maxres ) const
{

#ifdef DEBUG_CHECK_ASSERTIONS
   assert( red_or_black == 'r' || red_or_black == 'b' );
#endif
   const int offset = red_or_black == 'r' ? 0 : 1;
   tbox::Pointer< geom::CartesianPatchGeometry<DIM> > patch_geom
      = patch.getPatchGeometry();
   const hier::Box<DIM> &box=patch.getBox();
   const int *lower = box.lower();
   const int *upper = box.upper();
   const double *dx = patch_geom->getDx();

   tbox::Pointer<pdat::CellData<DIM,double> > scalar_field_data;
   double scalar_field_constant;
   tbox::Pointer<pdat::SideData<DIM,double> > diffcoef_data;
   double diffcoef_constant;

   if ( d_poisson_spec.cIsVariable() ) {
      scalar_field_data =
         patch.getPatchData( d_poisson_spec.getCPatchDataId() );
   }
   else if ( d_poisson_spec.cIsConstant() ) {
      scalar_field_constant = d_poisson_spec.getCConstant();
   }
   else {
      scalar_field_constant = 0.0;
   }
   if ( d_poisson_spec.dIsVariable() ) {
      diffcoef_data = patch.getPatchData( d_poisson_spec.getDPatchDataId() );
   }
   else {
      diffcoef_constant = d_poisson_spec.getDConstant();
   }

   double maxres=0.0;
   if ( d_poisson_spec.dIsVariable() && d_poisson_spec.cIsVariable() ) {
      if (DIM == 2) {
	 rbgswithfluxmaxvardcvarsf2d_(
	    flux_data.getPointer(0) ,
	    flux_data.getPointer(1) ,
	    &flux_data.getGhostCellWidth()[0],
	    &flux_data.getGhostCellWidth()[1] ,
	    diffcoef_data->getPointer(0) ,
	    diffcoef_data->getPointer(1) ,
	    &diffcoef_data->getGhostCellWidth()[0],
	    &diffcoef_data->getGhostCellWidth()[1] ,
	    rhs_data.getPointer() ,
	    &rhs_data.getGhostCellWidth()[0],
	    &rhs_data.getGhostCellWidth()[1] ,
	    scalar_field_data->getPointer() ,
	    &scalar_field_data->getGhostCellWidth()[0],
	    &scalar_field_data->getGhostCellWidth()[1] ,
	    soln_data.getPointer() ,
	    &soln_data.getGhostCellWidth()[0],
	    &soln_data.getGhostCellWidth()[1] ,
	    &lower[0], &upper[0],
	    &lower[1], &upper[1] ,
	    dx ,
	    &offset, &maxres );
      } else if (DIM == 3) {
	 rbgswithfluxmaxvardcvarsf3d_(
	    flux_data.getPointer(0) ,
	    flux_data.getPointer(1) ,
	    flux_data.getPointer(2) ,
	    &flux_data.getGhostCellWidth()[0],
	    &flux_data.getGhostCellWidth()[1] ,
	    &flux_data.getGhostCellWidth()[2] ,
	    diffcoef_data->getPointer(0) ,
	    diffcoef_data->getPointer(1) ,
	    diffcoef_data->getPointer(2) ,
	    &diffcoef_data->getGhostCellWidth()[0],
	    &diffcoef_data->getGhostCellWidth()[1] ,
	    &diffcoef_data->getGhostCellWidth()[2] ,
	    rhs_data.getPointer() ,
	    &rhs_data.getGhostCellWidth()[0],
	    &rhs_data.getGhostCellWidth()[1] ,
	    &rhs_data.getGhostCellWidth()[2] ,
	    scalar_field_data->getPointer() ,
	    &scalar_field_data->getGhostCellWidth()[0],
	    &scalar_field_data->getGhostCellWidth()[1] ,
	    &scalar_field_data->getGhostCellWidth()[2] ,
	    soln_data.getPointer() ,
	    &soln_data.getGhostCellWidth()[0],
	    &soln_data.getGhostCellWidth()[1] ,
	    &soln_data.getGhostCellWidth()[2] ,
	    &lower[0], &upper[0],
	    &lower[1], &upper[1] ,
	    &lower[2], &upper[2] ,
	    dx ,
	    &offset, &maxres );
      }
   }
   else if ( d_poisson_spec.dIsVariable() && d_poisson_spec.cIsConstant() ) {
      if (DIM == 2) {
	 rbgswithfluxmaxvardcconsf2d_(
	    flux_data.getPointer(0) ,
	    flux_data.getPointer(1) ,
	    &flux_data.getGhostCellWidth()[0],
	    &flux_data.getGhostCellWidth()[1] ,
	    diffcoef_data->getPointer(0) ,
	    diffcoef_data->getPointer(1) ,
	    &diffcoef_data->getGhostCellWidth()[0],
	    &diffcoef_data->getGhostCellWidth()[1] ,
	    rhs_data.getPointer() ,
	    &rhs_data.getGhostCellWidth()[0],
	    &rhs_data.getGhostCellWidth()[1] ,
	    scalar_field_constant ,
	    soln_data.getPointer() ,
	    &soln_data.getGhostCellWidth()[0],
	    &soln_data.getGhostCellWidth()[1] ,
	    &lower[0], &upper[0],
	    &lower[1], &upper[1] ,
	    dx ,
	    &offset, &maxres );
      } else if (DIM == 3) {
	 rbgswithfluxmaxvardcconsf3d_(
	    flux_data.getPointer(0) ,
	    flux_data.getPointer(1) ,
	    flux_data.getPointer(2) ,
	    &flux_data.getGhostCellWidth()[0],
	    &flux_data.getGhostCellWidth()[1] ,
	    &flux_data.getGhostCellWidth()[2] ,
	    diffcoef_data->getPointer(0) ,
	    diffcoef_data->getPointer(1) ,
	    diffcoef_data->getPointer(2) ,
	    &diffcoef_data->getGhostCellWidth()[0],
	    &diffcoef_data->getGhostCellWidth()[1] ,
	    &diffcoef_data->getGhostCellWidth()[2] ,
	    rhs_data.getPointer() ,
	    &rhs_data.getGhostCellWidth()[0],
	    &rhs_data.getGhostCellWidth()[1] ,
	    &rhs_data.getGhostCellWidth()[2] ,
	    scalar_field_constant ,
	    soln_data.getPointer() ,
	    &soln_data.getGhostCellWidth()[0],
	    &soln_data.getGhostCellWidth()[1] ,
	    &soln_data.getGhostCellWidth()[2] ,
	    &lower[0], &upper[0],
	    &lower[1], &upper[1] ,
	    &lower[2], &upper[2] ,
	    dx ,
	    &offset, &maxres );
      }
   }
   else if ( d_poisson_spec.dIsVariable() && d_poisson_spec.cIsZero() ) {
      if (DIM == 2) {
	 rbgswithfluxmaxvardcconsf2d_(
	    flux_data.getPointer(0) ,
	    flux_data.getPointer(1) ,
	    &flux_data.getGhostCellWidth()[0],
	    &flux_data.getGhostCellWidth()[1] ,
	    diffcoef_data->getPointer(0) ,
	    diffcoef_data->getPointer(1) ,
	    &diffcoef_data->getGhostCellWidth()[0],
	    &diffcoef_data->getGhostCellWidth()[1] ,
	    rhs_data.getPointer() ,
	    &rhs_data.getGhostCellWidth()[0],
	    &rhs_data.getGhostCellWidth()[1] ,
	    0.0 ,
	    soln_data.getPointer() ,
	    &soln_data.getGhostCellWidth()[0],
	    &soln_data.getGhostCellWidth()[1] ,
	    &lower[0], &upper[0],
	    &lower[1], &upper[1] ,
	    dx ,
	    &offset, &maxres );
      } else if (DIM == 3) {
	 rbgswithfluxmaxvardcconsf3d_(
	    flux_data.getPointer(0) ,
	    flux_data.getPointer(1) ,
	    flux_data.getPointer(2) ,
	    &flux_data.getGhostCellWidth()[0],
	    &flux_data.getGhostCellWidth()[1] ,
	    &flux_data.getGhostCellWidth()[2] ,
	    diffcoef_data->getPointer(0) ,
	    diffcoef_data->getPointer(1) ,
	    diffcoef_data->getPointer(2) ,
	    &diffcoef_data->getGhostCellWidth()[0],
	    &diffcoef_data->getGhostCellWidth()[1] ,
	    &diffcoef_data->getGhostCellWidth()[2] ,
	    rhs_data.getPointer() ,
	    &rhs_data.getGhostCellWidth()[0],
	    &rhs_data.getGhostCellWidth()[1] ,
	    &rhs_data.getGhostCellWidth()[2] ,
	    0.0 ,
	    soln_data.getPointer() ,
	    &soln_data.getGhostCellWidth()[0],
	    &soln_data.getGhostCellWidth()[1] ,
	    &soln_data.getGhostCellWidth()[2] ,
	    &lower[0], &upper[0],
	    &lower[1], &upper[1] ,
	    &lower[2], &upper[2] ,
	    dx ,
	    &offset, &maxres );
      }
   }
   else if ( !d_poisson_spec.dIsVariable() && d_poisson_spec.cIsVariable() ) {
      if (DIM == 2) {
	 rbgswithfluxmaxcondcvarsf2d_(
	    flux_data.getPointer(0) ,
	    flux_data.getPointer(1) ,
	    &flux_data.getGhostCellWidth()[0],
	    &flux_data.getGhostCellWidth()[1] ,
	    diffcoef_constant ,
	    rhs_data.getPointer() ,
	    &rhs_data.getGhostCellWidth()[0],
	    &rhs_data.getGhostCellWidth()[1] ,
	    scalar_field_data->getPointer() ,
	    &scalar_field_data->getGhostCellWidth()[0],
	    &scalar_field_data->getGhostCellWidth()[1] ,
	    soln_data.getPointer() ,
	    &soln_data.getGhostCellWidth()[0],
	    &soln_data.getGhostCellWidth()[1] ,
	    &lower[0], &upper[0],
	    &lower[1], &upper[1] ,
	    dx ,
	    &offset, &maxres );
      } else if (DIM == 3) {
	 rbgswithfluxmaxcondcvarsf3d_(
	    flux_data.getPointer(0) ,
	    flux_data.getPointer(1) ,
	    flux_data.getPointer(2) ,
	    &flux_data.getGhostCellWidth()[0],
	    &flux_data.getGhostCellWidth()[1] ,
	    &flux_data.getGhostCellWidth()[2] ,
	    diffcoef_constant ,
	    rhs_data.getPointer() ,
	    &rhs_data.getGhostCellWidth()[0],
	    &rhs_data.getGhostCellWidth()[1] ,
	    &rhs_data.getGhostCellWidth()[2] ,
	    scalar_field_data->getPointer() ,
	    &scalar_field_data->getGhostCellWidth()[0],
	    &scalar_field_data->getGhostCellWidth()[1] ,
	    &scalar_field_data->getGhostCellWidth()[2] ,
	    soln_data.getPointer() ,
	    &soln_data.getGhostCellWidth()[0],
	    &soln_data.getGhostCellWidth()[1] ,
	    &soln_data.getGhostCellWidth()[2] ,
	    &lower[0], &upper[0],
	    &lower[1], &upper[1] ,
	    &lower[2], &upper[2] ,
	    dx ,
	    &offset, &maxres );
      }
   }
   else if ( !d_poisson_spec.dIsVariable() && d_poisson_spec.cIsConstant() ) {
      if (DIM == 2) {
	 rbgswithfluxmaxcondcconsf2d_(
	    flux_data.getPointer(0) ,
	    flux_data.getPointer(1) ,
	    &flux_data.getGhostCellWidth()[0],
	    &flux_data.getGhostCellWidth()[1] ,
	    diffcoef_constant ,
	    rhs_data.getPointer() ,
	    &rhs_data.getGhostCellWidth()[0],
	    &rhs_data.getGhostCellWidth()[1] ,
	    scalar_field_constant ,
	    soln_data.getPointer() ,
	    &soln_data.getGhostCellWidth()[0],
	    &soln_data.getGhostCellWidth()[1] ,
	    &lower[0], &upper[0],
	    &lower[1], &upper[1] ,
	    dx ,
	    &offset, &maxres );
      } else if (DIM == 3) {
	 rbgswithfluxmaxcondcconsf3d_(
	    flux_data.getPointer(0) ,
	    flux_data.getPointer(1) ,
	    flux_data.getPointer(2) ,
	    &flux_data.getGhostCellWidth()[0],
	    &flux_data.getGhostCellWidth()[1] ,
	    &flux_data.getGhostCellWidth()[2] ,
	    diffcoef_constant ,
	    rhs_data.getPointer() ,
	    &rhs_data.getGhostCellWidth()[0],
	    &rhs_data.getGhostCellWidth()[1] ,
	    &rhs_data.getGhostCellWidth()[2] ,
	    scalar_field_constant ,
	    soln_data.getPointer() ,
	    &soln_data.getGhostCellWidth()[0],
	    &soln_data.getGhostCellWidth()[1] ,
	    &soln_data.getGhostCellWidth()[2] ,
	    &lower[0], &upper[0],
	    &lower[1], &upper[1] ,
	    &lower[2], &upper[2] ,
	    dx ,
	    &offset, &maxres );
      }
   }
   else if ( !d_poisson_spec.dIsVariable() && d_poisson_spec.cIsZero() ) {
      if (DIM == 2) {
	 rbgswithfluxmaxcondcconsf2d_(
	    flux_data.getPointer(0) ,
	    flux_data.getPointer(1) ,
	    &flux_data.getGhostCellWidth()[0],
	    &flux_data.getGhostCellWidth()[1] ,
	    diffcoef_constant ,
	    rhs_data.getPointer() ,
	    &rhs_data.getGhostCellWidth()[0],
	    &rhs_data.getGhostCellWidth()[1] ,
	    0.0 ,
	    soln_data.getPointer() ,
	    &soln_data.getGhostCellWidth()[0],
	    &soln_data.getGhostCellWidth()[1] ,
	    &lower[0], &upper[0],
	    &lower[1], &upper[1] ,
	    dx ,
	    &offset, &maxres );
      } else if (DIM == 3) {
	 rbgswithfluxmaxcondcconsf3d_(
	    flux_data.getPointer(0) ,
	    flux_data.getPointer(1) ,
	    flux_data.getPointer(2) ,
	    &flux_data.getGhostCellWidth()[0],
	    &flux_data.getGhostCellWidth()[1] ,
	    &flux_data.getGhostCellWidth()[2] ,
	    diffcoef_constant ,
	    rhs_data.getPointer() ,
	    &rhs_data.getGhostCellWidth()[0],
	    &rhs_data.getGhostCellWidth()[1] ,
	    &rhs_data.getGhostCellWidth()[2] ,
	    0.0 ,
	    soln_data.getPointer() ,
	    &soln_data.getGhostCellWidth()[0],
	    &soln_data.getGhostCellWidth()[1] ,
	    &soln_data.getGhostCellWidth()[2] ,
	    &lower[0], &upper[0],
	    &lower[1], &upper[1] ,
	    &lower[2], &upper[2] ,
	    dx ,
	    &offset, &maxres );
      }
   }

   *p_maxres = maxres;
   return;
}



template<int DIM> void
CellPoissonFACOps<DIM>::xeqScheduleProlongation(
   int dst_id,
   int src_id,
   int scr_id,
   int dest_ln
   )
{
   if ( ! d_prolongation_refine_schedules[dest_ln] ) {
      TBOX_ERROR("Expected schedule not found.");
   }
   xfer::RefineAlgorithm<DIM> refiner;
   refiner.
      registerRefine( dst_id ,
                      src_id ,
                      scr_id ,
                      d_prolongation_refine_operator );
   refiner.
      resetSchedule(d_prolongation_refine_schedules[dest_ln]);
   d_prolongation_refine_schedules[dest_ln]->fillData(0.0);
   d_prolongation_refine_algorithm->
      resetSchedule(d_prolongation_refine_schedules[dest_ln]);
   return;
}



template<int DIM> void
CellPoissonFACOps<DIM>::xeqScheduleURestriction(
   int dst_id,
   int src_id,
   int dest_ln
   )
{
   if ( ! d_urestriction_coarsen_schedules[dest_ln] ) {
      TBOX_ERROR("Expected schedule not found.");
   }
   xfer::CoarsenAlgorithm<DIM> coarsener;
   coarsener.
      registerCoarsen( dst_id ,
                       src_id ,
                       d_urestriction_coarsen_operator );
   coarsener.
      resetSchedule(d_urestriction_coarsen_schedules[dest_ln]);
   d_urestriction_coarsen_schedules[dest_ln]->coarsenData();
   d_urestriction_coarsen_algorithm->
      resetSchedule(d_urestriction_coarsen_schedules[dest_ln]);
   return;
}



template<int DIM> void
CellPoissonFACOps<DIM>::xeqScheduleRRestriction(
   int dst_id,
   int src_id,
   int dest_ln
   )
{
   if ( ! d_rrestriction_coarsen_schedules[dest_ln] ) {
      TBOX_ERROR("Expected schedule not found.");
   }
   xfer::CoarsenAlgorithm<DIM> coarsener;
   coarsener.
      registerCoarsen( dst_id ,
                       src_id ,
                       d_rrestriction_coarsen_operator );
   coarsener.
      resetSchedule(d_rrestriction_coarsen_schedules[dest_ln]);
   d_rrestriction_coarsen_schedules[dest_ln]->coarsenData();
   d_rrestriction_coarsen_algorithm->
      resetSchedule(d_rrestriction_coarsen_schedules[dest_ln]);
   return;
}



template<int DIM> void
CellPoissonFACOps<DIM>::xeqScheduleFluxCoarsen(
   int dst_id,
   int src_id,
   int dest_ln
   )
{
   if ( ! d_flux_coarsen_schedules[dest_ln] ) {
      TBOX_ERROR("Expected schedule not found.");
   }
   xfer::CoarsenAlgorithm<DIM> coarsener;
   coarsener.
      registerCoarsen( dst_id ,
                       src_id ,
                       d_flux_coarsen_operator );
   coarsener.
      resetSchedule(d_flux_coarsen_schedules[dest_ln]);
   d_flux_coarsen_schedules[dest_ln]->coarsenData();
   d_flux_coarsen_algorithm->
      resetSchedule(d_flux_coarsen_schedules[dest_ln]);
   return;
}



template<int DIM> void
CellPoissonFACOps<DIM>::xeqScheduleGhostFill(
   int dst_id,
   int dest_ln
   )
{
   if ( ! d_ghostfill_refine_schedules[dest_ln] ) {
      TBOX_ERROR("Expected schedule not found.");
   }
   xfer::RefineAlgorithm<DIM> refiner;
   refiner.
      registerRefine( dst_id ,
                      dst_id ,
                      dst_id ,
                      d_ghostfill_refine_operator );
   refiner.
      resetSchedule(d_ghostfill_refine_schedules[dest_ln]);
   d_ghostfill_refine_schedules[dest_ln]->fillData(0.0);
   d_ghostfill_refine_algorithm->
      resetSchedule(d_ghostfill_refine_schedules[dest_ln]);
   return;
}



template<int DIM> void
CellPoissonFACOps<DIM>::xeqScheduleGhostFillNoCoarse(
   int dst_id,
   int dest_ln
   )
{
   if ( ! d_ghostfill_nocoarse_refine_schedules[dest_ln] ) {
      TBOX_ERROR("Expected schedule not found.");
   }
   xfer::RefineAlgorithm<DIM> refiner;
   refiner.
      registerRefine( dst_id ,
                      dst_id ,
                      dst_id ,
                      d_ghostfill_nocoarse_refine_operator );
   refiner.
      resetSchedule(d_ghostfill_nocoarse_refine_schedules[dest_ln]);
   d_ghostfill_nocoarse_refine_schedules[dest_ln]->fillData(0.0);
   d_ghostfill_nocoarse_refine_algorithm->
      resetSchedule(d_ghostfill_nocoarse_refine_schedules[dest_ln]);
   return;
}



template<int DIM> void
CellPoissonFACOps<DIM>::freeVariables()
{
   s_cell_scratch_var.setNull();
   s_flux_scratch_var.setNull();
   s_oflux_scratch_var.setNull();
   return;
}


}
}
#endif
