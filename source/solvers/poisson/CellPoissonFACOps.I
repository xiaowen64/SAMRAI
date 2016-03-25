/*
 * File:        CellPoissonFACOps.I
 * Package:     SAMRAI solvers
 * Copyright:	(c) 1997-2005 The Regents of the University of California
 * Revision:	$Revision: 453 $
 * Modified:	$Date: 2005-06-16 10:19:28 -0700 (Thu, 16 Jun 2005) $
 * Description:	Operator class for solving scalar Poisson using FAC
 */

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
    namespace solv {




template<int DIM> inline
void CellPoissonFACOps<DIM>::setPreconditioner(
   const FACPreconditioner<DIM> *preconditioner ) {
   d_preconditioner = preconditioner;
   return;
}


#ifdef HAVE_HYPRE
template<int DIM> inline
void CellPoissonFACOps<DIM>::setUseSMG( bool use_smg )
{
   if ( d_hierarchy ) {
      TBOX_ERROR(d_object_name << ": setUseSMG(bool) may NOT be called\n"
                 << "while the solver state is initialized, as that\n"
                 << "would lead to a corrupted solver state.\n");
   }
   d_hypre_solver.setUseSMG(use_smg);
   return;
}
#endif



/*
********************************************************************
* Set the physical boundary condition object.                      *
********************************************************************
*/

template<int DIM> inline
void CellPoissonFACOps<DIM>::setPhysicalBcCoefObject(
   const RobinBcCoefStrategy<DIM> *physical_bc_coef )
{
   d_physical_bc_coef = physical_bc_coef;
   d_bc_helper.setCoefImplementation(physical_bc_coef);
#ifdef HAVE_HYPRE
   d_hypre_solver.setPhysicalBcCoefObject(d_physical_bc_coef);
#endif
   return;
}




/*
********************************************************************
* Set the object specifying the parameters of the Poisson equation *
********************************************************************
*/

template<int DIM> inline
void CellPoissonFACOps<DIM>::setPoissonSpecifications(
   const PoissonSpecifications &spec )
{
   d_poisson_spec = spec;
   return;
}




/*
********************************************************************
********************************************************************
*/

template<int DIM> inline
void CellPoissonFACOps<DIM>::enableLogging( bool enable_logging )
{
   d_enable_logging = enable_logging;
   return;
}




/*
********************************************************************
* Set the patch data id for the flux.                              *
********************************************************************
*/

template<int DIM> inline
void CellPoissonFACOps<DIM>::setFluxId( int flux_id ) {
   d_flux_id = flux_id;
#ifdef DEBUG_CHECK_ASSERTIONS
   checkInputPatchDataIndices();
#endif
}




/*
********************************************************************
* Set the choice for smoothing algorithm.                          *
********************************************************************
*/

template<int DIM> inline
void CellPoissonFACOps<DIM>::setSmoothingChoice(
   const string &smoothing_choice )
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if ( smoothing_choice != "redblack" ) {
      TBOX_ERROR(d_object_name << ": Bad smoothing choice '"
                 << smoothing_choice
                 << "' in CellPoissonFACOps<DIM>::setSmoothingChoice.");
   }
#endif
   d_smoothing_choice = smoothing_choice;
   return;
}


/*
********************************************************************
* Set the choice for the coarse level solver.                      *
********************************************************************
*/

template<int DIM> inline
void CellPoissonFACOps<DIM>::setCoarsestLevelSolverChoice( const string &choice ) {
#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef HAVE_HYPRE
   if ( choice == "hypre" ) {
      TBOX_ERROR(d_object_name << ": HYPRe library is not available.\n");
   }
#endif
#endif
   if (  choice == "redblack"
      || choice == "hypre" ) {
      d_coarse_solver_choice = choice;
   }
   else {
      TBOX_ERROR(d_object_name << ": Bad coarse level solver choice '"
                 << choice
                 << "' in scapCellPoissonOpsX::setCoarseLevelSolver.");
   }
   return;
}


/*
********************************************************************
* Set the tolerance for the coarse level solver.                   *
********************************************************************
*/

template<int DIM> inline
void CellPoissonFACOps<DIM>::setCoarsestLevelSolverTolerance( double tol ) {
   d_coarse_solver_tolerance = tol;
   return;
}


/*
********************************************************************
* Set the tolerance for the coarse level solver.                   *
********************************************************************
*/

template<int DIM> inline
void CellPoissonFACOps<DIM>::setCoarsestLevelSolverMaxIterations(
   int max_iterations ) {
#ifdef DEBUG_CHECK_ASSERTIONS
   if( max_iterations < 0 ) {
      TBOX_ERROR(d_object_name << ": Invalid number of max iterations\n");
   }
#endif
   d_coarse_solver_max_iterations = max_iterations;
   return;
}


/*
********************************************************************
* Set the coarse-fine discretization method.                       *
********************************************************************
*/

template<int DIM> inline
void CellPoissonFACOps<DIM>::setCoarseFineDiscretization(
   const string &coarsefine_method ) {
#ifdef DEBUG_CHECK_ASSERTIONS
   if ( d_hierarchy ) {
      TBOX_ERROR(d_object_name << ": Cannot change coarse-fine\n"
                 << "discretization method while operator state\n"
                 << "is initialized because that causes a\n"
                 << "corruption in the state.\n");
   }
#endif
   d_cf_discretization = coarsefine_method;
   return;
}


/*
********************************************************************
* Set the prolongation method                                      *
********************************************************************
*/

template<int DIM> inline
void CellPoissonFACOps<DIM>::setProlongationMethod(
   const string &prolongation_method ) {
#ifdef DEBUG_CHECK_ASSERTIONS
   if ( d_hierarchy ) {
      TBOX_ERROR(d_object_name << ": Cannot change prolongation method\n"
                 << "while operator state is initialized because that\n"
                 << "causes a corruption in the state.\n");
   }
#endif
   d_prolongation_method = prolongation_method;
   return;
}



}
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
