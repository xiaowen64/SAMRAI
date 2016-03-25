/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/examples/FAC/FACPoisson.C $
 * Package:     SAMRAI application
 * Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 2043 $
 * Modified:    $LastChangedDate: 2008-03-12 09:14:32 -0700 (Wed, 12 Mar 2008) $
 * Description: Numerical routines for example FAC Poisson solver
 */

#include "FACPoisson.h"

#include <iostream>

#include "IntVector.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "SimpleCellRobinBcCoefs.h"
#include "CellData.h"
#include "HierarchyCellDataOpsReal.h"
#include "SideData.h"
#include "PoissonSpecifications.h"
#include "tbox/Utilities.h"
#include "Variable.h"
#include "VariableDatabase.h"



extern "C" {
#if NDIM == 2
   void setexactandrhs_( const int &ifirst0,
                         const int &ilast0,
                         const int &ifirst1,
                         const int &ilast1,
                         double *exact,
                         double *rhs,
                         const double *dx,
                         const double *xlower );

#elif NDIM == 3
   void setexactandrhs_( const int &ifirst0,
                         const int &ilast0,
                         const int &ifirst1,
                         const int &ilast1,
                         const int &ifirst2,
                         const int &ilast2,
                         double *exact,
                         double *rhs,
                         const double *dx,
                         const double *xlower );
#endif
}

namespace SAMRAI {

/*
*************************************************************************
* Constructor creates a unique context for the object and register      *
* all its internal variables with the variable database.                *
*************************************************************************
*/
FACPoisson::FACPoisson(
  const std::string &object_name,
  tbox::Pointer<tbox::Database> database )
: d_object_name(object_name),
  d_hierarchy(NULL),
  d_poisson_fac_solver( object_name+"::poisson_hypre",
                        ( !database.isNull() &&
                          database->isDatabase("fac_solver") ) ?
                        database->getDatabase("fac_solver"):
                        tbox::Pointer<tbox::Database>(NULL) ),
  d_bc_coefs( object_name+"::bc_coefs",
              ( !database.isNull() &&
                database->isDatabase("bc_coefs") ) ?
              database->getDatabase("bc_coefs"):
              tbox::Pointer<tbox::Database>(NULL) ),
  d_context()
{

  hier::VariableDatabase<NDIM> *vdb =
     hier::VariableDatabase<NDIM>::getDatabase();

  /*
   * Get a unique context for variables owned by this object.
   */
  d_context = vdb->getContext( d_object_name + ":Context" );

  /*
   * Register variables with hier::VariableDatabase
   * and get the descriptor indices for those variables.
   */

  tbox::Pointer<pdat::CellVariable<NDIM,double> > comp_soln =
     new pdat::CellVariable<NDIM,double>(object_name+":computed solution",1);
  d_comp_soln_id
     = vdb->registerVariableAndContext(
     comp_soln ,
     d_context ,
     hier::IntVector<NDIM>(1) /* ghost cell width is 1 for stencil widths */ );

  tbox::Pointer<pdat::CellVariable<NDIM,double> > exact_solution =
     new pdat::CellVariable<NDIM,double>(object_name+":exact solution");
  d_exact_id
     = vdb->registerVariableAndContext(
     exact_solution ,
     d_context,
     hier::IntVector<NDIM>(1) /* ghost cell width is 1 in case needed */ );

  tbox::Pointer<pdat::CellVariable<NDIM,double> > rhs_variable =
     new pdat::CellVariable<NDIM,double>(object_name+":linear system right hand side");
  d_rhs_id
     = vdb->registerVariableAndContext(
     rhs_variable ,
     d_context,
     hier::IntVector<NDIM>(0) /* ghost cell width is 0 */ );

  /*
   * Specify an implementation of solv::RobinBcCoefStrategy<NDIM> for the solver to use.
   * We use the implementation solv::LocationIndexRobinBcCoefs<NDIM>, but other
   * implementations are possible, including user-implemented.
   */
  d_poisson_fac_solver.setBcObject( &d_bc_coefs );

  return;
}



/*
*************************************************************************
* Destructor does nothing interesting                                   *
*************************************************************************
*/
FACPoisson::~FACPoisson() 
{
} 





/*
*************************************************************************
* Initialize data on a level.                                           *
*                                                                       *
* Allocate the solution, exact solution and rhs memory.                 *
* Fill the rhs and exact solution.                                      *
*************************************************************************
*/
void FACPoisson::initializeLevelData (
   const tbox::Pointer<hier::BasePatchHierarchy<NDIM> > hierarchy_ , 
   const int level_number ,
   const double init_data_time ,
   const bool can_be_refined ,
   const bool initial_time ,
   const tbox::Pointer<hier::BasePatchLevel<NDIM> > old_level ,
   const bool allocate_data )
{

   tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy = hierarchy_;
   tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geom = hierarchy->getGridGeometry();

   tbox::Pointer<hier::PatchLevel<NDIM> > level
      = hierarchy->getPatchLevel(level_number);

   if ( allocate_data ) {
      level->allocatePatchData( d_comp_soln_id );
      level->allocatePatchData( d_rhs_id );
      level->allocatePatchData( d_exact_id );
   }

   /*
    * Initialize data in all patches in the level.
    */
   hier::PatchLevel<NDIM> ::Iterator pi(*level);
   for ( pi.initialize(*level); pi; pi++ ) {

      const int pn = *pi;
      tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(pn);
      if( patch.isNull() ) {
         TBOX_ERROR(d_object_name
                    << ": Cannot find patch.  Null patch pointer.");
      }
      hier::Box<NDIM> pbox = patch->getBox();
      tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > patch_geom
         = patch->getPatchGeometry();

      tbox::Pointer<pdat::CellData<NDIM,double> > exact_data
         = patch->getPatchData(d_exact_id);
      tbox::Pointer<pdat::CellData<NDIM,double> > rhs_data
         = patch->getPatchData(d_rhs_id);

      /*
       * Set source function and exact solution.
       */
#if NDIM == 2
      setexactandrhs_( pbox.lower()[0],
                       pbox.upper()[0],
                       pbox.lower()[1],
                       pbox.upper()[1],
                       exact_data->getPointer(),
                       rhs_data->getPointer(),
                       patch_geom->getDx(),
                       patch_geom->getXLower() );
#elif NDIM == 3
      setexactandrhs_( pbox.lower()[0],
                       pbox.upper()[0],
                       pbox.lower()[1],
                       pbox.upper()[1],
                       pbox.lower()[2],
                       pbox.upper()[2],
                       exact_data->getPointer(),
                       rhs_data->getPointer(),
                       patch_geom->getDx(),
                       patch_geom->getXLower() );
#endif

   }	// End patch loop.
}




/*
*************************************************************************
* Reset the hierarchy-dependent internal information.                   *
*************************************************************************
*/
void FACPoisson::resetHierarchyConfiguration (
   tbox::Pointer<hier::BasePatchHierarchy<NDIM> > new_hierarchy ,
   int coarsest_level ,
   int finest_level )
{
   d_hierarchy = new_hierarchy;
   return;
}





/*
*************************************************************************
* Set up the initial guess and problem parameters                       *
* and solve the Poisson problem.  We explicitly initialize and          *
* deallocate the solver state in this example.                          *
*************************************************************************
*/
int FACPoisson::solvePoisson()
{

  if ( d_hierarchy.isNull() ) {
    TBOX_ERROR(d_object_name
               << "Cannot solve using an uninitialized object.\n");
  }

  int ln;
  /*
   * Fill in the initial guess.
   */
  for ( ln=0; ln<=d_hierarchy->getFinestLevelNumber(); ++ln ) {
    tbox::Pointer<hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
    hier::PatchLevel<NDIM>::Iterator ip(*level);
    for ( ; ip; ip++ ) {
      tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(*ip);
      tbox::Pointer<pdat::CellData<NDIM,double> > data = patch->getPatchData(d_comp_soln_id);
      data->fill(0.0);
    }
  }

  /*
   * Set the parameters for the Poisson equation.
   * See classes solv::CellPoissonFACSolver<NDIM> or
   * solv::PoissonSpecifications.
   * (D is the diffusion coefficient.
   * C is the source term which is not used in this example.)
   */
  d_poisson_fac_solver.setDConstant(1.0);
  d_poisson_fac_solver.setCConstant(0.0);

  d_poisson_fac_solver.initializeSolverState(
    d_comp_soln_id,
    d_rhs_id,
    d_hierarchy,
    0,
    d_hierarchy->getFinestLevelNumber() );

  tbox::pout << "solving..." << std::endl;
  int solver_ret;
  solver_ret = d_poisson_fac_solver.solveSystem( d_comp_soln_id ,
                                                 d_rhs_id );
  /*
   * Present data on the solve.
   */
  double avg_factor, final_factor;
  d_poisson_fac_solver.getConvergenceFactors(avg_factor, final_factor);
  tbox::pout << "\t" << (solver_ret?"":"NOT ") << "converged " << "\n"
       << "	iterations: " << d_poisson_fac_solver.getNumberOfIterations() << "\n"
       << "	residual: " << d_poisson_fac_solver.getResidualNorm() << "\n"
       << "	average convergence: " << avg_factor << "\n"
       << "	final convergence: " << final_factor << "\n"
       << std::flush;


  d_poisson_fac_solver.deallocateSolverState();

  return 0;
}






/*
*************************************************************************
* Set up external plotter to plot internal data from this class.        *
* Tell the plotter about the refinement ratios.  Register variables     *
* appropriate for plotting.                                             *
*************************************************************************
*/
int FACPoisson::setupPlotter(
  appu::CartesianVizamraiDataWriter<NDIM> &plotter
) const {
   if ( d_hierarchy.isNull() ) {
      TBOX_ERROR(d_object_name << ": No hierarchy in\n"
                 << " FACPoisson::setupPlotter\n"
                 << "The hierarchy must be set before calling\n"
                 << "this function.\n");
   }

   /*
    * This must be done once (and again each time the data changes).
    */
   int ln;
   for ( ln=1; ln<d_hierarchy->getNumberOfLevels(); ln++ ) {
      tbox::Pointer<hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
      const hier::IntVector<NDIM> &lratio = level->getRatioToCoarserLevel();
      plotter.setRatioToCoarserLevel(ln, lratio);
   }
   plotter.registerPlotScalar("Computed solution", d_comp_soln_id);
   plotter.registerDerivedPlotScalar("Error" ,
                                        (appu::VisDerivedDataStrategy<NDIM> *)this );
   plotter.registerPlotScalar("Exact solution", d_exact_id);
   plotter.registerPlotScalar("Poisson source", d_rhs_id);

   return 0;
}






/*
*************************************************************************
* Set up external plotter to plot internal data from this class.        *
* Register variables appropriate for plotting.                          *
*************************************************************************
*/
int FACPoisson::setupPlotter(
  appu::VisItDataWriter<NDIM> &plotter
) const {
   if ( d_hierarchy.isNull() ) {
      TBOX_ERROR(d_object_name << ": No hierarchy in\n"
                 << " FACPoisson::setupPlotter\n"
                 << "The hierarchy must be set before calling\n"
                 << "this function.\n");
   }
   plotter.registerPlotQuantity("Computed solution",
                                "SCALAR",
                                d_comp_soln_id);
   plotter.registerDerivedPlotQuantity("Error" ,
                                       "SCALAR",
                                       (appu::VisDerivedDataStrategy<NDIM> *)this );
   plotter.registerPlotQuantity("Exact solution",
                                "SCALAR",
                                d_exact_id);
   plotter.registerPlotQuantity("Poisson source",
                                "SCALAR",
                                d_rhs_id);

   return 0;
}



/*
*************************************************************************
* Write derived data to the given stream.                               *
*************************************************************************
*/
bool FACPoisson::packDerivedDataIntoDoubleBuffer(
   double* buffer,
   const hier::Patch<NDIM> &patch ,
   const hier::Box<NDIM> &region ,
   const std::string &variable_name ,
   int depth_id) const
{


   pdat::CellData<NDIM,double>::Iterator icell( region );

   if ( variable_name == "Error" ) {
      tbox::Pointer<pdat::CellData<NDIM,double> > current_solution_
        = patch.getPatchData(d_comp_soln_id);
      tbox::Pointer<pdat::CellData<NDIM,double> > exact_solution_
        = patch.getPatchData(d_exact_id);
      pdat::CellData<NDIM,double> &current_solution = *current_solution_;
      pdat::CellData<NDIM,double> &exact_solution = *exact_solution_;
      for ( ; icell; icell++ ) {
         double diff = ( current_solution(*icell) - exact_solution(*icell) );
         *buffer = diff;
         buffer = buffer + 1;
      }
   }
   else {
      // Did not register this name.
      TBOX_ERROR("Unregistered variable name '" << variable_name << "' in\n"
                 << "FACPoissonX::packDerivedDataIntoDoubleBuffer");

   }
   // Return true if this patch has derived data on it. 
   // False otherwise.
   return(true);
}

}
