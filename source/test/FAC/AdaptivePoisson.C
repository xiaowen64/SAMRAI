/*
  File:		$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/test/FAC/AdaptivePoisson.C $
  Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
  Revision:	$LastChangedRevision: 1704 $
  Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
  Description:	AdaptivePoisson class implementation
*/

#include "SAMRAI_config.h"

#include "printObject.h"
#include "MDA_Access.h"
#include "arrayConversion.h"
#include "patchFcns.h"
#include "AdaptivePoisson.h"
#include "CellPoissonFACOps.h"

#include "tbox/Pointer.h"
#include "tbox/ConstPointer.h"
#include "tbox/Array.h"
#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"
#include "tbox/Database.h"
#include "tbox/InputDatabase.h"

#include "Patch.h"
#include "VariableDatabase.h"
#include "Variable.h"
#include "CartesianPatchGeometry.h"
#include "CartesianGridGeometry.h"
#include "HierarchyCellDataOpsReal.h"
#include "PatchCellDataOpsReal.h"
#include "CellData.h"
#include "CellIndex.h"
#include "CellVariable.h"
#include "SideData.h"
#include "OutersideData.h"
#include "Geometry.h"
#include "CoarsenAlgorithm.h"
#include "CoarsenOperator.h"
#include "CoarsenSchedule.h"
#include "RefineAlgorithm.h"
#include "RefineOperator.h"
#include "RefineSchedule.h"
#include "Index.h"
#include "CartesianCellDoubleLinearRefine.h"
#include "CartesianCellDoubleConservativeLinearRefine.h"
#include "CartesianCellDoubleWeightedAverage.h"
#include "CartesianSideDoubleWeightedAverage.h"
#include "CellDoubleConstantRefine.h"
#include "PatchCellDataOpsReal.h"
#include "HierarchyCellDataOpsReal.h"

#include IOMANIP_HEADER_FILE

using namespace SAMRAI;

AdaptivePoisson::AdaptivePoisson(
  const string &object_name
, tbox::Database &database
, /*! Standard output stream */ ostream *out_stream
, /*! Log output stream */ ostream *log_stream
)
: d_name(object_name) ,
  d_hierarchy() ,
  d_fac_ops( object_name+":scalar poisson operator" ,
	     database.getDatabase("ScalarPoissonOps") ) ,
  d_fac_preconditioner( "FAC preconditioner for Poisson's equation" ,
			d_fac_ops ) ,
  d_context_persistent("PERSISTENT") ,
  d_context_scratch("SCRATCH") ,
  d_diffcoef("solution:diffcoef",1) ,
  d_flux("flux",1) ,
  d_scalar("solution:scalar",1) ,
  d_constant_source("poisson source",1) ,
  d_ccoef("linear source coefficient",1) ,
  d_rhs("linear system rhs",1) ,
  d_exact("solution:exact",1) ,
  d_resid( new pdat::CellVariable<NDIM,double>(object_name+"residual") ) ,
  d_weight("vector weight",1) ,
  d_ostream(out_stream) ,
  d_lstream(log_stream) ,
  d_problem_name("sine") ,
  d_sps( object_name+"Poisson solver specifications" ),
  d_robin_refine_patch( object_name+"Refine patch implementation" ),
  d_physical_bc_coef(NULL),
  d_adaption_threshold(0.5),
  d_finest_plot_level(9999999) ,
  d_visit_writer() ,
  d_finest_dbg_plot_ln(database.getIntegerWithDefault("finest_dbg_plot_ln",99))
{

  /*
    Register variables with hier::VariableDatabase
    and get the descriptor indices for those variables.
    It is not necessary to save the indices obtained from
    the registration, because they can always be retrieved
    from the mapVariableAndContextToIndex, but we do it
    because we refer to the indices often.
  */
  {
    tbox::Pointer<hier::VariableContext>
      context_persistent_ptr(&d_context_persistent,false);
    tbox::Pointer<hier::VariableContext>
      context_scratch_ptr(&d_context_scratch,false);
    hier::VariableDatabase<NDIM> *variable_db = hier::VariableDatabase<NDIM>::getDatabase();

    /*
      Persistent data.
    */
    d_diffcoef_persistent
      = variable_db->registerVariableAndContext(
	tbox::Pointer< hier::Variable<NDIM> >(&d_diffcoef,false)
	, context_persistent_ptr
	);
    d_flux_persistent
      = variable_db->registerVariableAndContext(
	tbox::Pointer< hier::Variable<NDIM> >(&d_flux,false)
	, context_persistent_ptr
	);
    d_scalar_persistent
      = variable_db->registerVariableAndContext(
	tbox::Pointer< hier::Variable<NDIM> >(&d_scalar,false)
	, context_persistent_ptr
	, hier::IntVector<NDIM>(1) /* ghost cell width is 1 for stencil widths */
	);
    d_constant_source_persistent
      = variable_db->registerVariableAndContext(
	tbox::Pointer< hier::Variable<NDIM> >(&d_constant_source,false)
	, context_persistent_ptr
	);
    d_ccoef_persistent
      = variable_db->registerVariableAndContext(
	tbox::Pointer< hier::Variable<NDIM> >(&d_ccoef,false)
	, context_persistent_ptr
	);
    d_exact_persistent
      = variable_db->registerVariableAndContext(
	tbox::Pointer< hier::Variable<NDIM> >(&d_exact,false)
	, context_persistent_ptr
	);
    d_weight_persistent
      = variable_db->registerVariableAndContext(
	tbox::Pointer< hier::Variable<NDIM> >(&d_weight,false)
	, context_persistent_ptr
	);
    /*
      Scratch data.
    */
    d_rhs_scratch
      = variable_db->registerVariableAndContext(
	tbox::Pointer< hier::Variable<NDIM> >(&d_rhs,false)
	, context_scratch_ptr
	, hier::IntVector<NDIM>(0) /* ghost cell width is 0 */
	);
    d_resid_scratch
      = variable_db->registerVariableAndContext(
	d_resid ,
	context_scratch_ptr ,
	hier::IntVector<NDIM>(0) /* ghost cell width is 0 */
	);
  }

  d_finest_plot_level
    = database.getIntegerWithDefault( "finest_plot_level" ,
				      d_finest_plot_level );

  /*
    Experiment with algorithm choices in solv::FACPreconditioner<NDIM>.
  */
  string fac_algo = database.getStringWithDefault( "fac_algo", "default" );
  d_fac_preconditioner.setAlgorithmChoice(fac_algo);

  d_adaption_threshold
    = database.getDoubleWithDefault( "adaption_threshold" ,
				     d_adaption_threshold );

  /*
    Read in the possible solution-specific objects.
  */
  {
    if ( database.isDatabase( "sine_solution" ) ) {
      tbox::Pointer<tbox::Database> db =
	database.getDatabase( "sine_solution" );
      d_sine_solution.setFromDatabase( *db );
    }
    if ( database.isDatabase( "gaussian_solution" ) ) {
      tbox::Pointer<tbox::Database> db =
	database.getDatabase( "gaussian_solution" );
      d_gaussian_solution.setFromDatabase( *db );
    }
    if ( database.isDatabase( "multigaussian_solution" ) ) {
      tbox::Pointer<tbox::Database> db =
	database.getDatabase( "multigaussian_solution" );
      d_multigaussian_solution.setFromDatabase( *db );
    }
    if ( database.isDatabase( "polynomial_solution" ) ) {
      tbox::Pointer<tbox::Database> db =
	database.getDatabase( "polynomial_solution" );
      d_polynomial_solution.setFromDatabase( *db );
    }
    if ( database.isDatabase( "gaussian_diffcoef_solution" ) ) {
      tbox::Pointer<tbox::Database> db =
	database.getDatabase( "gaussian_diffcoef_solution" );
      d_gaussian_diffcoef_solution.setFromDatabase( *db );
    }
  }

  /*
    We are set up with a choice of problems to solve.
    The problem name identifies the specific one.
  */

  d_problem_name =
    database.getStringWithDefault("problem_name", d_problem_name);
  if ( d_problem_name != "sine"
       && d_problem_name != "sine-neumann"
       && d_problem_name != "gauss"
       && d_problem_name != "multigauss"
       && d_problem_name != "pernice"
       && d_problem_name != "poly"
       && d_problem_name != "gauss-coef"
       ) {
    TBOX_ERROR("Unrecognized problem name " << d_problem_name << "\n");
  }
  if ( d_problem_name == "sine" ) {
    d_sine_solution.setPoissonSpecifications(d_sps,
					     d_ccoef_persistent,
					     d_diffcoef_persistent);
    d_physical_bc_coef = &d_sine_solution;
  }
  else if ( d_problem_name == "gauss" ) {
    d_gaussian_solution.setPoissonSpecifications(d_sps,
						 d_ccoef_persistent,
						 d_diffcoef_persistent);
    d_physical_bc_coef = &d_gaussian_solution;
  }
  else if ( d_problem_name == "multigauss" ) {
    d_multigaussian_solution.setPoissonSpecifications(d_sps,
						      d_ccoef_persistent,
						      d_diffcoef_persistent);
    d_physical_bc_coef = &d_multigaussian_solution;
  }
  else if ( d_problem_name == "poly" ) {
    d_polynomial_solution.setPoissonSpecifications(d_sps,
						   d_ccoef_persistent,
						   d_diffcoef_persistent);
    d_physical_bc_coef = &d_polynomial_solution;
  }
  else if ( d_problem_name == "gauss-coef" ) {
    d_gaussian_diffcoef_solution.setPoissonSpecifications(d_sps,
							  d_ccoef_persistent,
							  d_diffcoef_persistent);
    d_physical_bc_coef = &d_gaussian_diffcoef_solution;
  }
  else {
    TBOX_ERROR("Unidentified problem name");
  }
  /*
    Tell ScalarPoissonOperator where to find some of the data
    we are providing it.
  */
  d_fac_ops.setPoissonSpecifications( d_sps );
  d_fac_ops.setFluxId( -1 );

  d_fac_ops.setPhysicalBcCoefObject(d_physical_bc_coef);

  d_fac_ops.setProlongationMethod("LINEAR_REFINE");

  tbox::plog << "Gaussian solution parameters:\n"
       << d_gaussian_solution << "\n\n" << endl;
#if 0
  tbox::plog << "Sine solution parameters:\n"
       << d_sine_solution << "\n\n" << endl;
  tbox::plog << "Polynomial solution parameters:\n"
       << d_polynomial_solution << "\n\n" << endl;
  tbox::plog << "Gaussian diffcoef solution parameters:\n"
       << d_gaussian_diffcoef_solution << "\n\n" << endl;
#endif
  tbox::plog << "Problem name is: " << d_problem_name << "\n\n" << endl;

  d_fac_ops.setPreconditioner( &d_fac_preconditioner );

  return;
}

void AdaptivePoisson::initializeLevelData (
  /*! Hierarchy to initialize */
  const tbox::Pointer< hier::BasePatchHierarchy<NDIM> > hierarchy ,
  /*! Level to initialize */
  const int ln ,
  const double init_data_time ,
  const bool can_be_refined ,
  /*! Whether level is being introduced for the first time */
  const bool initial_time ,
  /*! Level to copy data from */
  const tbox::Pointer< hier::BasePatchLevel<NDIM> > old_level ,
  const bool allocate_data )
{
  tbox::Pointer< hier::PatchHierarchy<NDIM> > patch_hierarchy = hierarchy;

  /*
    Reference the level object with the given index from the hierarchy.
  */
  tbox::Pointer< hier::PatchLevel<NDIM> > level
    = hierarchy->getPatchLevel(ln);

  /*
    If instructed, allocate all patch data on the level.
    Allocate only persistent data.  Scratch data will
    generally be allocated and deallocated as needed.
  */
  if ( allocate_data ) {
    if ( d_sps.dIsVariable() )
      level->allocatePatchData( d_diffcoef_persistent );
    level->allocatePatchData( d_flux_persistent );
    level->allocatePatchData( d_scalar_persistent );
    level->allocatePatchData( d_constant_source_persistent );
    if ( d_sps.cIsVariable() )
      level->allocatePatchData( d_ccoef_persistent );
    level->allocatePatchData( d_exact_persistent );
    level->allocatePatchData( d_weight_persistent );
  }

#if 0
  /*
    For debugging purpose, store the array of data for the
    patch 0 of levels 0 and 1
  */
  switch (ln) {
  case 0:
      hier::Patch<NDIM> &Patch0 = *level->getPatch(0);
      tbox::Pointer<pdat::CellData<NDIM,double> >
        Soln0_data = solution.getComponentPatchData ( 0 , Patch0 );
      MDA_Access<double,NDIM,MDA_OrderColMajor<NDIM> >
	d_array[0] = arrayData2ArrayAccess( Soln0_data->getArrayData() );
    break;
  case 1:
      hier::Patch<NDIM> &Patch1 = *level->getPatch(0);
      tbox::Pointer<pdat::CellData<NDIM,double> >
        Soln1_data = solution.getComponentPatchData ( 0 , Patch1 );
      MDA_Access<double,NDIM,MDA_OrderColMajor<NDIM> >
	d_array[1] = arrayData2ArrayAccess( Soln1_data->getArrayData() );
    break;
  }
#endif

  /*
    Initialize data in all patches in the level.
  */
  hier::PatchLevel<NDIM>::Iterator pi;
  for ( pi.initialize(*level); pi; pi++ ) {
    int pn = *pi;

    hier::Patch<NDIM> &patch = *(level->getPatch(pn));
    hier::Box<NDIM> pbox = patch.getBox();
    tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > patch_geom
      = patch.getPatchGeometry();

    tbox::Pointer<pdat::SideData<NDIM,double> > diffcoef_data
      = patch.getPatchData(d_diffcoef_persistent);
    tbox::Pointer<pdat::CellData<NDIM,double> > scalar_data
      = patch.getPatchData(d_scalar_persistent);
    tbox::Pointer<pdat::CellData<NDIM,double> > ccoef_data
      = patch.getPatchData(d_ccoef_persistent);
    tbox::Pointer<pdat::CellData<NDIM,double> > exact_data
      = patch.getPatchData(d_exact_persistent);
    tbox::Pointer<pdat::CellData<NDIM,double> > source_data
      = patch.getPatchData(d_constant_source_persistent);

    /* Set source function and exact solution. */
    if ( d_problem_name == "sine" ) {
      d_sine_solution.setGridData( patch ,
				   *diffcoef_data ,
				   *ccoef_data ,
				   *exact_data ,
				   *source_data );
    }
    else if ( d_problem_name == "gauss" ) {
      d_gaussian_solution.setGridData( patch ,
				    *diffcoef_data ,
				    *ccoef_data ,
				    *exact_data ,
				    *source_data );
    }
    else if ( d_problem_name == "multigauss" ) {
      d_multigaussian_solution.setGridData( patch ,
					    *diffcoef_data ,
					    *ccoef_data ,
					    *exact_data ,
					    *source_data );
    }
    else if ( d_problem_name == "poly" ) {
      d_polynomial_solution.setGridData( patch ,
					 *diffcoef_data ,
					 *ccoef_data ,
					 *exact_data ,
					 *source_data );
    }
    else if ( d_problem_name == "gauss-coef" ) {
      d_gaussian_diffcoef_solution.setGridData( patch ,
						*diffcoef_data ,
						*ccoef_data ,
						*exact_data ,
						*source_data );
    }
    else {
      TBOX_ERROR("Unidentified problem name");
    }

  }

  /*
    Refine solution data from coarser level and, if provided, old level.
  */
  int oldln=-1;
  if ( old_level ) oldln = old_level->getLevelNumber();
  {
    xfer::RefineAlgorithm<NDIM> refiner;
    tbox::Pointer<geom::CartesianGridGeometry<NDIM> >
      grid_geometry_ = patch_hierarchy->getGridGeometry();
    geom::CartesianGridGeometry<NDIM> &grid_geometry = *grid_geometry_;
    tbox::Pointer< xfer::RefineOperator<NDIM> > accurate_refine_op
      = grid_geometry.
      lookupRefineOperator( tbox::Pointer< hier::Variable<NDIM> >(&d_scalar,false),
			    "CONSERVATIVE_LINEAR_REFINE" );
    TBOX_ASSERT( accurate_refine_op );
    refiner.registerRefine( d_scalar_persistent ,
			    d_scalar_persistent ,
			    d_scalar_persistent ,
			    accurate_refine_op );
    tbox::Pointer< xfer::RefineSchedule<NDIM> > refine_schedule;
    if ( ln > 0 ) {
      /*
	Include coarser levels in setting data
      */
      refine_schedule =
	refiner.createSchedule( level ,
				old_level ,
				ln-1 ,
				hierarchy ,
				&d_robin_refine_patch );
    }
    else {
      /*
	There is no coarser level, and source data comes only
	from old_level, if any.
      */
      if ( old_level ) {
	refine_schedule =
	  refiner.createSchedule( level ,
				  old_level ,
				  &d_robin_refine_patch );
      }
    }
    if ( refine_schedule ) {
      d_robin_refine_patch.setCoefImplementation(d_physical_bc_coef);
      d_robin_refine_patch.setTargetDataId(d_scalar_persistent);
      d_robin_refine_patch.setHomogeneousBc(false);
      refine_schedule->fillData(0.0);
      // It is null if this is the bottom level.
    }
    else {
      math::HierarchyCellDataOpsReal<NDIM,double> hcellmath(hierarchy,ln,ln);
      hcellmath.setToScalar( d_scalar_persistent, 0.0, false );
    }
    if(0) {
      // begin debug code
      math::HierarchyCellDataOpsReal<NDIM,double> hcellmath(hierarchy);
      hcellmath.printData( d_scalar_persistent, tbox::pout, false );
      // end debug code
    }
  }

  /* Set vector weight. */
  d_fac_ops.computeVectorWeights( hierarchy, d_weight_persistent );
}

void AdaptivePoisson::resetHierarchyConfiguration (
  /*! New hierarchy */ tbox::Pointer< hier::BasePatchHierarchy<NDIM> > new_hierarchy ,
  /*! Coarsest level */ int coarsest_level ,
  /*! Finest level */ int finest_level )
{
  d_hierarchy = new_hierarchy;
  /*
    Recompute or reset internal data tied to the hierarchy,
    if any.  None at this time.
  */
  /*
    Log the new hierarchy.
  */
  if ( d_lstream ) {
    *d_lstream
      << "AdaptivePoisson<NDIM>::resetHierarchyConfiguration\n";
    d_hierarchy->recursivePrint( *d_lstream, "    ", 2 );
  }
  return;
}

void AdaptivePoisson::applyGradientDetector(
  const tbox::Pointer< hier::BasePatchHierarchy<NDIM> > hierarchy_,
  const int ln,
  const double error_data_time,
  const int tag_index,
  const bool initial_time,
  const bool uses_richardson_extrapolation )
{
  if ( d_lstream ) {
    *d_lstream
      << "AdaptivePoisson<NDIM>(" << d_name << ")::applyGradientDetector"
      << endl;
  }
  const tbox::Pointer< hier::PatchHierarchy<NDIM> > hierarchy__ = hierarchy_;
  hier::PatchHierarchy<NDIM> &hierarchy = *hierarchy__;
  tbox::Pointer<geom::CartesianGridGeometry<NDIM> >
    grid_geometry_ = hierarchy.getGridGeometry();
  hier::PatchLevel<NDIM> &level 
     = (hier::PatchLevel<NDIM>&) *hierarchy.getPatchLevel(ln);
  hier::PatchLevel<NDIM>::Iterator pi;
  int ntag=0, ntotal=0;
  double maxestimate=0;
  for ( pi.initialize(level); pi; pi++ ) {
    int pn = pi();
    hier::Patch<NDIM> &patch = *level.getPatch(pn);
    tbox::Pointer< hier::PatchData<NDIM> >
      tag_data = patch.getPatchData( tag_index );
    ntotal += patch.getBox().numberCells().getProduct();
    if ( tag_data.isNull() ) {
      TBOX_ERROR("Data index " << tag_index << " does not exist for patch.\n");
    }
    tbox::Pointer<pdat::CellData<NDIM,int> > tag_cell_data_ = tag_data;
    if ( tag_cell_data_.isNull() ) {
      TBOX_ERROR("Data index " << tag_index << " is not cell int data.\n");
    }
    tbox::Pointer< hier::PatchData<NDIM> >
      soln_data = patch.getPatchData( d_scalar_persistent );
    if ( soln_data.isNull() ) {
      TBOX_ERROR("Data index " << d_scalar_persistent
		 << " does not exist for patch.\n");
    }
    tbox::Pointer<pdat::CellData<NDIM,double> > soln_cell_data_ = soln_data;
    if ( soln_cell_data_.isNull() ) {
      TBOX_ERROR("Data index " << d_scalar_persistent
		 << " is not cell int data.\n");
    }
    pdat::CellData<NDIM,double> &soln_cell_data = *soln_cell_data_;
    pdat::CellData<NDIM,int> &tag_cell_data = *tag_cell_data_;
    pdat::CellData<NDIM,double> estimate_data( patch.getBox(),
					  1,
					  hier::IntVector<NDIM>(0) );
    computeAdaptionEstimate( estimate_data,
			     soln_cell_data);
    tag_cell_data.fill(0);
    hier::Box<NDIM>::Iterator i;
    for ( i.initialize(patch.getBox()); i; i++ ) {
      if ( maxestimate < estimate_data(*i) ) maxestimate = estimate_data(*i);
      if ( estimate_data(*i) > d_adaption_threshold ) {
	tag_cell_data(*i) = 1;
	++ntag;
      }
    }
  }
  tbox::plog << "Adaption threshold is " << d_adaption_threshold << "\n";
  tbox::plog << "Number of cells tagged on level " << ln << " is "
       << ntag << "/" << ntotal << "\n";
  tbox::plog << "Max estimate is " << maxestimate << "\n";

  return;
}

int AdaptivePoisson::registerVariablesWithPlotter(
  appu::CartesianVizamraiDataWriter<NDIM> &viz_writer
) const {

   /*
     Tell the plotter about the refinement ratios.
     This must be done once (and again each time the data changes).
   */
   if ( d_hierarchy.isNull() ) {
      TBOX_ERROR("No hierarchy in AdaptivePoisson::registerVariablesWithPlotter\n"
                 << "The hierarchy must be built before calling this function.\n");
   }
   int ln;
   for ( ln=1; ln<d_hierarchy->getNumberLevels(); ln++ ) {
     tbox::Pointer<hier::PatchLevel<NDIM> > level 
	= d_hierarchy->getPatchLevel(ln);
     const hier::IntVector<NDIM> &lratio = level->getRatioToCoarserLevel();
      viz_writer.setRatioToCoarserLevel(ln, lratio);
   }
    /*
      Set the finest level to plot (optional).
    */
   int finest_plot_level =
     tbox::MathUtilities<int>::Min( d_finest_plot_level,
			            d_hierarchy->getFinestLevelNumber() );
   viz_writer.setFinestLevelToPlot( finest_plot_level );
   /*
     Register variables with plotter.
   */
   viz_writer.registerPlotScalar("Computed solution", d_scalar_persistent);
   viz_writer.registerDerivedPlotScalar("Error" ,
     (appu::VisDerivedDataStrategy<NDIM> *)this );
   viz_writer.registerPlotScalar("Exact solution", d_exact_persistent);
   if(1) viz_writer.registerDerivedPlotScalar("Gradient Function" ,
     (appu::VisDerivedDataStrategy<NDIM> *)this );
   viz_writer.registerPlotScalar("Poisson source", d_constant_source_persistent);
   if(1) viz_writer.registerDerivedPlotScalar("Patch level number" ,
     (appu::VisDerivedDataStrategy<NDIM> *)this );

   return 0;
}

int AdaptivePoisson::registerVariablesWithPlotter(
  appu::VisItDataWriter<NDIM> &visit_writer
) {

   visit_writer.registerPlotQuantity("Computed solution",
                                     "SCALAR",
                                     d_scalar_persistent);
   visit_writer.registerPlotQuantity("Exact solution",
                                     "SCALAR",
                                     d_exact_persistent);
   visit_writer.registerPlotQuantity("Poisson source",
                                     "SCALAR",
                                     d_constant_source_persistent);
   visit_writer.registerDerivedPlotQuantity("Error" ,
                                            "SCALAR",
                                            this,
                                            1.0,
                                            "CELL");
   visit_writer.registerDerivedPlotQuantity("Gradient Function" ,
                                            "SCALAR",
                                            this,
                                            1.0,
                                            "CELL" );
   visit_writer.registerDerivedPlotQuantity("Patch level number" ,
                                          "SCALAR",
                                            this,
                                            1.0,
                                            "CELL" );

#if 0
   /*
     Tell the plotter about the refinement ratios.
     This must be done once (and again each time the data changes).
   */
   if ( d_hierarchy.isNull() ) {
      TBOX_ERROR("No hierarchy in AdaptivePoisson::registerVariablesWithPlotter\n"
                 << "The hierarchy must be built before calling this function.\n");
   }
   int ln;
   for ( ln=1; ln<d_hierarchy->getNumberLevels(); ln++ ) {
     tbox::Pointer<hier::PatchLevel<NDIM> > level 
	= d_hierarchy->getPatchLevel(ln);
     const hier::IntVector<NDIM> &lratio = level->getRatioToCoarserLevel();
      visit_writer.setRatioToCoarserLevel(ln, lratio);
   }
    /*
      Set the finest level to plot (optional).
    */
   int finest_plot_level =
     tbox::MathUtilities<int>::Min( d_finest_plot_level,
			            d_hierarchy->getFinestLevelNumber() );
   visit_writer.setFinestLevelToPlot( finest_plot_level );
   /*
     Register variables with plotter.
   */
   visit_writer.registerPlotScalar("Computed solution", d_scalar_persistent);
   visit_writer.registerDerivedPlotScalar("Error" ,
     (appu::VisDerivedDataStrategy<NDIM> *)this );
   visit_writer.registerPlotScalar("Exact solution", d_exact_persistent);
   if(1) visit_writer.registerDerivedPlotScalar("Gradient Function" ,
     (appu::VisDerivedDataStrategy<NDIM> *)this );
   visit_writer.registerPlotScalar("Poisson source", d_constant_source_persistent);
   if(1) visit_writer.registerDerivedPlotScalar("Patch level number" ,
     (appu::VisDerivedDataStrategy<NDIM> *)this );
#endif

   return 0;
}

bool AdaptivePoisson::packDerivedDataIntoDoubleBuffer(
  double* buffer ,
  const hier::Patch<NDIM> &patch ,
  const hier::Box<NDIM> &region ,
  const string &variable_name ,
  int depth_id) const
{

  // begin debug code
  // math::HierarchyCellDataOpsReal<NDIM,double> hcellmath(d_hierarchy);
  // hcellmath.printData( d_exact_persistent, pout, false );
  // end debug code

  const int *lower = region.lower();
  const int *upper = region.upper();

  if ( variable_name == "Error" ) {
    tbox::Pointer<pdat::CellData<NDIM,double> > current_solution
      = patch.getPatchData(d_scalar_persistent);
    tbox::Pointer<pdat::CellData<NDIM,double> > exact_solution
      = patch.getPatchData(d_exact_persistent);
    MDA_AccessConst<double,NDIM,MDA_OrderColMajor<NDIM> > ex =
      arrayData2ArrayAccess( exact_solution->getArrayData() );
    MDA_AccessConst<double,NDIM,MDA_OrderColMajor<NDIM> > co =
      arrayData2ArrayAccess( current_solution->getArrayData() );
#if NDIM == 2
    for ( int j=lower[1]; j<=upper[1]; ++j ) {
      for ( int i=lower[0]; i<=upper[0]; ++i ) {
	*(buffer++) = co(i,j) - ex(i,j);
      }
    }
#endif
#if NDIM == 3
    for ( int k=lower[2]; k<=upper[2]; ++k ) {
      for ( int j=lower[1]; j<=upper[1]; ++j ) {
	for ( int i=lower[0]; i<=upper[0]; ++i ) {
	  *(buffer++) = co(i,j,k) - ex(i,j,k);
	}
      }
    }
#endif
  }
  else if ( variable_name == "Gradient Function" ) {
    tbox::Pointer<pdat::CellData<NDIM,double> > soln_cell_data_ =
      patch.getPatchData(d_scalar_persistent);
    const pdat::CellData<NDIM,double> &soln_cell_data = *soln_cell_data_;
    pdat::CellData<NDIM,double> estimate_data( region,
					  1,
					  hier::IntVector<NDIM>(0) );
    computeAdaptionEstimate( estimate_data,
			     soln_cell_data);
    // tbox::plog << "estimate data: " << patch.getBox().size() << "\n";
    // estimate_data.print(region,0,tbox::plog);
    memcpy(buffer, estimate_data.getPointer(), sizeof(double)*region.size());
    // stream.pack( estimate_data.getPointer(), region.size() );
  }
  else if ( variable_name == "Patch level number" ) {
    double pln = patch.getPatchLevelNumber();
    int i, size=region.size();
    for ( i=0; i<size; ++i ) buffer[i] = pln;
  }
  else {
    // Did not register this name.
    TBOX_ERROR("Unregistered variable name '" << variable_name << "' in\n"
	       << "AdaptivePoisson<NDIM>::packDerivedPatchDataIntoDoubleBuffer");
  }

  // Return TRUE if this patch has derived data on it. 
  // FALSE otherwise.
  return(true);
}

void AdaptivePoisson::computeAdaptionEstimate(
  pdat::CellData<NDIM,double> &estimate_data,
  const pdat::CellData<NDIM,double> &soln_cell_data) const
{
  const int *lower = estimate_data.getBox().lower();
  const int *upper = estimate_data.getBox().upper();
  MDA_AccessConst<double,NDIM,MDA_OrderColMajor<NDIM> > co =
    arrayData2ArrayAccess( soln_cell_data.getArrayData() );
  MDA_Access<double,NDIM,MDA_OrderColMajor<NDIM> > es =
    arrayData2ArrayAccess( estimate_data.getArrayData() );
#if NDIM == 2
  int i, j;
  double estimate, est0, est1, est2, est3, est4, est5;
  for ( j=lower[1]; j<=upper[1]; ++j ) {
    for ( i=lower[0]; i<=upper[0]; ++i ) {
      est0 = 
         tbox::MathUtilities<double>::Abs( co(i+1,j) + co(i-1,j) - 2*co(i,j) );
      est1 = 
         tbox::MathUtilities<double>::Abs( co(i,j+1) + co(i,j-1) - 2*co(i,j) );
      est2 = 0.5 *
         tbox::MathUtilities<double>::Abs( co(i+1,j+1) + co(i-1,j-1) - 2*co(i,j) );
      est3 = 0.5 * 
         tbox::MathUtilities<double>::Abs( co(i+1,j-1) + co(i-1,j+1) - 2*co(i,j) );
      est4 = tbox::MathUtilities<double>::Max(est0, est1);
      est5 = tbox::MathUtilities<double>::Max(est2, est3);
      estimate = tbox::MathUtilities<double>::Max(est4, est5);
      es(i,j) = estimate;
    }
  }
#endif
#if NDIM == 3
  // math::PatchCellDataOpsReal<NDIM,double> cops;
  // cops.printData( soln_cell_data_, soln_cell_data_->getGhostBox(), tbox::plog );
  int i, j, k;
  double estimate, est0, est1, est2, est3, est4, est5, est6, est7, est8,
    esta, estb, estc, estd, este, estf, estg;
  for ( k=lower[2]; k<=upper[2]; ++k ) {
    for ( j=lower[1]; j<=upper[1]; ++j ) {
      for ( i=lower[0]; i<=upper[0]; ++i ) {
	est0 = tbox::MathUtilities<double>::Abs( co(i+1,j,k) + co(i-1,j,k) - 2*co(i,j,k) );
	est1 = tbox::MathUtilities<double>::Abs( co(i,j+1,k) + co(i,j-1,k) - 2*co(i,j,k) );
	est2 = tbox::MathUtilities<double>::Abs( co(i,j,k+1) + co(i,j,k-1) - 2*co(i,j,k) );
	est3 = 0.5*tbox::MathUtilities<double>::Abs( co(i,j+1,k+1) + co(i,j-1,k-1) - 2*co(i,j,k) );
	est4 = 0.5*tbox::MathUtilities<double>::Abs( co(i,j+1,k-1) + co(i,j-1,k+1) - 2*co(i,j,k) );
	est5 = 0.5*tbox::MathUtilities<double>::Abs( co(i+1,j,k+1) + co(i-1,j,k-1) - 2*co(i,j,k) );
	est6 = 0.5*tbox::MathUtilities<double>::Abs( co(i+1,j,k-1) + co(i-1,j,k+1) - 2*co(i,j,k) );
	est7 = 0.5*tbox::MathUtilities<double>::Abs( co(i+1,j+1,k) + co(i-1,j-1,k) - 2*co(i,j,k) );
	est8 = 0.5*tbox::MathUtilities<double>::Abs( co(i+1,j-1,k) + co(i-1,j+1,k) - 2*co(i,j,k) );
	esta = tbox::MathUtilities<double>::Max(est0, est1);
	estb = tbox::MathUtilities<double>::Max(est2, est3);
	estc = tbox::MathUtilities<double>::Max(est4, est5);
	estd = tbox::MathUtilities<double>::Max(est6, est7);
	este = tbox::MathUtilities<double>::Max(esta, estb);
	estf = tbox::MathUtilities<double>::Max(estc, estd);
	estg = tbox::MathUtilities<double>::Max(este, estf);
	estimate = tbox::MathUtilities<double>::Max(estg, est8);
	es(i,j,k) = estimate;
      }
    }
  }
#endif
  return;
}

int AdaptivePoisson::computeError(
  const hier::PatchHierarchy<NDIM> &hierarchy ,
  double *l2norm ,
  double *linorm ,
  tbox::Array<double> &l2norms ,
  tbox::Array<double> &linorms ) const
{

  int ln;

  /*
    Compute error on all levels, all patches.
  */
  double diff=0;
  double l2n=0, wts=0, lin=0;
  /*
    We give wtsum twice the space required so we can combine
    the l2norms during the sumReduction, saving a little
    parallel overhead.
  */
  const int nlevels = hierarchy.getNumberLevels();
  tbox::Array<double> wtsums(2*nlevels);
  for ( ln = nlevels-1; ln >= 0; --ln ) {
    tbox::Pointer< hier::PatchLevel<NDIM> > level = hierarchy.getPatchLevel(ln);

    double &levelwts(wtsums[ln]);
    double &levell2n(wtsums[ln+nlevels]); // l2n and wts combined in 1 array.
    double &levellin(linorms[ln]);

    levell2n = levellin = levelwts = 0.0;

    for ( hier::PatchLevel<NDIM>::Iterator p(level); p; p++ ) {
      tbox::Pointer< hier::Patch<NDIM> > patch = level->getPatch(p());

      tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > patch_geom =
	patch->getPatchGeometry();
      const double *dx( patch_geom->getDx() );
      double vol_weight = dx[0];
#if NDIM >= 2
      vol_weight *= dx[1];
#endif
#if NDIM >= 3
      vol_weight *= dx[2];
#endif

      /*
	Get the patch data.
      */
      tbox::Pointer<pdat::CellData<NDIM,double> > current_solution
	= patch->getPatchData(d_scalar_persistent);
      tbox::Pointer<pdat::CellData<NDIM,double> > exact_solution
	= patch->getPatchData(d_exact_persistent);
      tbox::Pointer<pdat::CellData<NDIM,double> > weight
	= patch->getPatchData(d_weight_persistent);

      {
	MDA_AccessConst<double,NDIM,MDA_OrderColMajor<NDIM> > ex =
          arrayData2ArrayAccess( exact_solution->getArrayData() );
	MDA_AccessConst<double,NDIM,MDA_OrderColMajor<NDIM> > co =
          arrayData2ArrayAccess( current_solution->getArrayData() );
	MDA_AccessConst<double,NDIM,MDA_OrderColMajor<NDIM> > wt =
          arrayData2ArrayAccess( weight->getArrayData() );
        const int *lower = current_solution->getBox().lower();
        const int *upper = current_solution->getBox().upper();
#if NDIM == 2
	for ( int j=lower[1]; j<=upper[1]; ++j ) {
	  for ( int i=lower[0]; i<=upper[0]; ++i ) {
	    /*
	      Disregard zero weights in error computations
	      because they are on coarse grids covered by finer grids.
	    */
	    if ( wt(i,j) > 0 ) {
	      diff = tbox::MathUtilities<double>::Abs( co(i,j) - ex(i,j) );
	      if ( levellin < diff ) levellin = diff;
	      levell2n += wt(i,j)*diff*diff;
	      levelwts += wt(i,j);
	    }
	  }
	}
#endif
#if NDIM == 3
	for ( int k=lower[2]; k<=upper[2]; ++k ) {
	  for ( int j=lower[1]; j<=upper[1]; ++j ) {
	    for ( int i=lower[0]; i<=upper[0]; ++i ) {
	      /*
		Disregard zero weights in error computations
		because they are on coarse grids covered by finer grids.
	      */
	      if ( wt(i,j,k) > 0 ) {
		diff = tbox::MathUtilities<double>::Abs( co(i,j,k) - ex(i,j,k) );
		if ( levellin < diff ) levellin = diff;
		levell2n += wt(i,j,k)*diff*diff;
		levelwts += wt(i,j,k);
	      }
	    }
	  }
	}
#endif
      }
    } // end patch loop

  } // end level loop

  if ( tbox::SAMRAI_MPI::getNodes() > 1 ) {
    /*
      Communicate global data if in parallel.
      We temporarily combine l2norms and wtsum so we can sumReduction
      in one shot, saving some parallel overhead.
    */
    tbox::SAMRAI_MPI::sumReduction( wtsums.getPointer(), 2*nlevels );
    tbox::SAMRAI_MPI::maxReduction( linorms.getPointer(), nlevels );
  }

  for ( ln=0; ln<nlevels; ++ln ) {
    /* Copy l2norm accumulated temporarily in wtsums */
    l2norms[ln] = wtsums[ln+nlevels];
    /* Data for whole hierarchy. */
    l2n += l2norms[ln];
    wts += wtsums[ln];
    lin = linorms[ln] > lin ? linorms[ln] : lin;
    /*
      Data for level ln.
      If a level is completely covered by a finer level,
      wtsums[ln] will legitimately be zero, so that protect
      it from a zero divide.
     */
    l2norms[ln] = wtsums[ln] == 0 ? 0 : sqrt(l2norms[ln]/wtsums[ln]);

  }

  if (wts != 0 ) *l2norm = sqrt( l2n/wts ); else *l2norm = 0.0;
  *linorm = lin;

  return 0;
}

int AdaptivePoisson::solvePoisson(
  tbox::Pointer< hier::PatchHierarchy<NDIM> > hierarchy ,
  int max_cycles ,
  double residual_tol ,
  int pre_sweeps ,
  int post_sweeps ,
  string initial_u )
{

  int ln;
  const int finest_ln = hierarchy->getFinestLevelNumber();
  const int coarsest_ln = 0;

  /*
    Allocate scratch data for use in the solve.
  */
  for ( ln=coarsest_ln; ln<=finest_ln; ++ln ) {
    tbox::Pointer< hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
    level->allocatePatchData(d_rhs_scratch);
  }

  /*
    Create vectors x and b for solving Ax=b.
  */
  solv::SAMRAIVectorReal<NDIM,double>
    x("solution",hierarchy,coarsest_ln,finest_ln) ,
    b("rhs",hierarchy,coarsest_ln,finest_ln) ;
  x.addComponent( tbox::Pointer< hier::Variable<NDIM> >(&d_scalar,false) ,
		  d_scalar_persistent, d_weight_persistent );
  b.addComponent( tbox::Pointer< hier::Variable<NDIM> >(&d_rhs,false) ,
		  d_rhs_scratch, d_weight_persistent );

  /*
    Fill the rhs vector (by filling the d_rhs_scratch data
    that forms its component).
    Fill the boundary condition coefficient data.
  */
  for ( int ln=coarsest_ln; ln <=finest_ln; ++ln ) {
    tbox::Pointer< hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);

    hier::PatchLevel<NDIM>::Iterator pi;
    for ( pi.initialize(*level); pi; pi++ ) {
      int pn = *pi;
      tbox::Pointer< hier::Patch<NDIM> > patch = level->getPatch(pn);

      tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > patch_geom
	= patch->getPatchGeometry();
      const hier::Box<NDIM> &box=patch->getBox();

      tbox::Pointer<pdat::CellData<NDIM,double> >
	source_data = patch->getPatchData( d_constant_source_persistent );
      tbox::Pointer<pdat::CellData<NDIM,double> >
	rhs_data = patch->getPatchData( d_rhs_scratch );
      tbox::Pointer<pdat::CellData<NDIM,double> >
	scalar_data = patch->getPatchData( d_scalar_persistent );
      math::PatchCellDataOpsReal<NDIM,double> cell_ops;
      cell_ops.scale( rhs_data, 1.0, source_data, box );

    }
  }

  // debug:
  if (0) {
    x.setToScalar( 1.0 );
    double weightsums =
      d_fac_ops.computeResidualNorm( x,finest_ln,coarsest_ln );
    tbox::pout << "weightsums: " << weightsums << endl;
  }

  /*
    Fill the vector x with the initial guess, if one is given.
    If not given, we assume the initial guess is in place.
  */
  if ( ! initial_u.empty() ) {
    if ( initial_u == "random" ) {
      x.setRandomValues( 1.0, 0.0 );
    }
    else {
      x.setToScalar( atof(initial_u.c_str()) );
    }
  }

  /*
    Create the viz data writer for use in debugging.
  */
#if 1
  d_visit_writer = new appu::VisItDataWriter<NDIM>("Internal VisIt Writer",
                                                   "ap-debug.visit");
  registerVariablesWithPlotter(*d_visit_writer);
#endif

  /*
    Set up FAC preconditioner object.
  */
  d_fac_preconditioner.setMaxCycles(max_cycles);
  d_fac_preconditioner.setResidualTolerance(residual_tol);
  d_fac_preconditioner.setPresmoothingSweeps(pre_sweeps);
  d_fac_preconditioner.setPostsmoothingSweeps(post_sweeps);
  if ( d_lstream ) d_fac_preconditioner.printClassData(*d_lstream);
  d_fac_preconditioner.initializeSolverState(x,b);

  /*
    Solve the system.
  */
  d_fac_preconditioner.solveSystem ( x , b );
  if ( d_lstream ) *d_lstream
    << "FAC solve completed with\n"
    << setw(30) << "number of iterations: " << d_fac_preconditioner.getNumberIterations() << "\n"
    << setw(30) << "residual norm: " << d_fac_preconditioner.getResidualNorm() << "\n"
    ;
  d_fac_preconditioner.deallocateSolverState();

  /*
    Get data on the solve.
  */
  double avg_convergence_factor, final_convergence_factor;
  d_fac_preconditioner.getConvergenceFactors( avg_convergence_factor ,
					      final_convergence_factor );
  if ( d_lstream ) *d_lstream
    << "Final result: \n"
    << setw(30) << "average convergence factor: "
    << avg_convergence_factor << "\n"
    << setw(30) << "final convergence factor: "
    << final_convergence_factor << "\n"
    ;

  /*
    Fill in boundary ghosts here to get the correct ghost cells
    values used to compute the gradient estimator when plotting.
    We are not sure what state ghost cell values are in after
    the solver finishes.
  */
  {
    xfer::RefineAlgorithm<NDIM> refiner;
    tbox::Pointer<geom::CartesianGridGeometry<NDIM> >
      grid_geometry_ = hierarchy->getGridGeometry();
    geom::CartesianGridGeometry<NDIM> &grid_geometry = *grid_geometry_;
    tbox::Pointer< xfer::RefineOperator<NDIM> > accurate_refine_op
      = grid_geometry.
      lookupRefineOperator( tbox::Pointer< hier::Variable<NDIM> >(&d_scalar,false),
			    "LINEAR_REFINE" );
    TBOX_ASSERT( accurate_refine_op );
    refiner.registerRefine( d_scalar_persistent ,
			    d_scalar_persistent ,
			    d_scalar_persistent ,
			    accurate_refine_op );
    d_robin_refine_patch.setTargetDataId(d_scalar_persistent);
    d_robin_refine_patch.setHomogeneousBc(false);
    d_robin_refine_patch.setCoefImplementation(d_physical_bc_coef);
    tbox::Pointer< xfer::RefineSchedule<NDIM> > refine_schedule;
    for ( ln=coarsest_ln; ln<=finest_ln; ++ln ) {
      tbox::Pointer< hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
      if ( ln > 0 ) {
        /* Include coarser levels in setting data */
	refine_schedule =
	  refiner.createSchedule( level ,
				  ln-1 ,
				  hierarchy ,
				  &d_robin_refine_patch );
      }
      else {
        /* Exclude coarser levels in setting data */
	refine_schedule =
	  refiner.createSchedule( level ,
				  &d_robin_refine_patch );
      }
      refine_schedule->fillData(0.0);
    }
  }
  // x.print(plog,false);

  /*
    Deallocate scratch data.
  */
  for ( int ln=coarsest_ln; ln<=finest_ln; ++ln ) {
    tbox::Pointer< hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
    level->deallocatePatchData(d_rhs_scratch);
  }

#if 1
  /*
    Destroy the viz data writer used for debugging.
  */
  d_visit_writer.setNull();
#endif

  return 0;
}

