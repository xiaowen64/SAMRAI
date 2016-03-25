/*
  File:		$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-1/source/test/FAC/main.C $
  Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
  Revision:	$LastChangedRevision: 1878 $
  Modified:	$LastChangedDate: 2008-01-22 08:25:20 -0800 (Tue, 22 Jan 2008) $
  Description:	Program for poisson solver on adaptive grid using FAC
*/

#include "SAMRAI_config.h"


#include IOMANIP_HEADER_FILE
#include <fstream>
using namespace std;

#include "printObject.h"
#include "AdaptivePoisson.h"
#include "get-input-filename.h"

/*
  Headers for basic SAMRAI objects used in this code.
*/
#include "tbox/IOStream.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/InputManager.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/SAMRAIManager.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"

/*
  Headers for major algorithm/data structure objects from SAMRAI
*/
#include "CartesianVizamraiDataWriter.h"
#include "VisItDataWriter.h"
#include "CartesianGridGeometry.h"
#include "GridGeometry.h"
#include "PatchHierarchy.h"
#include "VariableDatabase.h"
#include "BergerRigoutsos.h"
#include "GriddingAlgorithm.h"
#include "StandardTagAndInitialize.h"
#include "LoadBalancer.h"
#include "FACPreconditioner.h"

#ifndef LACKS_NAMESPACE
using namespace SAMRAI;
#endif

int main( int argc, char *argv[] )
{

   string input_filename;

   /*
     Initialize MPI, process argv, and initialize SAMRAI
   */
   tbox::SAMRAI_MPI::init(&argc, &argv);
   if ( get_input_filename(&argc, argv, input_filename) == 1 ) {
      tbox::pout << "Usage: " << argv[0] << " <input file>." << endl;
      tbox::SAMRAI_MPI::finalize();
      return 0;
   }
   tbox::SAMRAIManager::startup();

   bool error_ok=false;

   /*
    * Create block to force pointer deallocation.  If this is not done
    * then there will be memory leaks reported.
    */
   {


      /*
	These tests are set up to use hypre.
	Do not run them without hypre.
      */
#ifdef HAVE_HYPRE
      tbox::pout << "Input file is " << input_filename << endl;

      /*
	Create input database and parse all data in input file into it.
      */

      tbox::Pointer<tbox::Database> input_db = new tbox::InputDatabase("input_db");
      tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

      if ( input_db->isDatabase("TimerManager") ) {
	 tbox::TimerManager::createManager(input_db->getDatabase("TimerManager"));
      }

      /*
	Get the Main database part of the input database.
	This database contains information relevant to main.
      */

      tbox::Pointer<tbox::Database> main_db = input_db->getDatabase("Main");
      tbox::plog << "Main database:" << endl; main_db->printClassData(tbox::plog);

      /*
	Base filename info.
      */

      string base_name = main_db->getStringWithDefault("base_name", "noname");

      /*
	Log file info.
      */
      {
	 string log_filename
	    = main_db->getStringWithDefault("log_filename", base_name+".log");
	 bool log_all
	    = main_db->getBoolWithDefault("log_all", false);
	 if ( log_all )
	    tbox::PIO::logAllNodes(log_filename);
	 else 
	    tbox::PIO::logOnlyNodeZero(log_filename);
      }

      /*
	Create a patch hierarchy for use later.
	This object is a required input for these objects: adaptive_poisson.
      */
      tbox::Pointer< hier::PatchHierarchy<NDIM> > patch_hierarchy;
      {
	 /*
	   Create a grid geometry required for the patchHierarchy object.
	 */
	 tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geometry =
	    new geom::CartesianGridGeometry<NDIM>("CartesianGridGeometry",
						  input_db->getDatabase("CartesianGridGeometry"));
	 tbox::plog << "Grid Geometry:" << endl;
	 grid_geometry->printClassData(tbox::plog);
	 patch_hierarchy =
	    new hier::PatchHierarchy<NDIM>("Patch Hierarchy", grid_geometry);
      }

      /*
	Create the problem-specific object implementing the required
	SAMRAI virtual functions.
      */
      AdaptivePoisson
	 adaptive_poisson( "AdaptivePoisson" ,
			   *(input_db->getDatabase("AdaptivePoisson")) ,
			   &tbox::pout ,
			   &tbox::plog);

      tbox::Pointer< mesh::GriddingAlgorithm<NDIM> > gridding_algorithm;
      {
	 /*
	   Create the tag-and-initializer, box-generator and load-balancer
	   object references required by the gridding_algorithm object.
	 */
	 tbox::Pointer< mesh::StandardTagAndInitialize<NDIM> > tag_and_initializer =
	    new mesh::StandardTagAndInitialize<NDIM>(
	       "CellTaggingMethod",
	       tbox::Pointer< mesh::StandardTagAndInitStrategy<NDIM> >(&adaptive_poisson,false) ,
	       input_db->getDatabase("StandardTagAndInitialize")
	       );
	 tbox::Pointer< mesh::BergerRigoutsos<NDIM> > box_generator =
	    new mesh::BergerRigoutsos<NDIM>();
	 tbox::Pointer< mesh::LoadBalancer<NDIM> > load_balancer =
	    new mesh::LoadBalancer<NDIM>(input_db->getDatabase("LoadBalancer"));

	 /*
	   Create the gridding algorithm used to generate the SAMR grid
	   and create the grid.
	 */
	 gridding_algorithm =
	    new mesh::GriddingAlgorithm<NDIM>("Gridding Algorithm",
					      input_db->getDatabase("GriddingAlgorithm"),
					      tag_and_initializer,
					      box_generator,
					      load_balancer);
	 tbox::plog << "Gridding algorithm:" << endl;
	 gridding_algorithm->printClassData(tbox::plog);
	 /*
	   Make the coarse patch level.
	 */
	 gridding_algorithm->makeCoarsestLevel(patch_hierarchy,0.0);
      }

      int ln;

      /* Whether to plot */
      string vis_filename
	 = main_db->getStringWithDefault("vis_filename", base_name);
      bool do_plot
	 = main_db->getBoolWithDefault("do_plot", false );

      /*
	After creating all objects and initializing their state,
	we print the input database and variable database contents
	to the log file.
      */
      tbox::plog << "\nCheck input data and variables before simulation:" << endl;
      tbox::plog << "Input database..." << endl;
      input_db->printClassData(tbox::plog);
      tbox::plog << "\nVariable database..." << endl;
      hier::VariableDatabase<NDIM>::getDatabase()->printClassData(tbox::plog);

      tbox::plog << "\n\nFinal Hierarchy:\n";
      patch_hierarchy->recursivePrint( tbox::plog, "\t", 2 );

      double target_l2norm = 1e-6;
      target_l2norm = main_db->getDoubleWithDefault( "target_l2norm" ,
						     target_l2norm );
      double l2norm, linorm;
      int max_adaptions = 1;
      max_adaptions = main_db->getIntegerWithDefault( "max_adaptions" ,
						      max_adaptions );
      int adaption_number = 0;
      bool done=false;
      do {
	 /*
	   Solve.
	 */
	 string max_cycles_str = "max_cycles";
	 int max_cycles = main_db->getIntegerWithDefault(max_cycles_str, 10);
	 double residual_tol = main_db->getDoubleWithDefault("residual_tol", 1e-6);
	 tbox::pout.setf(ios::scientific);
	 int pre_sweeps = main_db->getIntegerWithDefault("pre_sweeps", 5);
	 int post_sweeps = main_db->getIntegerWithDefault("post_sweeps", 5);
	 string initial_u = main_db->getStringWithDefault("initial_u", "0.0");
	 adaptive_poisson.solvePoisson( patch_hierarchy ,
					max_cycles ,
					residual_tol ,
					pre_sweeps ,
					post_sweeps ,
					adaption_number ? string() : initial_u );
	 tbox::Array<double> l2norms(patch_hierarchy->getNumberLevels())
	    , linorms(patch_hierarchy->getNumberLevels());
	 adaptive_poisson.computeError( *patch_hierarchy, &l2norm, &linorm
					, l2norms, linorms );
	 error_ok = l2norm <= target_l2norm;
	 tbox::plog << "Err " << (error_ok ? "" : "NOT ") << "ok, err norm/target: "
		    << l2norm << '/' << target_l2norm << endl;
	 tbox::pout << "Err result after " << adaption_number << " adaptions: \n"
		    << setw(15) << "l2: " << setw(10) << l2norm
		    << setw(15) << "li: " << setw(10) << linorm
		    << "\n";
	 for ( ln=0; ln<patch_hierarchy->getNumberLevels(); ++ln ) {
	    tbox::pout << setw(10) << "l2[" << setw(2) << ln << "]: "
		       << setw(10) << l2norms[ln]
		       << setw(10) << "li[" << setw(2) << ln << "]: "
		       << setw(10) << linorms[ln]
		       << "\n";
	 }

	 /* Write the plot file. */
	 if ( do_plot ) {
	    tbox::Pointer<appu::VisItDataWriter<NDIM> > visit_writer =
	       new appu::VisItDataWriter<NDIM>("VisIt Writer", vis_filename+".visit");
	    adaptive_poisson.registerVariablesWithPlotter(*visit_writer);
	    visit_writer->writePlotData( patch_hierarchy ,
					 adaption_number );
	    tbox::pout << "Wrote viz file " << vis_filename << " for grid number "
		       << adaption_number << '\n';
	 }

	 /*
	   Done when max adaptions or convergence reached.
	 */
	 done = error_ok || ( adaption_number >= max_adaptions );

	 if ( !done ) {
	    /*
	      Adapt grid.
	    */
	    ++adaption_number;
	    tbox::pout << "Adaption number " << adaption_number << "\n";
	    tbox::Array<int> tag_buffer( gridding_algorithm->getMaxLevels() );
	    for ( ln=0; ln<tag_buffer.getSize(); ++ln ) {
	       tag_buffer[ln] = 1;
	    }
	    gridding_algorithm->regridAllFinerLevels(patch_hierarchy,
						     0,
						     0.0,
						     tag_buffer );
	    tbox::plog << "Newly adapted hierarchy\n";
	    patch_hierarchy->recursivePrint( tbox::plog, "    ", 1 );
	    if ( 0 ) {
	       /* Write post-adapt viz file for debugging */
	       tbox::Pointer<appu::VisItDataWriter<NDIM> > visit_writer =
		  new appu::VisItDataWriter<NDIM>("VisIt Writer", "postadapt.visit");
	       adaptive_poisson.registerVariablesWithPlotter(*visit_writer);
	       visit_writer->writePlotData( patch_hierarchy ,
					    adaption_number-1 );
	       tbox::plog << "Wrote viz file " << "postadapt.visit" << '\n';
	    }
	 }
      } while ( !done );

      tbox::pout << "After " << adaption_number << "/" << max_adaptions
		 << " adaptions, err is " << l2norm << "/" << target_l2norm
		 << endl;

      tbox::TimerManager::getManager()->print(tbox::plog);
#else
   error_ok=true;
#endif

      if ( error_ok ) {
	 tbox::pout << "\nPASSED:  FAC" << endl;
      }
      else {
	 TBOX_WARNING("Failed to meet accuracy specifications.");
      }
   }

   /*
     Exit properly by shutting down services in correct order.
   */
   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAI_MPI::finalize();

   return(0);
}
