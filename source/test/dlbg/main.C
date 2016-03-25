/*
  File:		$RCSfile$
  Copyright:	(c) 1997-2003 The Regents of the University of California
  Revision:	$Revision: 586 $
  Modified:	$Date: 2005-08-23 10:49:46 -0700 (Tue, 23 Aug 2005) $
  Description:	Test program for asynchronous BR implementation
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include "SAMRAI_config.h"
#include "tbox/SAMRAIManager.h"
#include "tbox/MPI.h"

/*
  Headers for basic SAMRAI objects used in this code.
*/
#include "tbox/SAMRAIManager.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/InputManager.h"
#include "tbox/MPI.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "tbox/Statistician.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"

/*
  Headers for major algorithm/data structure objects from SAMRAI
*/
#include "CartesianVizamraiDataWriter.h"
#include "VisItDataWriter.h"
#include "CartesianGridGeometry.h"
#include "BoxList.h"
#include "GridGeometry.h"
#include "LayerEdgeSet.h"
#include "PatchHierarchy.h"
#include "VariableDatabase.h"
#include "AsyncBergerRigoutsos.h"
#include "BergerRigoutsos.h"
#include "GriddingAlgorithm.h"
#include "StandardTagAndInitialize.h"
#include "LoadBalancer.h"
#include "FACPreconditioner.h"
#include "DLBGTest.h"

#include <assert.h>

#include "get-input-filename.h"
#include "intToString.h"

#ifndef LACKS_NAMESPACE
using namespace SAMRAI;
using namespace tbox;
using namespace hier;
using namespace mesh;
#endif


int createAndTestDLBG( tbox::Database &main_db,
                       PatchHierarchy<NDIM> &patch_hierarchy );


int main( int argc, char **argv )
{
  string input_filename;

  /*
    Initialize MPI, process argv, and initialize SAMRAI
  */
  tbox::MPI::init(&argc, &argv);
  if ( get_input_filename(&argc, argv, input_filename) == 1 ) {
    cout << "Usage: " << argv[0]
	 << " <input file>."
	 << endl;
    tbox::MPI::finalize();
    return 0;
  }
  tbox::SAMRAIManager::startup();


  tbox::pout << "Input file is " << input_filename << endl;


  string case_name;
  if ( argc >= 2 ) {
    case_name = argv[1];
  }


  /*
    Create input database and parse all data in input file into it.
  */

  tbox::Pointer<tbox::Database> input_db = new tbox::InputDatabase("input_db");
  tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

   /*
    * Retrieve "GlobalInputs" section of the input database and set
    * values accordingly.
    */

   if (input_db->keyExists("GlobalInputs")) {
      tbox::Pointer<tbox::Database> global_db =
         input_db->getDatabase("GlobalInputs");
      if (global_db->keyExists("tag_clustering_method")) {
         string tag_clustering_method =
            global_db->getString("tag_clustering_method");
         mesh::BergerRigoutsos<NDIM>::setClusteringOption(tag_clustering_method);
      }
      if (global_db->keyExists("call_abort_in_serial_instead_of_exit")) {
         bool flag = global_db->
            getBool("call_abort_in_serial_instead_of_exit");
         tbox::MPI::setCallAbortInSerialInsteadOfExit(flag);
      }
   }


  /*
    Get the Main database part of the input database.
    This database contains information relevant to main.
  */

  tbox::Pointer<tbox::Database> main_db = input_db->getDatabase("Main");
  tbox::plog << "Main database:" << endl; main_db->printClassData(tbox::plog);


  if ( input_db->isDatabase("TimerManager") ) {
    tbox::TimerManager::createManager(input_db->getDatabase("TimerManager"));
  }


  /*
    Base filename info.
  */

  string base_name = main_db->getStringWithDefault("base_name", "fp");


  /*
    Modify basename for this particular run.
    Add the number of processes and the case name.
  */
  if ( ! case_name.empty() ) {
    base_name = base_name + '-' + case_name;
  }
  if ( tbox::MPI::getNodes() > 1 ) {
    base_name = base_name + '-' + intToString( tbox::MPI::getNodes(), 5 );
  }
  tbox::pout << "Added case name (" << case_name << ") and nprocs ("
       << tbox::MPI::getNodes() << ") to base name -> '"
       << base_name << "'\n";

  /*
    Set the vis filename, defaults to base_name.
  */
  string vis_filename
     = main_db->getStringWithDefault("vis_filename", base_name);


  /*
    Log file info.
  */

  string log_filename
    = main_db->getStringWithDefault("log_filename", base_name+".log");
  bool log_all = false;
  log_all = main_db->getBoolWithDefault("log_all", log_all);
  if ( log_all && tbox::MPI::getNodes() > 1 ) {
     tbox::PIO::logAllNodes(log_filename);
  }
  else {
     tbox::PIO::logOnlyNodeZero(log_filename);
  }


  if ( ! case_name.empty() ) {
    tbox::plog << "Added case name (" << case_name << ") and nprocs ("
	 << tbox::MPI::getNodes() << ") to base name -> '"
	 << base_name << "'\n";
  }
  tbox::pout << "Running on " << tbox::MPI::getNodes() << " processes.\n";


  /*
    Choose which BR implementation to use.
  */
  char which_br = 'o';
  which_br = main_db->getCharWithDefault("which_br", which_br);
  tbox::pout << "which_br is " << which_br << endl;


  int ln;


  int plot_step = main_db->getIntegerWithDefault("plot_step", 0 );


  /*
    Create a patch hierarchy for use later.
    This object is a required input for these objects: abrtest.
  */
  tbox::Pointer<hier::PatchHierarchy<NDIM> > patch_hierarchy;
  {
     /*
       Create a grid geometry required for the
       hier::PatchHierarchy<NDIM> object.
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
  tbox::plog << "Creating abrtest.\n";
  DLBGTest<NDIM> abrtest ( "DLBGTest" ,
                          input_db->getDatabase("DLBGTest") );



  tbox::plog << "Creating box generator.\n";
  tbox::Pointer<mesh::BergerRigoutsos<NDIM> > old_br =
    new mesh::BergerRigoutsos<NDIM>();
  tbox::Pointer<mesh::AsyncBergerRigoutsos<NDIM> > new_br =
    new mesh::AsyncBergerRigoutsos<NDIM> ( input_db->isDatabase("AsyncBergerRigoutsos") ?
				    input_db->getDatabase("AsyncBergerRigoutsos") :
				    tbox::Pointer<tbox::Database>(NULL) );

  tbox::Pointer<mesh::BoxGeneratorStrategy<NDIM> > box_generator = which_br == 'o' ?
    tbox::Pointer<mesh::BoxGeneratorStrategy<NDIM> >(old_br) :
    tbox::Pointer<mesh::BoxGeneratorStrategy<NDIM> >(new_br) ;
  TBOX_ASSERT( !box_generator.isNull() );

  tbox::plog << "Creating grid algorithm.\n";
  tbox::Pointer<mesh::GriddingAlgorithm<NDIM> > gridding_algorithm;
  {
    /*
      Create the tag-and-initializer, box-generator and load-balancer
      object references required by the gridding_algorithm object.
    */
    tbox::Pointer<mesh::StandardTagAndInitialize<NDIM> > tag_and_initializer =
      new mesh::StandardTagAndInitialize<NDIM>(
          "CellTaggingMethod" ,
	  tbox::Pointer<mesh::StandardTagAndInitStrategy<NDIM> >(abrtest.getStandardTagAndInitObject(),false) ,
	  input_db->getDatabase("StandardTagAndInitialize")
        );
    tbox::Pointer<mesh::LoadBalancer<NDIM> > load_balancer =
      new mesh::LoadBalancer<NDIM>("load balancer",
			     input_db->getDatabase("LoadBalancer"));

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

  }

  bool log_hierarchy = false;
  log_hierarchy = main_db->getBoolWithDefault("log_hierarchy",
                                              log_hierarchy);
  int num_steps = num_steps = main_db->getIntegerWithDefault("num_steps",0);

  /*
    After setting up the problem and initializing the object states,
    we print the input database and variable database contents
    to the log file.
  */
  tbox::plog << "\nCheck input data:" << endl;
  tbox::plog << "Input database..." << endl;
  input_db->printClassData(tbox::plog);
  tbox::plog << "\nVariable database..." << endl;
  hier::VariableDatabase<NDIM>::getDatabase()->printClassData(tbox::plog);

  /*
    Make the patch levels.
  */

  static tbox::Pointer<tbox::Timer> t_generate_mesh =
     tbox::TimerManager::getManager()->
     getTimer("apps::main::generate_mesh");
  t_generate_mesh->start();
  gridding_algorithm->makeCoarsestLevel(patch_hierarchy,0.0);
  bool done=false;
  for (ln = 0; gridding_algorithm->levelCanBeRefined(ln) && !done; ln++) {
     tbox::pout << "Adding finer levels with ln = " << ln << endl;
     tbox::Pointer<hier::PatchLevel<NDIM> > level_ =
        patch_hierarchy->getPatchLevel(ln);
     gridding_algorithm->makeFinerLevel( patch_hierarchy ,
                                         /* simulation time */ 0.0 ,
                                         /* whether initial time */ true ,
                                         /* tag buffer size */ 0 );
     tbox::pout << "Just added finer level " << ln << " -> " << ln+1;
     tbox::Pointer<hier::PatchLevel<NDIM> > new_level =
        patch_hierarchy->getPatchLevel(ln+1);
     if ( patch_hierarchy->getNumberLevels() < ln+2 ) {
	tbox::pout << " (no new level!)" << endl;
     }
     else {
        tbox::Pointer<hier::PatchLevel<NDIM> > finer_level_ =
           patch_hierarchy->getPatchLevel(ln+1);
	tbox::pout
	   << " (" << level_->getNumberOfPatches()
	   << " -> " << finer_level_->getNumberOfPatches()
	   << " patches)"
	   << endl;
     }
     done = !(patch_hierarchy->finerLevelExists(ln));
     
  }
  t_generate_mesh->stop();


  if ( tbox::MPI::getRank() == 0 ) {
     tbox::pout << "Hierarchy generated:" << endl;
     patch_hierarchy->recursivePrint( tbox::pout, string("    "), 1 );
  }
  if ( log_hierarchy ) {
     tbox::plog << "Hierarchy generated:" << endl;
     patch_hierarchy->recursivePrint( tbox::plog, string("H-> "), 3 );
  }




  /*
    Write a plot file.
  */
  /* Get the output filename. */
  if ( plot_step > 0 ) {
    if(0) {
      /* Create the vizamrai data writer. */
      tbox::Pointer<appu::CartesianVizamraiDataWriter<NDIM> > viz_data_writer =
	new appu::CartesianVizamraiDataWriter<NDIM>("Viz Writer");
      /* Register variables with plotter. */
      abrtest.registerVariablesWithPlotter(viz_data_writer);
      /*
	Tell the plotter about the refinement ratios.
	This must be done once (and again each time the data changes).
      */
      for ( int ln=1; ln<patch_hierarchy->getNumberLevels(); ln++ ) {
         tbox::Pointer<hier::PatchLevel<NDIM> >
            level_ =patch_hierarchy->getPatchLevel(ln);
	const hier::IntVector<NDIM> &lratio =
	  level_->getRatioToCoarserLevel();
	viz_data_writer->setRatioToCoarserLevel(ln, lratio);
      }
      /* Write the plot file. */
      viz_data_writer->writePlotData( patch_hierarchy , vis_filename );
    }
    {
      const string visit_filename = vis_filename+".visit";
      /* Create the VisIt data writer. */
      tbox::Pointer<appu::VisItDataWriter<NDIM> > visit_data_writer =
	new appu::VisItDataWriter<NDIM>("VisIt Writer", visit_filename);
      /* Register variables with plotter. */
      abrtest.registerVariablesWithPlotter(visit_data_writer);
      /* Write the plot file. */
      visit_data_writer->writePlotData( patch_hierarchy , 0 );
    }
  }


  createAndTestDLBG( *main_db, *patch_hierarchy );

  /*
    Adapt the grid.
  */
  tbox::Array<int> tag_buffer(10);
  for ( int i=0; i<tag_buffer.size(); ++i ) tag_buffer[i] = 1;



  for ( int istep=0; istep<num_steps; ++istep ) {

     tbox::pout << "Adaption number " << istep << endl;

     // Recompute the front-dependent data at next time step.
     abrtest.computeHierarchyData( *patch_hierarchy,
                                   double(istep+1) );


     tbox::Array<double> regrid_start_time(0);
     for ( int i=0; i<regrid_start_time.size(); ++i )
        regrid_start_time[i] = istep;

     gridding_algorithm->regridAllFinerLevels( patch_hierarchy,
                                               0,
                                               double(istep+1),
                                               tag_buffer,
                                               regrid_start_time );


     if ( tbox::MPI::getRank() == 0 ) {
        patch_hierarchy->recursivePrint( tbox::pout, string("    "), 1 );
     }
     if ( log_hierarchy ) {
        tbox::plog << "Hierarchy adapted:" << endl;
        patch_hierarchy->recursivePrint( tbox::plog, string("H-> "), 3 );
     }


     if ( plot_step > 0 && (istep+1)%plot_step == 0  ) {
        if(0) {
           /* Create the vizamrai data writer. */
           tbox::Pointer<appu::CartesianVizamraiDataWriter<NDIM> > viz_data_writer =
              new appu::CartesianVizamraiDataWriter<NDIM>("Viz Writer");
           /* Register variables with plotter. */
           abrtest.registerVariablesWithPlotter(viz_data_writer);
           /*
             Tell the plotter about the refinement ratios.
             This must be done once (and again each time the data changes).
           */
           for ( int ln=1; ln<patch_hierarchy->getNumberLevels(); ln++ ) {
              tbox::Pointer<hier::PatchLevel<NDIM> >
                 level_ =patch_hierarchy->getPatchLevel(ln);
              const hier::IntVector<NDIM> &lratio =
                 level_->getRatioToCoarserLevel();
              viz_data_writer->setRatioToCoarserLevel(ln, lratio);
           }
           /* Write the plot file. */
           viz_data_writer->writePlotData( patch_hierarchy , vis_filename );
        }
        {
           const string visit_filename = vis_filename+".visit";
           /* Create the VisIt data writer. */
           tbox::Pointer<appu::VisItDataWriter<NDIM> > visit_data_writer =
              new appu::VisItDataWriter<NDIM>("VisIt Writer", visit_filename);
           /* Register variables with plotter. */
           abrtest.registerVariablesWithPlotter(visit_data_writer);
           /* Write the plot file. */
           visit_data_writer->writePlotData( patch_hierarchy , istep+1 );
        }
     }

     createAndTestDLBG( *main_db, *patch_hierarchy );

  }




  tbox::TimerManager::getManager()->print(tbox::plog);


  tbox::pout << "\nPASSED:  DLBG" << endl;

  /*
    Exit properly by shutting down services in correct order.
  */
  tbox::plog << "\nShutting down..." << endl;
  tbox::SAMRAIManager::shutdown();
  tbox::MPI::finalize();


  return(0);
}




int createAndTestDLBG( tbox::Database &main_db,
                       PatchHierarchy<NDIM> &patch_hierarchy )
{
  /*
    Generate the nested-level box-graph.
  */

  const bool build_cross_edge =
    main_db.getBoolWithDefault("build_cross_edge", false );
  const bool build_peer_edge =
    main_db.getBoolWithDefault("build_peer_edge", false );

  const bool globalize_node_sets =
    main_db.getBoolWithDefault("globalize_node_sets", false );

  const int node_log_detail =
    main_db.getIntegerWithDefault("node_log_detail", -1 );
  const int edge_log_detail =
    main_db.getIntegerWithDefault("edge_log_detail", -1 );

  static tbox::Pointer<tbox::Timer> t_generate_nlbg_in_parallel =
    tbox::TimerManager::getManager()->
    getTimer("apps::main::generate_nlbg_in_parallel");

  pout << "Generating boxgraph in parallel." << endl;

  //hier_NestedLevelBoxGraphX boxgraph("test boxgraph");

  LayerNodeSet<NDIM> *node_sets = NULL;
  LayerEdgeSet<NDIM> *crse_edge_sets = NULL;
  LayerEdgeSet<NDIM> *fine_edge_sets = NULL;
  LayerEdgeSet<NDIM> *peer_edge_sets = NULL;

  int ln;

  node_sets = new LayerNodeSet<NDIM>[patch_hierarchy.getNumberLevels()];
  if ( build_peer_edge ) {
     peer_edge_sets = new LayerEdgeSet<NDIM>[patch_hierarchy.getNumberLevels()];
  }
  if ( build_cross_edge ) {
     crse_edge_sets = new LayerEdgeSet<NDIM>[patch_hierarchy.getNumberLevels()];
     fine_edge_sets = new LayerEdgeSet<NDIM>[patch_hierarchy.getNumberLevels()];
  }

  /*
    Set the layer nodes.
  */
  for ( ln=0; ln<patch_hierarchy.getNumberLevels(); ++ln ) {
     Pointer<PatchLevel<NDIM> > level_ptr = patch_hierarchy.getPatchLevel(ln);
     PatchLevel<NDIM> &level = *level_ptr;
     node_sets[ln].setTo(level);
     plog << "****************************************\n";
     plog << "node_sets[" << ln << "]:\n";
     plog << "****************************************\n";
     node_sets[ln].printClassData( plog, node_log_detail );
     if ( globalize_node_sets ) {
        pout << "Globalizing NodeSet " << ln << ".\n";
        node_sets[ln].setParallelState(LayerNodeSet<NDIM>::GLOBALIZED);
     }
     pout << "NodeSet " << ln << " done.\n";
  }

  /*
    Compute the cross edges by searching globalized node layers.
  */
  if ( build_cross_edge ) {
     for ( ln=0; ln<patch_hierarchy.getNumberLevels(); ++ln ) {
        Pointer<PatchLevel<NDIM> > level_ptr = patch_hierarchy.getPatchLevel(ln);
        PatchLevel<NDIM> &level = *level_ptr;
        if ( ln < patch_hierarchy.getNumberLevels()-1 ) {
           fine_edge_sets[ln].initialize(LayerEdgeSet<NDIM>::DISTRIBUTED,
                                         node_sets[ln],
                                         node_sets[ln+1].getRefinementRatio(),
                                         IntVector<NDIM>(1) );
        }
        if ( ln > 0 ) {
           crse_edge_sets[ln].initialize(LayerEdgeSet<NDIM>::DISTRIBUTED,
                                         node_sets[ln],
                                         node_sets[ln-1].getRefinementRatio(),
                                         level.getRatioToCoarserLevel() );
           crse_edge_sets[ln].attachPartner(fine_edge_sets[ln-1]);
           crse_edge_sets[ln].findEdges();
           if ( edge_log_detail >= 0 ) {
              plog << "****************************************\n";
              plog << "fine_edge_sets[" << ln-1 << "]:\n";
              plog << "****************************************\n";
              fine_edge_sets[ln-1].printClassData( plog, edge_log_detail );
              plog << "****************************************\n";
              plog << "crse_edge_sets[" << ln << "]:\n";
              plog << "****************************************\n";
              crse_edge_sets[ln].printClassData( plog, edge_log_detail );
           }
        }
     }
  }


  /*
    Compute the peer edges by bridging centered on the coarse edge
    if available, or by searching globalized node layers.
  */
  if ( build_peer_edge ) {
     for ( ln=0; ln<patch_hierarchy.getNumberLevels(); ++ln ) {
        peer_edge_sets[ln].initialize(LayerEdgeSet<NDIM>::DISTRIBUTED,
                                      node_sets[ln],
                                      node_sets[ln].getRefinementRatio(),
                                      IntVector<NDIM>(1) );
        peer_edge_sets[ln].attachPartner(peer_edge_sets[ln]);
        if ( build_cross_edge && ln > 0 ) {
           peer_edge_sets[ln].bridge( fine_edge_sets[ln-1], fine_edge_sets[ln-1] );
           // peer_edge_sets[ln].findEdges();
        }
        else {
           peer_edge_sets[ln].findEdges();
        }
        if ( edge_log_detail >= 0 ) {
           plog << "****************************************\n";
           plog << "peer_edge_sets[" << ln << "]:\n";
           plog << "****************************************\n";
           peer_edge_sets[ln].printClassData( plog, edge_log_detail );
        }
     }
  }



  /*
    Check for accuracy.
  */

  if ( build_cross_edge ) {
     for ( ln=0; ln<patch_hierarchy.getNumberLevels()-1; ++ln ) {
        fine_edge_sets[ln].checkNodeConsistency();
        plog << "fine_edge_sets[" << ln   << "] passed checkNodeConsistency().\n";
        crse_edge_sets[ln+1].checkNodeConsistency();
        plog << "crse_edge_sets[" << ln+1 << "] passed checkNodeConsistency().\n";
        fine_edge_sets[ln].checkConnectivity(node_sets[ln+1]);
        plog << "fine_edge_sets[" << ln   << "] passed checkConnectivity().\n";
        crse_edge_sets[ln+1].checkConnectivity(node_sets[ln]);
        plog << "crse_edge_sets[" << ln+1 << "] passed checkConnectivity().\n";
     }
  }

  if ( build_peer_edge ) {
     for ( ln=0; ln<patch_hierarchy.getNumberLevels(); ++ln ) {
        peer_edge_sets[ln].checkNodeConsistency();
        plog << "peer_edge_sets[" << ln << "] passed checkNodeConsistency().\n";
        peer_edge_sets[ln].checkConnectivity(node_sets[ln]);
        plog << "peer_edge_sets[" << ln << "] passed checkConnectivity().\n";
     }
  }


   if ( node_sets      ) delete [] node_sets     ;
   if ( crse_edge_sets ) delete [] crse_edge_sets;
   if ( fine_edge_sets ) delete [] fine_edge_sets;
   if ( peer_edge_sets ) delete [] peer_edge_sets;

   return 0;
}
