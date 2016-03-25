//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/examples/Euler/Euler.C $
// Package:     SAMRAI application
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Numerical routines for Euler equations SAMRAI example
//

#include "Euler.h"

#include <iostream>
#include <iomanip>
#include <fstream>

#ifndef LACKS_SSTREAM
#ifndef included_sstream
#define included_sstream
#include <sstream>
#endif
#else
#ifndef included_strstream
#define included_strstream
#include <strstream>
#endif
#endif

using namespace std;

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include "BoxArray.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellIndex.h"
#include "CellIterator.h"
#include "FaceData.h"
#include "FaceIndex.h"
#include "FaceVariable.h"
#include "Index.h"
#include "LoadBalancer.h"
#include "tbox/PIO.h"
#include "tbox/RestartManager.h"
#include "VariableDatabase.h"
#include "tbox/TimerManager.h"
#include "tbox/Timer.h"
#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"

//integer constants for boundary conditions
#define CHECK_BDRY_DATA  (0)
#include "CartesianBoundaryDefines.h"

//integer constant for debugging improperly set boundary dat
#define BOGUS_BDRY_DATA   (-9999)

// routines for managing boundary data
#if (NDIM == 2) 
#include "CartesianBoundaryUtilities2.h"
#endif
#if (NDIM == 3)
#include "CartesianBoundaryUtilities3.h"
#endif

// External definitions for Fortran numerical routines
#include "EulerFort.h"

// Number of entries in state vector (NDIM velocity comps + pressure + density)
#define NEQU         (NDIM + 2)

// Number of ghosts cells used for each variable quantity.
#define CELLG           (4)
#define FACEG           (4)
#define FLUXG           (1)

// defines for initialization
#define PIECEWISE_CONSTANT_X        (10)
#define PIECEWISE_CONSTANT_Y        (11)
#define PIECEWISE_CONSTANT_Z        (12)
#define SPHERE                      (40)
#define STEP                        (80)

// defines for Riemann solver used in Godunov flux calculation
#define APPROX_RIEM_SOLVE   (20) // Colella-Glaz approx Riemann solver 
#define EXACT_RIEM_SOLVE    (21) // Exact Riemann solver
#define HLLC_RIEM_SOLVE     (22) // Harten, Lax, van Leer approx Riemann solver

// defines for cell tagging routines
#define RICHARDSON_NEWLY_TAGGED (-10)
#define RICHARDSON_ALREADY_TAGGED (-11)
#ifndef TRUE
#define TRUE (1)
#endif
#ifndef FALSE
#define FALSE (0)
#endif


// Version of Euler restart file data
#define EULER_VERSION (3)

tbox::Pointer<tbox::Timer> Euler::t_init;
tbox::Pointer<tbox::Timer> Euler::t_compute_dt;
tbox::Pointer<tbox::Timer> Euler::t_compute_fluxes;
tbox::Pointer<tbox::Timer> Euler::t_conservdiff;
tbox::Pointer<tbox::Timer> Euler::t_setphysbcs;
tbox::Pointer<tbox::Timer> Euler::t_taggradient;

/*
*************************************************************************
*                                                                       *
* The constructor for Euler class sets data members to defualt values,  *
* creates variables that define the solution state for the Euler        *
* equations.                                                            * 
*                                                                       *
* After default values are set, this routine calls getFromRestart()     *
* if execution from a restart file is specified.  Finally,              *
* getFromInput() is called to read values from the given input          *
* database (potentially overriding those found in the restart file).    *
*                                                                       *
*************************************************************************
*/

Euler::Euler(const string& object_name,
             tbox::Pointer<tbox::Database> input_db,
             tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geom)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(!input_db.isNull());
   TBOX_ASSERT(!grid_geom.isNull());
#endif

   d_object_name = object_name;
   tbox::RestartManager::getManager()->registerRestartItem(d_object_name, this);

   if ( t_init.isNull() ) {
      t_init = tbox::TimerManager::getManager()->
         getTimer("apps::Euler::initializeDataOnPatch()");
      t_compute_dt = tbox::TimerManager::getManager()->
         getTimer("apps::Euler::computeStableDtOnPatch()");
      t_compute_fluxes = tbox::TimerManager::getManager()->
         getTimer("apps::Euler::computeFluxesOnPatch()");
      t_conservdiff = tbox::TimerManager::getManager()->
         getTimer("apps::Euler::conservativeDifferenceOnPatch()");
      t_setphysbcs = tbox::TimerManager::getManager()->
         getTimer("apps::Euler::setPhysicalBoundaryConditions()");
      t_taggradient = tbox::TimerManager::getManager()->
         getTimer("apps::Euler::tagGradientDetectorCells()");
   }

   d_grid_geometry = grid_geom;
  
   d_use_nonuniform_workload = false;

   /*
    *hier::Variable<NDIM> quantities that define state of Euler problem.
    */
   d_density     = new pdat::CellVariable<NDIM,double>("density", 1);
   d_velocity    = new pdat::CellVariable<NDIM,double>("velocity", NDIM);
   d_pressure    = new pdat::CellVariable<NDIM,double>("pressure", 1);
   d_flux = new pdat::FaceVariable<NDIM,double>("flux", NEQU);

   /*
    * Default parameters for physical constants
    */

   d_gamma = 1.4;  // specific heat ratio for ideal diatomic gas (e.g., air)
 
   /*
    * Default parameters for numerical methods
    */

   d_riemann_solve = "APPROX_RIEM_SOLVE";
   d_godunov_order = 1;
   d_corner_transport = "CORNER_TRANSPORT_1";
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(CELLG == FACEG);
#endif
   d_nghosts = hier::IntVector<NDIM>(CELLG);
   d_fluxghosts = hier::IntVector<NDIM>(FLUXG); 

   /*
    * Defaults for problem type and initial data
    */

   d_radius = tbox::MathUtilities<double>::getSignalingNaN();
   tbox::MathUtilities<double>::setArrayToSignalingNaN(d_center, NDIM); 
   d_density_inside = tbox::MathUtilities<double>::getSignalingNaN();
   tbox::MathUtilities<double>::setArrayToSignalingNaN(d_velocity_inside, NDIM); 
   d_pressure_inside = tbox::MathUtilities<double>::getSignalingNaN();
   d_density_outside = tbox::MathUtilities<double>::getSignalingNaN();   
   tbox::MathUtilities<double>::setArrayToSignalingNaN(d_velocity_outside, NDIM); 
   d_pressure_outside = tbox::MathUtilities<double>::getSignalingNaN();   

   d_number_of_intervals = 0;
   d_front_position.resizeArray(0);
   d_interval_density.resizeArray(0);
   d_interval_velocity.resizeArray(0);
   d_interval_pressure.resizeArray(0);

   /*
    * Defaults for boundary conditions. Set to bogus values
    * for error checking.
    */

#if (NDIM == 2) 
   d_master_bdry_edge_conds.resizeArray(NUM_2D_EDGES);
   d_scalar_bdry_edge_conds.resizeArray(NUM_2D_EDGES);
   d_vector_bdry_edge_conds.resizeArray(NUM_2D_EDGES);
   for (int ei = 0; ei < NUM_2D_EDGES; ei++) {
      d_master_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
      d_scalar_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
      d_vector_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
   }

   d_master_bdry_node_conds.resizeArray(NUM_2D_NODES);
   d_scalar_bdry_node_conds.resizeArray(NUM_2D_NODES);
   d_vector_bdry_node_conds.resizeArray(NUM_2D_NODES);
   d_node_bdry_edge.resizeArray(NUM_2D_NODES);

   for (int ni = 0; ni < NUM_2D_NODES; ni++) {
      d_master_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
      d_scalar_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
      d_vector_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
      d_node_bdry_edge[ni] = BOGUS_BDRY_DATA;
   }

   d_bdry_edge_density.resizeArray(NUM_2D_EDGES);
   d_bdry_edge_velocity.resizeArray(NUM_2D_EDGES*NDIM);
   d_bdry_edge_pressure.resizeArray(NUM_2D_EDGES);
   tbox::MathUtilities<double>::setArrayToSignalingNaN(d_bdry_edge_density);
   tbox::MathUtilities<double>::setArrayToSignalingNaN(d_bdry_edge_velocity);
   tbox::MathUtilities<double>::setArrayToSignalingNaN(d_bdry_edge_pressure);
#endif
#if (NDIM == 3)
   d_master_bdry_face_conds.resizeArray(NUM_3D_FACES);
   d_scalar_bdry_face_conds.resizeArray(NUM_3D_FACES);
   d_vector_bdry_face_conds.resizeArray(NUM_3D_FACES);
   for (int fi = 0; fi < NUM_3D_FACES; fi++) {
      d_master_bdry_face_conds[fi] = BOGUS_BDRY_DATA;
      d_scalar_bdry_face_conds[fi] = BOGUS_BDRY_DATA;
      d_vector_bdry_face_conds[fi] = BOGUS_BDRY_DATA;
   }

   d_master_bdry_edge_conds.resizeArray(NUM_3D_EDGES);
   d_scalar_bdry_edge_conds.resizeArray(NUM_3D_EDGES);
   d_vector_bdry_edge_conds.resizeArray(NUM_3D_EDGES);
   d_edge_bdry_face.resizeArray(NUM_3D_EDGES);
   for (int ei = 0; ei < NUM_3D_EDGES; ei++) {
      d_master_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
      d_scalar_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
      d_vector_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
      d_edge_bdry_face[ei] = BOGUS_BDRY_DATA;
   }

   d_master_bdry_node_conds.resizeArray(NUM_3D_NODES);
   d_scalar_bdry_node_conds.resizeArray(NUM_3D_NODES);
   d_vector_bdry_node_conds.resizeArray(NUM_3D_NODES);
   d_node_bdry_face.resizeArray(NUM_3D_NODES);

   for (int ni = 0; ni < NUM_3D_NODES; ni++) {
      d_master_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
      d_scalar_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
      d_vector_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
      d_node_bdry_face[ni] = BOGUS_BDRY_DATA;
   }

   d_bdry_face_density.resizeArray(NUM_3D_FACES);
   d_bdry_face_velocity.resizeArray(NUM_3D_FACES*NDIM);
   d_bdry_face_pressure.resizeArray(NUM_3D_FACES);
   tbox::MathUtilities<double>::setArrayToSignalingNaN(d_bdry_face_density);
   tbox::MathUtilities<double>::setArrayToSignalingNaN(d_bdry_face_velocity);
   tbox::MathUtilities<double>::setArrayToSignalingNaN(d_bdry_face_pressure);
#endif

   /*
    * Initialize object with data read from given input/restart databases.
    */
   bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();
   if (is_from_restart) {
      getFromRestart();
   }
   getFromInput(input_db, is_from_restart);

   /*
    * Set problem data to values read from input/restart.
    */

   if (d_riemann_solve == "APPROX_RIEM_SOLVE") {
      d_riemann_solve_int = APPROX_RIEM_SOLVE;
   } else if (d_riemann_solve == "EXACT_RIEM_SOLVE") {
      d_riemann_solve_int = EXACT_RIEM_SOLVE;
   } else if (d_riemann_solve == "HLLC_RIEM_SOLVE") {
      d_riemann_solve_int = HLLC_RIEM_SOLVE;
   } else {
      TBOX_ERROR(d_object_name << ": "
         << "Unknown d_riemann_solve string = "
         << d_riemann_solve << " encountered in constructor" << endl);
   }

   if (d_data_problem == "PIECEWISE_CONSTANT_X") {
      d_data_problem_int = PIECEWISE_CONSTANT_X;
   } else if (d_data_problem == "PIECEWISE_CONSTANT_Y") {
      d_data_problem_int = PIECEWISE_CONSTANT_Y;
   } else if (d_data_problem == "PIECEWISE_CONSTANT_Z") {
      d_data_problem_int = PIECEWISE_CONSTANT_Z;
   } else if (d_data_problem == "SPHERE") {
      d_data_problem_int = SPHERE;
   } else if (d_data_problem == "STEP") {
      d_data_problem_int = STEP;
   } else {
      TBOX_ERROR(d_object_name << ": "
         << "Unknown d_data_problem string = "
         << d_data_problem << " encountered in constructor" << endl);
   }

   /*
    * Postprocess boundary data from input/restart values.
    */
#if (NDIM == 2) 
   for (int i = 0; i < NUM_2D_EDGES; i++) {
      d_scalar_bdry_edge_conds[i] = d_master_bdry_edge_conds[i];
      d_vector_bdry_edge_conds[i] = d_master_bdry_edge_conds[i];

      if (d_master_bdry_edge_conds[i] == REFLECT_BC) {
         d_scalar_bdry_edge_conds[i] = FLOW_BC;
      }
   }

   for (int i = 0; i < NUM_2D_NODES; i++) {
      d_scalar_bdry_node_conds[i] = d_master_bdry_node_conds[i];
      d_vector_bdry_node_conds[i] = d_master_bdry_node_conds[i];

      if (d_master_bdry_node_conds[i] == XREFLECT_BC) {
         d_scalar_bdry_node_conds[i] = XFLOW_BC;
      }
      if (d_master_bdry_node_conds[i] == YREFLECT_BC) {
         d_scalar_bdry_node_conds[i] = YFLOW_BC;
      }

      if (d_master_bdry_node_conds[i] != BOGUS_BDRY_DATA) {
         d_node_bdry_edge[i] =
            appu::CartesianBoundaryUtilities2::getEdgeLocationForNodeBdry(
                                            i, d_master_bdry_node_conds[i]);
      }
   }
#endif
#if (NDIM == 3)
   for (int i = 0; i < NUM_3D_FACES; i++) {
      d_scalar_bdry_face_conds[i] = d_master_bdry_face_conds[i];
      d_vector_bdry_face_conds[i] = d_master_bdry_face_conds[i];

      if (d_master_bdry_face_conds[i] == REFLECT_BC) {
         d_scalar_bdry_face_conds[i] = FLOW_BC;
      }
   }

   for (int i = 0; i < NUM_3D_EDGES; i++) {
      d_scalar_bdry_edge_conds[i] = d_master_bdry_edge_conds[i];
      d_vector_bdry_edge_conds[i] = d_master_bdry_edge_conds[i];

      if (d_master_bdry_edge_conds[i] == XREFLECT_BC) {
         d_scalar_bdry_edge_conds[i] = XFLOW_BC;
      }
      if (d_master_bdry_edge_conds[i] == YREFLECT_BC) {
         d_scalar_bdry_edge_conds[i] = YFLOW_BC;
      }
      if (d_master_bdry_edge_conds[i] == ZREFLECT_BC) {
         d_scalar_bdry_edge_conds[i] = ZFLOW_BC;
      }

      if (d_master_bdry_edge_conds[i] != BOGUS_BDRY_DATA) {
         d_edge_bdry_face[i] =
            appu::CartesianBoundaryUtilities3::getFaceLocationForEdgeBdry(
                                            i, d_master_bdry_edge_conds[i]);
      }
   }

   for (int i = 0; i < NUM_3D_NODES; i++) {
      d_scalar_bdry_node_conds[i] = d_master_bdry_node_conds[i];
      d_vector_bdry_node_conds[i] = d_master_bdry_node_conds[i];

      if (d_master_bdry_node_conds[i] == XREFLECT_BC) {
         d_scalar_bdry_node_conds[i] = XFLOW_BC;
      }
      if (d_master_bdry_node_conds[i] == YREFLECT_BC) {
         d_scalar_bdry_node_conds[i] = YFLOW_BC;
      }
      if (d_master_bdry_node_conds[i] == ZREFLECT_BC) {
         d_scalar_bdry_node_conds[i] = ZFLOW_BC;
      }

      if (d_master_bdry_node_conds[i] != BOGUS_BDRY_DATA) {
         d_node_bdry_face[i] =
            appu::CartesianBoundaryUtilities3::getFaceLocationForNodeBdry(
                                            i, d_master_bdry_node_conds[i]);
      }
   }

#endif

   stufprobc_(APPROX_RIEM_SOLVE,EXACT_RIEM_SOLVE,HLLC_RIEM_SOLVE,
              PIECEWISE_CONSTANT_X,PIECEWISE_CONSTANT_Y,PIECEWISE_CONSTANT_Z,
              SPHERE,STEP,
              CELLG,FACEG,FLUXG);
}


/*
*************************************************************************
*                                                                       *
* Empty destructor for Euler class.                                     *
*                                                                       *
*************************************************************************
*/

Euler::~Euler() 
{
   t_init = NULL;
   t_compute_dt = NULL;
   t_compute_fluxes = NULL;
   t_conservdiff = NULL;
   t_setphysbcs = NULL;
   t_taggradient = NULL;
} 

/*
*************************************************************************
*                                                                       *
* Register density, velocity, pressure (i.e., solution state variables),* 
* and flux variables with hyperbolic integrator that manages storage    *
* for those quantities.  Also, register plot data with the vis tool     *
* (Vizamrai or Visit).                                                  *
*                                                                       *
* Note that density coarsening/refining uses standard conservative      *
* operations provided in SAMRAI library.   Velocity and pressure        *
* are not conserved.  The Euler code provides operations to coarsen/    *
* refine momentum and total energy conservatively.  Velocity and        *
* pressure are calculated from the conserved quantities.                *
*                                                                       *
*************************************************************************
*/

void Euler::registerModelVariables(
   algs::HyperbolicLevelIntegrator<NDIM>* integrator)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(integrator != (algs::HyperbolicLevelIntegrator<NDIM>*)NULL);
   TBOX_ASSERT(CELLG == FACEG);
#endif

   integrator->registerVariable(d_density ,d_nghosts, 
                                algs::HyperbolicLevelIntegrator<NDIM>::TIME_DEP,
                                d_grid_geometry,
                                "CONSERVATIVE_COARSEN",
                                "CONSERVATIVE_LINEAR_REFINE");

   integrator->registerVariable(d_velocity, d_nghosts, 
                                algs::HyperbolicLevelIntegrator<NDIM>::TIME_DEP,
                                d_grid_geometry,
                                "USER_DEFINED_COARSEN",
                                "USER_DEFINED_REFINE");

   integrator->registerVariable(d_pressure, d_nghosts, 
                                algs::HyperbolicLevelIntegrator<NDIM>::TIME_DEP,
                                d_grid_geometry,
                                "USER_DEFINED_COARSEN",
                                "USER_DEFINED_REFINE");

   integrator->registerVariable(d_flux, d_fluxghosts, 
                                algs::HyperbolicLevelIntegrator<NDIM>::FLUX,
                                d_grid_geometry,
                                "CONSERVATIVE_COARSEN",
                                "NO_REFINE");

   hier::VariableDatabase<NDIM>* vardb = hier::VariableDatabase<NDIM>::getDatabase();

   d_plot_context = integrator->getPlotContext();

   if (!(d_vizamrai_writer.isNull())) {
      
      d_vizamrai_writer->registerPlotScalar("Density",
                                       vardb->mapVariableAndContextToIndex(
                                          d_density, d_plot_context));
      
      d_vizamrai_writer->registerPlotVector("Velocity",
                                       vardb->mapVariableAndContextToIndex(
                                          d_velocity, d_plot_context));

      d_vizamrai_writer->registerPlotScalar("Pressure",
                                       vardb->mapVariableAndContextToIndex(
                                          d_pressure, d_plot_context));
      
      d_vizamrai_writer->registerDerivedPlotScalar("Total Energy", this);
      d_vizamrai_writer->registerDerivedPlotVector("Momentum", NDIM, this);
   }

#ifdef HAVE_HDF5   
   if (!(d_visit_writer.isNull())) {
      d_visit_writer->registerPlotQuantity("Density",
                                           "SCALAR",
                                           vardb->mapVariableAndContextToIndex(
                                              d_density, d_plot_context));
      
      d_visit_writer->registerPlotQuantity("Velocity",
                                           "VECTOR",
                                           vardb->mapVariableAndContextToIndex(
                                              d_velocity, d_plot_context));

      d_visit_writer->registerPlotQuantity("Pressure",
                                           "SCALAR",
                                           vardb->mapVariableAndContextToIndex(
                                              d_pressure, d_plot_context));
      
      d_visit_writer->registerDerivedPlotQuantity("Total Energy", 
                                                  "SCALAR",
                                                  this);
      d_visit_writer->registerDerivedPlotQuantity("Momentum", 
                                                  "VECTOR",
                                                  this);
   }

   if (d_vizamrai_writer.isNull() && d_visit_writer.isNull()) {
      TBOX_WARNING(d_object_name << ": registerModelVariables()\n"
		   << "Neither a Vizamrai nor Visit data writer was\n"
		   << "registered.  Consequently, no plot data will\n"
		   << "be written." << endl);
   }
#else
   if (d_vizamrai_writer.isNull()) {
      TBOX_WARNING(d_object_name << ": registerModelVariables()\n"
		   << "A Vizamrai Visit data writer was\n"
		   << "registered.  Consequently, no plot data will\n"
		   << "be written." << endl);
   }      
#endif
	 
}

/*
*************************************************************************
*                                                                       *
* Set up parameters for nonuniform load balancing, if used.             *
*                                                                       *
*************************************************************************
*/

void Euler::setupLoadBalancer(algs::HyperbolicLevelIntegrator<NDIM>* integrator,
                              mesh::GriddingAlgorithm<NDIM>* gridding_algorithm)
{
   (void) integrator;

   hier::VariableDatabase<NDIM>* vardb = hier::VariableDatabase<NDIM>::getDatabase();

   if (d_use_nonuniform_workload && gridding_algorithm) {
       tbox::Pointer<mesh::LoadBalancer<NDIM> > load_balancer =  
          gridding_algorithm->getLoadBalanceStrategy();

       if (!load_balancer.isNull()) {
          d_workload_variable = new pdat::CellVariable<NDIM,double>("workload_variable", 1);
          d_workload_data_id = 
             vardb->registerVariableAndContext(d_workload_variable,
                                               vardb->getContext("WORKLOAD"));
          load_balancer->setWorkloadPatchDataIndex(d_workload_data_id);
          vardb->registerPatchDataForRestart(d_workload_data_id);
       } else {
          TBOX_WARNING(d_object_name << ": "
             << "  Unknown load balancer used in gridding algorithm." 
             << "  Ignoring request for nonuniform load balancing." << endl);
          d_use_nonuniform_workload = false;
       }
   } else {
      d_use_nonuniform_workload = false; 
   }

}

/*
*************************************************************************
*                                                                       *
* Set initial data for solution variables on patch interior.            *
* This routine is called whenever a new patch is introduced to the      *
* AMR patch hierarchy.  Note that the routine does nothing unless       *
* we are at the initial time.  In all other cases, conservative         *
* interpolation from coarser levels and copies from patches at the      *
* same mesh resolution are sufficient to set data.                      * 
*                                                                       *
*************************************************************************
*/

void Euler::initializeDataOnPatch(hier::Patch<NDIM>& patch,
                                  const double data_time,
                                  const bool initial_time)
{
   (void) data_time;

   t_init->start();

   if (initial_time) {

      const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
      const double* dx  = pgeom->getDx();
      const double* xlo = pgeom->getXLower();
      const double* xhi = pgeom->getXUpper();

      tbox::Pointer< pdat::CellData<NDIM,double> > density  = 
          patch.getPatchData(d_density, getDataContext());
      tbox::Pointer< pdat::CellData<NDIM,double> > velocity = 
          patch.getPatchData(d_velocity, getDataContext());
      tbox::Pointer< pdat::CellData<NDIM,double> > pressure = 
          patch.getPatchData(d_pressure, getDataContext());

#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(!density.isNull());
      TBOX_ASSERT(!velocity.isNull());
      TBOX_ASSERT(!pressure.isNull());
#endif
      hier::IntVector<NDIM> ghost_cells = density->getGhostCellWidth();
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(velocity->getGhostCellWidth() == ghost_cells);
      TBOX_ASSERT(pressure->getGhostCellWidth() == ghost_cells);
#endif

      const hier::Index<NDIM> ifirst=patch.getBox().lower();
      const hier::Index<NDIM> ilast =patch.getBox().upper();

      if (d_data_problem == "SPHERE") {

         eulerinitsphere_(d_data_problem_int,
                          dx,xlo,xhi,
                          ifirst(0),ilast(0),
                          ifirst(1),ilast(1),
#if (NDIM>2)
                          ifirst(2),ilast(2),
#endif
                          ghost_cells(0),
                          ghost_cells(1),
#if (NDIM>2)
                          ghost_cells(2),
#endif
                          d_gamma,
                          density->getPointer(),
                          velocity->getPointer(),
                          pressure->getPointer(),
                          d_density_inside,
                          d_velocity_inside,
                          d_pressure_inside,
                          d_density_outside,
                          d_velocity_outside,
                          d_pressure_outside,
                          d_center, d_radius);

      } else {

         eulerinit_(d_data_problem_int,
                    dx,xlo,xhi,
                    ifirst(0),ilast(0),
                    ifirst(1),ilast(1),
#if (NDIM>2)
                    ifirst(2),ilast(2),
#endif
                    ghost_cells(0),
                    ghost_cells(1),
#if (NDIM>2)
                    ghost_cells(2),
#endif
                    d_gamma,
                    density->getPointer(),
                    velocity->getPointer(),
                    pressure->getPointer(),
                    d_number_of_intervals,
                    d_front_position.getPointer(),
                    d_interval_density.getPointer(),
                    d_interval_velocity.getPointer(),
                    d_interval_pressure.getPointer());

      }

   }

   if (d_use_nonuniform_workload) {
      if (!patch.checkAllocated(d_workload_data_id)) {
         patch.allocatePatchData(d_workload_data_id);
      }
      tbox::Pointer< pdat::CellData<NDIM,double> > workload_data =
         patch.getPatchData(d_workload_data_id);
      workload_data->fillAll(1.0);
   }

   t_init->stop();

}

/*
*************************************************************************
*                                                                       *
* Compute stable time increment for patch.  Return this value.          *
*                                                                       *
*************************************************************************
*/

double Euler::computeStableDtOnPatch(
   hier::Patch<NDIM>& patch,
   const bool initial_time, 
   const double dt_time) 
{
   (void) initial_time;
   (void) dt_time;

   t_compute_dt->start();

   const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > patch_geom = patch.getPatchGeometry();
   const double* dx  = patch_geom->getDx();
 
   const hier::Index<NDIM> ifirst=patch.getBox().lower();
   const hier::Index<NDIM> ilast =patch.getBox().upper();

   const tbox::Pointer< pdat::CellData<NDIM,double> > density  = 
      patch.getPatchData(d_density, getDataContext());
   const tbox::Pointer< pdat::CellData<NDIM,double> > velocity = 
      patch.getPatchData(d_velocity, getDataContext());
   const tbox::Pointer< pdat::CellData<NDIM,double> > pressure = 
      patch.getPatchData(d_pressure, getDataContext());

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!density.isNull());
   TBOX_ASSERT(!velocity.isNull());
   TBOX_ASSERT(!pressure.isNull());
#endif 
   hier::IntVector<NDIM> ghost_cells = density->getGhostCellWidth();
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(velocity->getGhostCellWidth() == ghost_cells);
   TBOX_ASSERT(pressure->getGhostCellWidth() == ghost_cells);
#endif 

   double stabdt = 0.;
   stabledt_(dx,
             ifirst(0),ilast(0),
             ifirst(1),ilast(1),
#if (NDIM>2)
             ifirst(2),ilast(2),
#endif
             ghost_cells(0),
             ghost_cells(1),
#if (NDIM>2)
             ghost_cells(2),
#endif
             d_gamma,
             density->getPointer(),
             velocity->getPointer(),
             pressure->getPointer(),
             stabdt);

   t_compute_dt->stop();
   return stabdt;
}

/*
*************************************************************************
*                                                                       *
* Compute time integral of numerical fluxes for finite difference       *
* at each cell face on patch.  When NDIM == 3, there are two options    *
* for the transverse flux correction.  Otherwise, there is only one.    * 
*                                                                       *
*************************************************************************
*/

void Euler::computeFluxesOnPatch(hier::Patch<NDIM>& patch,
                                 const double time, 
                                 const double dt)
{
   (void) time;

   t_compute_fluxes->start();

#if (NDIM == 3)
   if (d_corner_transport == "CORNER_TRANSPORT_2") {
      compute3DFluxesWithCornerTransport2(patch, dt);
   } else {
      compute3DFluxesWithCornerTransport1(patch, dt);
   }
#endif  

#if (NDIM == 2)

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(CELLG == FACEG);
#endif

   const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > patch_geom = patch.getPatchGeometry();
   const double* dx  = patch_geom->getDx();

   hier::Box<NDIM> pbox = patch.getBox();
   const hier::Index<NDIM> ifirst=pbox.lower();
   const hier::Index<NDIM> ilast =pbox.upper();

   tbox::Pointer< pdat::CellData<NDIM,double> > density  = 
      patch.getPatchData(d_density, getDataContext());
   tbox::Pointer< pdat::CellData<NDIM,double> > velocity = 
      patch.getPatchData(d_velocity, getDataContext());
   tbox::Pointer< pdat::CellData<NDIM,double> > pressure = 
      patch.getPatchData(d_pressure, getDataContext());
   tbox::Pointer< pdat::FaceData<NDIM,double> > flux     = 
      patch.getPatchData(d_flux, getDataContext());

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!density.isNull());
   TBOX_ASSERT(!velocity.isNull());
   TBOX_ASSERT(!pressure.isNull());
   TBOX_ASSERT(!flux.isNull());
   TBOX_ASSERT(density->getGhostCellWidth() == d_nghosts);
   TBOX_ASSERT(velocity->getGhostCellWidth() == d_nghosts);
   TBOX_ASSERT(pressure->getGhostCellWidth() == d_nghosts);
   TBOX_ASSERT(flux->getGhostCellWidth() == d_fluxghosts);
#endif

   /*
    * Allocate patch data for temporaries local to this routine.
    */
   pdat::FaceData<NDIM,double> traced_left(pbox, NEQU, d_nghosts);
   pdat::FaceData<NDIM,double> traced_right(pbox, NEQU, d_nghosts);
   pdat::CellData<NDIM,double> sound_speed(pbox, 1, d_nghosts);

   /*
    *  Initialize traced states (w^R and w^L) with proper cell-centered values.
    */
   inittraceflux_(ifirst(0),ilast(0),
                  ifirst(1),ilast(1),
                  density->getPointer(),
                  velocity->getPointer(),
                  pressure->getPointer(),
                  traced_left.getPointer(0),
                  traced_left.getPointer(1),
                  traced_right.getPointer(0),
                  traced_right.getPointer(1),
                  flux->getPointer(0),
                  flux->getPointer(1));
   
   /*
    * If Godunov method requires slopes with order greater than one, perform
    * characteristic tracing to compute higher-order slopes.
    */
   if (d_godunov_order > 1) {

      /*
       * Prepare temporary data for characteristic tracing.
       */
      int Mcells = 0;
      for (int k=0;k<NDIM;k++) {
         Mcells = tbox::MathUtilities<int>::Max(Mcells, pbox.numberCells(k));
      }

      // Face-centered temporary arrays
      tbox::Array<double> ttedgslp((2*FACEG+1+Mcells)*NEQU);
      tbox::Array<double> ttraclft((2*FACEG+1+Mcells)*NEQU);
      tbox::Array<double> ttracrgt((2*FACEG+1+Mcells)*NEQU);

      // Cell-centered temporary arrays
      tbox::Array<double> ttsound((2*CELLG+Mcells));
      tbox::Array<double> ttcelslp((2*CELLG+Mcells)*NEQU);

      /*
       * Compute local sound speed in each computational cell.
       */
      computesound_(ifirst(0),ilast(0),
                    ifirst(1),ilast(1),
                    d_gamma,
                    density->getPointer(),
                    velocity->getPointer(),
                    pressure->getPointer(),
                    sound_speed.getPointer());

      /*
       *  Apply characteristic tracing to compute initial estimate of
       *  traces w^L and w^R at faces.
       *  Inputs: sound_speed, w^L, w^R (traced_left/right)
       *  Output: w^L, w^R
       */
      chartracing0_(dt,
                    ifirst(0),ilast(0),
                    ifirst(1),ilast(1),
                    Mcells, dx[0], d_gamma, d_godunov_order,
                    sound_speed.getPointer(),
                    traced_left.getPointer(0),
                    traced_right.getPointer(0),
                    ttcelslp.getPointer(),
                    ttedgslp.getPointer(),
                    ttsound.getPointer() ,
                    ttraclft.getPointer(),
                    ttracrgt.getPointer());

      chartracing1_(dt,
                    ifirst(0),ilast(0),
                    ifirst(1),ilast(1),
                    Mcells, dx[1], d_gamma, d_godunov_order,
                    sound_speed.getPointer(),
                    traced_left.getPointer(1),
                    traced_right.getPointer(1),
                    ttcelslp.getPointer(),
                    ttedgslp.getPointer(),
                    ttsound.getPointer() ,
                    ttraclft.getPointer(),
                    ttracrgt.getPointer());

   }  // if (d_godunov_order > 1) ...

// fluxcalculation_(dt,*,1,dx, to get artificial viscosity
// fluxcalculation_(dt,*,0,dx, to get NO artificial viscosity

   /*
    *  Compute preliminary fluxes at faces by solving approximate
    *  Riemann problem using the trace states computed so far.
    *  Inputs: P, rho, v, w^L, w^R (traced_left/right)
    *  Output: F (flux)
    */
   fluxcalculation_(dt,1,0,dx,
                    ifirst(0),ilast(0),ifirst(1),ilast(1),
                    d_gamma,
                    d_riemann_solve_int,
                    density->getPointer(),
                    velocity->getPointer(),
                    pressure->getPointer(),
                    flux->getPointer(0),
                    flux->getPointer(1),
                    traced_left.getPointer(0),
                    traced_left.getPointer(1),
                    traced_right.getPointer(0),
                    traced_right.getPointer(1));

   /*
    *  Update trace states at cell faces with transverse correction applied.
    *  Inputs: F (flux)
    *  Output: w^L, w^R (traced_left/right)
    */
   fluxcorrec_(dt,
               ifirst(0),ilast(0),ifirst(1),ilast(1),
               dx, d_gamma,
               density->getPointer(),
               velocity->getPointer(),
               pressure->getPointer(),
               flux->getPointer(0),
               flux->getPointer(1),
               traced_left.getPointer(0),
               traced_left.getPointer(1),
               traced_right.getPointer(0),
               traced_right.getPointer(1));

   boundaryReset(patch, traced_left, traced_right);

   /*
    *  Re-compute fluxes with updated trace states.
    *  Inputs: w^L, w^R (traced_left/right)
    *  Output: F (flux)
    */
   fluxcalculation_(dt,0,0,dx,
                    ifirst(0),ilast(0),ifirst(1),ilast(1),
                    d_gamma,
                    d_riemann_solve_int,
                    density->getPointer(),
                    velocity->getPointer(),
                    pressure->getPointer(),
                    flux->getPointer(0),
                    flux->getPointer(1),
                    traced_left.getPointer(0),
                    traced_left.getPointer(1),
                    traced_right.getPointer(0),
                    traced_right.getPointer(1));

#endif

   t_compute_fluxes->stop();
}

/*
*************************************************************************
*                                                                       *
* Compute numerical approximations to flux terms using an extension     *
* to three dimensions of Collella's corner transport upwind approach.   *
* I.E. input value corner_transport = "CORNER_TRANSPORT_1"              *
*                                                                       *
*************************************************************************
*/

#if (NDIM == 3)
void Euler::compute3DFluxesWithCornerTransport1(hier::Patch<NDIM>& patch,
                                                const double dt)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(CELLG == FACEG);
#endif

   const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > patch_geom = patch.getPatchGeometry();
   const double* dx  = patch_geom->getDx();

   hier::Box<NDIM> pbox = patch.getBox();
   const hier::Index<NDIM> ifirst=pbox.lower();
   const hier::Index<NDIM> ilast =pbox.upper();

   tbox::Pointer< pdat::CellData<NDIM,double> > density  = 
      patch.getPatchData(d_density, getDataContext());
   tbox::Pointer< pdat::CellData<NDIM,double> > velocity = 
      patch.getPatchData(d_velocity, getDataContext());
   tbox::Pointer< pdat::CellData<NDIM,double> > pressure = 
      patch.getPatchData(d_pressure, getDataContext());
   tbox::Pointer< pdat::FaceData<NDIM,double> > flux     = 
      patch.getPatchData(d_flux, getDataContext());

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!density.isNull());
   TBOX_ASSERT(!velocity.isNull());
   TBOX_ASSERT(!pressure.isNull());
   TBOX_ASSERT(!flux.isNull());
   TBOX_ASSERT(density->getGhostCellWidth() == d_nghosts);
   TBOX_ASSERT(velocity->getGhostCellWidth() == d_nghosts);
   TBOX_ASSERT(pressure->getGhostCellWidth() == d_nghosts);
   TBOX_ASSERT(flux->getGhostCellWidth() == d_fluxghosts);
#endif

   /*
    * Allocate patch data for temporaries local to this routine.
    */
   pdat::FaceData<NDIM,double> traced_left(pbox, NEQU, d_nghosts);
   pdat::FaceData<NDIM,double> traced_right(pbox, NEQU, d_nghosts);
   pdat::CellData<NDIM,double> sound_speed(pbox, 1, d_nghosts);
   pdat::FaceData<NDIM,double> temp_flux(pbox, NEQU, d_fluxghosts);
   pdat::FaceData<NDIM,double> temp_traced_left(pbox, NEQU, d_nghosts);
   pdat::FaceData<NDIM,double> temp_traced_right(pbox, NEQU, d_nghosts);

   /*
    *  Initialize traced states (w^R and w^L) with proper cell-centered values.
    */
   inittraceflux_(ifirst(0),ilast(0),
                  ifirst(1),ilast(1),
                  ifirst(2),ilast(2),
                  density->getPointer(),
                  velocity->getPointer(),
                  pressure->getPointer(),
                  traced_left.getPointer(0),
                  traced_left.getPointer(1),
                  traced_left.getPointer(2),
                  traced_right.getPointer(0),
                  traced_right.getPointer(1),
                  traced_right.getPointer(2),
                  flux->getPointer(0),
                  flux->getPointer(1),
                  flux->getPointer(2));

   /*
    * If Godunov method requires slopes with order greater than one, perform
    * characteristic tracing to compute higher-order slopes. 
    */
   if (d_godunov_order > 1) {

      /*
       * Prepare temporary data for characteristic tracing.
       */
      int Mcells = 0;
      for (int k=0;k<NDIM;k++) {
         Mcells = tbox::MathUtilities<int>::Max(Mcells, pbox.numberCells(k));
      }

      // Face-centered temporary arrays
      tbox::Array<double>  ttedgslp((2*FACEG+1+Mcells)*NEQU);
      tbox::Array<double>  ttraclft((2*FACEG+1+Mcells)*NEQU);
      tbox::Array<double>  ttracrgt((2*FACEG+1+Mcells)*NEQU);

      // Cell-centered temporary arrays
      tbox::Array<double>  ttsound((2*CELLG+Mcells));
      tbox::Array<double>  ttcelslp((2*CELLG+Mcells)*NEQU); 

      /*
       * Compute local sound speed in each computational cell.
       */
      computesound_(ifirst(0),ilast(0),
                    ifirst(1),ilast(1),
                    ifirst(2),ilast(2),
                    d_gamma, 
                    density->getPointer(), 
                    velocity->getPointer(),
                    pressure->getPointer(),
                    sound_speed.getPointer());

      /*
       *  Apply characteristic tracing to compute initial estimate of 
       *  traces w^L and w^R at faces.
       *  Inputs: sound_speed, w^L, w^R (traced_left/right)
       *  Output: w^L, w^R
       */
      chartracing0_(dt,
                    ifirst(0),ilast(0),
                    ifirst(1),ilast(1),
                    ifirst(2),ilast(2),
                    Mcells, dx[0], d_gamma, d_godunov_order,
                    sound_speed.getPointer(),
                    traced_left.getPointer(0),
                    traced_right.getPointer(0),
                    ttcelslp.getPointer(),
                    ttedgslp.getPointer(), 
                    ttsound.getPointer() ,
                    ttraclft.getPointer(),
                    ttracrgt.getPointer());

      chartracing1_(dt,
                    ifirst(0),ilast(0),
                    ifirst(1),ilast(1),
                    ifirst(2),ilast(2),
                    Mcells, dx[1], d_gamma, d_godunov_order,
                    sound_speed.getPointer(),
                    traced_left.getPointer(1),
                    traced_right.getPointer(1),
                    ttcelslp.getPointer(),
                    ttedgslp.getPointer(),
                    ttsound.getPointer() ,
                    ttraclft.getPointer(),
                    ttracrgt.getPointer());

      chartracing2_(dt,
                    ifirst(0),ilast(0),
                    ifirst(1),ilast(1),
                    ifirst(2),ilast(2),
                    Mcells, dx[2], d_gamma, d_godunov_order,
                    sound_speed.getPointer(),
                    traced_left.getPointer(2),
                    traced_right.getPointer(2),
                    ttcelslp.getPointer(),
                    ttedgslp.getPointer(),
                    ttsound.getPointer() ,
                    ttraclft.getPointer(),
                    ttracrgt.getPointer());

   }  // if (d_godunov_order > 1) ...

   /*
    *  Compute preliminary fluxes at faces by solving approximate 
    *  Riemann problem using the trace states computed so far.
    *  Inputs: P, rho, v, w^L, w^R (traced_left/right)
    *  Output: F (flux)
    */
//  fluxcalculation_(dt,*,*,1,dx,  to do artificial viscosity
//  fluxcalculation_(dt,*,*,0,dx,  to do NO artificial viscosity
   fluxcalculation_(dt,1,0,0,dx,
                    ifirst(0),ilast(0),ifirst(1),ilast(1),ifirst(2),ilast(2),
                    d_gamma,
                    d_riemann_solve_int,
                    density->getPointer(),
                    velocity->getPointer(),
                    pressure->getPointer(),
                    flux->getPointer(0),
                    flux->getPointer(1),
                    flux->getPointer(2),
                    traced_left.getPointer(0),
                    traced_left.getPointer(1),
                    traced_left.getPointer(2),
                    traced_right.getPointer(0),
                    traced_right.getPointer(1),
                    traced_right.getPointer(2));

   /*
    *  Re-compute face traces to include one set of correction terms with 
    *  transverse flux differences.  Store result in temporary vectors 
    *  (i.e. temp_traced_left/right).  
    *  Inputs: F (flux), P, rho, v, w^L, w^R (traced_left/right)
    *  Output: temp_traced_left/right
    */
   fluxcorrec2d_(dt,
                 ifirst(0),ilast(0),ifirst(1),ilast(1),ifirst(2),ilast(2),
                 dx,d_gamma,1,
                 density->getPointer(),
                 velocity->getPointer(),
                 pressure->getPointer(),
                 flux->getPointer(0),
                 flux->getPointer(1),
                 flux->getPointer(2),
                 traced_left.getPointer(0),
                 traced_left.getPointer(1),
                 traced_left.getPointer(2),
                 traced_right.getPointer(0),
                 traced_right.getPointer(1),
                 traced_right.getPointer(2),
                 temp_traced_left.getPointer(0),
                 temp_traced_left.getPointer(1),
                 temp_traced_left.getPointer(2),
                 temp_traced_right.getPointer(0),
                 temp_traced_right.getPointer(1),
                 temp_traced_right.getPointer(2));
 
   boundaryReset(patch, traced_left, traced_right);

   /*
    *  Compute fluxes with partially-corrected trace states.  Store 
    *  result in temporary flux vector. 
    *  Inputs: P, rho, v, temp_traced_left/right
    *  Output: temp_flux 
    */
   fluxcalculation_(dt,0,1,0,dx,
                    ifirst(0),ilast(0),ifirst(1),ilast(1),ifirst(2),ilast(2),
                    d_gamma,
                    d_riemann_solve_int,
                    density->getPointer(),
                    velocity->getPointer(),
                    pressure->getPointer(),
                    temp_flux.getPointer(0),
                    temp_flux.getPointer(1),
                    temp_flux.getPointer(2),
                    temp_traced_left.getPointer(0),
                    temp_traced_left.getPointer(1),
                    temp_traced_left.getPointer(2),
                    temp_traced_right.getPointer(0),
                    temp_traced_right.getPointer(1),
                    temp_traced_right.getPointer(2));
   /*
    *  Compute face traces with other transverse correction flux
    *  difference terms included.  Store result in temporary vectors 
    *  (i.e. temp_traced_left/right).
    *  Inputs: F (flux), P, rho, v, w^L, w^R (traced_left/right)
    *  Output: temp_traced_left/right
    */
   fluxcorrec2d_(dt,
                 ifirst(0),ilast(0),ifirst(1),ilast(1),ifirst(2),ilast(2),
                 dx,d_gamma,-1,
                 density->getPointer(),
                 velocity->getPointer(),
                 pressure->getPointer(),
                 flux->getPointer(0),
                 flux->getPointer(1),
                 flux->getPointer(2),
                 traced_left.getPointer(0),
                 traced_left.getPointer(1),
                 traced_left.getPointer(2),
                 traced_right.getPointer(0),
                 traced_right.getPointer(1),
                 traced_right.getPointer(2),
                 temp_traced_left.getPointer(0),
                 temp_traced_left.getPointer(1),
                 temp_traced_left.getPointer(2),
                 temp_traced_right.getPointer(0),
                 temp_traced_right.getPointer(1),
                 temp_traced_right.getPointer(2));
 
   boundaryReset(patch, traced_left, traced_right);

   /*
    *  Compute final predicted fluxes with both sets of transverse flux
    *  differences included.  Store the result in regular flux vector.  
    *  NOTE:  the fact that we store  these fluxes in the regular (i.e. 
    *  not temporary) flux vector does NOT indicate this is the final result.
    *  Rather, the flux vector is used as a convenient storage location.
    *  Inputs: P, rho, v, temp_traced_left/right
    *  Output: flux 
    */
   fluxcalculation_(dt,1,0,0,dx,
                    ifirst(0),ilast(0),ifirst(1),ilast(1),ifirst(2),ilast(2),
                    d_gamma,
                    d_riemann_solve_int,
                    density->getPointer(),
                    velocity->getPointer(),
                    pressure->getPointer(),
                    flux->getPointer(0),
                    flux->getPointer(1),
                    flux->getPointer(2),
                    temp_traced_left.getPointer(0),
                    temp_traced_left.getPointer(1),
                    temp_traced_left.getPointer(2),
                    temp_traced_right.getPointer(0),
                    temp_traced_right.getPointer(1),
                    temp_traced_right.getPointer(2));
   /*
    *  Compute the final trace state vectors at cell faces using transverse
    *  differences of final predicted fluxes.  Store the result in w^L 
    *  (traced_left) and w^R (traced_right) vectors.
    *  Inputs: temp_flux, flux
    *  Output: w^L, w^R (traced_left/right)
    */
   fluxcorrec3d_(dt,
                 ifirst(0),ilast(0),ifirst(1),ilast(1),ifirst(2),ilast(2),
                 dx,d_gamma,
                 density->getPointer(),
                 velocity->getPointer(),
                 pressure->getPointer(),
                 temp_flux.getPointer(0),
                 temp_flux.getPointer(1),
                 temp_flux.getPointer(2),
                 flux->getPointer(0),
                 flux->getPointer(1),
                 flux->getPointer(2),
                 traced_left.getPointer(0),
                 traced_left.getPointer(1),
                 traced_left.getPointer(2),
                 traced_right.getPointer(0),
                 traced_right.getPointer(1),
                 traced_right.getPointer(2));

   /*
    *  Final flux calculation using corrected trace states.  
    *  Inputs:  w^L, w^R (traced_left/right)
    *  Output:  F (flux)
    */
   fluxcalculation_(dt,0,0,0,dx,
                    ifirst(0),ilast(0),ifirst(1),ilast(1),ifirst(2),ilast(2),
                    d_gamma,
                    d_riemann_solve_int,
                    density->getPointer(),
                    velocity->getPointer(),
                    pressure->getPointer(),
                    flux->getPointer(0),
                    flux->getPointer(1),
                    flux->getPointer(2),
                    traced_left.getPointer(0),
                    traced_left.getPointer(1),
                    traced_left.getPointer(2),
                    traced_right.getPointer(0),
                    traced_right.getPointer(1),
                    traced_right.getPointer(2)); 
}
#endif

/*
*************************************************************************
*                                                                       *
* Compute numerical approximations to flux terms using John             *
* Trangenstein's interpretation of the three-dimensional version of     *
* Collella's corner transport upwind approach.                          *
* I.E. input value corner_transport = "CORNER_TRANSPORT_2"              *
*                                                                       *
*************************************************************************
*/

#if (NDIM==3)
void Euler::compute3DFluxesWithCornerTransport2(hier::Patch<NDIM>& patch,
                                                const double dt)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(CELLG == FACEG);
#endif

   const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > patch_geom = patch.getPatchGeometry();
   const double* dx  = patch_geom->getDx();

   hier::Box<NDIM> pbox = patch.getBox();
   const hier::Index<NDIM> ifirst=pbox.lower();
   const hier::Index<NDIM> ilast =pbox.upper();

   tbox::Pointer< pdat::CellData<NDIM,double> > density      = 
      patch.getPatchData(d_density, getDataContext());
   tbox::Pointer< pdat::CellData<NDIM,double> > velocity     = 
      patch.getPatchData(d_velocity, getDataContext());
   tbox::Pointer< pdat::CellData<NDIM,double> > pressure     = 
      patch.getPatchData(d_pressure, getDataContext());
   tbox::Pointer< pdat::FaceData<NDIM,double> > flux         = 
      patch.getPatchData(d_flux, getDataContext());

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!density.isNull());
   TBOX_ASSERT(!velocity.isNull());
   TBOX_ASSERT(!pressure.isNull());
   TBOX_ASSERT(!flux.isNull());
   TBOX_ASSERT(density->getGhostCellWidth() == d_nghosts);
   TBOX_ASSERT(velocity->getGhostCellWidth() == d_nghosts);
   TBOX_ASSERT(pressure->getGhostCellWidth() == d_nghosts);
   TBOX_ASSERT(flux->getGhostCellWidth() == d_fluxghosts);
#endif

   /*
    * Allocate patch data for temporaries local to this routine.
    */
   pdat::FaceData<NDIM,double> traced_left(pbox, NEQU, d_nghosts);
   pdat::FaceData<NDIM,double> traced_right(pbox, NEQU, d_nghosts);
   pdat::CellData<NDIM,double> sound_speed(pbox, 1, d_nghosts);
   pdat::FaceData<NDIM,double> temp_flux(pbox, NEQU, d_fluxghosts);
   pdat::CellData<NDIM,double> third_state(pbox, NEQU, d_nghosts); 

   /*
    *  Initialize traced states (w^R and w^L) with proper cell-centered values.
    */
   inittraceflux_(ifirst(0),ilast(0),
                  ifirst(1),ilast(1),
                  ifirst(2),ilast(2),
                  density->getPointer(),
                  velocity->getPointer(),
                  pressure->getPointer(),
                  traced_left.getPointer(0),
                  traced_left.getPointer(1),
                  traced_left.getPointer(2),
                  traced_right.getPointer(0),
                  traced_right.getPointer(1),
                  traced_right.getPointer(2),
                  flux->getPointer(0),
                  flux->getPointer(1),
                  flux->getPointer(2));

   /*
    *  Compute fluxes at faces by solving approximate Riemann problem
    *  using initial trace states.
    *  Inputs: P, rho, v, w^L, w^R (traced_left/right)
    *  Output: F (flux)
    */
//  fluxcalculation_(dt,*,*,1,dx,  to do artificial viscosity
//  fluxcalculation_(dt,*,*,0,dx,  to do NO artificial viscosity
   fluxcalculation_(dt,1,1,0,dx,
                    ifirst(0),ilast(0),ifirst(1),ilast(1),ifirst(2),ilast(2),
                    d_gamma,
                    d_riemann_solve_int,
                    density->getPointer(),
                    velocity->getPointer(),
                    pressure->getPointer(),
                    flux->getPointer(0),
                    flux->getPointer(1),
                    flux->getPointer(2),
                    traced_left.getPointer(0),
                    traced_left.getPointer(1),
                    traced_left.getPointer(2),
                    traced_right.getPointer(0),
                    traced_right.getPointer(1),
                    traced_right.getPointer(2));

   /*
    * If Godunov method requires slopes with order greater than one, perform
    * characteristic tracing to compute higher-order slopes.
    */
   if (d_godunov_order > 1) {

      /*
       * Prepare temporary data for characteristic tracing.
       */
      int Mcells = 0;
      for (int k=0;k<NDIM;k++) {
         Mcells = tbox::MathUtilities<int>::Max(Mcells, pbox.numberCells(k));
      }

      // Face-centered temporary arrays
      tbox::Array<double> ttedgslp((2*FACEG+1+Mcells)*NEQU);
      tbox::Array<double> ttraclft((2*FACEG+1+Mcells)*NEQU);
      tbox::Array<double> ttracrgt((2*FACEG+1+Mcells)*NEQU);

      // Cell-centered temporary arrays
      tbox::Array<double> ttsound((2*CELLG+Mcells));
      tbox::Array<double> ttcelslp((2*CELLG+Mcells)*NEQU);

      /*
       * Compute local sound speed in each computational cell.
       */
      computesound_(ifirst(0),ilast(0),
                    ifirst(1),ilast(1),
                    ifirst(2),ilast(2),
                    d_gamma,
                    density->getPointer(),
                    velocity->getPointer(),
                    pressure->getPointer(),
                    sound_speed.getPointer());

     /*
      *  Apply characteristic tracing to update traces w^L and w^R at faces.
      *  Inputs: sound_speed, w^L, w^R (traced_left/right)
      *  Output: w^L, w^R
      */
      chartracing0_(dt,
                    ifirst(0),ilast(0),
                    ifirst(1),ilast(1),
                    ifirst(2),ilast(2),
                    Mcells, dx[0], d_gamma, d_godunov_order,
                    sound_speed.getPointer(),
                    traced_left.getPointer(0),
                    traced_right.getPointer(0),
                    ttcelslp.getPointer(),
                    ttedgslp.getPointer(),
                    ttsound.getPointer() ,
                    ttraclft.getPointer(),
                    ttracrgt.getPointer());

      chartracing1_(dt,
                    ifirst(0),ilast(0),
                    ifirst(1),ilast(1),
                    ifirst(2),ilast(2),
                    Mcells, dx[1], d_gamma, d_godunov_order,
                    sound_speed.getPointer(),
                    traced_left.getPointer(1),
                    traced_right.getPointer(1),
                    ttcelslp.getPointer(),
                    ttedgslp.getPointer(),
                    ttsound.getPointer() ,
                    ttraclft.getPointer(),
                    ttracrgt.getPointer());

      chartracing2_(dt,
                    ifirst(0),ilast(0),
                    ifirst(1),ilast(1),
                    ifirst(2),ilast(2),
                    Mcells, dx[2], d_gamma, d_godunov_order,
                    sound_speed.getPointer(),
                    traced_left.getPointer(2),
                    traced_right.getPointer(2),
                    ttcelslp.getPointer(),
                    ttedgslp.getPointer(),
                    ttsound.getPointer() ,
                    ttraclft.getPointer(),
                    ttracrgt.getPointer());

   } // if (d_godunov_order > 1) ...

   for (int idir = 0; idir < NDIM; idir++ ) {

      /*
       *    Approximate traces at cell centers (in idir direction); 
       * i.e.,  "1/3 state".
       *    Inputs:  F (flux), rho, v, P
       *    Output:  third_state
       */
      onethirdstate_(dt,dx,idir, 
                     ifirst(0),ilast(0),ifirst(1),ilast(1),ifirst(2),ilast(2),
                     d_gamma,
                     density->getPointer(),
                     velocity->getPointer(),
                     pressure->getPointer(),
                     flux->getPointer(0),
                     flux->getPointer(1),
                     flux->getPointer(2),
                     third_state.getPointer());

      /*
       *    Compute fluxes using 1/3 state traces, in the two directions OTHER
       *    than idir.
       *    Inputs:  third_state, rho, v, P
       *    Output:  temp_flux (only directions other than idir are modified)  
       */
      fluxthird_(dt,dx,idir,
                 ifirst(0),ilast(0),ifirst(1),ilast(1),ifirst(2),ilast(2),
                 d_gamma,
                 d_riemann_solve_int,
                 density->getPointer(),
                 velocity->getPointer(),
                 pressure->getPointer(),
                 third_state.getPointer(),
                 temp_flux.getPointer(0),
                 temp_flux.getPointer(1),
                 temp_flux.getPointer(2));
      /*
       *    Compute transverse corrections for the traces in the two 
       *    directions (OTHER than idir) using the flux differences 
       *    computed in those directions.
       *    Inputs:  temp_flux, rho, v, P 
       *    Output:  w^L, w^R (traced_left/right)
       */
      fluxcorrecjt_(dt,dx,idir,
                    ifirst(0),ilast(0),ifirst(1),ilast(1),ifirst(2),ilast(2),
                    d_gamma,
                    density->getPointer(),
                    velocity->getPointer(),
                    pressure->getPointer(),
                    temp_flux.getPointer(0),
                    temp_flux.getPointer(1),
                    temp_flux.getPointer(2),
                    traced_left.getPointer(0),
                    traced_left.getPointer(1),
                    traced_left.getPointer(2),
                    traced_right.getPointer(0),
                    traced_right.getPointer(1),
                    traced_right.getPointer(2));

   }  // loop over directions...

   boundaryReset(patch, traced_left, traced_right);
 
   /*
    *  Final flux calculation using corrected trace states.  
    *  Inputs:  w^L, w^R (traced_left/right)
    *  Output:  F (flux)
    */
   fluxcalculation_(dt,0,0,0,dx,
                    ifirst(0),ilast(0),ifirst(1),ilast(1),ifirst(2),ilast(2),
                    d_gamma,
                    d_riemann_solve_int,
                    density->getPointer(),
                    velocity->getPointer(),
                    pressure->getPointer(),
                    flux->getPointer(0),
                    flux->getPointer(1),
                    flux->getPointer(2),
                    traced_left.getPointer(0),
                    traced_left.getPointer(1),
                    traced_left.getPointer(2),
                    traced_right.getPointer(0),
                    traced_right.getPointer(1),
                    traced_right.getPointer(2));

}
#endif

/*
*************************************************************************
*                                                                       *
* Update Euler solution variables by performing a conservative          *
* difference with the fluxes calculated in computeFluxesOnPatch().      *
* Although, "primitive" variables are maintained (i.e., density,        *
* velocity, pressure), "conserved" variables (i.e., density,            *
* momentum, total energy) are conserved.                                *
*                                                                       *
*************************************************************************
*/

void Euler::conservativeDifferenceOnPatch(hier::Patch<NDIM>& patch,
                                          const double time,
                                          const double dt,
                                          bool at_syncronization)
{
   (void) time;
   (void) dt;
   (void) at_syncronization;

   t_conservdiff->start();

   const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > patch_geom = patch.getPatchGeometry();
   const double* dx  = patch_geom->getDx();

   hier::Box<NDIM> pbox = patch.getBox();
   const hier::Index<NDIM> ifirst=pbox.lower();
   const hier::Index<NDIM> ilast =pbox.upper();
 
   tbox::Pointer< pdat::FaceData<NDIM,double> > flux        = 
      patch.getPatchData(d_flux, getDataContext());
   tbox::Pointer< pdat::CellData<NDIM,double> > density     = 
      patch.getPatchData(d_density, getDataContext());
   tbox::Pointer< pdat::CellData<NDIM,double> > velocity    = 
      patch.getPatchData(d_velocity, getDataContext());
   tbox::Pointer< pdat::CellData<NDIM,double> > pressure    = 
      patch.getPatchData(d_pressure, getDataContext());

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!density.isNull());
   TBOX_ASSERT(!velocity.isNull());
   TBOX_ASSERT(!pressure.isNull());
   TBOX_ASSERT(!flux.isNull());
   TBOX_ASSERT(density->getGhostCellWidth() == d_nghosts);
   TBOX_ASSERT(velocity->getGhostCellWidth() == d_nghosts);
   TBOX_ASSERT(pressure->getGhostCellWidth() == d_nghosts);
   TBOX_ASSERT(flux->getGhostCellWidth() == d_fluxghosts);
#endif

#if (NDIM==2)
   consdiff_(ifirst(0),ilast(0),ifirst(1),ilast(1),dx,
             flux->getPointer(0),
             flux->getPointer(1),
             d_gamma,
             density->getPointer(),
             velocity->getPointer(), 
             pressure->getPointer());
#endif
#if (NDIM==3)
   consdiff_(ifirst(0),ilast(0),ifirst(1),ilast(1),ifirst(2),ilast(2),dx,
             flux->getPointer(0),
             flux->getPointer(1),
             flux->getPointer(2),
             d_gamma,
             density->getPointer(),
             velocity->getPointer(), 
             pressure->getPointer());
#endif

   t_conservdiff->stop();

}

/*
*************************************************************************
*                                                                       *
* Reset physical boundary values for special cases, such as those       *
* involving reflective boundary conditions and when the "STEP"          *
* problem is run.                                                       *
*                                                                       *
*************************************************************************
*/

void Euler::boundaryReset(hier::Patch<NDIM>& patch,
                          pdat::FaceData<NDIM,double>& traced_left,  
                          pdat::FaceData<NDIM,double>& traced_right) const
{
   const hier::Index<NDIM> ifirst  =patch.getBox().lower();
   const hier::Index<NDIM> ilast   =patch.getBox().upper();
   int i,idir;
   bool bdry_cell = true;

   const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > patch_geom = patch.getPatchGeometry();
   hier::BoxArray<NDIM> domain_boxes;
   d_grid_geometry->computePhysicalDomain(domain_boxes, patch_geom->getRatio());
   int num_domain_boxes = domain_boxes.getNumberOfBoxes();
   const double* dx = patch_geom->getDx();
   const double* xpatchhi = patch_geom->getXUpper();
   const double* xdomainhi = d_grid_geometry->getXUpper();

   pdat::CellIndex<NDIM> icell = ifirst;
   hier::BoxArray<NDIM> bdrybox(2*NDIM);
   hier::Index<NDIM> ibfirst = ifirst;
   hier::Index<NDIM> iblast  = ilast;
   int bdry_case;

   for (idir=0;idir<NDIM; idir++) {
      ibfirst(idir) = ifirst(idir) -1;
      iblast(idir)  = ifirst(idir) -1;
      bdrybox.getBox(2*idir) = hier::Box<NDIM>(ibfirst,iblast);

      ibfirst(idir) = ilast(idir) +1;
      iblast(idir)  = ilast(idir) +1;
      bdrybox.getBox(2*idir+1) = hier::Box<NDIM>(ibfirst,iblast);
   }

   for (idir=0;idir<NDIM; idir++) {
      int bside = 2*idir;
#if (NDIM == 2) 
      bdry_case = d_master_bdry_edge_conds[bside];
#endif
#if (NDIM == 3)
      bdry_case = d_master_bdry_face_conds[bside];
#endif
      if (bdry_case == REFLECT_BC) {
         for (pdat::CellIterator<NDIM> ic(bdrybox.getBox(bside)); ic; ic++) {
            for ( i = 0; i < num_domain_boxes; i++ ) {
               if (domain_boxes.getBox(i).contains(ic()))
                  bdry_cell = false;
            }
            if (bdry_cell) {
              pdat::FaceIndex<NDIM> sidein = pdat::FaceIndex<NDIM>(ic(),idir,1);
              (traced_left)(sidein,0) = (traced_right)(sidein,0);
            }
         }
      }

      int bnode = 2*idir+1;
#if (NDIM == 2) 
      bdry_case = d_master_bdry_edge_conds[bnode];
#endif
#if (NDIM == 3)
      bdry_case = d_master_bdry_face_conds[bnode];
#endif
// BEGIN SIMPLE-MINDED FIX FOR STEP PROBLEM
      if ((d_data_problem == "STEP") && (bnode == 1) && 
          (tbox::MathUtilities<double>::Abs(xpatchhi[0]-xdomainhi[0]) < dx[0]) )
      {
          bdry_case = FLOW_BC;
      }
// END SIMPLE-MINDED FIX FOR STEP PROBLEM
      if (bdry_case == REFLECT_BC) {
         for (pdat::CellIterator<NDIM> ic(bdrybox.getBox(bnode)); ic; ic++) {
            for ( i = 0; i < num_domain_boxes; i++ ) {
               if (domain_boxes.getBox(i).contains(ic()))
                  bdry_cell = false;
            }
            if (bdry_cell) {
              pdat::FaceIndex<NDIM> sidein = pdat::FaceIndex<NDIM>(ic(),idir,0);
              (traced_right)(sidein,0) = (traced_left)(sidein,0);
            }
         }
      }
   }
}

/*
*************************************************************************
*                                                                       *
* Refine velocity and pressure by conservatively refining               *
* momentum and total energy.                                            *
*                                                                       *
*************************************************************************
*/

hier::IntVector<NDIM> Euler::getRefineOpStencilWidth() const
{
   return(hier::IntVector<NDIM>(1));
}

void Euler::postprocessRefine(hier::Patch<NDIM>& fine,
                              const hier::Patch<NDIM>& coarse,
                              const hier::Box<NDIM>& fine_box,
                              const hier::IntVector<NDIM>& ratio)
{

   tbox::Pointer< pdat::CellData<NDIM,double> > cdensity     = 
      coarse.getPatchData(d_density, getDataContext());
   tbox::Pointer< pdat::CellData<NDIM,double> > cvelocity    = 
      coarse.getPatchData(d_velocity, getDataContext());
   tbox::Pointer< pdat::CellData<NDIM,double> > cpressure    = 
      coarse.getPatchData(d_pressure, getDataContext());

   tbox::Pointer< pdat::CellData<NDIM,double> > fdensity     = 
      fine.getPatchData(d_density, getDataContext());
   tbox::Pointer< pdat::CellData<NDIM,double> > fvelocity    = 
      fine.getPatchData(d_velocity, getDataContext());
   tbox::Pointer< pdat::CellData<NDIM,double> > fpressure    = 
      fine.getPatchData(d_pressure, getDataContext());

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!cdensity.isNull());
   TBOX_ASSERT(!cvelocity.isNull());
   TBOX_ASSERT(!cpressure.isNull());
   TBOX_ASSERT(!fdensity.isNull());
   TBOX_ASSERT(!fvelocity.isNull());
   TBOX_ASSERT(!fpressure.isNull());

   hier::IntVector<NDIM> gccheck = cdensity->getGhostCellWidth();
   TBOX_ASSERT(cvelocity->getGhostCellWidth() == gccheck);
   TBOX_ASSERT(cpressure->getGhostCellWidth() == gccheck);

   gccheck = fdensity->getGhostCellWidth();
   TBOX_ASSERT(fvelocity->getGhostCellWidth() == gccheck);
   TBOX_ASSERT(fpressure->getGhostCellWidth() == gccheck);
#endif

   const hier::Box<NDIM> cgbox(cdensity->getGhostBox());

   const hier::Index<NDIM> cilo = cgbox.lower();
   const hier::Index<NDIM> cihi = cgbox.upper();
   const hier::Index<NDIM> filo = fdensity->getGhostBox().lower();
   const hier::Index<NDIM> fihi = fdensity->getGhostBox().upper();

   const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > cgeom = coarse.getPatchGeometry();
   const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > fgeom = fine.getPatchGeometry();

   const hier::Box<NDIM> coarse_box = hier::Box<NDIM>::coarsen(fine_box, ratio);
   const hier::Index<NDIM> ifirstc = coarse_box.lower();
   const hier::Index<NDIM> ilastc = coarse_box.upper();
   const hier::Index<NDIM> ifirstf = fine_box.lower();
   const hier::Index<NDIM> ilastf = fine_box.upper();

   const hier::IntVector<NDIM> cons_ghosts(1);
   pdat::CellData<NDIM,double> conserved(coarse_box, 1, cons_ghosts);

   const hier::IntVector<NDIM> tmp_ghosts(0);

   double* diff0 = new double[coarse_box.numberCells(0)+1];
   pdat::CellData<NDIM,double> slope0(coarse_box, 1, tmp_ghosts);

   double* diff1 = new double[coarse_box.numberCells(1)+1];
   pdat::CellData<NDIM,double> slope1(coarse_box, 1, tmp_ghosts);

#if (NDIM == 3)
   double* diff2 = new double[coarse_box.numberCells(2)+1];
   pdat::CellData<NDIM,double> slope2(coarse_box, 1, tmp_ghosts);
#endif

#if (NDIM == 2)
   pdat::CellData<NDIM,double> flat0(coarse_box, 1, tmp_ghosts);
   pdat::CellData<NDIM,double> flat1(coarse_box, 1, tmp_ghosts);
   int mc = cihi(0)-cilo(0) + 1;
   mc = tbox::MathUtilities<int>::Max(mc,cihi(1)-cilo(1) + 1);
   double* tflat  = new double[mc];
   double* tflat2 = new double[mc];
   double* tsound = new double[mc];
   double* tdensc = new double[mc];
   double* tpresc = new double[mc];
   double* tvelc  = new double[mc];
   conservlinint2d_(ifirstc(0),ifirstc(1),ilastc(0),ilastc(1),  /* input */
                    ifirstf(0),ifirstf(1),ilastf(0),ilastf(1),
                    cilo(0),cilo(1),cihi(0),cihi(1),
                    filo(0),filo(1),fihi(0),fihi(1),
                    ratio,
                    cgeom->getDx(),
                    fgeom->getDx(),
                    d_gamma,
                    cdensity->getPointer(),
                    fdensity->getPointer(),
                    cvelocity->getPointer(),
                    cpressure->getPointer(),
                    fvelocity->getPointer(),                   /* output */
                    fpressure->getPointer(),
                    conserved.getPointer(),                   /* temporaries */
                    tflat,tflat2,tsound,mc,
                    tdensc,tpresc,tvelc,
                    flat0.getPointer(),
                    flat1.getPointer(),
                    diff0,slope0.getPointer(),
                    diff1,slope1.getPointer());
#endif
#if (NDIM == 3)
   pdat::CellData<NDIM,double> flat0(coarse_box, 1, tmp_ghosts);
   pdat::CellData<NDIM,double> flat1(coarse_box, 1, tmp_ghosts);
   pdat::CellData<NDIM,double> flat2(coarse_box, 1, tmp_ghosts);
   int mc = cihi(0)-cilo(0) + 1;
   mc = tbox::MathUtilities<int>::Max(mc,cihi(1)-cilo(1) + 1);
   mc = tbox::MathUtilities<int>::Max(mc,cihi(2)-cilo(2) + 1);
   double* tflat  = new double[mc];
   double* tflat2 = new double[mc];
   double* tsound = new double[mc];
   double* tdensc = new double[mc];
   double* tpresc = new double[mc];
   double* tvelc  = new double[mc];
   conservlinint3d_(ifirstc(0),ifirstc(1),ifirstc(2),         /* input */
                    ilastc(0),ilastc(1),ilastc(2),
                    ifirstf(0),ifirstf(1),ifirstf(2),
                    ilastf(0),ilastf(1),ilastf(2),
                    cilo(0),cilo(1),cilo(2),cihi(0),cihi(1),cihi(2),
                    filo(0),filo(1),filo(2),fihi(0),fihi(1),fihi(2),
                    ratio,
                    cgeom->getDx(),
                    fgeom->getDx(),
                    d_gamma,
                    cdensity->getPointer(),
                    fdensity->getPointer(),
                    cvelocity->getPointer(),
                    cpressure->getPointer(),
                    fvelocity->getPointer(),                   /* output */
                    fpressure->getPointer(),
                    conserved.getPointer(),                   /* temporaries */
                    tflat,tflat2,tsound,mc,
                    tdensc,tpresc,tvelc,
                    flat0.getPointer(),
                    flat1.getPointer(),
                    flat2.getPointer(),
                    diff0,slope0.getPointer(),
                    diff1,slope1.getPointer(),
                    diff2,slope2.getPointer());
#endif
   delete [] tflat;
   delete [] tflat2;
   delete [] tsound;
   delete [] tdensc;
   delete [] tpresc;
   delete [] tvelc;

   delete [] diff0;
   delete [] diff1;
#if (NDIM == 3)
   delete [] diff2;
#endif

}

/*
*************************************************************************
*                                                                       *
* Coarsen velocity and pressure by conservatively coarsening            *
* momentum and total energy.                                            *
*                                                                       *
*************************************************************************
*/

hier::IntVector<NDIM> Euler::getCoarsenOpStencilWidth() const
{
   return(hier::IntVector<NDIM>(0));
}

void Euler::postprocessCoarsen(hier::Patch<NDIM>& coarse,
                               const hier::Patch<NDIM>& fine,
                               const hier::Box<NDIM>& coarse_box,
                               const hier::IntVector<NDIM>& ratio)
{

   tbox::Pointer< pdat::CellData<NDIM,double> > fdensity     = fine.getPatchData( 
                                              d_density, getDataContext());
   tbox::Pointer< pdat::CellData<NDIM,double> > fvelocity    = fine.getPatchData( 
                                              d_velocity, getDataContext());
   tbox::Pointer< pdat::CellData<NDIM,double> > fpressure    = fine.getPatchData( 
                                              d_pressure, getDataContext());
   tbox::Pointer< pdat::CellData<NDIM,double> > cdensity     = coarse.getPatchData( 
                                              d_density, getDataContext());
   tbox::Pointer< pdat::CellData<NDIM,double> > cvelocity    = coarse.getPatchData( 
                                              d_velocity, getDataContext());
   tbox::Pointer< pdat::CellData<NDIM,double> > cpressure    = coarse.getPatchData( 
                                              d_pressure, getDataContext());

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!cdensity.isNull());
   TBOX_ASSERT(!cvelocity.isNull());
   TBOX_ASSERT(!cpressure.isNull());
   TBOX_ASSERT(!fdensity.isNull());
   TBOX_ASSERT(!fvelocity.isNull());
   TBOX_ASSERT(!fpressure.isNull());

   hier::IntVector<NDIM> gccheck = cdensity->getGhostCellWidth();
   TBOX_ASSERT(cvelocity->getGhostCellWidth() == gccheck);
   TBOX_ASSERT(cpressure->getGhostCellWidth() == gccheck);

   gccheck = fdensity->getGhostCellWidth();
   TBOX_ASSERT(fvelocity->getGhostCellWidth() == gccheck);
   TBOX_ASSERT(fpressure->getGhostCellWidth() == gccheck);
#endif

   const hier::Index<NDIM> filo = fdensity->getGhostBox().lower();
   const hier::Index<NDIM> fihi = fdensity->getGhostBox().upper();
   const hier::Index<NDIM> cilo = cdensity->getGhostBox().lower();
   const hier::Index<NDIM> cihi = cdensity->getGhostBox().upper();

   const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > fgeom = fine.getPatchGeometry();
   const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > cgeom = coarse.getPatchGeometry();

   const hier::Box<NDIM> fine_box = hier::Box<NDIM>::refine(coarse_box, ratio);
   const hier::Index<NDIM> ifirstc = coarse_box.lower();
   const hier::Index<NDIM> ilastc = coarse_box.upper();
   const hier::Index<NDIM> ifirstf = fine_box.lower();
   const hier::Index<NDIM> ilastf = fine_box.upper();

   const hier::IntVector<NDIM> cons_ghosts(0);
   pdat::CellData<NDIM,double> conserved(fine_box, 1, cons_ghosts);

#if (NDIM == 2)
   conservavg2d_(ifirstf(0),ifirstf(1),ilastf(0),ilastf(1),  /* input */
                 ifirstc(0),ifirstc(1),ilastc(0),ilastc(1),
                 filo(0),filo(1),fihi(0),fihi(1),
                 cilo(0),cilo(1),cihi(0),cihi(1),
                 ratio,
                 fgeom->getDx(),
                 cgeom->getDx(),
                 d_gamma,
                 fdensity->getPointer(),
                 cdensity->getPointer(),
                 fvelocity->getPointer(),
                 fpressure->getPointer(),
                 cvelocity->getPointer(),                    /* output */
                 cpressure->getPointer(),
                 conserved.getPointer());                   /* temporary */
#endif
#if (NDIM == 3)
   conservavg3d_(ifirstf(0),ifirstf(1),ifirstf(2),           /* input */
                 ilastf(0),ilastf(1),ilastf(2),
                 ifirstc(0),ifirstc(1),ifirstc(2),
                 ilastc(0),ilastc(1),ilastc(2),
                 filo(0),filo(1),filo(2),fihi(0),fihi(1),fihi(2),
                 cilo(0),cilo(1),cilo(2),cihi(0),cihi(1),cihi(2),
                 ratio,
                 fgeom->getDx(),
                 cgeom->getDx(),
                 d_gamma,
                 fdensity->getPointer(),
                 cdensity->getPointer(),
                 fvelocity->getPointer(),
                 fpressure->getPointer(),
                 cvelocity->getPointer(),                    /* output */
                 cpressure->getPointer(),
                 conserved.getPointer());
#endif

}

/*
*************************************************************************
*                                                                       *
* Set the data in ghost cells corresponding to physical boundary        *
* conditions.  Note that boundary geometry configuration information    *
* (i.e., faces, edges, and nodes) is obtained from the patch geometry   *
* object owned by the patch.                                            *
*                                                                       *
*************************************************************************
*/

void Euler::setPhysicalBoundaryConditions(
   hier::Patch<NDIM>& patch, 
   const double fill_time,
   const hier::IntVector<NDIM>& ghost_width_to_fill)
{
   (void) fill_time;
   t_setphysbcs->start();

   tbox::Pointer< pdat::CellData<NDIM,double> > density  =
      patch.getPatchData(d_density, getDataContext());
   tbox::Pointer< pdat::CellData<NDIM,double> > velocity =
      patch.getPatchData(d_velocity, getDataContext());
   tbox::Pointer< pdat::CellData<NDIM,double> > pressure =
      patch.getPatchData(d_pressure, getDataContext());

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!density.isNull());
   TBOX_ASSERT(!velocity.isNull());
   TBOX_ASSERT(!pressure.isNull());
#endif
   hier::IntVector<NDIM> ghost_cells = density->getGhostCellWidth();
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(velocity->getGhostCellWidth() == ghost_cells);
   TBOX_ASSERT(pressure->getGhostCellWidth() == ghost_cells);
#endif

#if (NDIM == 2) 

   /*
    * Set boundary conditions for cells corresponding to patch edges.
    *
    * Note: We apply a simple-minded adjustment for the "STEP" problem
    *       so that the right edge of the domain gets (out)FLOW conditions
    *       whereas the right edge at the step gets REFLECT condtions (from input),
    */
   tbox::Array<int> tmp_edge_scalar_bcond(NUM_2D_EDGES);
   tbox::Array<int> tmp_edge_vector_bcond(NUM_2D_EDGES);
   for (int i = 0; i < NUM_2D_EDGES; i++) {
      tmp_edge_scalar_bcond[i] = d_scalar_bdry_edge_conds[i];
      tmp_edge_vector_bcond[i] = d_vector_bdry_edge_conds[i];
   }

   if (d_data_problem == "STEP") {

      const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > patch_geom = 
         patch.getPatchGeometry();
      const double* dx = patch_geom->getDx();
      const double* xpatchhi = patch_geom->getXUpper();
      const double* xdomainhi = d_grid_geometry->getXUpper();

      if (tbox::MathUtilities<double>::Abs(xpatchhi[0]-xdomainhi[0]) < dx[0]) {
         tmp_edge_scalar_bcond[XHI] = FLOW_BC;
         tmp_edge_vector_bcond[XHI] = FLOW_BC;
      }

   }

   appu::CartesianBoundaryUtilities2::
      fillEdgeBoundaryData("density", density,
                           patch,
                           ghost_width_to_fill,
                           tmp_edge_scalar_bcond,
                           d_bdry_edge_density);
   appu::CartesianBoundaryUtilities2::
      fillEdgeBoundaryData("velocity", velocity,
                           patch,
                           ghost_width_to_fill,
                           tmp_edge_vector_bcond,
                           d_bdry_edge_velocity);
   appu::CartesianBoundaryUtilities2::
      fillEdgeBoundaryData("pressure", pressure,
                           patch,
                           ghost_width_to_fill,
                           tmp_edge_scalar_bcond,
                           d_bdry_edge_pressure);

#ifdef DEBUG_CHECK_ASSERTIONS
#if CHECK_BDRY_DATA
   checkBoundaryData(EDGE2D_BDRY_TYPE, patch, ghost_width_to_fill,
                     tmp_edge_scalar_bcond, tmp_edge_vector_bcond);
#endif
#endif

   /*
    *  Set boundary conditions for cells corresponding to patch nodes.
    */

   appu::CartesianBoundaryUtilities2::
      fillNodeBoundaryData("density", density,
                           patch,
                           ghost_width_to_fill,
                           d_scalar_bdry_node_conds,
                           d_bdry_edge_density);
   appu::CartesianBoundaryUtilities2::
      fillNodeBoundaryData("velocity", velocity,
                           patch,
                           ghost_width_to_fill,
                           d_vector_bdry_node_conds,
                           d_bdry_edge_velocity);
   appu::CartesianBoundaryUtilities2::
      fillNodeBoundaryData("pressure", pressure,
                           patch,
                           ghost_width_to_fill,
                           d_scalar_bdry_node_conds,
                           d_bdry_edge_pressure);


#ifdef DEBUG_CHECK_ASSERTIONS
#if CHECK_BDRY_DATA
   checkBoundaryData(NODE2D_BDRY_TYPE, patch, ghost_width_to_fill,
                     d_scalar_bdry_node_conds, d_vector_bdry_node_conds);
#endif
#endif

#endif // NDIM == 2

#if (NDIM == 3)

   /*
    *  Set boundary conditions for cells corresponding to patch faces.
    */

   appu::CartesianBoundaryUtilities3::
      fillFaceBoundaryData("density", density,
                           patch,
                           ghost_width_to_fill,
                           d_scalar_bdry_face_conds,
                           d_bdry_face_density);
   appu::CartesianBoundaryUtilities3::
      fillFaceBoundaryData("velocity", velocity,
                           patch,
                           ghost_width_to_fill,
                           d_vector_bdry_face_conds,
                           d_bdry_face_velocity);
   appu::CartesianBoundaryUtilities3::
      fillFaceBoundaryData("pressure", pressure,
                           patch,
                           ghost_width_to_fill,
                           d_scalar_bdry_face_conds,
                           d_bdry_face_pressure);

#ifdef DEBUG_CHECK_ASSERTIONS
#if CHECK_BDRY_DATA
   checkBoundaryData(FACE3D_BDRY_TYPE, patch, ghost_width_to_fill,
                     d_scalar_bdry_face_conds, d_vector_bdry_face_conds);
#endif
#endif

   /*
    *  Set boundary conditions for cells corresponding to patch edges.
    */

   appu::CartesianBoundaryUtilities3::
      fillEdgeBoundaryData("density", density,
                           patch,
                           ghost_width_to_fill,
                           d_scalar_bdry_edge_conds,
                           d_bdry_face_density);
   appu::CartesianBoundaryUtilities3::
      fillEdgeBoundaryData("velocity", velocity,
                           patch,
                           ghost_width_to_fill,
                           d_vector_bdry_edge_conds,
                           d_bdry_face_velocity);
   appu::CartesianBoundaryUtilities3::
      fillEdgeBoundaryData("pressure", pressure,
                           patch,
                           ghost_width_to_fill,
                           d_scalar_bdry_edge_conds,
                           d_bdry_face_pressure);

#ifdef DEBUG_CHECK_ASSERTIONS
#if CHECK_BDRY_DATA
   checkBoundaryData(EDGE3D_BDRY_TYPE, patch, ghost_width_to_fill,
                     d_scalar_bdry_edge_conds, d_vector_bdry_edge_conds);
#endif
#endif

   /*
    *  Set boundary conditions for cells corresponding to patch nodes.
    */

   appu::CartesianBoundaryUtilities3::
      fillNodeBoundaryData("density", density,
                           patch,
                           ghost_width_to_fill,
                           d_scalar_bdry_node_conds,
                           d_bdry_face_density);
   appu::CartesianBoundaryUtilities3::
      fillNodeBoundaryData("velocity", velocity,
                           patch,
                           ghost_width_to_fill,
                           d_vector_bdry_node_conds,
                           d_bdry_face_velocity);
   appu::CartesianBoundaryUtilities3::
      fillNodeBoundaryData("pressure", pressure,
                           patch,
                           ghost_width_to_fill,
                           d_scalar_bdry_node_conds,
                           d_bdry_face_pressure);

#ifdef DEBUG_CHECK_ASSERTIONS
#if CHECK_BDRY_DATA
   checkBoundaryData(NODE3D_BDRY_TYPE, patch, ghost_width_to_fill,
                     d_scalar_bdry_node_conds, d_scalar_bdry_node_conds);
#endif
#endif

#endif // NDIM == 3

   t_setphysbcs->stop();
}

/*
*************************************************************************
*                                                                       *
* Tag cells for refinement using gradient detector.  Tagging criteria   *
* defined in input.                                                     *
*                                                                       *
*************************************************************************
*/

void Euler::tagGradientDetectorCells(hier::Patch<NDIM>& patch,
   const double regrid_time,
   const bool initial_error,
   const int tag_indx,
   const bool uses_richardson_extrapolation_too)
{
   t_taggradient->start();

   const int error_level_number = patch.getPatchLevelNumber();

   const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > patch_geom = patch.getPatchGeometry();
   const double* dx  = patch_geom->getDx();

   tbox::Pointer< pdat::CellData<NDIM,int> > tags        = patch.getPatchData(tag_indx);

   hier::Box<NDIM> pbox = patch.getBox(); 
   hier::Box<NDIM> pboxm1 = pbox.grow(pbox,-1);
   hier::BoxArray<NDIM> domain_boxes;
   d_grid_geometry->computePhysicalDomain(domain_boxes, patch_geom->getRatio());
   /*
    * Construct domain bounding box
    */
   hier::Box<NDIM> domain;
   for( int i=0; i < domain_boxes.getNumberOfBoxes(); i++ ) {
     domain += domain_boxes(i);
   }

   const hier::Index<NDIM> domfirst=domain.lower();
   const hier::Index<NDIM> domlast =domain.upper();
   const hier::Index<NDIM> ifirst=patch.getBox().lower();
   const hier::Index<NDIM> ilast =patch.getBox().upper();

   /*
    * Create a set of temporary tags and set to untagged value.
    */
 
   tbox::Pointer< pdat::CellData<NDIM,int> > temp_tags  = new pdat::CellData<NDIM,int>(pbox, 1, d_nghosts);
   temp_tags->fillAll(FALSE);

#if (NDIM==2)
   /*
    * Problem specific criteria for step case.
    */
   if (initial_error && d_data_problem == "STEP") { 
      if (error_level_number < 2) { 
         hier::Box<NDIM> tagbox(hier::Index<NDIM>(9,0), hier::Index<NDIM>(9,3));
         if (error_level_number == 1) {
            tagbox.refine(hier::IntVector<NDIM>(2));
         }
         hier::Box<NDIM> ibox = pbox * tagbox;
     
         for (pdat::CellIterator<NDIM> itc(ibox); itc; itc++) {
            (*temp_tags)(itc(),0) = TRUE; 
         }
      }
   }
#endif

   /*
    * Possible tagging criteria includes 
    *    DENSITY_DEVIATION, DENSITY_GRADIENT, DENSITY_SHOCK
    *    PRESSURE_DEVIATION, PRESSURE_GRADIENT, PRESSURE_SHOCK
    * The criteria is specified over a time interval.
    *
    * Loop over criteria provided and check to make sure we are in the
    * specified time interval.  If so, apply appropriate tagging for
    * the level.
    */
   for (int ncrit = 0; ncrit < d_refinement_criteria.getSize(); ncrit++) {

      string ref = d_refinement_criteria[ncrit];
      tbox::Pointer< pdat::CellData<NDIM, double > > var;     
      int size = 0;
      double tol = 0.;
      double onset = 0.;
      double dev = 0.;
      bool time_allowed = false;

      if (ref == "DENSITY_DEVIATION") {
         var = patch.getPatchData(d_density, getDataContext());
         size = d_density_dev_tol.getSize();
         tol = ( ( error_level_number < size) 
                 ? d_density_dev_tol[error_level_number] 
                 : d_density_dev_tol[size-1] );
         size = d_density_dev.getSize();
         dev = ( ( error_level_number < size) 
                 ? d_density_dev[error_level_number] 
                 : d_density_dev[size-1] );        
         size = d_density_dev_time_min.getSize();
         double time_min = ( ( error_level_number < size) 
                 ? d_density_dev_time_min[error_level_number] 
                 : d_density_dev_time_min[size-1] );
         size = d_density_dev_time_max.getSize();
         double time_max = ( ( error_level_number < size) 
                 ? d_density_dev_time_max[error_level_number] 
                 : d_density_dev_time_max[size-1] );
         time_allowed = (time_min <= regrid_time) && (time_max > regrid_time);
      }

      if (ref == "DENSITY_GRADIENT") {
         var = patch.getPatchData(d_density, getDataContext());
         size = d_density_grad_tol.getSize();
         tol = ( ( error_level_number < size) 
                 ? d_density_grad_tol[error_level_number] 
                 : d_density_grad_tol[size-1] );
         size = d_density_grad_time_min.getSize();
         double time_min = ( ( error_level_number < size) 
                 ? d_density_grad_time_min[error_level_number] 
                 : d_density_grad_time_min[size-1] );
         size = d_density_grad_time_max.getSize();
         double time_max = ( ( error_level_number < size) 
                 ? d_density_grad_time_max[error_level_number] 
                 : d_density_grad_time_max[size-1] );
         time_allowed = (time_min <= regrid_time) && (time_max > regrid_time);
      }

      if (ref == "DENSITY_SHOCK") {
         var = patch.getPatchData(d_density, getDataContext());
         size = d_density_shock_tol.getSize();
         tol = ( ( error_level_number < size) 
                 ? d_density_shock_tol[error_level_number] 
                 : d_density_shock_tol[size-1] );
         size = d_density_shock_onset.getSize();
         onset = ( ( error_level_number < size) 
                 ? d_density_shock_onset[error_level_number] 
                 : d_density_shock_onset[size-1] );
         size = d_density_shock_time_min.getSize();
         double time_min = ( ( error_level_number < size) 
                 ? d_density_shock_time_min[error_level_number] 
                 : d_density_shock_time_min[size-1] );
         size = d_density_shock_time_max.getSize();
         double time_max = ( ( error_level_number < size) 
                 ? d_density_shock_time_max[error_level_number] 
                 : d_density_shock_time_max[size-1] );
         time_allowed = (time_min <= regrid_time) && (time_max > regrid_time);
      }

      if (ref == "PRESSURE_DEVIATION") {
         var = patch.getPatchData(d_pressure, getDataContext());
         size = d_pressure_dev_tol.getSize();
         tol = ( ( error_level_number < size) 
                 ? d_pressure_dev_tol[error_level_number] 
                 : d_pressure_dev_tol[size-1] );
         size = d_pressure_dev.getSize();
         dev = ( ( error_level_number < size) 
                 ? d_pressure_dev[error_level_number] 
                 : d_pressure_dev[size-1] );        
         size = d_pressure_dev_time_min.getSize();
         double time_min = ( ( error_level_number < size) 
                 ? d_pressure_dev_time_min[error_level_number] 
                 : d_pressure_dev_time_min[size-1] );
         size = d_pressure_dev_time_max.getSize();
         double time_max = ( ( error_level_number < size) 
                 ? d_pressure_dev_time_max[error_level_number] 
                 : d_pressure_dev_time_max[size-1] );
         time_allowed = (time_min <= regrid_time) && (time_max > regrid_time);
      }

      if (ref == "PRESSURE_GRADIENT") {
         var = patch.getPatchData(d_pressure, getDataContext());
         size = d_pressure_grad_tol.getSize();
         tol = ( ( error_level_number < size) 
                 ? d_pressure_grad_tol[error_level_number] 
                 : d_pressure_grad_tol[size-1] );
         size = d_pressure_grad_time_min.getSize();
         double time_min = ( ( error_level_number < size) 
                 ? d_pressure_grad_time_min[error_level_number] 
                 : d_pressure_grad_time_min[size-1] );
         size = d_pressure_grad_time_max.getSize();
         double time_max = ( ( error_level_number < size) 
                 ? d_pressure_grad_time_max[error_level_number] 
                 : d_pressure_grad_time_max[size-1] );
         time_allowed = (time_min <= regrid_time) && (time_max > regrid_time);
      }

      if (ref == "PRESSURE_SHOCK") {
         var = patch.getPatchData(d_pressure, getDataContext());
         size = d_pressure_shock_tol.getSize();
         tol = ( ( error_level_number < size) 
                 ? d_pressure_shock_tol[error_level_number] 
                 : d_pressure_shock_tol[size-1] );
         size = d_pressure_shock_onset.getSize();
         onset = ( ( error_level_number < size) 
                 ? d_pressure_shock_onset[error_level_number] 
                 : d_pressure_shock_onset[size-1] );
         size = d_pressure_shock_time_min.getSize();
         double time_min = ( ( error_level_number < size) 
                 ? d_pressure_shock_time_min[error_level_number] 
                 : d_pressure_shock_time_min[size-1] );
         size = d_pressure_shock_time_max.getSize();
         double time_max = ( ( error_level_number < size) 
                 ? d_pressure_shock_time_max[error_level_number] 
                 : d_pressure_shock_time_max[size-1] );
         time_allowed = (time_min <= regrid_time) && (time_max > regrid_time);
      }

      if (time_allowed) {

#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!var.isNull());
#endif

         hier::IntVector<NDIM> vghost = var->getGhostCellWidth();
         hier::IntVector<NDIM> tagghost = tags->getGhostCellWidth();

         if (ref == "DENSITY_DEVIATION" || ref == "PRESSURE_DEVIATION") {

            /*
             * Check for tags that have already been set in a previous
             * step.  Do NOT consider values tagged with value
             * RICHARDSON_NEWLY_TAGGED since these were set most recently
             * by Richardson extrapolation.
             */ 
            for (pdat::CellIterator<NDIM> ic(pbox); ic; ic++) {
               double locden = tol;
               int tag_val = (*tags)(ic(),0);
               if (tag_val) {
                  if (tag_val != RICHARDSON_NEWLY_TAGGED) {
                     locden *= 0.75;
                  }
               }
               if (tbox::MathUtilities<double>::Abs((*var)(ic())-dev) > locden)
               {
                  (*temp_tags)(ic(),0) = TRUE; 
               }
            }
         }

         if (ref == "DENSITY_GRADIENT" || ref == "PRESSURE_GRADIENT") {
            detectgrad_(ifirst(0),ilast(0),
                        ifirst(1),ilast(1),
#if (NDIM>2)
                        ifirst(2),ilast(2),
#endif
                        vghost(0),tagghost(0),d_nghosts(0),
                        vghost(1),tagghost(1),d_nghosts(1),
#if (NDIM>2)
                        vghost(2),tagghost(2),d_nghosts(2),
#endif
                        dx,
                        tol,
                        TRUE, FALSE,
                        var->getPointer(),
                        tags->getPointer(),temp_tags->getPointer());
         }

         if (ref == "DENSITY_SHOCK" || ref == "PRESSURE_SHOCK") {  
            detectshock_(ifirst(0),ilast(0),
                         ifirst(1),ilast(1),
#if (NDIM>2)
                         ifirst(2),ilast(2),
#endif
                         vghost(0),tagghost(0),d_nghosts(0),
                         vghost(1),tagghost(1),d_nghosts(1),
#if (NDIM>2)
                         vghost(2),tagghost(2),d_nghosts(2),
#endif
                         dx,
                         tol,
                         onset,
                         TRUE, FALSE,
                         var->getPointer(),
                         tags->getPointer(),temp_tags->getPointer());
         }
      
      }  // if time_allowed 

   }  // loop over criteria

   /*
    * Adjust temp_tags from those tags set in Richardson extrapolation.
    */
   if (uses_richardson_extrapolation_too) {
      for (pdat::CellIterator<NDIM> ic(pbox); ic; ic++) {              
         if ((*tags)(ic(),0) == RICHARDSON_ALREADY_TAGGED || 
             (*tags)(ic(),0) == RICHARDSON_NEWLY_TAGGED) {
            (*temp_tags)(ic(),0) = TRUE; 
         }
      }
   }

   /*
    * Update tags
    */
   for (pdat::CellIterator<NDIM> ic(pbox); ic; ic++) {
      (*tags)(ic(),0) = (*temp_tags)(ic(),0);
   }

   t_taggradient->stop();
}

/*
*************************************************************************
*                                                                       *
* Tag cells for refinement using Richardson extrapolation.  Criteria    *
* defined in input.                                                     *
*                                                                       *
*************************************************************************
*/

void Euler::tagRichardsonExtrapolationCells(
   hier::Patch<NDIM>& patch, 
   const int error_level_number,
   const tbox::Pointer<hier::VariableContext> coarsened_fine, 
   const tbox::Pointer<hier::VariableContext> advanced_coarse,
   const double regrid_time,
   const double deltat,
   const int error_coarsen_ratio,
   const bool initial_error,
   const int tag_index,
   const bool uses_gradient_detector_too)
{
   (void) initial_error;

   const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > patch_geom = 
      patch.getPatchGeometry();
   const double* xdomainlo = d_grid_geometry->getXLower();
   const double* xdomainhi = d_grid_geometry->getXUpper();

   hier::Box<NDIM> pbox = patch.getBox();

   tbox::Pointer< pdat::CellData<NDIM,int> > tags     = patch.getPatchData(tag_index);

   /*
    * Possible tagging criteria includes 
    *    DENSITY_RICHARDSON, PRESSURE_RICHARDSON
    * The criteria is specified over a time interval.
    *
    * Loop over criteria provided and check to make sure we are in the
    * specified time interval.  If so, apply appropriate tagging for
    * the level.
    */
   for (int ncrit = 0; ncrit < d_refinement_criteria.getSize(); ncrit++) {

      string ref = d_refinement_criteria[ncrit];
      tbox::Pointer< pdat::CellData<NDIM, double > > coarsened_fine_var; 
      tbox::Pointer< pdat::CellData<NDIM, double > > advanced_coarse_var;          
      int size = 0;
      double tol = 0.;
      bool time_allowed = false;

      if (ref == "DENSITY_RICHARDSON") {
         coarsened_fine_var = patch.getPatchData(d_density, coarsened_fine);
         advanced_coarse_var = patch.getPatchData(d_density, advanced_coarse);
         size = d_density_rich_tol.getSize();
         tol = ( ( error_level_number < size) 
                 ? d_density_rich_tol[error_level_number] 
                 : d_density_rich_tol[size-1] );
         size = d_density_rich_time_min.getSize();
         double time_min = ( ( error_level_number < size) 
                 ? d_density_rich_time_min[error_level_number] 
                 : d_density_rich_time_min[size-1] );
         size = d_density_rich_time_max.getSize();
         double time_max = ( ( error_level_number < size) 
                 ? d_density_rich_time_max[error_level_number] 
                 : d_density_rich_time_max[size-1] );
         time_allowed = (time_min <= regrid_time) && (time_max > regrid_time);
      }

      if (ref == "PRESSURE_RICHARDSON") {
         coarsened_fine_var = patch.getPatchData(d_pressure, coarsened_fine);
         advanced_coarse_var = patch.getPatchData(d_pressure, advanced_coarse);
         size = d_pressure_rich_tol.getSize();
         tol = ( ( error_level_number < size) 
                 ? d_pressure_rich_tol[error_level_number] 
                 : d_pressure_rich_tol[size-1] );
         size = d_pressure_rich_time_min.getSize();
         double time_min = ( ( error_level_number < size) 
                 ? d_pressure_rich_time_min[error_level_number] 
                 : d_pressure_rich_time_min[size-1] );
         size = d_pressure_rich_time_max.getSize();
         double time_max = ( ( error_level_number < size) 
                 ? d_pressure_rich_time_max[error_level_number] 
                 : d_pressure_rich_time_max[size-1] );
         time_allowed = (time_min <= regrid_time) && (time_max > regrid_time);
      }

      if (time_allowed) {

#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(!coarsened_fine_var.isNull());
         TBOX_ASSERT(!advanced_coarse_var.isNull());
#endif

         if (ref == "DENSITY_RICHARDSON" || ref == "PRESSURE_RICHARDSON") {

            /*
             * We tag wherever the global error > specified tolerance.
             * The estimated global error is the
             * local truncation error * the approximate number of steps
             * used in the simulation.  Approximate the number of steps as:
             *
             *       steps = L / (s*deltat)
             * where
             *       L = length of problem domain
             *       s = wave speed
             *       delta t = timestep on current level
             *
             * Compute max wave speed from delta t.  This presumes that
             * deltat was computed as deltat = dx/s_max.  We have deltat
             * and dx, so back out s_max from this.
             */

             const double* dx  = patch_geom->getDx();
             double max_dx = 0.;
             double max_length = 0.;
             for (int idir = 0; idir < NDIM; idir++) {
                max_dx = tbox::MathUtilities<double>::Max(max_dx, dx[idir]);
                double length = xdomainhi[idir] - xdomainlo[idir];
                max_length = 
                   tbox::MathUtilities<double>::Max(max_length, length);
             }
             double max_wave_speed = max_dx / deltat;
             double steps = max_length / (max_wave_speed * deltat);


             /*
              * Tag cells where |w_c - w_f| * (r^n -1) * steps
              *
              * where
              *       w_c = soln on coarse level (pressure_crse)
              *       w_f = soln on fine level (pressure_fine)
              *       r   = error coarsen ratio
              *       n   = spatial order of scheme (1st or 2nd depending 
              *             on whether Godunov order is 1st or 2nd/4th)
              */
             int order = 1;
             if (d_godunov_order > 1) order = 2;
             double r = error_coarsen_ratio;
             double rnminus1 = pow(r,order) - 1;

             double diff = 0.;
             double error = 0.;

  
             for (pdat::CellIterator<NDIM> ic(pbox); ic; ic++) {

                /*
                 * Compute error norm
                 */
                diff = (*advanced_coarse_var)(ic(),0) - 
                       (*coarsened_fine_var)(ic(),0);
                error = 
                   tbox::MathUtilities<double>::Abs(diff) * rnminus1 * steps;

                /*
                 * Tag cell if error > prescribed threshold. Since we are
                 * operating on the actual tag values (not temporary ones)
                 * distinguish here tags that were previously set before
                 * coming into this routine and those that are set here.
                 *     RICHARDSON_ALREADY_TAGGED - tagged before coming
                 *                                 into this method.
                 *     RICHARDSON_NEWLY_TAGGED - newly tagged in this method
                 */
                if (error > tol) {
                   if ((*tags)(ic(),0)) {
                      (*tags)(ic(),0) = RICHARDSON_ALREADY_TAGGED;
                   } else {
                      (*tags)(ic(),0) = RICHARDSON_NEWLY_TAGGED;
                   }
                }
            }

         }

      } // time_allowed

   } // loop over refinement criteria

   /*
    * If we are NOT performing gradient detector (i.e. only
    * doing Richardson extrapolation) set tags marked in this method
    * to TRUE and all others false.  Otherwise, leave tags set to the 
    * RICHARDSON_ALREADY_TAGGED and RICHARDSON_NEWLY_TAGGED as we may 
    * use this information in the gradient detector.
    */
   if (!uses_gradient_detector_too) {
      for (pdat::CellIterator<NDIM> ic(pbox); ic; ic++) {
         if ( (*tags)(ic(),0) == RICHARDSON_ALREADY_TAGGED ||
             (*tags)(ic(),0) == RICHARDSON_NEWLY_TAGGED ) {
            (*tags)(ic(),0) = TRUE;
         } else {
            (*tags)(ic(),0) = FALSE;
         }
      }
   }

}

/*
*************************************************************************
*                                                                       *
* Register Vizamrai data writer to write data to plot files that may    *
* be postprocessed by the Vizamrai tool.                                *
*                                                                       *
*************************************************************************
*/

void Euler::registerVizamraiDataWriter(
   tbox::Pointer<appu::CartesianVizamraiDataWriter<NDIM> > viz_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!(viz_writer.isNull()));
#endif
   d_vizamrai_writer = viz_writer;
}

/*
*************************************************************************
*                                                                       *
* Register VisIt data writer to write data to plot files that may       *
* be postprocessed by the VisIt tool.                                   *
*                                                                       *
*************************************************************************
*/

#ifdef HAVE_HDF5
void Euler::registerVisItDataWriter(
   tbox::Pointer<appu::VisItDataWriter<NDIM> > viz_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!(viz_writer.isNull()));
#endif
   d_visit_writer = viz_writer;
}
#endif

/*
*************************************************************************
*                                                                       *
* Pack "total energy" and "momentum" (derived Vis plot quantities)      *
* for the patch into a double precision buffer.                         *
*                                                                       *
*************************************************************************
*/

bool Euler::packDerivedDataIntoDoubleBuffer(
   double* dbuffer,
   const hier::Patch<NDIM>& patch,
   const hier::Box<NDIM>& region,
   const string& variable_name,
   int depth_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((region * patch.getBox()) == region);
#endif

   bool data_on_patch = FALSE;
   
   tbox::Pointer< pdat::CellData<NDIM,double> > density  =
      patch.getPatchData(d_density, d_plot_context);
   tbox::Pointer< pdat::CellData<NDIM,double> > velocity =
      patch.getPatchData(d_velocity, d_plot_context);
   tbox::Pointer< pdat::CellData<NDIM,double> > pressure =
      patch.getPatchData(d_pressure, d_plot_context);

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!density .isNull());
   TBOX_ASSERT(!velocity.isNull());
   TBOX_ASSERT(!pressure.isNull());
   TBOX_ASSERT(density->getGhostBox() == patch.getBox());
   TBOX_ASSERT(velocity->getGhostBox() == patch.getBox());
   TBOX_ASSERT(pressure->getGhostBox() == patch.getBox());
#endif

   const hier::Box<NDIM>& data_box = density->getGhostBox();
   const int box_w0 = region.numberCells(0);
   const int dat_w0 = data_box.numberCells(0);
   const int box_w1 = region.numberCells(1);
#if (NDIM > 2)
   const int dat_w1 = data_box.numberCells(1);
   const int box_w2 = region.numberCells(2);
#endif

   if (variable_name == "Total Energy") {
      const double *const dens = density->getPointer(); 
      const double *const xvel = velocity->getPointer(0); 
      const double *const yvel = velocity->getPointer(1);
#if (NDIM > 2)
      const double *const zvel = velocity->getPointer(2);
#endif
      const double *const pres = pressure->getPointer(); 

      double valinv = 1.0/(d_gamma-1.0);
      int buf_b1 = 0;
      int dat_b2 = data_box.offset(region.lower());

#if (NDIM > 2)
      for (int i2 = 0; i2 < box_w2; i2++) {
#endif
         int dat_b1 = dat_b2;
         for (int i1 = 0; i1 < box_w1; i1++) {
            for (int i0 = 0; i0 < box_w0; i0++) {
               int dat_indx = dat_b1+i0;
               double v2norm = pow(xvel[dat_indx], 2.0)
                             + pow(yvel[dat_indx], 2.0)
#if (NDIM > 2)
                             + pow(zvel[dat_indx], 2.0)
#endif
               ;
               double rho = dens[dat_indx];
               double int_energy = 0.0;
               if (rho > 0.0) {
                  int_energy = valinv * pres[dat_indx] / dens[dat_indx];
               }
               dbuffer[buf_b1+i0] = 
                  dens[dat_indx] * (0.5 * v2norm + int_energy);
            } 
            dat_b1 += dat_w0;
            buf_b1 += box_w0;
         }
#if (NDIM > 2)
         dat_b2 += dat_w1 * dat_w0;
      }
#endif

      data_on_patch = TRUE;

   } else if (variable_name == "Momentum") {
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(depth_id < NDIM);
#endif

      const double *const dens = density->getPointer();
      const double *const vel = velocity->getPointer(depth_id);
      int buf_b1 = 0;
      int dat_b2 = data_box.offset(region.lower());

#if (NDIM > 2)
      for (int i2 = 0; i2 < box_w2; i2++) {
#endif
         int dat_b1 = dat_b2;
         for (int i1 = 0; i1 < box_w1; i1++) {
            for (int i0 = 0; i0 < box_w0; i0++) {
               int dat_indx = dat_b1+i0;
               dbuffer[buf_b1+i0] = dens[dat_indx] * vel[dat_indx];
            }
            dat_b1 += dat_w0;
            buf_b1 += box_w0;
         }
#if (NDIM > 2)
         dat_b2 += dat_w1 * dat_w0;
      }
#endif

      data_on_patch = TRUE;

   } else {
      TBOX_ERROR("Euler::packDerivedDataIntoDoubleBuffer()"
             << "\n    unknown variable_name " << variable_name << "\n");
   }

   return(data_on_patch);
   
}

/*
*************************************************************************
*                                                                       *
* Write 1d data intersection of patch and pencil box to file            *
* with given name for plotting with Matlab.                             *
*                                                                       *
*************************************************************************
*/

void Euler::writeData1dPencil(const tbox::Pointer<hier::Patch<NDIM> > patch,
                              const hier::Box<NDIM>& pencil_box,
                              const int idir,
                              ostream& file)
{

   const hier::Box<NDIM>& patchbox = patch->getBox();
   const hier::Box<NDIM> box = pencil_box * patchbox;

   if (!box.empty()) {

      tbox::Pointer< pdat::CellData<NDIM,double> > density  =
         patch->getPatchData(d_density, getDataContext());
      tbox::Pointer< pdat::CellData<NDIM,double> > velocity =
         patch->getPatchData(d_velocity, getDataContext());
      tbox::Pointer< pdat::CellData<NDIM,double> > pressure =
         patch->getPatchData(d_pressure, getDataContext());

#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(!density .isNull());
      TBOX_ASSERT(!velocity.isNull());
      TBOX_ASSERT(!pressure.isNull());
#endif

      const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
      const double* dx  = pgeom->getDx();
      const double* xlo = pgeom->getXLower();

      const double cell_center = xlo[idir] 
                                 + (double(box.lower(idir)
                                     - patchbox.lower(idir)) 
                                    + 0.5) * dx[idir];

      double valinv = 1.0/(d_gamma-1.0);

      int ccount = 0;
      for (pdat::CellIterator<NDIM> ic(box); ic; ic++) {
         file << cell_center + ccount*dx[idir] << " ";
         ccount++;

         double rho = (*density)(ic(),0);
         double vel = (*velocity)(ic(),idir);
         double p   = (*pressure)(ic(),0);

         double mom = rho * vel;
         double eint = 0.0;
         if (rho > 0.0) {
            eint = valinv*(p/rho); 
         }
         double etot = rho*(0.5*vel*vel + eint);

         /*
          * Write out conserved quantities.
          */ 
         file << rho << " ";
         file << mom << " ";
         file << etot << " ";

         /* 
          * Write out "primitive" quantities and internal energy.
          */
         file << p << " ";
         file << vel << " ";
         file << eint << " ";

         file << endl;
      }

   }

}

/*
*************************************************************************
*                                                                       *
* Write all class data members to specified output stream.              *
*                                                                       *
*************************************************************************
*/

void Euler::printClassData(ostream &os) const 
{
   int j,k;

   os << "\nEuler::printClassData..." << endl;
   os << "Euler: this = " << (Euler*)this << endl;
   os << "d_object_name = " << d_object_name << endl;
   os << "d_grid_geometry = "
      << (geom::CartesianGridGeometry<NDIM>*)d_grid_geometry << endl;

   os << "Parameters for physical problem ..." << endl;
   os << "   d_gamma = " << d_gamma << endl;

   os << "Numerical method description and ghost sizes..." << endl;
   os << "   d_riemann_solve = " << d_riemann_solve << endl;
   os << "   d_riemann_solve_int = " << d_riemann_solve_int << endl;
   os << "   d_godunov_order = " << d_godunov_order << endl;
   os << "   d_corner_transport = " << d_corner_transport << endl;
   os << "   d_nghosts = " << d_nghosts << endl;
   os << "   d_fluxghosts = " << d_fluxghosts << endl;

   os << "Problem description and initial data..." << endl;
   os << "   d_data_problem = " << d_data_problem << endl;
   os << "   d_data_problem_int = " << d_data_problem_int << endl;

   os << "       d_radius = " << d_radius << endl;
   os << "       d_center = " ;
      for (j=0;j<NDIM;j++) os << d_center[j]<<" ";
   os << endl;
   os << "       d_density_inside = " << d_density_inside << endl;
   os << "       d_velocity_inside = " ;
     for (j=0;j<NDIM;j++) os << d_velocity_inside[j]<<" ";
   os << endl;
   os << "       d_pressure_inside = " << d_pressure_inside << endl;
   os << "       d_density_outside = " << d_density_outside << endl;
   os << "       d_velocity_outside = " ;
     for (j=0;j<NDIM;j++) os << d_velocity_outside[j]<<" ";
   os << endl;
   os << "       d_pressure_outside = " << d_pressure_outside << endl;

   os << "       d_number_of_intervals = " << d_number_of_intervals << endl;
   os << "       d_front_position = ";
   for (k = 0; k < d_number_of_intervals-1; k++) {
      os << d_front_position[k] << "  "; 
   } 
   os << endl;
   os << "       d_interval_density = " << endl;
   for (k = 0; k < d_number_of_intervals; k++) {
      os << "            " << d_interval_density[k] << endl;
   } 
   os << "       d_interval_velocity = " << endl;
   for (k = 0; k < d_number_of_intervals; k++) {
      os << "            ";
      for (j = 0; j < NDIM; j++) {
         os << d_interval_velocity[k*NDIM+j] << "  "; 
      }
      os << endl;
   } 
   os << "       d_interval_pressure = " << endl;
   for (k = 0; k < d_number_of_intervals; k++) {
      os << "            " << d_interval_pressure[k] << endl;
   }

   os << "   Boundary condition data " << endl;

#if (NDIM == 2) 
   for (j = 0; j < d_master_bdry_edge_conds.getSize(); j++) {
      os << "\n       d_master_bdry_edge_conds[" << j << "] = "
         << d_master_bdry_edge_conds[j] << endl;
      os << "       d_scalar_bdry_edge_conds[" << j << "] = "
         << d_scalar_bdry_edge_conds[j] << endl;
      os << "       d_vector_bdry_edge_conds[" << j << "] = "
         << d_vector_bdry_edge_conds[j] << endl;
      if (d_master_bdry_edge_conds[j] == DIRICHLET_BC) {
         os << "         d_bdry_edge_density[" << j << "] = "
            << d_bdry_edge_density[j] << endl;
         os << "         d_bdry_edge_velocity[" << j << "] = "
            << d_bdry_edge_velocity[j*NDIM+0] << " , "
            << d_bdry_edge_velocity[j*NDIM+1] << endl;
         os << "         d_bdry_edge_pressure[" << j << "] = "
            << d_bdry_edge_pressure[j] << endl;
      }
   }
   os << endl;
   for (j = 0; j < d_master_bdry_node_conds.getSize(); j++) {
      os << "\n       d_master_bdry_node_conds[" << j << "] = "
         << d_master_bdry_node_conds[j] << endl;
      os << "       d_scalar_bdry_node_conds[" << j << "] = "
         << d_scalar_bdry_node_conds[j] << endl;
      os << "       d_vector_bdry_node_conds[" << j << "] = "
         << d_vector_bdry_node_conds[j] << endl;
      os << "       d_node_bdry_edge[" << j << "] = "
         << d_node_bdry_edge[j] << endl;
   }
#endif
#if (NDIM == 3)
   for (j = 0; j < d_master_bdry_face_conds.getSize(); j++) {
      os << "\n       d_master_bdry_face_conds[" << j << "] = "
         << d_master_bdry_face_conds[j] << endl;
      os << "       d_scalar_bdry_face_conds[" << j << "] = "
         << d_scalar_bdry_face_conds[j] << endl;
      os << "       d_vector_bdry_face_conds[" << j << "] = "
         << d_vector_bdry_face_conds[j] << endl;
      if (d_master_bdry_face_conds[j] == DIRICHLET_BC) {
         os << "         d_bdry_face_density[" << j << "] = "
            << d_bdry_face_density[j] << endl;
         os << "         d_bdry_face_velocity[" << j << "] = "
            << d_bdry_face_velocity[j*NDIM+0] << " , "
            << d_bdry_face_velocity[j*NDIM+1] << " , "
            << d_bdry_face_velocity[j*NDIM+2] << endl;
         os << "         d_bdry_face_pressure[" << j << "] = "
            << d_bdry_face_pressure[j] << endl;
      }
   }
   os << endl;
   for (j = 0; j < d_master_bdry_edge_conds.getSize(); j++) {
      os << "\n       d_master_bdry_edge_conds[" << j << "] = "
         << d_master_bdry_edge_conds[j] << endl;
      os << "       d_scalar_bdry_edge_conds[" << j << "] = "
         << d_scalar_bdry_edge_conds[j] << endl;
      os << "       d_vector_bdry_edge_conds[" << j << "] = "
         << d_vector_bdry_edge_conds[j] << endl;
      os << "       d_edge_bdry_face[" << j << "] = "
         << d_edge_bdry_face[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_master_bdry_node_conds.getSize(); j++) {
      os << "\n       d_master_bdry_node_conds[" << j << "] = "
         << d_master_bdry_node_conds[j] << endl;
      os << "       d_scalar_bdry_node_conds[" << j << "] = "
         << d_scalar_bdry_node_conds[j] << endl;
      os << "       d_vector_bdry_node_conds[" << j << "] = "
         << d_vector_bdry_node_conds[j] << endl;
      os << "       d_node_bdry_face[" << j << "] = "
         << d_node_bdry_face[j] << endl;
   }
#endif

   os << "   Refinement criteria parameters " << endl;
 
   for (j = 0; j < d_refinement_criteria.getSize(); j++) {
      os << "       d_refinement_criteria[" << j << "] = "
         << d_refinement_criteria[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_density_dev_tol.getSize(); j++) {
      os << "       d_density_dev_tol[" << j << "] = "
         << d_density_dev_tol[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_density_dev.getSize(); j++) {
      os << "       d_density_dev[" << j << "] = "
         << d_density_dev[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_density_dev_time_max.getSize(); j++) {
      os << "       d_density_dev_time_max[" << j << "] = "
         << d_density_dev_time_max[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_density_dev_time_min.getSize(); j++) {
      os << "       d_density_dev_time_min[" << j << "] = "
         << d_density_dev_time_min[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_density_grad_tol.getSize(); j++) {
      os << "       d_density_grad_tol[" << j << "] = "
         << d_density_grad_tol[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_density_grad_time_max.getSize(); j++) {
      os << "       d_density_grad_time_max[" << j << "] = "
         << d_density_grad_time_max[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_density_grad_time_min.getSize(); j++) {
      os << "       d_density_grad_time_min[" << j << "] = "
         << d_density_grad_time_min[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_density_shock_onset.getSize(); j++) {
      os << "       d_density_shock_onset[" << j << "] = "
         << d_density_shock_onset[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_density_shock_tol.getSize(); j++) {
      os << "       d_density_shock_tol[" << j << "] = "
         << d_density_shock_tol[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_density_shock_time_max.getSize(); j++) {
      os << "       d_density_shock_time_max[" << j << "] = "
         << d_density_shock_time_max[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_density_shock_time_min.getSize(); j++) {
      os << "       d_density_shock_time_min[" << j << "] = "
         << d_density_shock_time_min[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_density_rich_tol.getSize(); j++) {
      os << "       d_density_rich_tol[" << j << "] = "
         << d_density_rich_tol[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_density_rich_time_max.getSize(); j++) {
      os << "       d_density_rich_time_max[" << j << "] = "
         << d_density_rich_time_max[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_density_rich_time_min.getSize(); j++) {
      os << "       d_density_rich_time_min[" << j << "] = "
         << d_density_rich_time_min[j] << endl;
   }
   os << endl;

   for (j = 0; j < d_pressure_dev_tol.getSize(); j++) {
      os << "       d_pressure_dev_tol[" << j << "] = "
         << d_pressure_dev_tol[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_pressure_dev.getSize(); j++) {
      os << "       d_pressure_dev[" << j << "] = "
         << d_pressure_dev[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_pressure_dev_time_max.getSize(); j++) {
      os << "       d_pressure_dev_time_max[" << j << "] = "
         << d_pressure_dev_time_max[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_pressure_dev_time_min.getSize(); j++) {
      os << "       d_pressure_dev_time_min[" << j << "] = "
         << d_pressure_dev_time_min[j] << endl;
   }
   os << endl; 
   for (j = 0; j < d_pressure_grad_tol.getSize(); j++) {
      os << "       d_pressure_grad_tol[" << j << "] = "
         << d_pressure_grad_tol[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_pressure_grad_time_max.getSize(); j++) {
      os << "       d_pressure_grad_time_max[" << j << "] = "
         << d_pressure_grad_time_max[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_pressure_grad_time_min.getSize(); j++) {
      os << "       d_pressure_grad_time_min[" << j << "] = "
         << d_pressure_grad_time_min[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_pressure_shock_onset.getSize(); j++) {
      os << "       d_pressure_shock_onset[" << j << "] = "
         << d_pressure_shock_onset[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_pressure_shock_tol.getSize(); j++) {
      os << "       d_pressure_shock_tol[" << j << "] = "
         << d_pressure_shock_tol[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_pressure_shock_time_max.getSize(); j++) {
      os << "       d_pressure_shock_time_max[" << j << "] = "
         << d_pressure_shock_time_max[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_pressure_shock_time_min.getSize(); j++) {
      os << "       d_pressure_shock_time_min[" << j << "] = "
         << d_pressure_shock_time_min[j] << endl;
   }
   os << endl;
  for (j = 0; j < d_pressure_rich_tol.getSize(); j++) {
      os << "       d_pressure_rich_tol[" << j << "] = "
         << d_pressure_rich_tol[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_pressure_rich_time_max.getSize(); j++) {
      os << "       d_pressure_rich_time_max[" << j << "] = "
         << d_pressure_rich_time_max[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_pressure_rich_time_min.getSize(); j++) {
      os << "       d_pressure_rich_time_min[" << j << "] = "
         << d_pressure_rich_time_min[j] << endl;
   }
   os << endl;

}

/*
*************************************************************************
*                                                                       *
* Read data members from input.  Note all values set from restart       *
* can be overridden by values in the input database.                    *
*                                                                       *
*************************************************************************
*/

void Euler::getFromInput(
   tbox::Pointer<tbox::Database> db,
   bool is_from_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!db.isNull());
#endif

   /*
    * Note: if we are restarting, then we only allow nonuniform
    * workload to be used if nonuniform workload was used originally.
    */
   if (!is_from_restart) { 
      d_use_nonuniform_workload = 
         db->getBoolWithDefault("use_nonuniform_workload",
                                d_use_nonuniform_workload);
   } else {
      if (d_use_nonuniform_workload) {
         d_use_nonuniform_workload = 
            db->getBool("use_nonuniform_workload");
      }
   }

   if (!is_from_restart) {
      d_gamma = db->getDoubleWithDefault("gamma", d_gamma);
   }

   if (db->keyExists("riemann_solve")) {
      d_riemann_solve = db->getString("riemann_solve");
      if ( (d_riemann_solve != "APPROX_RIEM_SOLVE") &&
           (d_riemann_solve != "EXACT_RIEM_SOLVE") &&
           (d_riemann_solve != "HLLC_RIEM_SOLVE") ) {
         TBOX_ERROR(d_object_name << ": "
                    << "`riemann_solve' in input must be either string "
                    << "'APPROX_RIEM_SOLVE', 'EXACT_RIEM_SOLVE', "
                    << "'HLLC_RIEM_SOLVE'." << endl);

      }
   } else {
      d_riemann_solve = db->getStringWithDefault("d_riemann_solve",
                                                  d_riemann_solve);
   }

   if (db->keyExists("godunov_order")) {
      d_godunov_order = db->getInteger("godunov_order");
      if ( (d_godunov_order != 1) &&
           (d_godunov_order != 2) &&
           (d_godunov_order != 4) ) {
         TBOX_ERROR(d_object_name << ": "
            << "`godunov_order' in input must be 1, 2, or 4." << endl);

      }
   } else {
      d_godunov_order = db->getIntegerWithDefault("d_godunov_order",
                                                   d_godunov_order);
   }

   if (db->keyExists("corner_transport")) {
      d_corner_transport = db->getString("corner_transport");
      if ( (d_corner_transport != "CORNER_TRANSPORT_1") &&
           (d_corner_transport != "CORNER_TRANSPORT_2") ) {
         TBOX_ERROR(d_object_name << ": "
            << "`corner_transport' in input must be either string"
            << " 'CORNER_TRANSPORT_1' or 'CORNER_TRANSPORT_2'." << endl);
      }
   } else {
      d_corner_transport = db->getStringWithDefault("corner_transport",
                                                     d_corner_transport);
   }


   if (db->keyExists("Refinement_data")) {
      tbox::Pointer<tbox::Database> refine_db = db->getDatabase("Refinement_data");
      tbox::Array<string> refinement_keys = refine_db->getAllKeys();
      int num_keys = refinement_keys.getSize();

      if (refine_db->keyExists("refine_criteria")) {
         d_refinement_criteria = 
            refine_db->getStringArray("refine_criteria");
      } else {
         TBOX_WARNING(d_object_name << ": "
                   << "No key `refine_criteria' found in data for"
                   << " RefinementData. No refinement will occur." << endl);
      }
       
      tbox::Array<string> ref_keys_defined(num_keys);
      int def_key_cnt = 0;
      tbox::Pointer<tbox::Database> error_db;
      for (int i = 0; i < refinement_keys.getSize(); i++) {
           
         string error_key = refinement_keys[i];
         error_db.setNull();

         if ( !(error_key == "refine_criteria") ) {

            if ( !(error_key == "DENSITY_DEVIATION" ||
                   error_key == "DENSITY_GRADIENT" ||
                   error_key == "DENSITY_SHOCK" ||
                   error_key == "DENSITY_RICHARDSON" ||
                   error_key == "PRESSURE_DEVIATION" ||
                   error_key == "PRESSURE_GRADIENT" ||
                   error_key == "PRESSURE_SHOCK" ||
                   error_key == "PRESSURE_RICHARDSON") ) { 
               TBOX_ERROR(d_object_name << ": "
                         << "Unknown refinement criteria: " << error_key
                         << "\nin input." << endl);
            } else {
               error_db = refine_db->getDatabase(error_key);
               ref_keys_defined[def_key_cnt] = error_key;
               def_key_cnt++;
            }
               
            if (!error_db.isNull() && error_key == "DENSITY_DEVIATION") {

               if (error_db->keyExists("dev_tol")) {
                  d_density_dev_tol = 
                  error_db->getDoubleArray("dev_tol");
               } else {
                  TBOX_ERROR(d_object_name << ": "
                           << "No key `dev_tol' found in data for "
                           << error_key << endl);
               }

               if (error_db->keyExists("density_dev")) {
                  d_density_dev = 
                  error_db->getDoubleArray("density_dev");
               } else {
                  TBOX_ERROR(d_object_name << ": "
                           << "No key `density_dev' found in data for "
                           << error_key << endl);
               }

               if (error_db->keyExists("time_max")){
                  d_density_dev_time_max = 
                  error_db->getDoubleArray("time_max");
               } else {
                  d_density_dev_time_max.resizeArray(1);
                  d_density_dev_time_max[0] = 
                     tbox::MathUtilities<double>::getMax();
               }

               if (error_db->keyExists("time_min")){
                  d_density_dev_time_min = 
                  error_db->getDoubleArray("time_min");
               } else {
                  d_density_dev_time_min.resizeArray(1);
                  d_density_dev_time_min[0] = 0.;
               } 

            }
              
            if (!error_db.isNull() && error_key == "DENSITY_GRADIENT") {

               if (error_db->keyExists("grad_tol")) {
                  d_density_grad_tol = 
                  error_db->getDoubleArray("grad_tol");
               } else {
                  TBOX_ERROR(d_object_name << ": "
                           << "No key `grad_tol' found in data for "
                           << error_key << endl);
               }

               if (error_db->keyExists("time_max")){
                  d_density_grad_time_max = 
                  error_db->getDoubleArray("time_max");
               } else {
                  d_density_grad_time_max.resizeArray(1);
                  d_density_grad_time_max[0] = 
                     tbox::MathUtilities<double>::getMax();
               }

               if (error_db->keyExists("time_min")){
                  d_density_grad_time_min = 
                  error_db->getDoubleArray("time_min");
               } else {
                  d_density_grad_time_min.resizeArray(1);
                  d_density_grad_time_min[0] = 0.;
               } 

            }
              
            if (!error_db.isNull() && error_key == "DENSITY_SHOCK") {

               if (error_db->keyExists("shock_onset")) {
                  d_density_shock_onset = 
                  error_db->getDoubleArray("shock_onset");
               } else {
                  TBOX_ERROR(d_object_name << ": "
                           << "No key `shock_onset' found in data for "
                           << error_key << endl);
               }

               if (error_db->keyExists("shock_tol")) {
                  d_density_shock_tol = 
                  error_db->getDoubleArray("shock_tol");
               } else {
                  TBOX_ERROR(d_object_name << ": "
                           << "No key `shock_tol' found in data for "
                           << error_key << endl);
               }

               if (error_db->keyExists("time_max")){
                  d_density_shock_time_max = 
                  error_db->getDoubleArray("time_max");
               } else {
                  d_density_shock_time_max.resizeArray(1);
                  d_density_shock_time_max[0] = 
                     tbox::MathUtilities<double>::getMax();
               }

               if (error_db->keyExists("time_min")){
                  d_density_shock_time_min = 
                  error_db->getDoubleArray("time_min");
               } else {
                  d_density_shock_time_min.resizeArray(1);
                  d_density_shock_time_min[0] = 0.;
               }

            }

            if (!error_db.isNull() && error_key == "DENSITY_RICHARDSON") {

               if (error_db->keyExists("rich_tol")) {
                  d_density_rich_tol = 
                  error_db->getDoubleArray("rich_tol");
               } else {
                  TBOX_ERROR(d_object_name << ": "
                           << "No key `rich_tol' found in data for "
                           << error_key << endl);
               }

               if (error_db->keyExists("time_max")){
                  d_density_rich_time_max = 
                  error_db->getDoubleArray("time_max");
               } else {
                  d_density_rich_time_max.resizeArray(1);
                  d_density_rich_time_max[0] = 
                     tbox::MathUtilities<double>::getMax();
               }

               if (error_db->keyExists("time_min")){
                  d_density_rich_time_min = 
                  error_db->getDoubleArray("time_min");
               } else {
                  d_density_rich_time_min.resizeArray(1);
                  d_density_rich_time_min[0] = 0.;
               } 

            }

            if (!error_db.isNull() && error_key == "PRESSURE_DEVIATION") {

               if (error_db->keyExists("dev_tol")) {
                  d_pressure_dev_tol = 
                  error_db->getDoubleArray("dev_tol");
               } else {
                  TBOX_ERROR(d_object_name << ": "
                           << "No key `dev_tol' found in data for "
                           << error_key << endl);
               }

               if (error_db->keyExists("pressure_dev")) {
                  d_pressure_dev = 
                  error_db->getDoubleArray("pressure_dev");
               } else {
                  TBOX_ERROR(d_object_name << ": "
                           << "No key `pressure_dev' found in data for "
                           << error_key << endl);
               }

               if (error_db->keyExists("time_max")){
                  d_pressure_dev_time_max = 
                  error_db->getDoubleArray("time_max");
               } else {
                  d_pressure_dev_time_max.resizeArray(1);
                  d_pressure_dev_time_max[0] = 
                     tbox::MathUtilities<double>::getMax();
               }

               if (error_db->keyExists("time_min")){
                  d_pressure_dev_time_min = 
                  error_db->getDoubleArray("time_min");
               } else {
                  d_pressure_dev_time_min.resizeArray(1);
                  d_pressure_dev_time_min[0] = 0.;
               } 

            }

            if (!error_db.isNull() && error_key == "PRESSURE_GRADIENT") {

               if (error_db->keyExists("grad_tol")) {
                  d_pressure_grad_tol = 
                  error_db->getDoubleArray("grad_tol");
               } else {
                  TBOX_ERROR(d_object_name << ": "
                           << "No key `grad_tol' found in data for "
                           << error_key << endl);
               }

               if (error_db->keyExists("time_max")){
                  d_pressure_grad_time_max = 
                  error_db->getDoubleArray("time_max");
               } else {
                  d_pressure_grad_time_max.resizeArray(1);
                  d_pressure_grad_time_max[0] = 
                     tbox::MathUtilities<double>::getMax();
               }

               if (error_db->keyExists("time_min")){
                  d_pressure_grad_time_min = 
                  error_db->getDoubleArray("time_min");
               } else {
                  d_pressure_grad_time_min.resizeArray(1);
                  d_pressure_grad_time_min[0] = 0.;
               }

            }
      
            if (!error_db.isNull() && error_key == "PRESSURE_SHOCK") {

               if (error_db->keyExists("shock_onset")) {
                  d_pressure_shock_onset = 
                  error_db->getDoubleArray("shock_onset");
               } else {
                  TBOX_ERROR(d_object_name << ": "
                           << "No key `shock_onset' found in data for "
                           << error_key << endl);
               }

               if (error_db->keyExists("shock_tol")) {
                  d_pressure_shock_tol = 
                  error_db->getDoubleArray("shock_tol");
               } else {
                  TBOX_ERROR(d_object_name << ": "
                           << "No key `shock_tol' found in data for "
                           << error_key << endl);
               }

               if (error_db->keyExists("time_max")){
                  d_pressure_shock_time_max = 
                  error_db->getDoubleArray("time_max");
               } else {
                  d_pressure_shock_time_max.resizeArray(1);
                  d_pressure_shock_time_max[0] = 
                     tbox::MathUtilities<double>::getMax();
               }

               if (error_db->keyExists("time_min")){
                  d_pressure_shock_time_min = 
                  error_db->getDoubleArray("time_min");
               } else {
                  d_pressure_shock_time_min.resizeArray(1);
                  d_pressure_shock_time_min[0] = 0.;
               }

            }

            if (!error_db.isNull() && error_key == "PRESSURE_RICHARDSON") {

               if (error_db->keyExists("rich_tol")) {
                  d_pressure_rich_tol = 
                  error_db->getDoubleArray("rich_tol");
               } else {
                  TBOX_ERROR(d_object_name << ": "
                          << "No key `rich_tol' found in data for "
                          << error_key << endl);
               }

               if (error_db->keyExists("time_max")){
                  d_pressure_rich_time_max = 
                  error_db->getDoubleArray("time_max");
               } else {
                  d_pressure_rich_time_max.resizeArray(1);
                  d_pressure_rich_time_max[0] = 
                     tbox::MathUtilities<double>::getMax();
               }

               if (error_db->keyExists("time_min")){
                  d_pressure_rich_time_min = 
                  error_db->getDoubleArray("time_min");
               } else {
                  d_pressure_rich_time_min.resizeArray(1);
                  d_pressure_rich_time_min[0] = 0.;
               } 

            }

         } 

      } // loop over refine criteria     


      /* 
       * Check that input is found for each string identifier in key list.
       */
      for (int k0 = 0; k0 < d_refinement_criteria.getSize(); k0++) {
         string use_key = d_refinement_criteria[k0];
         bool key_found = false;
         for (int k1 = 0; k1 < def_key_cnt; k1++) {
             string def_key = ref_keys_defined[k1];
             if (def_key == use_key) key_found = true;
         }

         if (!key_found) {
             TBOX_ERROR(d_object_name << ": "
                       << "No input found for specified refine criteria: "
                       << d_refinement_criteria[k0] << endl);
         }
      }

   } // if "Refinement_data" db entry exists

   if (!is_from_restart) { 

      if (db->keyExists("data_problem")) {
         d_data_problem = db->getString("data_problem");
      } else {
         TBOX_ERROR(d_object_name << ": "
            << "`data_problem' value not found in input." << endl);
      }

      tbox::Pointer<tbox::Database> init_data_db;
      if (db->keyExists("Initial_data")) {
         init_data_db = db->getDatabase("Initial_data");
      } else {
         TBOX_ERROR(d_object_name << ": "
                  << "No `Initial_data' database found in input." << endl);
      }

      bool found_problem_data = false;

      if (d_data_problem == "SPHERE") {

         if (init_data_db->keyExists("radius")) {
            d_radius = init_data_db->getDouble("radius");
         } else {
            TBOX_ERROR(d_object_name << ": "
               << "`radius' input required for SPHERE problem." << endl);
         }
         if (init_data_db->keyExists("center")) {
            init_data_db->getDoubleArray("center", d_center, NDIM);
         } else {
            TBOX_ERROR(d_object_name << ": "
               << "`center' input required for SPHERE problem." << endl);
         }
         if (init_data_db->keyExists("density_inside")) {
            d_density_inside = init_data_db->getDouble("density_inside");
         } else {
            TBOX_ERROR(d_object_name << ": "
               << "`density_inside' input required for "
               << "SPHERE problem." << endl);
         }
         if (init_data_db->keyExists("velocity_inside")) {
            init_data_db->getDoubleArray("velocity_inside",
                                         d_velocity_inside, NDIM);
         } else {
            TBOX_ERROR(d_object_name << ": "
               << "`velocity_inside' input required for "
               << "SPHERE problem." << endl);
         }
         if (init_data_db->keyExists("pressure_inside")) {
            d_pressure_inside = init_data_db->getDouble("pressure_inside");
         } else {
            TBOX_ERROR(d_object_name << ": "
               << "`pressure_inside' input required for "
               << "SPHERE problem." << endl);
         }
         if (init_data_db->keyExists("density_outside")) {
            d_density_outside = init_data_db->getDouble("density_outside");
         } else {
            TBOX_ERROR(d_object_name << ": "
               << "`density_outside' input required for "
               << "SPHERE problem." << endl);
         }
         if (init_data_db->keyExists("velocity_outside")) {
            init_data_db->getDoubleArray("velocity_outside",
                                         d_velocity_outside, NDIM);
         } else {
            TBOX_ERROR(d_object_name << ": "
               << "`velocity_outside' input required for "
               << "SPHERE problem." << endl);
         }
         if (init_data_db->keyExists("pressure_outside")) {
            d_pressure_outside = init_data_db->getDouble("pressure_outside");
         } else {
            TBOX_ERROR(d_object_name << ": "
               << "`pressure_outside' input required for "
               << "SPHERE problem." << endl);
         }
   
         found_problem_data = true;

      }

      if (!found_problem_data &&
          (d_data_problem == "PIECEWISE_CONSTANT_X") ||
          (d_data_problem == "PIECEWISE_CONSTANT_Y") ||
          (d_data_problem == "PIECEWISE_CONSTANT_Z") ||
          (d_data_problem == "STEP")) {

          int idir = 0;
          if (d_data_problem == "PIECEWISE_CONSTANT_Y") {
             idir = 1;
          }

          if (d_data_problem == "PIECEWISE_CONSTANT_Z") {
             if (NDIM < 3) {
                TBOX_ERROR(d_object_name << ": `PIECEWISE_CONSTANT_Z' "
                           << "problem invalid in 2 dimensions." << endl);
             }
             idir = 2;
          }

          tbox::Array<string> init_data_keys = init_data_db->getAllKeys();

          if (init_data_db->keyExists("front_position")) {
             d_front_position = init_data_db->getDoubleArray("front_position");
          } else {
             TBOX_ERROR(d_object_name << ": "
                << "`front_position' input required for "
                << "PIECEWISE_CONSTANT_* problem." << endl);
          }

          d_number_of_intervals = 
             tbox::MathUtilities<int>::Min( d_front_position.getSize()+1,
                                            init_data_keys.getSize()-1 );

          d_front_position.resizeArray(d_front_position.getSize()+1);
          d_front_position[d_front_position.getSize()-1] = 
             d_grid_geometry->getXUpper()[idir]; 

          d_interval_density.resizeArray(d_number_of_intervals);
          d_interval_velocity.resizeArray(d_number_of_intervals*NDIM);
          d_interval_pressure.resizeArray(d_number_of_intervals);

          int i = 0;
          int nkey = 0;
          bool found_interval_data = false;

          while (   !found_interval_data 
                 && (i < d_number_of_intervals)
                 && (nkey < init_data_keys.getSize()) ) {

             if ( !(init_data_keys[nkey] == "front_position") ) {

                tbox::Pointer<tbox::Database> interval_db =
                      init_data_db->getDatabase(init_data_keys[nkey]);

                   readStateDataEntry(interval_db,
                                      init_data_keys[nkey],
                                      i,
                                      d_interval_density,
                                      d_interval_velocity,
                                      d_interval_pressure);

                i++;
   
                found_interval_data = (i == d_number_of_intervals);
   
             }
   
             nkey++;
   
          }

          if (!found_interval_data) {
             TBOX_ERROR(d_object_name << ": "
                        << "Insufficient interval data given in input"
                        << " for PIECEWISE_CONSTANT_* or STEP problem." << endl);
          }
      
          found_problem_data = true;

      }

      if (!found_problem_data) {
         TBOX_ERROR(d_object_name << ": "
            << "`Initial_data' database found in input." 
            << " But bad data supplied." << endl);
      } 
 
   } // if !is_from_restart read in problem data

   hier::IntVector<NDIM> periodic = d_grid_geometry->getPeriodicShift();
   int num_per_dirs = 0;
   for (int id = 0; id < NDIM; id++) {
      if (periodic(id)) num_per_dirs++;
   }

   if (num_per_dirs < NDIM) {

      if (db->keyExists("Boundary_data")) {

         tbox::Pointer<tbox::Database> bdry_db = db->getDatabase("Boundary_data");

#if (NDIM == 2)
         appu::CartesianBoundaryUtilities2::readBoundaryInput(this,
                                                        bdry_db,
                                                        d_master_bdry_edge_conds,
                                                        d_master_bdry_node_conds,
                                                        periodic);
#endif
#if (NDIM == 3)
         appu::CartesianBoundaryUtilities3::readBoundaryInput(this,
                                                        bdry_db,
                                                        d_master_bdry_face_conds,
                                                        d_master_bdry_edge_conds,
                                                        d_master_bdry_node_conds,
                                                        periodic);
#endif

      } else {
         TBOX_ERROR(d_object_name << ": "
            << "Key data `Boundary_data' not found in input. " << endl);
      }

   }

}

/*
*************************************************************************
*                                                                       *
* Routines to put/get data members to/from from restart database.       *
*                                                                       *
*************************************************************************
*/

void Euler::putToDatabase(tbox::Pointer<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!db.isNull());
#endif

   db->putInteger("EULER_VERSION", EULER_VERSION);

   db->putDouble("d_gamma", d_gamma);

   db->putString("d_riemann_solve", d_riemann_solve);
   db->putInteger("d_godunov_order", d_godunov_order);
   db->putString("d_corner_transport", d_corner_transport);
   db->putIntegerArray("d_nghosts", d_nghosts, NDIM);
   db->putIntegerArray("d_fluxghosts", d_fluxghosts, NDIM);

   db->putString("d_data_problem", d_data_problem);

   if (d_data_problem == "SPHERE") {
      db->putDouble("d_radius", d_radius);
      db->putDoubleArray("d_center", d_center, NDIM);
      db->putDouble("d_density_inside", d_density_inside);
      db->putDoubleArray("d_velocity_inside", d_velocity_inside, NDIM);
      db->putDouble("d_pressure_inside", d_pressure_inside);
      db->putDouble("d_density_outside", d_density_outside);
      db->putDoubleArray("d_velocity_outside", d_velocity_outside, NDIM);
      db->putDouble("d_pressure_outside", d_pressure_outside);
   }

   if ( (d_data_problem == "PIECEWISE_CONSTANT_X") ||
        (d_data_problem == "PIECEWISE_CONSTANT_Y") ||
        (d_data_problem == "PIECEWISE_CONSTANT_Z") ||
        (d_data_problem == "STEP") ) {
      db->putInteger("d_number_of_intervals", d_number_of_intervals);
      if (d_number_of_intervals > 0) {
         db->putDoubleArray("d_front_position", d_front_position);   
         db->putDoubleArray("d_interval_density", d_interval_density);   
         db->putDoubleArray("d_interval_velocity", d_interval_velocity);   
         db->putDoubleArray("d_interval_pressure", d_interval_pressure);   
      }
   }

   db->putIntegerArray("d_master_bdry_edge_conds", d_master_bdry_edge_conds);
   db->putIntegerArray("d_master_bdry_node_conds", d_master_bdry_node_conds);

#if (NDIM == 2) 
   db->putDoubleArray("d_bdry_edge_density", d_bdry_edge_density);
   db->putDoubleArray("d_bdry_edge_velocity", d_bdry_edge_velocity);
   db->putDoubleArray("d_bdry_edge_pressure", d_bdry_edge_pressure);
#endif
#if (NDIM == 3)
   db->putIntegerArray("d_master_bdry_face_conds", d_master_bdry_face_conds);

   db->putDoubleArray("d_bdry_face_density", d_bdry_face_density);
   db->putDoubleArray("d_bdry_face_velocity", d_bdry_face_velocity);
   db->putDoubleArray("d_bdry_face_pressure", d_bdry_face_pressure);
#endif

   if (d_refinement_criteria.getSize() > 0) {
      db->putStringArray("d_refinement_criteria", d_refinement_criteria);
   }
   for (int i = 0; i < d_refinement_criteria.getSize(); i++) {

      if (d_refinement_criteria[i] == "DENSITY_DEVIATION") {

         db->putDoubleArray("d_density_dev_tol", 
                            d_density_dev_tol);
         db->putDoubleArray("d_density_dev", 
                            d_density_dev);
         db->putDoubleArray("d_density_dev_time_max", 
                            d_density_dev_time_max);
         db->putDoubleArray("d_density_dev_time_min", 
                            d_density_dev_time_min);

      } else if (d_refinement_criteria[i] == "DENSITY_GRADIENT") {

         db->putDoubleArray("d_density_grad_tol", 
                            d_density_grad_tol);
         db->putDoubleArray("d_density_grad_time_max", 
                            d_density_grad_time_max);
         db->putDoubleArray("d_density_grad_time_min", 
                            d_density_grad_time_min);

      } else if  (d_refinement_criteria[i] == "DENSITY_SHOCK") {

         db->putDoubleArray("d_density_shock_onset", 
                            d_density_shock_onset);
         db->putDoubleArray("d_density_shock_tol", 
                            d_density_shock_tol);
         db->putDoubleArray("d_density_shock_time_max", 
                            d_density_shock_time_max);
         db->putDoubleArray("d_density_shock_time_min", 
                            d_density_shock_time_min);

      } else if  (d_refinement_criteria[i] == "DENSITY_RICHARDSON") {

         db->putDoubleArray("d_density_rich_tol", 
                            d_density_rich_tol);
         db->putDoubleArray("d_density_rich_time_max", 
                            d_density_rich_time_max);
         db->putDoubleArray("d_density_rich_time_min", 
                            d_density_rich_time_min);
 
      } else if (d_refinement_criteria[i] == "PRESSURE_DEVIATION") {

         db->putDoubleArray("d_pressure_dev_tol", 
                            d_pressure_dev_tol);
         db->putDoubleArray("d_pressure_dev", 
                            d_pressure_dev);
         db->putDoubleArray("d_pressure_dev_time_max", 
                            d_pressure_dev_time_max);
         db->putDoubleArray("d_pressure_dev_time_min", 
                            d_pressure_dev_time_min);

      } else if  (d_refinement_criteria[i] == "PRESSURE_GRADIENT") {

         db->putDoubleArray("d_pressure_grad_tol", 
                            d_pressure_grad_tol);
         db->putDoubleArray("d_pressure_grad_time_max", 
                            d_pressure_grad_time_max);
         db->putDoubleArray("d_pressure_grad_time_min", 
                            d_pressure_grad_time_min);

      } else if  (d_refinement_criteria[i] == "PRESSURE_SHOCK") {

         db->putDoubleArray("d_pressure_shock_onset", 
                            d_pressure_shock_onset);
         db->putDoubleArray("d_pressure_shock_tol", 
                            d_pressure_shock_tol);
         db->putDoubleArray("d_pressure_shock_time_max", 
                            d_pressure_shock_time_max);
         db->putDoubleArray("d_pressure_shock_time_min", 
                            d_pressure_shock_time_min);
 
      } else if  (d_refinement_criteria[i] == "PRESSURE_RICHARDSON") {

         db->putDoubleArray("d_pressure_rich_tol", 
                            d_pressure_rich_tol);
         db->putDoubleArray("d_pressure_rich_time_max", 
                            d_pressure_rich_time_max);
         db->putDoubleArray("d_pressure_rich_time_min", 
                            d_pressure_rich_time_min);
      }

   }

}

void Euler::getFromRestart()
{

   tbox::Pointer<tbox::Database> root_db = 
      tbox::RestartManager::getManager()->getRootDatabase();

   tbox::Pointer<tbox::Database> db;
   if ( root_db->isDatabase(d_object_name) ) {
      db = root_db->getDatabase(d_object_name);
   } else {
      TBOX_ERROR("Restart database corresponding to "
                 << d_object_name << " not found in restart file." << endl);
   }

   int ver = db->getInteger("EULER_VERSION");
   if (ver != EULER_VERSION) {
      TBOX_ERROR(d_object_name << ": "
          << "Restart file version different than class version." << endl);
   }

   d_gamma = db->getDouble("d_gamma");

   d_riemann_solve = db->getString("d_riemann_solve");
   d_godunov_order = db->getInteger("d_godunov_order");
   d_corner_transport = db->getString("d_corner_transport");

   int* tmp_nghosts = d_nghosts;
   db->getIntegerArray("d_nghosts", tmp_nghosts, NDIM);
   for (int i = 0; i < NDIM; i++) {
      if (d_nghosts(i) != CELLG) {
         TBOX_ERROR(d_object_name << ": "
            << "Key data `d_nghosts' in restart file != CELLG." << endl);
      }
   }
   int* tmp_fluxghosts = d_fluxghosts;
   db->getIntegerArray("d_fluxghosts", tmp_fluxghosts, NDIM);
   for (int i = 0; i < NDIM; i++) {
      if (d_fluxghosts(i) != FLUXG) {
         TBOX_ERROR(d_object_name << ": "
            << "Key data `d_fluxghosts' in restart file != FLUXG." << endl);
      }
   }

   d_data_problem = db->getString("d_data_problem");

   if (d_data_problem == "SPHERE") {
      d_radius = db->getDouble("d_radius");
      db->getDoubleArray("d_center", d_center, NDIM);
      d_density_inside = db->getDouble("d_density_inside");
      db->getDoubleArray("d_velocity_inside", d_velocity_inside, NDIM);
      d_pressure_inside = db->getDouble("d_pressure_inside");
      d_density_outside = db->getDouble("d_density_outside");
      db->getDoubleArray("d_velocity_outside", d_velocity_outside, NDIM);
      d_pressure_outside = db->getDouble("d_pressure_outside");
   } 

   if ( (d_data_problem == "PIECEWISE_CONSTANT_X") ||
        (d_data_problem == "PIECEWISE_CONSTANT_Y") ||
        (d_data_problem == "PIECEWISE_CONSTANT_Z") ||
        (d_data_problem == "STEP") ) {
      d_number_of_intervals = db->getInteger("d_number_of_intervals");
      if (d_number_of_intervals > 0) {
         d_front_position = db->getDoubleArray("d_front_position");
         d_interval_density = db->getDoubleArray("d_interval_density");
         d_interval_velocity = db->getDoubleArray("d_interval_velocity");
         d_interval_pressure = db->getDoubleArray("d_interval_pressure");
      }
   }

   d_master_bdry_edge_conds = db->getIntegerArray("d_master_bdry_edge_conds");
   d_master_bdry_node_conds = db->getIntegerArray("d_master_bdry_node_conds");

#if (NDIM == 2) 
   d_bdry_edge_density = db->getDoubleArray("d_bdry_edge_density");
   d_bdry_edge_velocity = db->getDoubleArray("d_bdry_edge_velocity");
   d_bdry_edge_pressure = db->getDoubleArray("d_bdry_edge_pressure");
#endif
#if (NDIM == 3)
   d_master_bdry_face_conds = db->getIntegerArray("d_master_bdry_face_conds");

   d_bdry_face_density = db->getDoubleArray("d_bdry_face_density");
   d_bdry_face_velocity = db->getDoubleArray("d_bdry_face_velocity");
   d_bdry_face_pressure = db->getDoubleArray("d_bdry_face_pressure");
#endif

   if (db->keyExists("d_refinement_criteria")) {
      d_refinement_criteria = db->getStringArray("d_refinement_criteria");
   }

   for (int i = 0; i < d_refinement_criteria.getSize(); i++) {

      if (d_refinement_criteria[i] == "DENSITY_DEVIATION") {

         d_density_dev_tol = 
            db->getDoubleArray("d_density_dev_tol");
         d_density_dev = 
            db->getDoubleArray("d_density_dev");
         d_density_dev_time_max = 
            db->getDoubleArray("d_density_dev_time_max");
         d_density_dev_time_min = 
            db->getDoubleArray("d_density_dev_time_min");

      } else if (d_refinement_criteria[i] == "DENSITY_GRADIENT") {

         d_density_grad_tol = 
            db->getDoubleArray("d_density_grad_tol");
         d_density_grad_time_max = 
            db->getDoubleArray("d_density_grad_time_max");
         d_density_grad_time_min = 
            db->getDoubleArray("d_density_grad_time_min");

      } else if  (d_refinement_criteria[i] == "DENSITY_SHOCK") {

         d_density_shock_onset = 
            db->getDoubleArray("d_density_shock_onset");
         d_density_shock_tol = 
            db->getDoubleArray("d_density_shock_tol");
         d_density_shock_time_max = 
            db->getDoubleArray("d_density_shock_time_max");
         d_density_shock_time_min = 
            db->getDoubleArray("d_density_shock_time_min");

      } else if  (d_refinement_criteria[i] == "DENSITY_RICHARDSON") {

         d_density_rich_tol = 
            db->getDoubleArray("d_density_rich_tol");
         d_density_rich_time_max = 
            db->getDoubleArray("d_density_rich_time_max");
         d_density_rich_time_min = 
            db->getDoubleArray("d_density_rich_time_min");

      } else if (d_refinement_criteria[i] == "PRESSURE_DEVIATION") {

         d_pressure_dev_tol = 
            db->getDoubleArray("d_pressure_dev_tol");
         d_pressure_dev = 
            db->getDoubleArray("d_pressure_dev");
         d_pressure_dev_time_max = 
            db->getDoubleArray("d_pressure_dev_time_max");
         d_pressure_dev_time_min = 
            db->getDoubleArray("d_pressure_dev_time_min");

      } else if  (d_refinement_criteria[i] == "PRESSURE_GRADIENT") {

         d_pressure_grad_tol = 
            db->getDoubleArray("d_pressure_grad_tol");
         d_pressure_grad_time_max = 
            db->getDoubleArray("d_pressure_grad_time_max");
         d_pressure_grad_time_min = 
            db->getDoubleArray("d_pressure_grad_time_min");

      } else if  (d_refinement_criteria[i] == "PRESSURE_SHOCK") {

         d_pressure_shock_onset = 
            db->getDoubleArray("d_pressure_shock_onset");
         d_pressure_shock_tol = 
            db->getDoubleArray("d_pressure_shock_tol");
         d_pressure_shock_time_max = 
            db->getDoubleArray("d_pressure_shock_time_max");
         d_pressure_shock_time_min = 
            db->getDoubleArray("d_pressure_shock_time_min");
 
      } else if  (d_refinement_criteria[i] == "PRESSURE_RICHARDSON") {

         d_pressure_rich_tol = 
            db->getDoubleArray("d_pressure_rich_tol");
         d_pressure_rich_time_max = 
            db->getDoubleArray("d_pressure_rich_time_max");
         d_pressure_rich_time_min = 
            db->getDoubleArray("d_pressure_rich_time_min");

      }

   }

}

/*
*************************************************************************
*                                                                       *
* Routines to read boundary data from input database.                   *
*                                                                       *
*************************************************************************
*/

void Euler::readDirichletBoundaryDataEntry(tbox::Pointer<tbox::Database> db,
                                           string& db_name,
                                           int bdry_location_index)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!db.isNull());
   TBOX_ASSERT(!db_name.empty());
#endif
#if (NDIM == 2)
   readStateDataEntry(db,
                      db_name,
                      bdry_location_index,
                      d_bdry_edge_density,
                      d_bdry_edge_velocity,
                      d_bdry_edge_pressure);
#endif
#if (NDIM == 3)
   readStateDataEntry(db,
                      db_name,
                      bdry_location_index,
                      d_bdry_face_density,
                      d_bdry_face_velocity,
                      d_bdry_face_pressure);
#endif
}

void Euler::readStateDataEntry(tbox::Pointer<tbox::Database> db,
                               const string& db_name,
                               int array_indx,
                               tbox::Array<double>& density,
                               tbox::Array<double>& velocity,
                               tbox::Array<double>& pressure)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!db.isNull());
   TBOX_ASSERT(!db_name.empty());
   TBOX_ASSERT(array_indx >= 0);
   TBOX_ASSERT(density.getSize() > array_indx);
   TBOX_ASSERT(velocity.getSize() > array_indx*NDIM);
   TBOX_ASSERT(pressure.getSize() > array_indx);
#endif

   if (db->keyExists("density")) {
      density[array_indx] = db->getDouble("density");
   } else {
      TBOX_ERROR(d_object_name << ": "
         << "`density' entry missing from " << db_name
         << " input database. " << endl);
   }
   if (db->keyExists("velocity")) {
      tbox::Array<double> tmp_vel(0);
      tmp_vel = db->getDoubleArray("velocity");
      if (tmp_vel.getSize() < NDIM) {
         TBOX_ERROR(d_object_name << ": "
                    << "Insufficient number `velocity' values"
                    << " given in " << db_name
                    << " input database." << endl);
      }
      for (int iv = 0; iv < NDIM; iv++) {
         velocity[array_indx*NDIM+iv] = tmp_vel[iv];
      }
   } else {
      TBOX_ERROR(d_object_name << ": "
         << "`velocity' entry missing from " << db_name
         << " input database. " << endl);
   }
   if (db->keyExists("pressure")) {
      pressure[array_indx] = db->getDouble("pressure");
   } else {
      TBOX_ERROR(d_object_name << ": "
         << "`pressure' entry missing from " << db_name 
         << " input database. " << endl);
   }

}

/*
*************************************************************************
*                                                                       *
* Routine to check boundary data when debugging.                        *
*                                                                       *
*************************************************************************
*/

void Euler::checkBoundaryData(int btype,
                              const hier::Patch<NDIM>& patch,
                              const hier::IntVector<NDIM>& ghost_width_to_check,
                              const tbox::Array<int>& scalar_bconds,
                              const tbox::Array<int>& vector_bconds) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
#if (NDIM == 2) 
   TBOX_ASSERT(btype == EDGE2D_BDRY_TYPE || 
          btype == NODE2D_BDRY_TYPE);
#endif
#if (NDIM == 3)
   TBOX_ASSERT(btype == FACE3D_BDRY_TYPE ||
          btype == EDGE3D_BDRY_TYPE ||
          btype == NODE3D_BDRY_TYPE);
#endif
#endif

   const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
   const tbox::Array<hier::BoundaryBox<NDIM> > bdry_boxes =
      pgeom->getCodimensionBoundary(btype);

   hier::VariableDatabase<NDIM>* vdb = hier::VariableDatabase<NDIM>::getDatabase();

   for (int i = 0; i < bdry_boxes.getSize(); i++ ) {
      hier::BoundaryBox<NDIM> bbox = bdry_boxes[i];
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(bbox.getBoundaryType() == btype);
#endif
      int bloc = bbox.getLocationIndex();

      int bscalarcase, bvelocitycase, refbdryloc;
#if (NDIM == 2)
      if (btype == EDGE2D_BDRY_TYPE) {
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(scalar_bconds.getSize() == NUM_2D_EDGES);
         TBOX_ASSERT(vector_bconds.getSize() == NUM_2D_EDGES);
#endif
         bscalarcase = scalar_bconds[bloc];
         bvelocitycase = vector_bconds[bloc];
         refbdryloc = bloc;
      } else { // btype == NODE2D_BDRY_TYPE
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(scalar_bconds.getSize() == NUM_2D_NODES);
         TBOX_ASSERT(vector_bconds.getSize() == NUM_2D_NODES);
#endif
         bscalarcase = scalar_bconds[bloc];
         bvelocitycase = vector_bconds[bloc];
         refbdryloc = d_node_bdry_edge[bloc];
      }
#endif
#if (NDIM == 3)
      if (btype == FACE3D_BDRY_TYPE) {
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(scalar_bconds.getSize() == NUM_3D_FACES);
         TBOX_ASSERT(vector_bconds.getSize() == NUM_3D_FACES);
#endif
         bscalarcase = scalar_bconds[bloc];
         bvelocitycase = vector_bconds[bloc];
         refbdryloc = bloc;
      } else if (btype == EDGE3D_BDRY_TYPE) {
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(scalar_bconds.getSize() == NUM_3D_EDGES);
         TBOX_ASSERT(vector_bconds.getSize() == NUM_3D_EDGES);
#endif
         bscalarcase = scalar_bconds[bloc];
         bvelocitycase = vector_bconds[bloc];
         refbdryloc = d_edge_bdry_face[bloc];
      } else { // btype == NODE3D_BDRY_TYPE
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(scalar_bconds.getSize() == NUM_3D_NODES);
         TBOX_ASSERT(vector_bconds.getSize() == NUM_3D_NODES);
#endif
         bscalarcase = scalar_bconds[bloc];
         bvelocitycase = vector_bconds[bloc];
         refbdryloc = d_node_bdry_face[bloc];
      }
#endif

      int num_bad_values = 0;

#if (NDIM == 2)
      num_bad_values = 
      appu::CartesianBoundaryUtilities2::checkBdryData(
         d_density->getName(),
         patch,
         vdb->mapVariableAndContextToIndex(d_density, getDataContext()), 0, 
         ghost_width_to_check,
         bbox,
         bscalarcase,
         d_bdry_edge_density[refbdryloc]);
#endif
#if (NDIM == 3)
      num_bad_values = 
      appu::CartesianBoundaryUtilities3::checkBdryData(
         d_density->getName(),
         patch,
         vdb->mapVariableAndContextToIndex(d_density, getDataContext()), 0, 
         ghost_width_to_check,
         bbox,
         bscalarcase,
         d_bdry_face_density[refbdryloc]);
#endif
#if (TESTING == 1)
      if (num_bad_values > 0) {
         tbox::perr << "\nEuler Boundary Test FAILED: \n" 
           << "     " << num_bad_values << " bad DENSITY values found for\n"
           << "     boundary type " << btype << " at location " << bloc << endl;
      }
#endif

#if (NDIM == 2)
      num_bad_values =
      appu::CartesianBoundaryUtilities2::checkBdryData(
         d_pressure->getName(),
         patch,
         vdb->mapVariableAndContextToIndex(d_pressure, getDataContext()), 
         0,
         ghost_width_to_check,
         bbox,
         bscalarcase,
         d_bdry_edge_density[refbdryloc]);
#endif
#if (NDIM == 3)
      num_bad_values =
      appu::CartesianBoundaryUtilities3::checkBdryData(
         d_pressure->getName(),
         patch,
         vdb->mapVariableAndContextToIndex(d_pressure, getDataContext()), 
         0,
         ghost_width_to_check,
         bbox,
         bscalarcase,
         d_bdry_face_density[refbdryloc]);
#endif
#if (TESTING == 1)
      if (num_bad_values > 0) {
         tbox::perr << "\nEuler Boundary Test FAILED: \n" 
           << "     " << num_bad_values << " bad PRESSURE values found for\n"
           << "     boundary type " << btype << " at location " << bloc << endl;
      }
#endif

      for (int idir = 0; idir < NDIM; idir++) {

         int vbcase = bscalarcase;
#if (NDIM == 2)
         if (btype == EDGE2D_BDRY_TYPE) {
            if ( (idir == 0 && (bloc == XLO || bloc == XHI)) ||
                 (idir == 1 && (bloc == YLO || bloc == YHI)) ) {
               vbcase = bvelocitycase;
            }
         } else if (btype == NODE2D_BDRY_TYPE) {            
            if ( (idir == 0 && bvelocitycase == XREFLECT_BC) ||
                 (idir == 1 && bvelocitycase == YREFLECT_BC) ) {
               vbcase = bvelocitycase;
            }
         }
#endif
#if (NDIM == 3)
         if (btype == FACE3D_BDRY_TYPE) {
            if ( (idir == 0 && (bloc == XLO || bloc == XHI)) ||
                 (idir == 1 && (bloc == YLO || bloc == YHI)) ||
                 (idir == 2 && (bloc == ZLO || bloc == ZHI)) ) {
               vbcase = bvelocitycase;
            }
         } else if (btype == EDGE3D_BDRY_TYPE || btype == NODE3D_BDRY_TYPE) {
            if ( (idir == 0 && bvelocitycase == XREFLECT_BC) ||
                 (idir == 1 && bvelocitycase == YREFLECT_BC) ||
                 (idir == 2 && bvelocitycase == ZREFLECT_BC) ) {
               vbcase = bvelocitycase;
            }
         }
#endif

#if (NDIM == 2)
         num_bad_values =
         appu::CartesianBoundaryUtilities2::checkBdryData(
               d_velocity->getName(),
               patch,
               vdb->mapVariableAndContextToIndex(d_velocity, getDataContext()), 
               idir,
               ghost_width_to_check,
               bbox,
               vbcase,
               d_bdry_edge_velocity[refbdryloc*NDIM+idir]);
#endif
#if (NDIM == 3)
         num_bad_values =
         appu::CartesianBoundaryUtilities3::checkBdryData(
               d_velocity->getName(),
               patch,
               vdb->mapVariableAndContextToIndex(d_velocity, getDataContext()), 
               idir,
               ghost_width_to_check,
               bbox,
               vbcase,
               d_bdry_face_velocity[refbdryloc*NDIM+idir]);
#endif
#if (TESTING == 1)
         if (num_bad_values > 0) {
            tbox::perr << "\nEuler Boundary Test FAILED: \n" 
              << "     " << num_bad_values 
              << " bad VELOCITY values found in direction " << idir << " for\n"
              << "     boundary type " << btype << " at location " << bloc << endl;
         }
#endif
      }

   }

}
