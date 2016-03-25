//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/examples/LinAdv/LinAdv.C $
// Package:     SAMRAI application
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Numerical routines for single patch in linear advection ex.
//
#include "LinAdv.h" 

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
#include <strstream.h>
#endif
#endif

using namespace std;

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include "tbox/Array.h"
#include "BoundaryBox.h"
#include "BoxArray.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellIndex.h"
#include "CellIterator.h"
#include "CellVariable.h"
#include "FaceData.h"
#include "FaceIndex.h"
#include "FaceVariable.h"
#include "Index.h"
#include "LoadBalancer.h"
#include "tbox/PIO.h"
#include "tbox/RestartManager.h"
#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"
#include "VariableDatabase.h"

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
#include "LinAdvFort.h"
 
// Number of ghosts cells used for each variable quantity
#define CELLG           (4)
#define FACEG           (4)
#define FLUXG           (1)

// defines for initialization
#define PIECEWISE_CONSTANT_X        (10)
#define PIECEWISE_CONSTANT_Y        (11)
#define PIECEWISE_CONSTANT_Z        (12)
#define SINE_CONSTANT_X             (20)
#define SINE_CONSTANT_Y             (21)
#define SINE_CONSTANT_Z             (22)
#define SPHERE                      (40)

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

// Version of LinAdv restart file data
#define LINADV_VERSION (3)

/*
*************************************************************************
*                                                                       *
* The constructor for LinAdv class sets data members to defualt values, *
* creates variables that define the solution state for the linear       *
* advection equation. 
*                                                                       *
* After default values are set, this routine calls getFromRestart()     *
* if execution from a restart file is specified.  Finally,              *
* getFromInput() is called to read values from the given input          *
* database (potentially overriding those found in the restart file).    *
*                                                                       *
*************************************************************************
*/

LinAdv::LinAdv(
   const string& object_name,
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

   d_grid_geometry = grid_geom;

   d_use_nonuniform_workload = false;

   /*
    * hier::Variable quantities that define state of linear advection problem.
    */
   d_uval         = new pdat::CellVariable<NDIM,double>("uval",1);
   d_flux         = new pdat::FaceVariable<NDIM,double>("flux",1);

   /*
    * Default parameters for the numerical method.
    */
   d_godunov_order = 1;
   d_corner_transport = "CORNER_TRANSPORT_1";
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(CELLG == FACEG);
#endif
   d_nghosts = hier::IntVector<NDIM>(CELLG);
   d_fluxghosts = hier::IntVector<NDIM>(FLUXG);

   /*
    * Defaults for problem type and initial data.
    */
   d_data_problem_int = tbox::MathUtilities<int>::getMax();

   int k;

   // SPHERE problem...
   d_radius = tbox::MathUtilities<double>::getSignalingNaN();
   tbox::MathUtilities<double>::setArrayToSignalingNaN(d_center, NDIM);
   d_uval_inside = tbox::MathUtilities<double>::getSignalingNaN();
   d_uval_outside = tbox::MathUtilities<double>::getSignalingNaN();

   d_number_of_intervals = 0;
   d_front_position.resizeArray(0);
   d_interval_uval.resizeArray(0);

   // SINE problem
   d_amplitude = 0.;
   for (k = 0; k < NDIM; k++) d_frequency[k] = 0.;

   /*
    * Defaults for boundary conditions. Set to bogus values
    * for error checking.
    */

#if (NDIM == 2)
   d_scalar_bdry_edge_conds.resizeArray(NUM_2D_EDGES);
   for (int ei = 0; ei < NUM_2D_EDGES; ei++) {
      d_scalar_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
   }

   d_scalar_bdry_node_conds.resizeArray(NUM_2D_NODES);
   d_node_bdry_edge.resizeArray(NUM_2D_NODES);

   for (int ni = 0; ni < NUM_2D_NODES; ni++) {
      d_scalar_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
      d_node_bdry_edge[ni] = BOGUS_BDRY_DATA;
   }

   d_bdry_edge_uval.resizeArray(NUM_2D_EDGES);
   tbox::MathUtilities<double>::setArrayToSignalingNaN(d_bdry_edge_uval);
#endif
#if (NDIM == 3)
   d_scalar_bdry_face_conds.resizeArray(NUM_3D_FACES);
   for (int fi = 0; fi < NUM_3D_FACES; fi++) {
      d_scalar_bdry_face_conds[fi] = BOGUS_BDRY_DATA;
   }

   d_scalar_bdry_edge_conds.resizeArray(NUM_3D_EDGES);
   d_edge_bdry_face.resizeArray(NUM_3D_EDGES);
   for (int ei = 0; ei < NUM_3D_EDGES; ei++) {
      d_scalar_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
      d_edge_bdry_face[ei] = BOGUS_BDRY_DATA;
   }

   d_scalar_bdry_node_conds.resizeArray(NUM_3D_NODES);
   d_node_bdry_face.resizeArray(NUM_3D_NODES);

   for (int ni = 0; ni < NUM_3D_NODES; ni++) {
      d_scalar_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
      d_node_bdry_face[ni] = BOGUS_BDRY_DATA;
   }

   d_bdry_face_uval.resizeArray(NUM_3D_FACES);
   tbox::MathUtilities<double>::setArrayToSignalingNaN(d_bdry_face_uval);
#endif

   /*
    * Initialize object with data read from given input/restart databases.
    */
   bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();
   if (is_from_restart){
      getFromRestart();
   }
   getFromInput(input_db, is_from_restart);

   /*
    * Set problem data to values read from input/restart.
    */

   if (d_data_problem == "PIECEWISE_CONSTANT_X") {
      d_data_problem_int = PIECEWISE_CONSTANT_X;
   } else if (d_data_problem == "PIECEWISE_CONSTANT_Y") {
      d_data_problem_int = PIECEWISE_CONSTANT_Y;
   } else if (d_data_problem == "PIECEWISE_CONSTANT_Z") {
      d_data_problem_int = PIECEWISE_CONSTANT_Z;
   } else if (d_data_problem == "SINE_CONSTANT_X") {
      d_data_problem_int = SINE_CONSTANT_X;
   } else if (d_data_problem == "SINE_CONSTANT_Y") {
      d_data_problem_int = SINE_CONSTANT_Y;
   } else if (d_data_problem == "SINE_CONSTANT_Z") {
      d_data_problem_int = SINE_CONSTANT_Z;
   } else if (d_data_problem == "SPHERE") {
      d_data_problem_int = SPHERE;
   } else {
      TBOX_ERROR(d_object_name << ": "
         << "Unknown d_data_problem string = "
         << d_data_problem << " encountered in constructor" << endl);
   }

   /*
    * Postprocess boundary data from input/restart values.  Note: scalar
    * quantity in this problem cannot have reflective boundary conditions
    * so we reset them to FLOW.
    */
#if (NDIM == 2) 
   for (int i = 0; i < NUM_2D_EDGES; i++) {
      if (d_scalar_bdry_edge_conds[i] == REFLECT_BC) {
         d_scalar_bdry_edge_conds[i] = FLOW_BC;
      }
   }

   for (int i = 0; i < NUM_2D_NODES; i++) {
      if (d_scalar_bdry_node_conds[i] == XREFLECT_BC) {
         d_scalar_bdry_node_conds[i] = XFLOW_BC;
      }
      if (d_scalar_bdry_node_conds[i] == YREFLECT_BC) {
         d_scalar_bdry_node_conds[i] = YFLOW_BC;
      }

      if (d_scalar_bdry_node_conds[i] != BOGUS_BDRY_DATA) {
         d_node_bdry_edge[i] =
	    appu::CartesianBoundaryUtilities2::getEdgeLocationForNodeBdry(
                                            i, d_scalar_bdry_node_conds[i]);
      }
   }
#endif
#if (NDIM == 3)
   for (int i = 0; i < NUM_3D_FACES; i++) {
      if (d_scalar_bdry_face_conds[i] == REFLECT_BC) {
         d_scalar_bdry_face_conds[i] = FLOW_BC;
      }
   }

   for (int i = 0; i < NUM_3D_EDGES; i++) {
      if (d_scalar_bdry_edge_conds[i] == XREFLECT_BC) {
         d_scalar_bdry_edge_conds[i] = XFLOW_BC;
      }
      if (d_scalar_bdry_edge_conds[i] == YREFLECT_BC) {
         d_scalar_bdry_edge_conds[i] = YFLOW_BC;
      }
      if (d_scalar_bdry_edge_conds[i] == ZREFLECT_BC) {
         d_scalar_bdry_edge_conds[i] = ZFLOW_BC;
      }

      if (d_scalar_bdry_edge_conds[i] != BOGUS_BDRY_DATA) {
         d_edge_bdry_face[i] =
            appu::CartesianBoundaryUtilities3::getFaceLocationForEdgeBdry(
                                            i, d_scalar_bdry_edge_conds[i]);
      }
   }

   for (int i = 0; i < NUM_3D_NODES; i++) {
      if (d_scalar_bdry_node_conds[i] == XREFLECT_BC) {
         d_scalar_bdry_node_conds[i] = XFLOW_BC;
      }
      if (d_scalar_bdry_node_conds[i] == YREFLECT_BC) {
         d_scalar_bdry_node_conds[i] = YFLOW_BC;
      }
      if (d_scalar_bdry_node_conds[i] == ZREFLECT_BC) {
         d_scalar_bdry_node_conds[i] = ZFLOW_BC;
      }

      if (d_scalar_bdry_node_conds[i] != BOGUS_BDRY_DATA) {
         d_node_bdry_face[i] =
            appu::CartesianBoundaryUtilities3::getFaceLocationForNodeBdry(
                                            i, d_scalar_bdry_node_conds[i]);
      }
   }

#endif

   stufprobc_(PIECEWISE_CONSTANT_X,PIECEWISE_CONSTANT_Y,PIECEWISE_CONSTANT_Z,
              SINE_CONSTANT_X,SINE_CONSTANT_Y,SINE_CONSTANT_Z,SPHERE,
              CELLG,FACEG,FLUXG);
   
}

/*
*************************************************************************
*                                                                       *
* Empty destructor for LinAdv class.                                    *
*                                                                       *
*************************************************************************
*/

LinAdv::~LinAdv() {
}

/*
*************************************************************************
*                                                                       *
* Register conserved variable (u) (i.e., solution state variable) and   *
* flux variable with hyperbolic integrator that manages storage for     *
* those quantities.  Also, register plot data with Vizamrai or VisIt.   *
*                                                                       *
*************************************************************************
*/

void LinAdv::registerModelVariables(
   algs::HyperbolicLevelIntegrator<NDIM>* integrator)
{

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(integrator != (algs::HyperbolicLevelIntegrator<NDIM>*)NULL);
   TBOX_ASSERT(CELLG == FACEG);
#endif

   integrator->registerVariable(d_uval, d_nghosts,
                                algs::HyperbolicLevelIntegrator<NDIM>::TIME_DEP,
                                d_grid_geometry,
                                "CONSERVATIVE_COARSEN",
                                "CONSERVATIVE_LINEAR_REFINE");

   integrator->registerVariable(d_flux, d_fluxghosts,
                                algs::HyperbolicLevelIntegrator<NDIM>::FLUX,
                                d_grid_geometry,
                                "CONSERVATIVE_COARSEN",
                                "NO_REFINE");

   hier::VariableDatabase<NDIM>* vardb = hier::VariableDatabase<NDIM>::getDatabase();

   if (!(d_vizamrai_writer.isNull())) {
      d_vizamrai_writer->
         registerPlotScalar("U",
                            vardb->mapVariableAndContextToIndex(
                               d_uval, integrator->getPlotContext()));
   }
#ifdef HAVE_HDF5
   if (!(d_visit_writer.isNull())) {
      d_visit_writer->
         registerPlotQuantity("U",
                              "SCALAR",
                            vardb->mapVariableAndContextToIndex(
                               d_uval, integrator->getPlotContext()));
   }
#endif

#ifdef HAVE_HDF5
   if (d_vizamrai_writer.isNull() && d_visit_writer.isNull()) {
      TBOX_WARNING(d_object_name << ": registerModelVariables()"
                   << "\nNeither a Vizamrai nor Visit data writer was"
                   << "\nregistered.  Consequently, no plot data will"
                   << "\nbe written." << endl);
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

void LinAdv::setupLoadBalancer(algs::HyperbolicLevelIntegrator<NDIM>* integrator,
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
void LinAdv::initializeDataOnPatch(hier::Patch<NDIM>& patch,
                                   const double data_time,
                                   const bool initial_time)
{
   (void) data_time;

   if (initial_time) {

      const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
      const double* dx  = pgeom->getDx();
      const double* xlo = pgeom->getXLower();
      const double* xhi = pgeom->getXUpper();
 
      tbox::Pointer< pdat::CellData<NDIM,double> > uval = 
         patch.getPatchData(d_uval, getDataContext());

#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(!uval.isNull());
#endif
      hier::IntVector<NDIM> ghost_cells = uval->getGhostCellWidth();

      const hier::Index<NDIM> ifirst=patch.getBox().lower();
      const hier::Index<NDIM> ilast =patch.getBox().upper();

      if ( (d_data_problem_int == SPHERE) ) {

         initsphere_(d_data_problem_int,dx,xlo,xhi, 
                     ifirst(0),ilast(0),
#if (NDIM>1)
                     ifirst(1),ilast(1),
#endif
#if (NDIM>2)
                     ifirst(2),ilast(2),
#endif
                     ghost_cells(0),
#if (NDIM>1)
                     ghost_cells(1),
#endif
#if (NDIM>2)
                     ghost_cells(2),
#endif

                     uval->getPointer(), 
                     d_uval_inside,
                     d_uval_outside,
                     d_center, d_radius);

      } else if ( d_data_problem_int == SINE_CONSTANT_X ||
                  d_data_problem_int == SINE_CONSTANT_Y ||
                  d_data_problem_int == SINE_CONSTANT_Z ) {

         const double* domain_xlo = d_grid_geometry->getXLower();
         const double* domain_xhi = d_grid_geometry->getXUpper();
         double domain_length[NDIM];
         for (int i = 0; i < NDIM; i++) {
            domain_length[i] = domain_xhi[i] - domain_xlo[i];
         }

         linadvinitsine_(d_data_problem_int,dx,xlo,
                        domain_xlo,domain_length,
                        ifirst(0),ilast(0),
#if (NDIM>1)
                        ifirst(1),ilast(1),
#endif
#if (NDIM>2)
                        ifirst(2),ilast(2),
#endif
                        ghost_cells(0),
#if (NDIM>1)
                        ghost_cells(1),
#endif
#if (NDIM>2)
                        ghost_cells(2),
#endif
                        uval->getPointer(),
                        d_number_of_intervals,
                        d_front_position.getPointer(),
                        d_interval_uval.getPointer(),
                        d_amplitude,
                        d_frequency);

      } else {

         linadvinit_(d_data_problem_int,dx,xlo,xhi,
                        ifirst(0),ilast(0),    
#if (NDIM>1)
                        ifirst(1),ilast(1),
#endif
#if (NDIM>2)
                        ifirst(2),ilast(2),
#endif
                        ghost_cells(0),
#if (NDIM>1)
                        ghost_cells(1),
#endif
#if (NDIM>2)
                        ghost_cells(2),
#endif                    
                        uval->getPointer(), 
                        d_number_of_intervals,
                        d_front_position.getPointer(),
                        d_interval_uval.getPointer());
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

}

/*
*************************************************************************
*                                                                       *
* Compute stable time increment for patch.  Return this value.          *
*                                                                       *
*************************************************************************
*/

double LinAdv::computeStableDtOnPatch(
   hier::Patch<NDIM>& patch,
   const bool initial_time,
   const double dt_time)
{
   (void) initial_time;
   (void) dt_time;

   const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > patch_geom = patch.getPatchGeometry();
   const double* dx  = patch_geom->getDx();
 
   const hier::Index<NDIM> ifirst=patch.getBox().lower();
   const hier::Index<NDIM> ilast =patch.getBox().upper();

   tbox::Pointer< pdat::CellData<NDIM,double> > uval = 
      patch.getPatchData(d_uval, getDataContext());

#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(!uval.isNull());
#endif 
   hier::IntVector<NDIM> ghost_cells = uval->getGhostCellWidth();

   double stabdt;
   stabledt_(dx,
             ifirst(0),ilast(0),
#if (NDIM>1)
             ifirst(1),ilast(1),
#endif
#if (NDIM>2)
             ifirst(2),ilast(2),
#endif
             ghost_cells(0),
#if (NDIM>1)
             ghost_cells(1),
#endif
#if (NDIM>2)
             ghost_cells(2),
#endif 
             d_advection_velocity, 
             uval->getPointer(),
             stabdt);

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

void LinAdv::computeFluxesOnPatch(hier::Patch<NDIM>& patch,
                                  const double time,
                                  const double dt)
{
   (void) time;

#if (NDIM == 3)

   if (d_corner_transport == "CORNER_TRANSPORT_2") {
      compute3DFluxesWithCornerTransport2(patch, dt);
   } else {
      compute3DFluxesWithCornerTransport1(patch, dt);
   }

#endif

#if (NDIM < 3)

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(CELLG == FACEG);
#endif

   const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > patch_geom = patch.getPatchGeometry();
   const double* dx  = patch_geom->getDx();

   hier::Box<NDIM> pbox = patch.getBox();
   const hier::Index<NDIM> ifirst=patch.getBox().lower();
   const hier::Index<NDIM> ilast =patch.getBox().upper();

   tbox::Pointer< pdat::CellData<NDIM,double> > uval =
      patch.getPatchData(d_uval, getDataContext());
   tbox::Pointer< pdat::FaceData<NDIM,double> > flux =
      patch.getPatchData(d_flux, getDataContext());

   /*
    * Verify that the integrator providing the context correctly
    * created it, and that the ghost cell width associated with the
    * context matches the ghosts defined in this class...
    */
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!uval.isNull());
   TBOX_ASSERT(!flux.isNull());
   TBOX_ASSERT(uval->getGhostCellWidth() == d_nghosts);
   TBOX_ASSERT(flux->getGhostCellWidth() == d_fluxghosts);
#endif 

   /*
    * Allocate patch data for temporaries local to this routine.
    */
   pdat::FaceData<NDIM,double> traced_left(pbox, 1, d_nghosts);
   pdat::FaceData<NDIM,double> traced_right(pbox, 1, d_nghosts);

   inittraceflux_(ifirst(0),ilast(0),
#if (NDIM==2)
                  ifirst(1),ilast(1),
#endif
                  uval->getPointer(),
                  traced_left.getPointer(0),
#if (NDIM==2)
                  traced_left.getPointer(1),
#endif
                  traced_right.getPointer(0),
#if (NDIM==2)
                  traced_right.getPointer(1),
#endif
                  flux->getPointer(0)
#if (NDIM==2)
                  ,flux->getPointer(1)
#endif
                  );
 
   if (d_godunov_order >1) {

      /*
       * Prepare temporary data for characteristic tracing.
       */
      int Mcells = 0;
      for (int k=0;k<NDIM;k++) {
         Mcells = tbox::MathUtilities<int>::Max(Mcells, pbox.numberCells(k));
      }

// Face-centered temporary arrays
      tbox::Array<double>  ttedgslp(2*FACEG+1+Mcells);
      tbox::Array<double>  ttraclft(2*FACEG+1+Mcells);
      tbox::Array<double>  ttracrgt(2*FACEG+1+Mcells);

// Cell-centered temporary arrays
      tbox::Array<double>  ttcelslp(2*CELLG+Mcells);

/*
 *  Apply characteristic tracing to compute initial estimate of
 *  traces w^L and w^R at faces.
 *  Inputs: w^L, w^R (traced_left/right)
 *  Output: w^L, w^R
 */
      chartracing0_(dt,
                    ifirst(0),ilast(0),
#if (NDIM==2)
                    ifirst(1),ilast(1),
#endif
                    Mcells,dx[0],d_advection_velocity[0],d_godunov_order,
                    uval->getPointer(),
                    traced_left.getPointer(0),
                    traced_right.getPointer(0),
                    ttcelslp.getPointer(),
                    ttedgslp.getPointer(),
                    ttraclft.getPointer(),
                    ttracrgt.getPointer());

#if (NDIM==2)
      chartracing1_(dt,
                    ifirst(0),ilast(0),ifirst(1),ilast(1),
                    Mcells,dx[1],d_advection_velocity[1],d_godunov_order,
                    uval->getPointer(),
                    traced_left.getPointer(1),
                    traced_right.getPointer(1),
                    ttcelslp.getPointer(),
                    ttedgslp.getPointer(),
                    ttraclft.getPointer(),
                    ttracrgt.getPointer());
#endif

   }  // if (d_godunov_order > 1) ...
 
#if (NDIM==2)
/*
 *  Compute fluxes at faces using the face states computed so far.
 *  Inputs: w^L, w^R (traced_left/right)
 *  Output: F (flux)
 */
// fluxcalculation_(dt,*,1,dx, to get artificial viscosity
// fluxcalculation_(dt,*,0,dx, to get NO artificial viscosity

   fluxcalculation_(dt,1,0,dx,
                    ifirst(0),ilast(0),ifirst(1),ilast(1),
                    d_advection_velocity,
                    uval->getPointer(),
                    flux->getPointer(0),
                    flux->getPointer(1),
                    traced_left.getPointer(0),
                    traced_left.getPointer(1),
                    traced_right.getPointer(0),
                    traced_right.getPointer(1));

/*
 *  Re-compute traces at cell faces with transverse correction applied.
 *  Inputs: F (flux)
 *  Output: w^L, w^R (traced_left/right)
 */
   fluxcorrec_(dt,ifirst(0),ilast(0),ifirst(1),ilast(1),
               dx,d_advection_velocity,
               uval->getPointer(),
               flux->getPointer(0),
               flux->getPointer(1),
               traced_left.getPointer(0),
               traced_left.getPointer(1),
               traced_right.getPointer(0),
               traced_right.getPointer(1));

   boundaryReset(patch, traced_left, traced_right);

/*
 *  Re-compute fluxes with updated traces.
 *  Inputs: w^L, w^R (traced_left/right)
 *  Output: F (flux)
 */
   fluxcalculation_(dt,0,0,dx,
                    ifirst(0),ilast(0),ifirst(1),ilast(1),
                    d_advection_velocity,
                    uval->getPointer(),
                    flux->getPointer(0),
                    flux->getPointer(1),
                    traced_left.getPointer(0),
                    traced_left.getPointer(1),
                    traced_right.getPointer(0),
                    traced_right.getPointer(1));

#endif

//     tbox::plog << "flux values: option1...." << endl;
//     flux->print(pbox, tbox::plog);
#endif
}


/*
*************************************************************************
*                                                                       *
* Compute numerical approximations to flux terms using an extension     *
* to three dimensions of Collella's corner transport upwind approach.   *
* I.E. input value corner_transport = CORNER_TRANSPORT_1                *
*                                                                       *
*************************************************************************
*/
#if (NDIM==3)
void LinAdv::compute3DFluxesWithCornerTransport1(hier::Patch<NDIM>& patch,
                                                 const double dt)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(CELLG == FACEG);
#endif

   const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > patch_geom = patch.getPatchGeometry();
   const double* dx  = patch_geom->getDx();

   hier::Box<NDIM> pbox = patch.getBox();
   const hier::Index<NDIM> ifirst=patch.getBox().lower();
   const hier::Index<NDIM> ilast =patch.getBox().upper();

   tbox::Pointer< pdat::CellData<NDIM,double> > uval =
      patch.getPatchData(d_uval, getDataContext());
   tbox::Pointer< pdat::FaceData<NDIM,double> > flux =
      patch.getPatchData(d_flux, getDataContext());

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!uval.isNull());
   TBOX_ASSERT(!flux.isNull());
   TBOX_ASSERT(uval->getGhostCellWidth() == d_nghosts);
   TBOX_ASSERT(flux->getGhostCellWidth() == d_fluxghosts);
#endif


   /*
    * Allocate patch data for temporaries local to this routine.
    */
   pdat::FaceData<NDIM,double> traced_left(pbox, 1, d_nghosts);
   pdat::FaceData<NDIM,double> traced_right(pbox, 1, d_nghosts);
   pdat::FaceData<NDIM,double> temp_flux(pbox, 1, d_fluxghosts);
   pdat::FaceData<NDIM,double> temp_traced_left(pbox, 1, d_nghosts);
   pdat::FaceData<NDIM,double> temp_traced_right(pbox, 1, d_nghosts);

   inittraceflux_(ifirst(0),ilast(0),
                  ifirst(1),ilast(1),
                  ifirst(2),ilast(2),
                  uval->getPointer(),
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
    if (d_godunov_order >1) {
 
       /*
        * Prepare temporary data for characteristic tracing.
        */
       int Mcells = 0;
      for (int k=0;k<NDIM;k++) {
         Mcells = tbox::MathUtilities<int>::Max(Mcells, pbox.numberCells(k));
      }

      // Face-centered temporary arrays
      tbox::Array<double>  ttedgslp(2*FACEG+1+Mcells);
      tbox::Array<double>  ttraclft(2*FACEG+1+Mcells);
      tbox::Array<double>  ttracrgt(2*FACEG+1+Mcells);

      // Cell-centered temporary arrays
      tbox::Array<double>  ttcelslp(2*CELLG+Mcells);

       /*
        *  Apply characteristic tracing to compute initial estimate of
        *  traces w^L and w^R at faces.
        *  Inputs: w^L, w^R (traced_left/right)
        *  Output: w^L, w^R
        */
      chartracing0_(dt,
                    ifirst(0),ilast(0),
                    ifirst(1),ilast(1),
                    ifirst(2),ilast(2),
                    Mcells,dx[0],d_advection_velocity[0],d_godunov_order,
                    uval->getPointer(),
                    traced_left.getPointer(0),
                    traced_right.getPointer(0),
                    ttcelslp.getPointer(),
                    ttedgslp.getPointer(),
                    ttraclft.getPointer(),
                    ttracrgt.getPointer());

      chartracing1_(dt,
                    ifirst(0),ilast(0),
		    ifirst(1),ilast(1),
                    ifirst(2),ilast(2),
                    Mcells,dx[1],d_advection_velocity[1],d_godunov_order,
                    uval->getPointer(),
                    traced_left.getPointer(1),
                    traced_right.getPointer(1),
                    ttcelslp.getPointer(),
                    ttedgslp.getPointer(),
                    ttraclft.getPointer(),
                    ttracrgt.getPointer());

      chartracing2_(dt,
                    ifirst(0),ilast(0),
		    ifirst(1),ilast(1),
		    ifirst(2),ilast(2),
                    Mcells,dx[2],d_advection_velocity[2],d_godunov_order,
                    uval->getPointer(),
                    traced_left.getPointer(2),
                    traced_right.getPointer(2),
                    ttcelslp.getPointer(),
                    ttedgslp.getPointer(),
                    ttraclft.getPointer(),
                    ttracrgt.getPointer());
   }

   /*
    *  Compute preliminary fluxes at faces using the face states computed 
    *  so far.
    *  Inputs: w^L, w^R (traced_left/right)
    *  Output: F (flux)
    */

//  fluxcalculation_(dt,*,*,1,dx,  to do artificial viscosity
//  fluxcalculation_(dt,*,*,0,dx,  to do NO artificial viscosity
   fluxcalculation_(dt,1,0,0,dx,
                    ifirst(0),ilast(0),ifirst(1),ilast(1),ifirst(2),ilast(2),
                    d_advection_velocity,
                    uval->getPointer(),
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
    *  Inputs: F (flux), w^L, w^R (traced_left/right)
    *  Output: temp_traced_left/right
    */
   fluxcorrec2d_(dt,ifirst(0),ilast(0),ifirst(1),ilast(1),ifirst(2),ilast(2),
                 dx,d_advection_velocity,1,
                 uval->getPointer(),
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
    *  Compute fluxes with partially-corrected trace states.  Store result in
    *  temporary flux vector.
    *  Inputs: temp_traced_left/right
    *  Output: temp_flux
    */
   fluxcalculation_(dt,0,1,0,dx,
                    ifirst(0),ilast(0),ifirst(1),ilast(1),ifirst(2),ilast(2),
                    d_advection_velocity,
                    uval->getPointer(),
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
    *  Inputs: F (flux), w^L, w^R (traced_left/right)
    *  Output: temp_traced_left/right
    */
   fluxcorrec2d_(dt,ifirst(0),ilast(0),ifirst(1),ilast(1),ifirst(2),ilast(2),
                 dx,d_advection_velocity,-1,
                 uval->getPointer(),
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
    *  Inputs: temp_traced_left/right
    *  Output: flux
    */
   fluxcalculation_(dt,1,0,0,dx,
                    ifirst(0),ilast(0),ifirst(1),ilast(1),ifirst(2),ilast(2),
                    d_advection_velocity,
                    uval->getPointer(),
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
    *  Compute the final trace state vectors at cell faces, using transverse
    *  differences of final predicted fluxes.  Store result w^L 
    *  (traced_left) and w^R (traced_right) vectors.
    *  Inputs: temp_flux, flux
    *  Output: w^L, w^R (traced_left/right)
    */
   fluxcorrec3d_(dt,ifirst(0),ilast(0),ifirst(1),ilast(1),ifirst(2),ilast(2),
                dx,d_advection_velocity,
                uval->getPointer(),
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
                    d_advection_velocity,
                    uval->getPointer(),
                    flux->getPointer(0),
                    flux->getPointer(1),
                    flux->getPointer(2),
                    traced_left.getPointer(0),
                    traced_left.getPointer(1),
                    traced_left.getPointer(2),
                    traced_right.getPointer(0),
                    traced_right.getPointer(1),
                    traced_right.getPointer(2));

//     tbox::plog << "flux values: option1...." << endl;
//     flux->print(pbox, tbox::plog);

}
#endif

/*
*************************************************************************
*                                                                       *
* Compute numerical approximations to flux terms using John             *
* Trangenstein's interpretation of the three-dimensional version of     *
* Collella's corner transport upwind approach.                          *
* I.E. input value corner_transport = CORNER_TRANSPORT_2                *
*                                                                       *
*************************************************************************
*/
#if (NDIM==3)
void LinAdv::compute3DFluxesWithCornerTransport2(hier::Patch<NDIM>& patch,
                                                 const double dt)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(CELLG == FACEG);
#endif

   const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > patch_geom = patch.getPatchGeometry();
   const double* dx  = patch_geom->getDx();

   hier::Box<NDIM> pbox = patch.getBox();
   const hier::Index<NDIM> ifirst=patch.getBox().lower();
   const hier::Index<NDIM> ilast =patch.getBox().upper();

   tbox::Pointer< pdat::CellData<NDIM,double> > uval = 
      patch.getPatchData(d_uval, getDataContext());
   tbox::Pointer< pdat::FaceData<NDIM,double> > flux = 
      patch.getPatchData(d_flux, getDataContext());

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!uval.isNull());
   TBOX_ASSERT(!flux.isNull());
   TBOX_ASSERT(uval->getGhostCellWidth() == d_nghosts);
   TBOX_ASSERT(flux->getGhostCellWidth() == d_fluxghosts);
#endif

   /*
    * Allocate patch data for temporaries local to this routine.
    */
   pdat::FaceData<NDIM,double> traced_left(pbox, 1, d_nghosts);
   pdat::FaceData<NDIM,double> traced_right(pbox, 1, d_nghosts);
   pdat::FaceData<NDIM,double> temp_flux(pbox, 1, d_fluxghosts);
   pdat::CellData<NDIM,double> third_state(pbox, 1, d_nghosts);

   /*
    *  Initialize trace fluxes (w^R and w^L) with cell-centered values.
    */
   inittraceflux_(ifirst(0),ilast(0),
		  ifirst(1),ilast(1),
		  ifirst(2),ilast(2),
                  uval->getPointer(),
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
    *  Compute preliminary fluxes at faces using the face states computed 
    *  so far.
    *  Inputs: w^L, w^R (traced_left/right)
    *  Output: F (flux)
    */
   fluxcalculation_(dt,1,1,0,dx,
                    ifirst(0),ilast(0),ifirst(1),ilast(1),ifirst(2),ilast(2),
                    d_advection_velocity,
                    uval->getPointer(),
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
   if (d_godunov_order >1) {

      /*
       * Prepare temporary data for characteristic tracing.
       */
      int Mcells = 0;
      for (int k=0;k<NDIM;k++) {
         Mcells = tbox::MathUtilities<int>::Max(Mcells, pbox.numberCells(k));
      }

      // Face-centered temporary arrays
      tbox::Array<double> ttedgslp(2*FACEG+1+Mcells);
      tbox::Array<double> ttraclft(2*FACEG+1+Mcells);
      tbox::Array<double> ttracrgt(2*FACEG+1+Mcells);

      // Cell-centered temporary arrays
      tbox::Array<double> ttcelslp(2*CELLG+Mcells);

      /*
       *  Apply characteristic tracing to update traces w^L and
       *  w^R at faces.
       *  Inputs: w^L, w^R (traced_left/right)
       *  Output: w^L, w^R
       */
      chartracing0_(dt,
                    ifirst(0),ilast(0),
                    ifirst(1),ilast(1),
                    ifirst(2),ilast(2),
                    Mcells,dx[0],d_advection_velocity[0],d_godunov_order,
                    uval->getPointer(),
                    traced_left.getPointer(0),
                    traced_right.getPointer(0),
                    ttcelslp.getPointer(),
                    ttedgslp.getPointer(),
                    ttraclft.getPointer(),
                    ttracrgt.getPointer());

      chartracing1_(dt,
                    ifirst(0),ilast(0),ifirst(1),ilast(1),
                    ifirst(2),ilast(2),
                    Mcells,dx[1],d_advection_velocity[1],d_godunov_order,
                    uval->getPointer(),
                    traced_left.getPointer(1),
                    traced_right.getPointer(1),
                    ttcelslp.getPointer(),
                    ttedgslp.getPointer(),
                    ttraclft.getPointer(),
                    ttracrgt.getPointer());

      chartracing2_(dt,
                    ifirst(0),ilast(0),ifirst(1),ilast(1),ifirst(2),ilast(2),
                    Mcells,dx[2],d_advection_velocity[2],d_godunov_order,
                    uval->getPointer(),
                    traced_left.getPointer(2),
                    traced_right.getPointer(2),
                    ttcelslp.getPointer(),
                    ttedgslp.getPointer(),
                    ttraclft.getPointer(),
                    ttracrgt.getPointer());

   } //  if (d_godunov_order > 1) ...


   for (int idir = 0; idir < NDIM; idir++) {

      /*
       *    Approximate traces at cell centers (in idir direction) - denoted
       *    1/3 state.
       *    Inputs:  F (flux)
       *    Output:  third_state
       */
      onethirdstate_(dt,dx,idir,
                     ifirst(0),ilast(0),ifirst(1),ilast(1),ifirst(2),ilast(2),
                     d_advection_velocity,
                     uval->getPointer(),
                     flux->getPointer(0),
                     flux->getPointer(1),
                     flux->getPointer(2),
                     third_state.getPointer());
      /*
       *    Compute fluxes using 1/3 state traces, in the two directions OTHER
       *    than idir.
       *    Inputs:  third_state
       *    Output:  temp_flux (only two directions (i.e. those other than idir)
       *             are modified)
       */
      fluxthird_(dt,dx,idir,
                 ifirst(0),ilast(0),ifirst(1),ilast(1),ifirst(2),ilast(2),
                 d_advection_velocity,
                 uval->getPointer(),
                 third_state.getPointer(),
                 temp_flux.getPointer(0),
                 temp_flux.getPointer(1),
                 temp_flux.getPointer(2));

      /*
       *    Compute transverse corrections for the traces in the two directions
       *    (OTHER than idir) using the differenced fluxes computed in those
       *    directions.
       *    Inputs:  temp_flux
       *    Output:  w^L, w^R (traced_left/right)
       */
      fluxcorrecjt_(dt,dx,idir,
                    ifirst(0),ilast(0),ifirst(1),ilast(1),ifirst(2),ilast(2),
                    d_advection_velocity,
                    uval->getPointer(),
                    temp_flux.getPointer(0),
                    temp_flux.getPointer(1),
                    temp_flux.getPointer(2),
                    traced_left.getPointer(0),
                    traced_left.getPointer(1),
                    traced_left.getPointer(2),
                    traced_right.getPointer(0),
                    traced_right.getPointer(1),
                    traced_right.getPointer(2));

   } // loop over directions...
 
   boundaryReset(patch, traced_left, traced_right);
 
   /*
    *  Final flux calculation using corrected trace states.
    *  Inputs:  w^L, w^R (traced_left/right)
    *  Output:  F (flux)
    */
   fluxcalculation_(dt,0,0,0,dx,
                    ifirst(0),ilast(0),ifirst(1),ilast(1),ifirst(2),ilast(2),
                    d_advection_velocity,
                    uval->getPointer(),
                    flux->getPointer(0),
                    flux->getPointer(1),
                    flux->getPointer(2),
                    traced_left.getPointer(0),
                    traced_left.getPointer(1),
                    traced_left.getPointer(2),
                    traced_right.getPointer(0),
                    traced_right.getPointer(1),
                    traced_right.getPointer(2));

//     tbox::plog << "flux values: option2...." << endl;
//     flux->print(pbox, tbox::plog);
}
#endif

/*
*************************************************************************
*                                                                       *
* Update solution variables by performing a conservative                *
* difference with the fluxes calculated in computeFluxesOnPatch().      *
*                                                                       *
*************************************************************************
*/

void LinAdv::conservativeDifferenceOnPatch(hier::Patch<NDIM>& patch,
                                           const double time,
                                           const double dt,
                                           bool at_syncronization)
{
   (void) time;
   (void) dt;
   (void) at_syncronization;

   const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > patch_geom = patch.getPatchGeometry();
   const double* dx  = patch_geom->getDx();

   const hier::Index<NDIM> ifirst=patch.getBox().lower();
   const hier::Index<NDIM> ilast =patch.getBox().upper();

   tbox::Pointer< pdat::CellData<NDIM,double> > uval =
      patch.getPatchData(d_uval, getDataContext());
   tbox::Pointer< pdat::FaceData<NDIM,double> > flux =
      patch.getPatchData(d_flux, getDataContext()); 

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!uval.isNull());
   TBOX_ASSERT(!flux.isNull());
   TBOX_ASSERT(uval->getGhostCellWidth() == d_nghosts);
   TBOX_ASSERT(flux->getGhostCellWidth() == d_fluxghosts);
#endif

#if (NDIM==1)
   consdiff_(ifirst(0),ilast(0),dx,
             flux->getPointer(0), 
             d_advection_velocity,
             uval->getPointer());
#endif
#if (NDIM==2)
   consdiff_(ifirst(0),ilast(0),ifirst(1),ilast(1),dx,
             flux->getPointer(0),
             flux->getPointer(1),
             d_advection_velocity,
             uval->getPointer());
#endif
#if (NDIM==3)
   consdiff_(ifirst(0),ilast(0),ifirst(1),ilast(1),ifirst(2),ilast(2),dx,
             flux->getPointer(0),
             flux->getPointer(1),
             flux->getPointer(2),
             d_advection_velocity,
             uval->getPointer());
#endif

}

/*
*************************************************************************
*                                                                       *
* Reset physical boundary values for special cases, such as those       *
* involving symmetric (i.e., reflective) boundary conditions and        *
* when the "STEP" problem is run.                                       *
*                                                                       *
*************************************************************************
*/
void LinAdv::boundaryReset(hier::Patch<NDIM>& patch,
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

   pdat::CellIndex<NDIM> icell = ifirst;
   hier::BoxArray<NDIM> bdrybox(2*NDIM);
   hier::Index<NDIM> ibfirst = ifirst;
   hier::Index<NDIM> iblast  = ilast;
   int bdry_case,bside;

   for (idir=0;idir<NDIM; idir++) {
      ibfirst(idir) = ifirst(idir) -1; 
      iblast(idir)  = ifirst(idir) -1; 
      bdrybox.getBox(2*idir) = hier::Box<NDIM>(ibfirst,iblast);

      ibfirst(idir) = ilast(idir) +1; 
      iblast(idir)  = ilast(idir) +1; 
      bdrybox.getBox(2*idir+1) = hier::Box<NDIM>(ibfirst,iblast);
   }

   for (idir=0;idir<NDIM; idir++) {
      bside = 2*idir;
#if (NDIM == 2)
      bdry_case = d_scalar_bdry_edge_conds[bside];
#endif
#if (NDIM == 3)
      bdry_case = d_scalar_bdry_face_conds[bside];
#endif
      if (bdry_case == REFLECT_BC) {
         for (pdat::CellIterator<NDIM> ic(bdrybox.getBox(bside)); ic; ic++){
            for ( i = 0; i < num_domain_boxes; i++ ){
               if (domain_boxes.getBox(i).contains(ic()))
                  bdry_cell = false;
            } 
            if (bdry_cell){
              pdat::FaceIndex<NDIM> sidein = pdat::FaceIndex<NDIM>(ic(),idir,1);
              (traced_left)(sidein,0) = (traced_right)(sidein,0);
            }
         } 
      } 
   
      int bnode = 2*idir+1;
#if (NDIM == 2)
      bdry_case = d_scalar_bdry_edge_conds[bnode];
#endif
#if (NDIM == 3)
      bdry_case = d_scalar_bdry_face_conds[bnode];
#endif
      if (bdry_case == REFLECT_BC) {
         for (pdat::CellIterator<NDIM> ic(bdrybox.getBox(bside)); ic; ic++){
            for ( i = 0; i < num_domain_boxes; i++ ){
               if (domain_boxes.getBox(i).contains(ic()))
                  bdry_cell = false;
            } 
            if (bdry_cell){
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
* Set the data in ghost cells corresponding to physical boundary        *
* conditions.  Note that boundary geometry configuration information    *
* (i.e., faces, edges, and nodes) is obtained from the patch geometry   *
* object owned by the patch.                                            *
*                                                                       *
*************************************************************************
*/

void LinAdv::setPhysicalBoundaryConditions(
   hier::Patch<NDIM>& patch, 
   const double fill_time,
   const hier::IntVector<NDIM>& ghost_width_to_fill)
{
   (void) fill_time;

   tbox::Pointer< pdat::CellData<NDIM,double> > uval =
      patch.getPatchData(d_uval, getDataContext());

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!uval.isNull());
#endif
   hier::IntVector<NDIM> ghost_cells = uval->getGhostCellWidth();
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(uval->getGhostCellWidth() == d_nghosts);
#endif

#if (NDIM == 2) 

   /*
    * Set boundary conditions for cells corresponding to patch edges.
    */
   appu::CartesianBoundaryUtilities2::
      fillEdgeBoundaryData("uval", uval,
                           patch,
                           ghost_width_to_fill,
                           d_scalar_bdry_edge_conds,
                           d_bdry_edge_uval);
 
#ifdef DEBUG_CHECK_ASSERTIONS
#if CHECK_BDRY_DATA
   checkBoundaryData(EDGE2D_BDRY_TYPE, patch, ghost_width_to_fill,
                     d_scalar_bdry_edge_conds);
#endif
#endif

   /*
    *  Set boundary conditions for cells corresponding to patch nodes.
    */

   appu::CartesianBoundaryUtilities2::
      fillNodeBoundaryData("uval", uval,
                           patch,
                           ghost_width_to_fill,
                           d_scalar_bdry_node_conds,
                           d_bdry_edge_uval);

#ifdef DEBUG_CHECK_ASSERTIONS
#if CHECK_BDRY_DATA
   checkBoundaryData(NODE2D_BDRY_TYPE, patch, ghost_width_to_fill,
                     d_scalar_bdry_node_conds);
#endif
#endif

#endif // NDIM == 2

#if (NDIM == 3)

   /*
    *  Set boundary conditions for cells corresponding to patch faces.
    */

   appu::CartesianBoundaryUtilities3::
      fillFaceBoundaryData("uval", uval,
                           patch,
                           ghost_width_to_fill,
                           d_scalar_bdry_face_conds,
                           d_bdry_face_uval);
#ifdef DEBUG_CHECK_ASSERTIONS
#if CHECK_BDRY_DATA
   checkBoundaryData(FACE3D_BDRY_TYPE, patch, ghost_width_to_fill,
                     d_scalar_bdry_face_conds);
#endif
#endif

   /*
    *  Set boundary conditions for cells corresponding to patch edges.
    */

   appu::CartesianBoundaryUtilities3::
      fillEdgeBoundaryData("uval", uval,
                           patch,
                           ghost_width_to_fill,
                           d_scalar_bdry_edge_conds,
                           d_bdry_face_uval);
#ifdef DEBUG_CHECK_ASSERTIONS
#if CHECK_BDRY_DATA
   checkBoundaryData(EDGE3D_BDRY_TYPE, patch, ghost_width_to_fill,
                     d_scalar_bdry_edge_conds);
#endif
#endif

   /*
    *  Set boundary conditions for cells corresponding to patch nodes.
    */

   appu::CartesianBoundaryUtilities3::
      fillNodeBoundaryData("uval", uval,
                           patch,
                           ghost_width_to_fill,
                           d_scalar_bdry_node_conds,
                           d_bdry_face_uval);
#ifdef DEBUG_CHECK_ASSERTIONS
#if CHECK_BDRY_DATA
   checkBoundaryData(NODE3D_BDRY_TYPE, patch, ghost_width_to_fill,
                     d_scalar_bdry_node_conds);
#endif
#endif

#endif // NDIM == 3

}

/*
*************************************************************************
*                                                                       *
* Tag cells for refinement using Richardson extrapolation.  Criteria    *
* defined in input.                                                     *
*                                                                       *
*************************************************************************
*/
void LinAdv::tagRichardsonExtrapolationCells(
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
   hier::Box<NDIM> pbox = patch.getBox();

   tbox::Pointer< pdat::CellData<NDIM,int> > tags     = patch.getPatchData(tag_index);

   /*
    * Possible tagging criteria includes 
    *    UVAL_RICHARDSON
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
      int size;
      double tol;
      bool time_allowed;

      if (ref == "UVAL_RICHARDSON") {
         coarsened_fine_var = patch.getPatchData(d_uval, coarsened_fine);
         advanced_coarse_var = patch.getPatchData(d_uval, advanced_coarse);
         size = d_rich_tol.getSize();
         tol = ( ( error_level_number < size) 
                 ? d_rich_tol[error_level_number] 
                 : d_rich_tol[size-1] );
         size = d_rich_time_min.getSize();
         double time_min = ( ( error_level_number < size) 
                 ? d_rich_time_min[error_level_number] 
                 : d_rich_time_min[size-1] );
         size = d_rich_time_max.getSize();
         double time_max = ( ( error_level_number < size) 
                 ? d_rich_time_max[error_level_number] 
                 : d_rich_time_max[size-1] );
         time_allowed = (time_min <= regrid_time) && (time_max > regrid_time);

         if (time_allowed) {

#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!coarsened_fine_var.isNull());
            TBOX_ASSERT(!advanced_coarse_var.isNull());
#endif
            /*
             * We tag wherever the global error > specified tolerance
             * (i.e. d_rich_tol).  The estimated global error is the
             * local truncation error * the approximate number of steps
             * used in the simulation.  Approximate the number of steps as:
             *
             *       steps = L / (s*deltat)
             * where
             *       L = length of problem domain
             *       s = wave speed
             *       delta t = timestep on current level
             *
             */
             const double* xdomainlo = d_grid_geometry->getXLower();
             const double* xdomainhi = d_grid_geometry->getXUpper();
             double max_length = 0.;
             double max_wave_speed = 0.;
             for (int idir = 0; idir < NDIM; idir++) {
                double length = xdomainhi[idir] - xdomainlo[idir];
                if (length > max_length) max_length = length;

                double wave_speed = d_advection_velocity[idir];
                if (wave_speed > max_wave_speed) max_wave_speed = wave_speed;
             }

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
                 *    
                 */
                if (error > tol) {
                   if ((*tags)(ic(),0)) {
                      (*tags)(ic(),0) = RICHARDSON_ALREADY_TAGGED;
                   } else {
                      (*tags)(ic(),0) = RICHARDSON_NEWLY_TAGGED;
                   }
                }
 
            }

         } // time_allowed

      } // if UVAL_RICHARDSON

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
* Tag cells for refinement using gradient detector.  Tagging criteria   *
* defined in input.                                                     *
*                                                                       *
*************************************************************************
*/

void LinAdv::tagGradientDetectorCells(hier::Patch<NDIM>& patch,
   const double regrid_time,
   const bool initial_error,
   const int tag_indx,
   const bool uses_richardson_extrapolation_too)
{

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

   hier::Index<NDIM> ict;

   int not_refine_tag_val = FALSE;
   int refine_tag_val = TRUE;
 
   /*
    * Create a set of temporary tags and set to untagged value.
    */
   tbox::Pointer< pdat::CellData<NDIM,int> > temp_tags  =
                               new pdat::CellData<NDIM,int>(pbox, 1, d_nghosts);
   temp_tags->fillAll(not_refine_tag_val);

   /*
    * Possible tagging criteria includes 
    *    UVAL_DEVIATION, UVAL_GRADIENT, UVAL_SHOCK
    * The criteria is specified over a time interval.
    *
    * Loop over criteria provided and check to make sure we are in the
    * specified time interval.  If so, apply appropriate tagging for
    * the level.
    */
   for (int ncrit = 0; ncrit < d_refinement_criteria.getSize(); ncrit++) {

      string ref = d_refinement_criteria[ncrit];
      tbox::Pointer< pdat::CellData<NDIM, double > > var = 
         patch.getPatchData(d_uval, getDataContext());;  
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(!var.isNull());
#endif   
      hier::IntVector<NDIM> vghost = var->getGhostCellWidth();
      hier::IntVector<NDIM> tagghost = tags->getGhostCellWidth();

      int size = 0;
      double tol = 0.;
      double onset = 0.;
      bool time_allowed = false;

      if (ref == "UVAL_DEVIATION") {
         size = d_dev_tol.getSize();
         tol = ( ( error_level_number < size) 
                 ? d_dev_tol[error_level_number] 
                 : d_dev_tol[size-1] );
         size = d_dev.getSize();
         double dev = ( ( error_level_number < size) 
                 ? d_dev[error_level_number] 
                 : d_dev[size-1] );        
         size = d_dev_time_min.getSize();
         double time_min = ( ( error_level_number < size) 
                 ? d_dev_time_min[error_level_number] 
                 : d_dev_time_min[size-1] );
         size = d_dev_time_max.getSize();
         double time_max = ( ( error_level_number < size) 
                 ? d_dev_time_max[error_level_number] 
                 : d_dev_time_max[size-1] );
         time_allowed = (time_min <= regrid_time) && (time_max > regrid_time);

         if (time_allowed) {

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
                  (*temp_tags)(ic(),0) = refine_tag_val; 
               }
            }
         }
      }

      if (ref == "UVAL_GRADIENT") {
         size = d_grad_tol.getSize();
         tol = ( ( error_level_number < size) 
                 ? d_grad_tol[error_level_number] 
                 : d_grad_tol[size-1] );
         size = d_grad_time_min.getSize();
         double time_min = ( ( error_level_number < size) 
                 ? d_grad_time_min[error_level_number] 
                 : d_grad_time_min[size-1] );
         size = d_grad_time_max.getSize();
         double time_max = ( ( error_level_number < size) 
                 ? d_grad_time_max[error_level_number] 
                 : d_grad_time_max[size-1] );
         time_allowed = (time_min <= regrid_time) && (time_max > regrid_time);

         if (time_allowed) {

            detectgrad_(
#if (NDIM==2)
                  ifirst(0),ilast(0), ifirst(1),ilast(1),
                  vghost(0),tagghost(0),d_nghosts(0),
                  vghost(1),tagghost(1),d_nghosts(1),
#endif
#if (NDIM==3)
                  ifirst(0),ilast(0), ifirst(1),ilast(1),ifirst(2),ilast(2),
                  vghost(0),tagghost(0),d_nghosts(0),
                  vghost(1),tagghost(1),d_nghosts(1),
                  vghost(2),tagghost(2),d_nghosts(2),
#endif
                  dx,
                  tol,
                  refine_tag_val,not_refine_tag_val,
                  var->getPointer(),
                  tags->getPointer(),temp_tags->getPointer());
         }

      }

      if (ref == "UVAL_SHOCK") {
         size = d_shock_tol.getSize();
         tol = ( ( error_level_number < size) 
                 ? d_shock_tol[error_level_number] 
                 : d_shock_tol[size-1] );
         size = d_shock_onset.getSize();
         onset = ( ( error_level_number < size) 
                 ? d_shock_onset[error_level_number] 
                 : d_shock_onset[size-1] );
         size = d_shock_time_min.getSize();
         double time_min = ( ( error_level_number < size) 
                 ? d_shock_time_min[error_level_number] 
                 : d_shock_time_min[size-1] );
         size = d_shock_time_max.getSize();
         double time_max = ( ( error_level_number < size) 
                 ? d_shock_time_max[error_level_number] 
                 : d_shock_time_max[size-1] );
         time_allowed = (time_min <= regrid_time) && (time_max > regrid_time);

         if (time_allowed) {

            detectshock_(
#if (NDIM==2)
                  ifirst(0),ilast(0), ifirst(1),ilast(1),
                  vghost(0),tagghost(0),d_nghosts(0),
                  vghost(1),tagghost(1),d_nghosts(1),
#endif
#if (NDIM==3)
                  ifirst(0),ilast(0), ifirst(1),ilast(1),ifirst(2),ilast(2),
                  vghost(0),tagghost(0),d_nghosts(0),
                  vghost(1),tagghost(1),d_nghosts(1),
                  vghost(2),tagghost(2),d_nghosts(2),
#endif
                  dx,
                  tol,
                  onset,
                  refine_tag_val,not_refine_tag_val,
                  var->getPointer(),
                  tags->getPointer(),temp_tags->getPointer());
         }
      
      }  

   }  // loop over criteria

   /*
    * Adjust temp_tags from those tags set in Richardson extrapolation.
    * Here, we just reset any tags that were set in Richardson extrapolation
    * to be the designated "refine_tag_val".
    */
   if (uses_richardson_extrapolation_too) {
      for (pdat::CellIterator<NDIM> ic(pbox); ic; ic++) {
         if ((*tags)(ic(),0) == RICHARDSON_ALREADY_TAGGED ||
             (*tags)(ic(),0) == RICHARDSON_NEWLY_TAGGED) {
            (*temp_tags)(ic(),0) = refine_tag_val;
         }
      }
   }

   /*
    * Update tags.
    */
   for (pdat::CellIterator<NDIM> ic(pbox); ic; ic++) {
      (*tags)(ic(),0) = (*temp_tags)(ic(),0);
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

void LinAdv::registerVizamraiDataWriter(
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
void LinAdv::registerVisItDataWriter(
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
* Write LinAdv object state to specified stream.                        *
*                                                                       *
*************************************************************************
*/

void LinAdv::printClassData(ostream &os) const 
{
   int j,k;

   os << "\nLinAdv::printClassData..." << endl;
   os << "LinAdv: this = " << (LinAdv*)this << endl;
   os << "d_object_name = " << d_object_name << endl;
   os << "d_grid_geometry = "
      << (geom::CartesianGridGeometry<NDIM>*)d_grid_geometry << endl;


   os << "Parameters for numerical method ..." << endl;
   os << "   d_advection_velocity = ";
   for (j=0;j<NDIM;j++) os << d_advection_velocity[j] << " ";
   os << endl;
   os << "   d_godunov_order = " << d_godunov_order << endl;
   os << "   d_corner_transport = " << d_corner_transport << endl;
   os << "   d_nghosts = " << d_nghosts << endl;
   os << "   d_fluxghosts = " << d_fluxghosts << endl;

   os << "Problem description and initial data..." << endl;
   os << "   d_data_problem = " << d_data_problem << endl;
   os << "   d_data_problem_int = " << d_data_problem << endl;

   os << "       d_radius = " << d_radius << endl;
   os << "       d_center = " ;
      for (j=0;j<NDIM;j++) os << d_center[j]<<" ";
   os << endl;
   os << "       d_uval_inside = " << d_uval_inside << endl;
   os << "       d_uval_outside = " << d_uval_outside << endl;

   os << "       d_number_of_intervals = " << d_number_of_intervals << endl;
   os << "       d_front_position = ";
   for (k = 0; k < d_number_of_intervals-1; k++) {
      os << d_front_position[k] << "  ";
   }
   os << endl;
   os << "       d_interval_uval = " << endl;
   for (k = 0; k < d_number_of_intervals; k++) {
      os << "            " << d_interval_uval[k] << endl;
   }
   os << "   Boundary condition data " << endl;

#if (NDIM == 2) 
   for (j = 0; j < d_scalar_bdry_edge_conds.getSize(); j++) {
      os << "       d_scalar_bdry_edge_conds[" << j << "] = "
         << d_scalar_bdry_edge_conds[j] << endl;
      if (d_scalar_bdry_edge_conds[j] == DIRICHLET_BC) {
         os << "         d_bdry_edge_uval[" << j << "] = "
            << d_bdry_edge_uval[j] << endl;
      }
   }
   os << endl;
   for (j = 0; j < d_scalar_bdry_node_conds.getSize(); j++) {
      os << "       d_scalar_bdry_node_conds[" << j << "] = "
         << d_scalar_bdry_node_conds[j] << endl;
      os << "       d_node_bdry_edge[" << j << "] = "
         << d_node_bdry_edge[j] << endl;
   }
#endif
#if (NDIM == 3)
   for (j = 0; j < d_scalar_bdry_face_conds.getSize(); j++) {
      os << "       d_scalar_bdry_face_conds[" << j << "] = "
         << d_scalar_bdry_face_conds[j] << endl;
      if (d_scalar_bdry_face_conds[j] == DIRICHLET_BC) {
         os << "         d_bdry_face_uval[" << j << "] = "
            << d_bdry_face_uval[j] << endl;
      }
   }
   os << endl;
   for (j = 0; j < d_scalar_bdry_edge_conds.getSize(); j++) {
      os << "       d_scalar_bdry_edge_conds[" << j << "] = "
         << d_scalar_bdry_edge_conds[j] << endl;
      os << "       d_edge_bdry_face[" << j << "] = "
         << d_edge_bdry_face[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_scalar_bdry_node_conds.getSize(); j++) {
      os << "       d_scalar_bdry_node_conds[" << j << "] = "
         << d_scalar_bdry_node_conds[j] << endl;
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
   for (j = 0; j < d_dev_tol.getSize(); j++) {
      os << "       d_dev_tol[" << j << "] = "
         << d_dev_tol[j] << endl;
   }
  for (j = 0; j < d_dev.getSize(); j++) {
      os << "       d_dev[" << j << "] = "
         << d_dev[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_dev_time_max.getSize(); j++) {
      os << "       d_dev_time_max[" << j << "] = "
         << d_dev_time_max[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_dev_time_min.getSize(); j++) {
      os << "       d_dev_time_min[" << j << "] = "
         << d_dev_time_min[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_grad_tol.getSize(); j++) {
      os << "       d_grad_tol[" << j << "] = "
         << d_grad_tol[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_grad_time_max.getSize(); j++) {
      os << "       d_grad_time_max[" << j << "] = "
         << d_grad_time_max[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_grad_time_min.getSize(); j++) {
      os << "       d_grad_time_min[" << j << "] = "
         << d_grad_time_min[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_shock_onset.getSize(); j++) {
      os << "       d_shock_onset[" << j << "] = "
         << d_shock_onset[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_shock_tol.getSize(); j++) {
      os << "       d_shock_tol[" << j << "] = "
         << d_shock_tol[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_shock_time_max.getSize(); j++) {
      os << "       d_shock_time_max[" << j << "] = "
         << d_shock_time_max[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_shock_time_min.getSize(); j++) {
      os << "       d_shock_time_min[" << j << "] = "
         << d_shock_time_min[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_rich_tol.getSize(); j++) {
      os << "       d_rich_tol[" << j << "] = "
         << d_rich_tol[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_rich_time_max.getSize(); j++) {
      os << "       d_rich_time_max[" << j << "] = "
         << d_rich_time_max[j] << endl;
   }
   os << endl;
   for (j = 0; j < d_rich_time_min.getSize(); j++) {
      os << "       d_rich_time_min[" << j << "] = "
         << d_rich_time_min[j] << endl;
   }
   os << endl;

}


/*
*************************************************************************
*                                                                       *
* Read data members from input.  All values set from restart can be	*
* overridden by values in the input database.
*                                                                       *
*************************************************************************
*/
void LinAdv::getFromInput( 
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

   if (db->keyExists("advection_velocity")){
      db->getDoubleArray("advection_velocity",
                         d_advection_velocity, NDIM);
   } else {
      TBOX_ERROR(d_object_name << ":  "
         << "Key data `advection_velocity' not found in input.");
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

   if (db->keyExists("Refinement_data")){
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

            if ( !(error_key == "UVAL_DEVIATION" ||
                   error_key == "UVAL_GRADIENT" ||
                   error_key == "UVAL_SHOCK" ||
                   error_key == "UVAL_RICHARDSON" ) ) { 
                TBOX_ERROR(d_object_name << ": "
                         << "Unknown refinement criteria: " << error_key
                         << "\nin input." << endl);
            } else {
               error_db = refine_db->getDatabase(error_key);
               ref_keys_defined[def_key_cnt] = error_key;
               def_key_cnt++;
            }

            if (!error_db.isNull() && error_key == "UVAL_DEVIATION") {

               if (error_db->keyExists("dev_tol")) {
                  d_dev_tol = 
                  error_db->getDoubleArray("dev_tol");
               } else {
                  TBOX_ERROR(d_object_name << ": "
                           << "No key `dev_tol' found in data for "
                           << error_key << endl);
               }

               if (error_db->keyExists("uval_dev")) {
                  d_dev = 
                  error_db->getDoubleArray("uval_dev");
               } else {
                  TBOX_ERROR(d_object_name << ": "
                           << "No key `uval_dev' found in data for "
                           << error_key << endl);
               }

               if (error_db->keyExists("time_max")){
                  d_dev_time_max = 
                  error_db->getDoubleArray("time_max");
               } else {
                  d_dev_time_max.resizeArray(1);
                  d_dev_time_max[0] = tbox::MathUtilities<double>::getMax();
               }

               if (error_db->keyExists("time_min")){
                  d_dev_time_min = 
                  error_db->getDoubleArray("time_min");
               } else {
                  d_dev_time_min.resizeArray(1);
                  d_dev_time_min[0] = 0.;
               } 

            }
      
            if (!error_db.isNull() && error_key == "UVAL_GRADIENT") {

               if (error_db->keyExists("grad_tol")) {
                  d_grad_tol = 
                  error_db->getDoubleArray("grad_tol");
               } else {
                  TBOX_ERROR(d_object_name << ": "
                           << "No key `grad_tol' found in data for "
                           << error_key << endl);
               }

               if (error_db->keyExists("time_max")){
                  d_grad_time_max = 
                  error_db->getDoubleArray("time_max");
               } else {
                  d_grad_time_max.resizeArray(1);
                  d_grad_time_max[0] = tbox::MathUtilities<double>::getMax();
               }

               if (error_db->keyExists("time_min")){
                  d_grad_time_min = 
                  error_db->getDoubleArray("time_min");
               } else {
                  d_grad_time_min.resizeArray(1);
                  d_grad_time_min[0] = 0.;
               } 

            }
              
            if (!error_db.isNull() && error_key == "UVAL_SHOCK") {

               if (error_db->keyExists("shock_onset")) {
                  d_shock_onset = 
                  error_db->getDoubleArray("shock_onset");
               } else {
                  TBOX_ERROR(d_object_name << ": "
                           << "No key `shock_onset' found in data for "
                           << error_key << endl);
               }

               if (error_db->keyExists("shock_tol")) {
                  d_shock_tol = 
                  error_db->getDoubleArray("shock_tol");
               } else {
                  TBOX_ERROR(d_object_name << ": "
                           << "No key `shock_tol' found in data for "
                           << error_key << endl);
               }

               if (error_db->keyExists("time_max")){
                  d_shock_time_max = 
                  error_db->getDoubleArray("time_max");
               } else {
                  d_shock_time_max.resizeArray(1);
                  d_shock_time_max[0] = tbox::MathUtilities<double>::getMax();
               }

               if (error_db->keyExists("time_min")){
                  d_shock_time_min = 
                  error_db->getDoubleArray("time_min");
               } else {
                  d_shock_time_min.resizeArray(1);
                  d_shock_time_min[0] = 0.;
               }

            }

            if (!error_db.isNull() && error_key == "UVAL_RICHARDSON") {

               if (error_db->keyExists("rich_tol")) {
                  d_rich_tol = 
                  error_db->getDoubleArray("rich_tol");
               } else {
                  TBOX_ERROR(d_object_name << ": "
                           << "No key `rich_tol' found in data for "
                           << error_key << endl);
               }

               if (error_db->keyExists("time_max")){
                  d_rich_time_max = 
                  error_db->getDoubleArray("time_max");
               } else {
                  d_rich_time_max.resizeArray(1);
                  d_rich_time_max[0] = tbox::MathUtilities<double>::getMax();
               }

               if (error_db->keyExists("time_min")){
                  d_rich_time_min = 
                  error_db->getDoubleArray("time_min");
               } else {
                  d_rich_time_min.resizeArray(1);
                  d_rich_time_min[0] = 0.;
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

   } // refine db entry exists 

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
         if (init_data_db->keyExists("uval_inside")) {
            d_uval_inside = init_data_db->getDouble("uval_inside");
         } else {
            TBOX_ERROR(d_object_name << ": "
               << "`uval_inside' input required for "
               << "SPHERE problem." << endl);
         }
         if (init_data_db->keyExists("uval_outside")) {
            d_uval_outside = init_data_db->getDouble("uval_outside");
         } else {
            TBOX_ERROR(d_object_name << ": "
               << "`uval_outside' input required for "
               << "SPHERE problem." << endl);
         }

         found_problem_data = true;

      } 

      if (!found_problem_data && 
          (d_data_problem == "PIECEWISE_CONSTANT_X") ||
          (d_data_problem == "PIECEWISE_CONSTANT_Y") ||
          (d_data_problem == "PIECEWISE_CONSTANT_Z") ||
          (d_data_problem == "SINE_CONSTANT_X") ||
          (d_data_problem == "SINE_CONSTANT_Y") ||
          (d_data_problem == "SINE_CONSTANT_Z")) {

          int idir = 0;
          if (d_data_problem == "PIECEWISE_CONSTANT_Y") {
             if (NDIM < 2) {
                TBOX_ERROR(d_object_name << ": `PIECEWISE_CONSTANT_Y' "
                           << "problem invalid in 1 dimension." << endl);
             }
             idir = 1;
          }

          if (d_data_problem == "PIECEWISE_CONSTANT_Z") {
             if (NDIM < 3) {
                TBOX_ERROR(d_object_name << ": `PIECEWISE_CONSTANT_Z' "
                           << "problem invalid in 1 or 2 dimensions." << endl);
             }
             idir = 2;
          }

          tbox::Array<string> init_data_keys = init_data_db->getAllKeys();

          if (init_data_db->keyExists("front_position")) {
             d_front_position = init_data_db->getDoubleArray("front_position");
          } else {
             TBOX_ERROR(d_object_name << ": "
                << "`front_position' input required for "
                << d_data_problem << " problem." << endl);
          }

          d_number_of_intervals = 
             tbox::MathUtilities<int>::Min( d_front_position.getSize()+1,
                                            init_data_keys.getSize()-1);

          d_front_position.resizeArray(d_front_position.getSize()+1);
          d_front_position[d_front_position.getSize()-1] = 
             d_grid_geometry->getXUpper()[idir]; 

          d_interval_uval.resizeArray(d_number_of_intervals);

          int i = 0;
          int nkey = 0;
          bool found_interval_data = false;

          while (   !found_interval_data 
                 && (i < d_number_of_intervals)
                 && (nkey < init_data_keys.getSize()) ) {

             if ( !(init_data_keys[nkey] == "front_position") ) {
   
                tbox::Pointer<tbox::Database> interval_db = 
                   init_data_db->getDatabase(init_data_keys[nkey]);
   
                if (interval_db->keyExists("uval")) {
                   d_interval_uval[i] = interval_db->getDouble("uval");
                } else {
                   TBOX_ERROR(d_object_name << ": "
                          << "`uval' data missing in input for key = "
                          << init_data_keys[nkey] << endl);
                }
                i++;
   
                found_interval_data = (i == d_number_of_intervals);
   
             }
   
             nkey++;
   
          }

          if ((d_data_problem == "SINE_CONSTANT_X") ||
              (d_data_problem == "SINE_CONSTANT_Y") ||
              (d_data_problem == "SINE_CONSTANT_Z")) {
               if (init_data_db->keyExists("amplitude")) {
                  d_amplitude = init_data_db->getDouble("amplitude");
               }
               if (init_data_db->keyExists("frequency")) {
                  init_data_db->getDoubleArray("frequency",d_frequency,NDIM);
               } else {
                  TBOX_ERROR(d_object_name << ": "
                   << "`frequency' input required for SINE problem." << endl);
               }
          }

          if (!found_interval_data) {
             TBOX_ERROR(d_object_name << ": "
                        << "Insufficient interval data given in input"
                        << " for PIECEWISE_CONSTANT_*problem." << endl);
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

   if (db->keyExists("Boundary_data")) {
      
      tbox::Pointer<tbox::Database> bdry_db = db->getDatabase("Boundary_data");
      
#if (NDIM == 2)
      appu::CartesianBoundaryUtilities2::readBoundaryInput(this,
                                                     bdry_db,
                                                     d_scalar_bdry_edge_conds,
                                                     d_scalar_bdry_node_conds,
                                                     periodic);
#endif
#if (NDIM == 3)
      appu::CartesianBoundaryUtilities3::readBoundaryInput(this,
                                                     bdry_db,
                                                     d_scalar_bdry_face_conds,
                                                     d_scalar_bdry_edge_conds,
                                                     d_scalar_bdry_node_conds,
                                                     periodic);
#endif
      
   } else {
      TBOX_ERROR(d_object_name << ": "
                 << "Key data `Boundary_data' not found in input. " << endl);
   }
   
}


/*
*************************************************************************
*                                                                       *
* Routines to put/get data members to/from restart database.            *
*                                                                       *
*************************************************************************
*/

void LinAdv::putToDatabase( 
   tbox::Pointer<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!db.isNull());
#endif

   db->putInteger("LINADV_VERSION",LINADV_VERSION);

   db->putDoubleArray("d_advection_velocity",d_advection_velocity,NDIM);

   db->putInteger("d_godunov_order", d_godunov_order);
   db->putString("d_corner_transport", d_corner_transport);
   db->putIntegerArray("d_nghosts", (int*)d_nghosts, NDIM);
   db->putIntegerArray("d_fluxghosts", (int*)d_fluxghosts, NDIM);

   db->putString("d_data_problem", d_data_problem);
   
   if (d_data_problem == "SPHERE") {
      db->putDouble("d_radius", d_radius);
      db->putDoubleArray("d_center", d_center, NDIM);
      db->putDouble("d_uval_inside", d_uval_inside);
      db->putDouble("d_uval_outside", d_uval_outside);
   }
   
   if ( (d_data_problem == "PIECEWISE_CONSTANT_X") ||
        (d_data_problem == "PIECEWISE_CONSTANT_Y") ||
        (d_data_problem == "PIECEWISE_CONSTANT_Z") ||
        (d_data_problem == "SINE_CONSTANT_X") ||
        (d_data_problem == "SINE_CONSTANT_Y") ||
        (d_data_problem == "SINE_CONSTANT_Z") ) {
      db->putInteger("d_number_of_intervals", d_number_of_intervals);
      if (d_number_of_intervals > 0) {
         db->putDoubleArray("d_front_position", d_front_position);   
         db->putDoubleArray("d_interval_uval", d_interval_uval);   
      }
   }

   db->putIntegerArray("d_scalar_bdry_edge_conds", d_scalar_bdry_edge_conds);
   db->putIntegerArray("d_scalar_bdry_node_conds", d_scalar_bdry_node_conds);

#if (NDIM == 2) 
   db->putDoubleArray("d_bdry_edge_uval", d_bdry_edge_uval);
#endif
#if (NDIM == 3)
   db->putIntegerArray("d_scalar_bdry_face_conds", d_scalar_bdry_face_conds);
   db->putDoubleArray("d_bdry_face_uval", d_bdry_face_uval);
#endif

   if (d_refinement_criteria.getSize() > 0) {
      db->putStringArray("d_refinement_criteria", d_refinement_criteria);
   }
   for (int i = 0; i < d_refinement_criteria.getSize(); i++) {

      if (d_refinement_criteria[i] == "UVAL_DEVIATION") {
         db->putDoubleArray("d_dev_tol", d_dev_tol);
         db->putDoubleArray("d_dev", d_dev);
         db->putDoubleArray("d_dev_time_max", d_dev_time_max);
         db->putDoubleArray("d_dev_time_min", d_dev_time_min);
      } else if (d_refinement_criteria[i] == "UVAL_GRADIENT") {
         db->putDoubleArray("d_grad_tol", d_grad_tol);
         db->putDoubleArray("d_grad_time_max", d_grad_time_max);
         db->putDoubleArray("d_grad_time_min", d_grad_time_min);
      } else if  (d_refinement_criteria[i] == "UVAL_SHOCK") {
         db->putDoubleArray("d_shock_onset", d_shock_onset);
         db->putDoubleArray("d_shock_tol", d_shock_tol);
         db->putDoubleArray("d_shock_time_max", d_shock_time_max);
         db->putDoubleArray("d_shock_time_min", d_shock_time_min);
      } else if  (d_refinement_criteria[i] == "UVAL_RICHARDSON") {
         db->putDoubleArray("d_rich_tol", d_rich_tol);
         db->putDoubleArray("d_rich_time_max", d_rich_time_max);
         db->putDoubleArray("d_rich_time_min", d_rich_time_min);
      }

   }

}


/*
*************************************************************************
*                                                                       *
*    Access class information from restart database.                    *
*                                                                       *
*************************************************************************
*/
void LinAdv::getFromRestart() 
{
   tbox::Pointer<tbox::Database> root_db = 
      tbox::RestartManager::getManager()->getRootDatabase();

   tbox::Pointer<tbox::Database> db;
   if ( root_db->isDatabase(d_object_name) ) {
      db = root_db->getDatabase(d_object_name);
   } else {
      TBOX_ERROR("Restart database corresponding to "
              << d_object_name << " not found in restart file.");
   }

   int ver = db->getInteger("LINADV_VERSION");
   if (ver != LINADV_VERSION){
      TBOX_ERROR(d_object_name << ":  "
              << "Restart file version different than class version.");
   }

   db->getDoubleArray("d_advection_velocity",d_advection_velocity,NDIM);

   d_godunov_order = db->getInteger("d_godunov_order");
   d_corner_transport = db->getString("d_corner_transport");

   int* tmp_nghosts = d_nghosts;
   db->getIntegerArray("d_nghosts", tmp_nghosts, NDIM);
   if ( !(d_nghosts == CELLG) ) {
      TBOX_ERROR(d_object_name << ": "
                 << "Key data `d_nghosts' in restart file != CELLG." << endl);
   }
   int* tmp_fluxghosts = d_fluxghosts;
   db->getIntegerArray("d_fluxghosts", tmp_fluxghosts, NDIM);
   if ( !(d_fluxghosts == FLUXG) ) {
      TBOX_ERROR(d_object_name << ": "
         << "Key data `d_fluxghosts' in restart file != FLUXG." << endl);
   }

   d_data_problem = db->getString("d_data_problem");

   if (d_data_problem == "SPHERE") {
      d_data_problem_int = SPHERE;
      d_radius = db->getDouble("d_radius");
      db->getDoubleArray("d_center", d_center, NDIM);
      d_uval_inside = db->getDouble("d_uval_inside");
      d_uval_outside = db->getDouble("d_uval_outside");
   }
   
   if ( (d_data_problem == "PIECEWISE_CONSTANT_X") ||
        (d_data_problem == "PIECEWISE_CONSTANT_Y") ||
        (d_data_problem == "PIECEWISE_CONSTANT_Z") ||
        (d_data_problem == "SINE_CONSTANT_X") ||
        (d_data_problem == "SINE_CONSTANT_Y") ||
        (d_data_problem == "SINE_CONSTANT_Z") ) {
      d_number_of_intervals = db->getInteger("d_number_of_intervals");
      if (d_number_of_intervals > 0) {
         d_front_position = db->getDoubleArray("d_front_position");
         d_interval_uval = db->getDoubleArray("d_interval_uval");
      }
   }

   d_scalar_bdry_edge_conds = db->getIntegerArray("d_scalar_bdry_edge_conds");
   d_scalar_bdry_node_conds = db->getIntegerArray("d_scalar_bdry_node_conds");

#if (NDIM == 2) 
   d_bdry_edge_uval = db->getDoubleArray("d_bdry_edge_uval");
#endif
#if (NDIM == 3)
   d_scalar_bdry_face_conds = db->getIntegerArray("d_scalar_bdry_face_conds");

   d_bdry_face_uval = db->getDoubleArray("d_bdry_face_uval");
#endif

   if (db->keyExists("d_refinement_criteria")) {
      d_refinement_criteria = db->getStringArray("d_refinement_criteria");
   }
   for (int i = 0; i < d_refinement_criteria.getSize(); i++) {

      if (d_refinement_criteria[i] == "UVAL_DEVIATION") {
         d_dev_tol = db->getDoubleArray("d_dev_tol");
         d_dev_time_max = db->getDoubleArray("d_dev_time_max");
         d_dev_time_min = db->getDoubleArray("d_dev_time_min");
      } else if (d_refinement_criteria[i] == "UVAL_GRADIENT") {
         d_grad_tol = db->getDoubleArray("d_grad_tol");
         d_grad_time_max = db->getDoubleArray("d_grad_time_max");
         d_grad_time_min = db->getDoubleArray("d_grad_time_min");
      } else if  (d_refinement_criteria[i] == "UVAL_SHOCK") {
         d_shock_onset = db->getDoubleArray("d_shock_onset");
         d_shock_tol = db->getDoubleArray("d_shock_tol");
         d_shock_time_max = db->getDoubleArray("d_shock_time_max");
         d_shock_time_min = db->getDoubleArray("d_shock_time_min");
      } else if  (d_refinement_criteria[i] == "UVAL_RICHARDSON") {
         d_rich_tol = db->getDoubleArray("d_rich_tol");
         d_rich_time_max = db->getDoubleArray("d_rich_time_max");
         d_rich_time_min = db->getDoubleArray("d_rich_time_min");
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

void LinAdv::readDirichletBoundaryDataEntry(tbox::Pointer<tbox::Database> db,
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
                      d_bdry_edge_uval);
#endif
#if (NDIM == 3)
   readStateDataEntry(db,
                      db_name,
                      bdry_location_index,
                      d_bdry_face_uval);
#endif
}

void LinAdv::readStateDataEntry(tbox::Pointer<tbox::Database> db,
                               const string& db_name,
                               int array_indx,
                               tbox::Array<double>& uval)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!db.isNull());
   TBOX_ASSERT(!db_name.empty());
   TBOX_ASSERT(array_indx >= 0);
   TBOX_ASSERT(uval.getSize() > array_indx);
#endif

   if (db->keyExists("uval")) {
      uval[array_indx] = db->getDouble("uval");
   } else {
      TBOX_ERROR(d_object_name << ": "
         << "`uval' entry missing from " << db_name
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

void LinAdv::checkBoundaryData(int btype,
                              const hier::Patch<NDIM>& patch,
                              const hier::IntVector<NDIM>& ghost_width_to_check,
                              const tbox::Array<int>& scalar_bconds) const
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

      int bscalarcase, refbdryloc;
#if (NDIM == 2)
      if (btype == EDGE2D_BDRY_TYPE) {
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(scalar_bconds.getSize() == NUM_2D_EDGES);
#endif
         bscalarcase = scalar_bconds[bloc];
         refbdryloc = bloc;
      } else { // btype == NODE2D_BDRY_TYPE
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(scalar_bconds.getSize() == NUM_2D_NODES);
#endif
         bscalarcase = scalar_bconds[bloc];
         refbdryloc = d_node_bdry_edge[bloc];
      }
#endif
#if (NDIM == 3)
      if (btype == FACE3D_BDRY_TYPE) {
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(scalar_bconds.getSize() == NUM_3D_FACES);
#endif
         bscalarcase = scalar_bconds[bloc];
         refbdryloc = bloc;
      } else if (btype == EDGE3D_BDRY_TYPE) {
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(scalar_bconds.getSize() == NUM_3D_EDGES);
#endif
         bscalarcase = scalar_bconds[bloc];
         refbdryloc = d_edge_bdry_face[bloc];
      } else { // btype == NODE3D_BDRY_TYPE
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(scalar_bconds.getSize() == NUM_3D_NODES);
#endif
         bscalarcase = scalar_bconds[bloc];
         refbdryloc = d_node_bdry_face[bloc];
      }
#endif

      int num_bad_values = 0;

#if (NDIM == 2)
      num_bad_values = 
      appu::CartesianBoundaryUtilities2::checkBdryData(
         d_uval->getName(),
         patch,
         vdb->mapVariableAndContextToIndex(d_uval, getDataContext()), 0, 
         ghost_width_to_check,
         bbox,
         bscalarcase,
         d_bdry_edge_uval[refbdryloc]);
#endif
#if (NDIM == 3)
      num_bad_values = 
      appu::CartesianBoundaryUtilities3::checkBdryData(
         d_uval->getName(),
         patch,
         vdb->mapVariableAndContextToIndex(d_uval, getDataContext()), 0, 
         ghost_width_to_check,
         bbox,
         bscalarcase,
         d_bdry_face_uval[refbdryloc]);
#endif
#if (TESTING == 1)
      if (num_bad_values > 0) {
         tbox::perr << "\nLinAdv Boundary Test FAILED: \n" 
           << "     " << num_bad_values << " bad UVAL values found for\n"
           << "     boundary type " << btype << " at location " << bloc << endl;
      }
#endif

   }

}

