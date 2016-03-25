//
// File:        ConvDiff.C
// Package:     SAMRAI application
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 712 $
// Modified:    $Date: 2005-11-08 13:45:37 -0800 (Tue, 08 Nov 2005) $
// Description: Numerical routines for single patch in convection
//              diffusion example.
//

#include "ConvDiff.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

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

#include "tbox/Array.h"
#include "BoundaryBox.h"
#include "BoxArray.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellIndex.h"
#include "CellIterator.h"
#include "CellVariable.h"
#include "tbox/IEEE.h"
#include "Index.h"
#include "tbox/PIO.h"
#include "tbox/RestartManager.h"
#include "tbox/Utilities.h"
#include "VariableDatabase.h"

//integer constants for boundary conditions
#define CHECK_BDRY_DATA  (1)
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


// external definitions for Fortran numerical routines
#include "ConvDiffFort.h"

// defines for initialization
#define SPHERE          (40)

// Number of ghosts cells used for each variable quantity.
#define CELLG           (1)

// Define class version number
#define CONV_DIFF_VERSION (2)

/*
*************************************************************************
*                                                                       *
* The constructor for ConvDiff class sets data members to defualt       *
* values, creates variables that define the solution state for the      *
* convection diffusion equation. 
*                                                                       *
* After default values are set, this routine calls getFromRestart()     *
* if execution from a restart file is specified.  Finally,              *
* getFromInput() is called to read values from the given input          *
* database (potentially overriding those found in the restart file).    *
*                                                                       *
*************************************************************************
*/

ConvDiff::ConvDiff(
   const string& object_name,
   tbox::Pointer<tbox::Database> input_db,
   tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geom) 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!object_name.empty());
   assert(!input_db.isNull());
   assert(!grid_geom.isNull());
#endif

   d_object_name = object_name;
   tbox::RestartManager::getManager()->registerRestartItem(d_object_name, this);

   d_grid_geometry = grid_geom;

   int k;

   /*
    *hier::Variable<NDIM> quantities that define state of convection diffusion
    * problem.
    *
    *     dU/dt + alpha dU/dx = mu d^2U/dx^2 + gamma
    *
    *     U     = primitive variable(s)
    *     F(U)  = function evaluation
    *     alpha = convection coefficient
    *     mu    = diffusion coefficient
    *     gamma = source coefficient
    */
   d_primitive_vars = new pdat::CellVariable<NDIM,double>("primitive_vars",1);
   d_function_eval  = new pdat::CellVariable<NDIM,double>("function_eval",1);
   d_diffusion_coeff = 1.;
   for (k = 0; k < NDIM; k++) d_convection_coeff[k] = 0.;
   d_source_coeff = 0.;

   //Set ghost width 
   d_nghosts = 1;
   d_zero_ghosts = 0;

   // Physics parameters
   d_cfl = 0.9;
   for (k = 0; k < NEQU; k++) d_tolerance[k] = 0.;

   /*
    * Defaults for problem type and initial data.  Set initial 
    * data to NaNs so we make sure input has set it to appropriate
    * problem.
    */
   d_data_problem = tbox::IEEE::getINT_MAX();

   // SPHERE problem...
   d_radius = tbox::IEEE::getSignalingNaN();
   tbox::IEEE::initializeArrayToSignalingNaN(d_center, NDIM);
   tbox::IEEE::initializeArrayToSignalingNaN(d_val_inside, NEQU);
   tbox::IEEE::initializeArrayToSignalingNaN(d_val_inside, NEQU);

   /*
    * Boundary condition initialization.
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

   d_bdry_edge_val.resizeArray(NUM_2D_EDGES);
   tbox::IEEE::initializeArrayToSignalingNaN(d_bdry_edge_val);
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

   d_bdry_face_val.resizeArray(NUM_3D_FACES);
   tbox::IEEE::initializeArrayToSignalingNaN(d_bdry_face_val);
#endif

   /*
    * Initialize object with data read from given input/restart databases.
    */
   bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();
   if (is_from_restart){
      getFromRestart();
   }
   getFromInput(input_db,is_from_restart);

   /*
    * Set problem data to values read from input/restart.
    */

   if (d_data_problem == "SPHERE") {
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

}

/*************************************************/
ConvDiff::~ConvDiff() 
{
} 

/*
*************************************************************************
*                                                                       *
* Register variables with the MOL integrator.  Since the integrator     *
* that manages storage for these quantities, this is how it finds       *
* out about them.                                                       *
*                                                                       *
* Register variables also with the vizamrai data writer so they may be  *
* written out.                                                          *
*                                                                       *
*************************************************************************
*/
void ConvDiff::registerModelVariables(algs::MethodOfLinesIntegrator<NDIM>* integrator) 
{
   /*
    * Only two types of variables are used by the integrator - SOLN
    * and RHS.  The primitive variables are of type SOLN while 
    * the function evaluation is of type RHS.  
    *
    * SOLN needs two contexts - current (without ghosts) and 
    * scratch (with ghosts). Current is maintained between timesteps, 
    * while scratch is created and destroyed within a timestep.  
    */

   integrator->registerVariable(d_primitive_vars,d_nghosts,
                                algs::MethodOfLinesIntegrator<NDIM>::SOLN,
                                d_grid_geometry,
                                "CONSERVATIVE_COARSEN", 
                                "LINEAR_REFINE");

   /*
    * RHS needs only one context - scratch (with ghosts).  It is used
    * to store communicated ghost information within the timestep.
    * The function evaluation is not communicated across levels, 
    * it is recomputed using the interpolated primitive variables.
    * Hence, we don't need to define a coarsen or refine operator.
    */
   integrator->registerVariable(d_function_eval, d_nghosts,
                                algs::MethodOfLinesIntegrator<NDIM>::RHS,
                                d_grid_geometry,
                                "NO_COARSEN", 
                                "NO_REFINE");

   /*
    * Loop over primitive variables and register each with the 
    * vizamrai data writer.
    */
   hier::VariableDatabase<NDIM>* vardb = hier::VariableDatabase<NDIM>::getDatabase();
   
   int prim_var_id = vardb->mapVariableAndContextToIndex(
                     d_primitive_vars, getInteriorContext());

   string dump_name = "Primitive Var #";
   const int size = dump_name.length() + 16;
   char *buffer = new char[size];

   for (int n = 0; n < NEQU; n++) {
      sprintf(buffer, "%s%01d", dump_name.c_str(), n);
      string variable_name(buffer);
      if (!(d_vizamrai_writer.isNull())) {
         d_vizamrai_writer->
            registerPlotScalar(variable_name, prim_var_id, n);
      }
#ifdef HAVE_HDF5
      if (!(d_visit_writer.isNull())) {
         d_visit_writer->
            registerPlotQuantity(variable_name, "SCALAR",
                                 prim_var_id, n);
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

   delete [] buffer;
      
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
void ConvDiff::initializeDataOnPatch(hier::Patch<NDIM>& patch, 
                                    const double time,
                                    const bool initial_time) const
{
   (void) time;

   if (initial_time) {

     const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > patch_geom = 
                                 patch.getPatchGeometry();

     const double* dx  = patch_geom->getDx();
     const double* xlo = patch_geom->getXLower();
     const double* xhi = patch_geom->getXUpper();

     tbox::Pointer< pdat::CellData<NDIM,double> > primitive_vars = 
        patch.getPatchData(d_primitive_vars, getInteriorContext());

#ifdef DEBUG_CHECK_ASSERTIONS
     assert(!primitive_vars.isNull());
#endif
     hier::IntVector<NDIM> ghost_cells = primitive_vars->getGhostCellWidth();

     hier::Box<NDIM> pbox = patch.getBox();
     tbox::Pointer< pdat::CellData<NDIM,double> > temp_var;
     temp_var = new pdat::CellData<NDIM,double>(pbox, 1, ghost_cells);

   
     const hier::Index<NDIM> ifirst = patch.getBox().lower();
     const hier::Index<NDIM> ilast  = patch.getBox().upper();

     FORT_SPHERE_INIT(dx, xlo, xhi,
         ifirst(0),ilast(0), ifirst(1),ilast(1),
#if (NDIM==3)
         ifirst(2),ilast(2), 
#endif
         ghost_cells(0), ghost_cells(1),
#if (NDIM==3)
         ghost_cells(2),
#endif
         primitive_vars->getPointer(),
         d_val_inside,
         d_val_outside,
         d_center, d_radius,
         NEQU); 
 
      tbox::plog << "Level:" << patch.getPatchLevelNumber() << "\n" << endl;
      tbox::plog << "Patch:" << endl;   
      primitive_vars->print(patch.getBox(),tbox::plog);
   }
}

/*
*************************************************************************
*                                                                       *
* Compute stable time increment for patch.  Return this value.          *
* (See Hirsch, Vol 1, pp 448 for description of stability analysis)     *
*                                                                       *
*************************************************************************
*/
double ConvDiff::computeStableDtOnPatch( hier::Patch<NDIM>& patch,
                                         const double time) const
{
   (void) time;

   const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > patch_geom = patch.getPatchGeometry();
   const double* dx  = patch_geom->getDx();

   const hier::Index<NDIM> ifirst=patch.getBox().lower();
   const hier::Index<NDIM> ilast =patch.getBox().upper();

   tbox::Pointer< pdat::CellData<NDIM,double> > primitive_vas   = 
      patch.getPatchData(d_primitive_vars, getInteriorContext());

   double stabdt;

//  Use the condition defined on pg. 448 of Hirsch, Vol. 1.
//     for system du/dt = beta*u_xx
//           0 <= beta*dt/(dx**2) <= 1   
//     assume Beta=d_cfl for this case.

#if (NDIM==2)
   stabdt = d_cfl*((*dx)*(*dx) );
#endif
#if (NDIM==3)
   stabdt = d_cfl*((*dx)*(*dx)*(*dx) );
#endif


//   Alternatively, one could use a fortran function here if you want
//   something more complex. 
//
//   FORT_STABLE_DT(stabdt
//                  ifirst(0),ilast(0),
//                  ifirst(1),ilast(1),
//                  ifirst(2),ilast(2),
//                  d_params,
//                  density->getPointer(),
//                  velocity->getPointer(),
//                  pressure->getPointer(),
//                  stabdt);
   return stabdt;
  
}

/*
*************************************************************************
*                                                                       *
* Perform a single Runge-Kutta sub-iteration using the passed-in        *
* alpha.                                                                *
*                                                                       *
*************************************************************************
*/
void ConvDiff::singleStep(hier::Patch<NDIM>& patch, 
                         const double dt, 
                         const double alpha_1,
                         const double alpha_2,
                         const double beta) const
{

  tbox::Pointer< pdat::CellData<NDIM,double> > prim_var_updated  = 
      patch.getPatchData(d_primitive_vars, getInteriorWithGhostsContext());
  
  tbox::Pointer< pdat::CellData<NDIM,double> > prim_var_fixed  = 
      patch.getPatchData(d_primitive_vars, getInteriorContext());
  
  tbox::Pointer< pdat::CellData<NDIM,double> > function_eval     = 
      patch.getPatchData(d_function_eval, getInteriorContext());

  const hier::Index<NDIM> ifirst = patch.getBox().lower();
  const hier::Index<NDIM> ilast  = patch.getBox().upper();

  const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > patch_geom = 
     patch.getPatchGeometry();
  const double* dx  = patch_geom->getDx();

//      tbox::plog << "----primitive_var_current" << endl;
//      prim_var_current->print(prim_var_current->getGhostBox(),tbox::plog);
//      tbox::plog << "----primitive_var_scratch" << endl;
//      prim_var_scratch->print(prim_var_scratch->getGhostBox(),tbox::plog);
//
// Evaluate Right hand side F(prim_var_scratch) 
//
   FORT_COMP_RHS(ifirst(0),ilast(0),ifirst(1),ilast(1),
#if (NDIM==3)
                ifirst(2),ilast(2),
#endif
                d_nghosts(0),d_nghosts(1),
#if (NDIM==3)
                d_nghosts(2),
#endif
                dx,
                d_convection_coeff,
                d_diffusion_coeff,
                d_source_coeff,
                prim_var_updated->getPointer(), 
                function_eval->getPointer(),
                NEQU ); 

//    tbox::plog << "Function Evaluation" << endl;
//    function_eval->print(function_eval->getBox());
//
// Take RK step
//
   FORT_RK_STEP(ifirst(0),ilast(0),ifirst(1),ilast(1),
#if (NDIM==3)
                ifirst(2),ilast(2),
#endif
                d_nghosts(0),d_nghosts(1),
#if (NDIM==3)
                d_nghosts(2),
#endif
                dt, alpha_1, alpha_2, beta,
                d_convection_coeff,
                d_diffusion_coeff,
                d_source_coeff, 
                prim_var_updated->getPointer(),
                prim_var_fixed->getPointer(),
                function_eval->getPointer(),
                NEQU );
//        tbox::plog << "----prim_var_scratch after RK step" << endl;
//        prim_var_scratch->print(prim_var_scratch->getGhostBox(),tbox::plog);


}

/*
*************************************************************************
*                                                                       *
*  Cell tagging routine - tag cells that require refinement based on    *
*  a provided condition.                                                *
*                                                                       *
*************************************************************************
*/
void ConvDiff::tagGradientDetectorCells(
   hier::Patch<NDIM>& patch, 
   const double regrid_time, 
   const bool initial_error, 
   const int tag_index,
   const bool uses_richardson_extrapolation_too)
{
   (void) regrid_time; 
   (void) initial_error; 
   (void) uses_richardson_extrapolation_too; 

   tbox::Pointer< pdat::CellData<NDIM,int> > tags = patch.getPatchData(tag_index);
   tbox::Pointer< pdat::CellData<NDIM,double> > primitive_vars   = 
      patch.getPatchData(d_primitive_vars, getInteriorContext());

   const hier::Index<NDIM> ifirst = patch.getBox().lower();
   const hier::Index<NDIM> ilast  = patch.getBox().upper();
   
   FORT_TAG_CELLS(ifirst(0),ilast(0),ifirst(1),ilast(1),           
#if (NDIM==3)
                  ifirst(2),ilast(2),
#endif
                  tags->getPointer(),
                  primitive_vars->getPointer(),
                  true,
                  d_tolerance,
                  NEQU );

//      tbox::plog << "/n----Tagged Error Cells" << endl;
//      tags->print(tags->getBox(),tbox::plog);
//

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
void ConvDiff::setPhysicalBoundaryConditions(
   hier::Patch<NDIM>& patch,
   const double fill_time,
   const hier::IntVector<NDIM>& ghost_width_to_fill)
{
   (void) fill_time;

   tbox::Pointer< pdat::CellData<NDIM,double> > primitive_vars =
      patch.getPatchData(d_primitive_vars, getInteriorWithGhostsContext());

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!primitive_vars.isNull());
#endif
   hier::IntVector<NDIM> ghost_cells = primitive_vars->getGhostCellWidth();
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(primitive_vars->getGhostCellWidth() == d_nghosts);
#endif

#if (NDIM == 2) 

   /*
    * Set boundary conditions for cells corresponding to patch edges.
    */
   appu::CartesianBoundaryUtilities2::
      fillEdgeBoundaryData("primitive_vars", primitive_vars,
                           patch,
                           ghost_width_to_fill,
                           d_scalar_bdry_edge_conds,
                           d_bdry_edge_val);
 
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
      fillNodeBoundaryData("primitive_vars", primitive_vars,
                           patch,
                           ghost_width_to_fill,
                           d_scalar_bdry_node_conds,
                           d_bdry_edge_val);

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
      fillFaceBoundaryData("primitive_vars", primitive_vars,
                           patch,
                           ghost_width_to_fill,
                           d_scalar_bdry_face_conds,
                           d_bdry_face_val);
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
      fillEdgeBoundaryData("primitive_vars", primitive_vars,
                           patch,
                           ghost_width_to_fill,
                           d_scalar_bdry_edge_conds,
                           d_bdry_face_val);
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
      fillNodeBoundaryData("primitive_vars", primitive_vars,
                           patch,
                           ghost_width_to_fill,
                           d_scalar_bdry_node_conds,
                           d_bdry_face_val);
#ifdef DEBUG_CHECK_ASSERTIONS
#if CHECK_BDRY_DATA
   checkBoundaryData(NODE3D_BDRY_TYPE, patch, ghost_width_to_fill,
                     d_scalar_bdry_node_conds);
#endif
#endif

#endif // NDIM == 3

//    tbox::plog << "----Boundary Conditions "  << endl;
//    primitive_vars->print(primitive_vars->getGhostBox());

}


/*
*************************************************************************
*                                                                       *
* Register Vizamrai data writer to write data to plot files that may    *
* be postprocessed by the Vizamrai tool.                                *
*                                                                       *
*************************************************************************
*/

void ConvDiff::registerVizamraiDataWriter(
   tbox::Pointer<appu::CartesianVizamraiDataWriter<NDIM> > viz_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(viz_writer.isNull()));
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
void ConvDiff::registerVisItDataWriter(
   tbox::Pointer<appu::VisItDataWriter<NDIM> > viz_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(viz_writer.isNull()));
#endif
   d_visit_writer = viz_writer;
}
#endif

/*
*************************************************************************
*                                                                       *
* Prints class data - writes out info in class if exception is thrown   *
*                                                                       *
*************************************************************************
*/

void ConvDiff::printClassData(ostream &os) const
{
   fflush(stdout);
   int j;

   os << "ptr ConvDiff = " << (ConvDiff*) this << endl;
   os << "ptr grid geometry = "
      << (geom::CartesianGridGeometry<NDIM>*) d_grid_geometry << endl;

   os << "Coefficients..." << endl;
   for (j = 0; j < NDIM; j++) os << "d_convection_coeff[" << j << "] = " 
       << d_convection_coeff[j] << endl;
   os << "d_diffusion_coeff = " << d_diffusion_coeff << endl;
   os << "d_source_coeff = " << d_source_coeff << endl;


   os << "Problem description and initial data..." << endl;
   os << "   d_data_problem = " << d_data_problem << endl;

   os << "       d_radius = " << d_radius << endl;
   os << "       d_center = " ;
      for (j=0;j<NDIM;j++) os << d_center[j]<<" ";
   os << endl;
   os << "       d_val_inside = " ;
      for (j=0;j<NEQU;j++) os << d_val_inside[j]<<" ";
   os << endl;
   os << "       d_val_outside = " ;
      for (j=0;j<NEQU;j++) os << d_val_outside[j]<<" ";
   os << endl;


   os << "Boundary Condition data..." << endl;
#if (NDIM == 2) 
   for (j = 0; j < d_scalar_bdry_edge_conds.getSize(); j++) {
      os << "       d_scalar_bdry_edge_conds[" << j << "] = "
         << d_scalar_bdry_edge_conds[j] << endl;
      if (d_scalar_bdry_edge_conds[j] == DIRICHLET_BC) {
         os << "         d_bdry_edge_val[" << j << "] = "
            << d_bdry_edge_val[j] << endl;
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
         os << "         d_bdry_face_val[" << j << "] = "
            << d_bdry_face_val[j] << endl;
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

   os << "         d_nghosts = " << d_nghosts << endl;
   os << "         d_zero_ghosts = " << d_zero_ghosts << endl;
   os << "         d_cfl = " << d_cfl << endl;

}



/*
*************************************************************************
*                                                                       *
*                                                                       *
*************************************************************************
*/
void ConvDiff::getFromInput(tbox::Pointer<tbox::Database> db,
                            bool is_from_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!db.isNull());
#endif

   if (db->keyExists("convection_coeff")){
       db->getDoubleArray("convection_coeff",
                          d_convection_coeff, NDIM);
   } else {
       TBOX_ERROR(d_object_name << ":  "
          << "Key data `convection_coeff' not found in input.");
   }
   if (db->keyExists("diffusion_coeff")){
       d_diffusion_coeff = db->getDouble("diffusion_coeff");
   } else {
       TBOX_ERROR(d_object_name << ":  "
          << "Key data `diffusion_coeff' not found in input.");
   }
   if (db->keyExists("source_coeff")){
       d_source_coeff = db->getDouble("source_coeff");
   } else {
       TBOX_ERROR(d_object_name << ":  "
          << "Key data `source_coeff' not found in input.");
   }

   d_cfl = db->getDoubleWithDefault("cfl", d_cfl);

   if (db->keyExists("cell_tagging_tolerance")){
       db->getDoubleArray("cell_tagging_tolerance",
                          d_tolerance, NEQU);
   } else {
       TBOX_ERROR(d_object_name << ":  "
          << "Key data `cell_tagging_tolerance' not found in input.");
   }

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
         if (init_data_db->keyExists("val_inside")) {
            init_data_db->getDoubleArray("val_inside", d_val_inside, NEQU);
         } else {
            TBOX_ERROR(d_object_name << ": "
               << "val_inside' input required for "
               << "SPHERE problem." << endl);
         }
         if (init_data_db->keyExists("val_outside")) {
            init_data_db->getDoubleArray("val_outside", d_val_outside, NEQU);
         } else {
            TBOX_ERROR(d_object_name << ": "
               << "`val_outside' input required for "
               << "SPHERE problem." << endl);
         }

         found_problem_data = true;

      }

      if (!found_problem_data) {
         TBOX_ERROR(d_object_name << ": "
            << "Bad data in `Initial_data' database." << endl);
      }
   
      hier::IntVector<NDIM> periodic = d_grid_geometry->getPeriodicShift();
      int num_per_dirs = 0;
      for (int id = 0; id < NDIM; id++) {
         if (periodic(id)) num_per_dirs++;
      }

      if (db->keyExists("Boundary_data")) {
         tbox::Pointer<tbox::Database> boundary_db = db->getDatabase("Boundary_data");
         
#if (NDIM == 2)
         appu::CartesianBoundaryUtilities2::readBoundaryInput(this,
                                                        boundary_db,
                                                        d_scalar_bdry_edge_conds,
                                                        d_scalar_bdry_node_conds,
                                                        periodic);
#endif
#if (NDIM == 3)
         appu::CartesianBoundaryUtilities3::readBoundaryInput(this,
                                                        boundary_db,
                                                        d_scalar_bdry_face_conds,
                                                        d_scalar_bdry_edge_conds,
                                                        d_scalar_bdry_node_conds,
                                                        periodic);
#endif
         
      } else {
         TBOX_WARNING(d_object_name << ": "
                      << "Key data `Boundary_data' not found in input. "
                      << "Using default FLOW boundary conditions." << endl);
      }
   }
}


/*
*************************************************************************
*                                                                       *
* Routines to put/get data members to/from restart database.            *
*                                                                       *
*************************************************************************
*/

void ConvDiff::putToDatabase( tbox::Pointer<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!db.isNull());
#endif

   db->putInteger("CONV_DIFF_VERSION",CONV_DIFF_VERSION);

   db->putDouble("d_diffusion_coeff",d_diffusion_coeff);
   db->putDoubleArray("d_convection_coeff",d_convection_coeff,NDIM);
   db->putDouble("d_source_coeff",d_source_coeff);
   db->putIntegerArray("d_nghosts", (int*)d_nghosts, NDIM);

   db->putString("d_data_problem", d_data_problem);
   
   db->putDouble("d_radius", d_radius);
   db->putDoubleArray("d_center", d_center, NDIM);
   db->putDoubleArray("d_val_inside", d_val_inside, NEQU);
   db->putDoubleArray("d_val_outside", d_val_outside, NEQU);

   db->putDouble("d_cfl",d_cfl);

   db->putIntegerArray("d_scalar_bdry_edge_conds", d_scalar_bdry_edge_conds);
   db->putIntegerArray("d_scalar_bdry_node_conds", d_scalar_bdry_node_conds);

#if (NDIM == 2) 
   db->putDoubleArray("d_bdry_edge_val", d_bdry_edge_val);
#endif
#if (NDIM == 3)
   db->putIntegerArray("d_scalar_bdry_face_conds", d_scalar_bdry_face_conds);
   db->putDoubleArray("d_bdry_face_val", d_bdry_face_val);
#endif

}

/*
*************************************************************************
*                                                                       *
*                                                                       *
*************************************************************************
*/
void ConvDiff::getFromRestart()
{

   tbox::Pointer<tbox::Database> root_db = 
      tbox::RestartManager::getManager()->getRootDatabase();

   tbox::Pointer<tbox::Database> db;
   if ( root_db->isDatabase(d_object_name) ) {
      db = root_db->getDatabase(d_object_name);
   } else {
      TBOX_ERROR("Restart database corresponding to "
              << d_object_name << " not found in the restart file.");
   } 
   
   int ver  = db->getInteger("CONV_DIFF_VERSION");
   if (ver != CONV_DIFF_VERSION) {
      TBOX_ERROR(d_object_name << ":  "
              << "Restart file version different than class version.");
   }

   d_diffusion_coeff = db->getDouble("d_diffusion_coeff");
   db->getDoubleArray("d_convection_coeff", d_convection_coeff, NDIM);
   d_source_coeff = db->getDouble("d_source_coeff");
   db->getIntegerArray("d_nghosts", d_nghosts, NDIM);

   d_data_problem = db->getInteger("d_data_problem");

   d_radius = db->getDouble("d_radius");
   db->getDoubleArray("d_center", d_center, NDIM);
   db->getDoubleArray("d_val_inside", d_val_inside, NEQU);
   db->getDoubleArray("d_val_outside", d_val_outside, NEQU);

   d_cfl = db->getDouble("d_cfl");

   d_scalar_bdry_edge_conds = db->getIntegerArray("d_scalar_bdry_edge_conds");
   d_scalar_bdry_node_conds = db->getIntegerArray("d_scalar_bdry_node_conds");

#if (NDIM == 2) 
   d_bdry_edge_val = db->getDoubleArray("d_bdry_edge_val");
#endif
#if (NDIM == 3)
   d_scalar_bdry_face_conds = db->getIntegerArray("d_scalar_bdry_face_conds");

   d_bdry_face_val = db->getDoubleArray("d_bdry_face_val");
#endif


}

/*
*************************************************************************
*                                                                       *
* Routines to read boundary data from input database.                   *
*                                                                       *
*************************************************************************
*/

void ConvDiff::readDirichletBoundaryDataEntry(tbox::Pointer<tbox::Database> db,
                                           string& db_name,
                                           int bdry_location_index)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!db.isNull());
   assert(!db_name.empty());
#endif
#if (NDIM == 2)
   readStateDataEntry(db,
                      db_name,
                      bdry_location_index,
                      d_bdry_edge_val);
#endif
#if (NDIM == 3)
   readStateDataEntry(db,
                      db_name,
                      bdry_location_index,
                      d_bdry_face_val);
#endif
}

void ConvDiff::readNeumannBoundaryDataEntry(tbox::Pointer<tbox::Database> db,
                                           string& db_name,
                                           int bdry_location_index)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!db.isNull());
   assert(!db_name.empty());
#endif
#if (NDIM == 2)
   readStateDataEntry(db,
                      db_name,
                      bdry_location_index,
                      d_bdry_edge_val);
#endif
#if (NDIM == 3)
   readStateDataEntry(db,
                      db_name,
                      bdry_location_index,
                      d_bdry_face_val);
#endif
}

void ConvDiff::readStateDataEntry(tbox::Pointer<tbox::Database> db,
                               const string& db_name,
                               int array_indx,
                               tbox::Array<double>& val)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!db.isNull());
   assert(!db_name.empty());
   assert(array_indx >= 0);
   assert(val.getSize() > array_indx);
#endif

   if (db->keyExists("val")) {
      val[array_indx] = db->getDouble("val");
   } else {
      TBOX_ERROR(d_object_name << ": "
         << "`val' entry missing from " << db_name
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

void ConvDiff::checkBoundaryData(int btype,
                              const hier::Patch<NDIM>& patch,
                              const hier::IntVector<NDIM>& ghost_width_to_check,
                              const tbox::Array<int>& scalar_bconds) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
#if (NDIM == 2) 
   assert(btype == EDGE2D_BDRY_TYPE || 
          btype == NODE2D_BDRY_TYPE);
#endif
#if (NDIM == 3)
   assert(btype == FACE3D_BDRY_TYPE ||
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
      assert(bbox.getBoundaryType() == btype);
#endif
      int bloc = bbox.getLocationIndex();

      int bscalarcase, refbdryloc;
#if (NDIM == 2)
      if (btype == EDGE2D_BDRY_TYPE) {
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(scalar_bconds.getSize() == NUM_2D_EDGES);
#endif
         bscalarcase = scalar_bconds[bloc];
         refbdryloc = bloc;
      } else { // btype == NODE2D_BDRY_TYPE
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(scalar_bconds.getSize() == NUM_2D_NODES);
#endif
         bscalarcase = scalar_bconds[bloc];
         refbdryloc = d_node_bdry_edge[bloc];
      }
#endif
#if (NDIM == 3)
      if (btype == FACE3D_BDRY_TYPE) {
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(scalar_bconds.getSize() == NUM_3D_FACES);
#endif
         bscalarcase = scalar_bconds[bloc];
         refbdryloc = bloc;
      } else if (btype == EDGE3D_BDRY_TYPE) {
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(scalar_bconds.getSize() == NUM_3D_EDGES);
#endif
         bscalarcase = scalar_bconds[bloc];
         refbdryloc = d_edge_bdry_face[bloc];
      } else { // btype == NODE3D_BDRY_TYPE
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(scalar_bconds.getSize() == NUM_3D_NODES);
#endif
         bscalarcase = scalar_bconds[bloc];
         refbdryloc = d_node_bdry_face[bloc];
      }
#endif

      int num_bad_values = 0;

#if (NDIM == 2)
      num_bad_values = 
      appu::CartesianBoundaryUtilities2::checkBdryData(
         d_primitive_vars->getName(),
         patch,
         vdb->mapVariableAndContextToIndex(d_primitive_vars, 
                                           getInteriorWithGhostsContext()), 
         0, 
         ghost_width_to_check,
         bbox,
         bscalarcase,
         d_bdry_edge_val[refbdryloc]);
#endif
#if (NDIM == 3)
      num_bad_values = 
      appu::CartesianBoundaryUtilities3::checkBdryData(
         d_primitive_vars->getName(),
         patch,
         vdb->mapVariableAndContextToIndex(d_primitive_vars, 
                                           getInteriorWithGhostsContext()),
         0, 
         ghost_width_to_check,
         bbox,
         bscalarcase,
         d_bdry_face_val[refbdryloc]);
#endif
#if (TESTING == 1)
      if (num_bad_values > 0) {
         tbox::perr << "\nConvDiff Boundary Test FAILED: \n" 
           << "     " << num_bad_values << " bad VAL values found for\n"
           << "     boundary type " << btype << " at location " << bloc << endl;
      }
#endif

   }

}

