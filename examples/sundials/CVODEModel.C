//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/examples/sundials/CVODEModel.C $
// Package:     SAMRAI mesh
// Copyright:   (c) 1997-2002 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2043 $
// Modified:    $LastChangedDate: 2008-03-12 09:14:32 -0700 (Wed, 12 Mar 2008) $
// Description: 
//

#include "CVODEModel.h"

#if defined(HAVE_SUNDIALS) && defined(HAVE_HYPRE)


#include "tbox/Array.h"
#include "BoundaryBox.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellIndex.h"
#include "CellIterator.h"
#include "CoarsenAlgorithm.h"
#include "CoarsenOperator.h"
#include "CoarsenSchedule.h"
#include "FaceData.h"
#include "HierarchyDataOpsReal.h"
#include "HierarchyCellDataOpsReal.h"
#include "Index.h"
#include "PatchCellDataOpsReal.h"
#include "OuterfaceData.h"
#include "Patch.h"
#include "PatchData.h"
// #include "PoissonHYPRESolver<NDIM>.h"
#include "RefineOperator.h"
#include "tbox/RestartManager.h"
#include "SAMRAIVectorReal.h"
#include "SideData.h"
#include "tbox/MathUtilities.h"
#include "tbox/Utilities.h"
#include "VariableDatabase.h"


//integer constants for boundary conditions
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

// Define class version number
#define CVODE_MODEL_VERSION (1)

// This is used in the cell tagging routine.
#ifndef TRUE
#define TRUE (1)
#endif
#ifndef FALSE
#define FALSE (0)
#endif


extern "C" {
#if (NDIM == 2)
void comprhs_(const int&, const int&,
              const int&, const int&,
              const int&, const int&,
              const double*,
              const double*,
              const double*, const double*,
              double*);
#ifdef USE_FAC_PRECONDITIONER
void setneufluxvalues_(const int&, const int&,
                       const int&, const int&,
                       const int*, const double*,
                       int*, int*, int*, int*,
                       double*, double*, double*, double*);
#endif
#endif
#if (NDIM == 3)
void comprhs_(const int&, const int&, 
              const int&, const int&, 
              const int&, const int&,
              const int&, const int&, const int&,
              const double*,
              const double*,
              const double*, const double*, const double*,
              double*);
#ifdef USE_FAC_PRECONDITIONER
void setneufluxvalues_(const int&, const int&,
                       const int&, const int&,
                       const int&, const int&,
                       const int*, const double*,
                       int*, int*, int*, int*, int*, int*,
                       double*, double*, double*, double*, double*, double*);
#endif
#endif
}


/*************************************************************************
 *
 * Constructor and Destructor for CVODEModel class.
 *
 ************************************************************************/

CVODEModel::CVODEModel(
   const string& object_name,
   Pointer<Database> input_db,
   Pointer<CartesianGridGeometry<NDIM> > grid_geom)
   : d_FAC_solver(object_name+":FAC solver")
{

   d_object_name = object_name;

   d_grid_geometry = grid_geom;
   
   /* 
    * set up variables and contexts 
    */
   VariableDatabase<NDIM>* variable_db = VariableDatabase<NDIM>::getDatabase();

   d_soln_var = new CellVariable<NDIM,double>("soln",1);
   
   d_cur_cxt = variable_db->getContext("CURRENT");
   d_scr_cxt = variable_db->getContext("SCRATCH");
 
   d_soln_cur_id = variable_db->registerVariableAndContext(d_soln_var,
                                                           d_cur_cxt,
                                                           IntVector<NDIM>(0));
   d_soln_scr_id = variable_db->registerVariableAndContext(d_soln_var,
                                                           d_scr_cxt,
                                                           IntVector<NDIM>(1));
#ifdef USE_FAC_PRECONDITIONER
   d_diff_var = new SideVariable<NDIM,double>("diffusion",1);
   
   d_diff_id = variable_db->registerVariableAndContext(d_diff_var,
                                                       d_cur_cxt,
                                                       IntVector<NDIM>(0));

   /*
    * Set default values for preconditioner.
    */
   d_max_fac_its   = 1;
   d_fac_tol       = 1.0e-40;
   d_max_hypre_its = 1;
   d_hypre_tol     = 1.0e-40;
   d_use_neumann_bcs      = false;
   // d_FAC_solver    = (PoissonHierarchySolver<NDIM>*)NULL;
   // d_FAC_solver_allocated = false;
   // d_level_solver  = (PoissonLevelStrategy<NDIM>*)NULL;
   // d_level_solver_allocated = false;

   d_current_soln_time = 0.;
   
#endif 
 
   /*
    * Print solver data.
    */
   d_print_solver_info = false; 

   /*
    * Counters.
    */
   d_number_rhs_eval = 0;
   d_number_precond_setup = 0;
   d_number_precond_solve = 0;
   
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
   tbox::MathUtilities<double>::setArrayToSignalingNaN(d_bdry_edge_val);
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
   tbox::MathUtilities<double>::setArrayToSignalingNaN(d_bdry_face_val);
#endif

   /*
    * Initialize object with data read from given input/restart databases.
    */
   bool is_from_restart = RestartManager::getManager()->isFromRestart();
   if (is_from_restart){
      getFromRestart();
   }
   getFromInput(input_db,is_from_restart);

#ifdef USE_FAC_PRECONDITIONER
   /*
    * Construct outerface variable to hold boundary flags and Neumann fluxes.
    */
   if (d_use_neumann_bcs) {
      d_flag_var = new OuterfaceVariable<NDIM,int>("bdryflag",1);
      d_flag_id  = variable_db->registerVariableAndContext(d_flag_var,
                                                           d_cur_cxt,
                                                           IntVector<NDIM>(0));
      d_neuf_var = new OuterfaceVariable<NDIM,double>("neuflux",1);
      d_neuf_id  = variable_db->registerVariableAndContext(d_neuf_var,
                                                           d_cur_cxt,
                                                           IntVector<NDIM>(0));
   } else {
      d_flag_id = -1;
      d_neuf_id = -1;
   }

   /*
    * Set boundary types for FAC preconditioner.
    *  bdry_types holds a flag where 0 = dirichlet, 1 = neumann
    */
#if (NDIM == 2)
   for (int i = 0; i < NUM_2D_EDGES; i++) {
      d_bdry_types[i] = 0;
      if (d_scalar_bdry_edge_conds[i] == DIRICHLET_BC) d_bdry_types[i] = 0;
      if (d_scalar_bdry_edge_conds[i] == NEUMANN_BC)   d_bdry_types[i] = 1;
   }
#endif
#if (NDIM == 3)
   for (int i = 0; i < NUM_3D_FACES; i++) {
      d_bdry_types[i] = 0;
      if (d_scalar_bdry_face_conds[i] == DIRICHLET_BC) d_bdry_types[i] = 0;
      if (d_scalar_bdry_face_conds[i] == NEUMANN_BC)   d_bdry_types[i] = 1;
   }
#endif
#endif
      

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
            CartesianBoundaryUtilities2::getEdgeLocationForNodeBdry(
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
            CartesianBoundaryUtilities3::getFaceLocationForEdgeBdry(
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
            CartesianBoundaryUtilities3::getFaceLocationForNodeBdry(
                                            i, d_scalar_bdry_node_conds[i]);
      }
   }

#endif


}

CVODEModel::~CVODEModel()
{
   Pointer< SAMRAIVectorReal<NDIM,double> > soln_samvect = 
      Sundials_SAMRAIVector<NDIM>::getSAMRAIVector(d_solution_vector);
   Sundials_SAMRAIVector<NDIM>::destroySundialsVector(d_solution_vector);

   soln_samvect->freeVectorComponents();
   soln_samvect.setNull();

   // if (d_level_solver_allocated) delete d_level_solver;
   // d_level_solver_allocated = false;
   // if (d_FAC_solver_allocated) delete d_FAC_solver;
   // d_FAC_solver_allocated = false;
   
}

/*************************************************************************
 *
 * Methods inherited from mesh::StandardTagAndInitStrategy<NDIM>.
 *
 ************************************************************************/

void
CVODEModel::initializeLevelData(
   const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
   const int level_number,
   const double time,
   const bool can_be_refined,
   const bool initial_time,
   const Pointer<BasePatchLevel<NDIM> > old_level,
   const bool allocate_data)
{
   (void) hierarchy;
   (void) level_number;
   (void) time;
   (void) can_be_refined;
   (void) time;
   (void) old_level;
   (void) allocate_data;

   // This method is empty because initialization is taken care of 
   // by the setInitialConditions() method below.  If there is any
   // data that is not managed inside the SAMRAI CVODESolver class
   // but that must be set on the level, do it here.

}

void 
CVODEModel::resetHierarchyConfiguration(
   const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
   const int coarsest_level,
   const int finest_level)
{
   (void) hierarchy;
   (void) coarsest_level;
   (void) finest_level;

   // This method is empty because this example does not exercise the
   // situation when the grid changes, so it effectively is never called. 
   // This is a subject for future work...
}

/*
*************************************************************************
*                                                                       *
* Cell tagging and patch level data initialization routines declared    *
* in the GradientDetectorStrategy interface.  They are used to          *
* construct the hierarchy initially.                                    *
*                                                                       *
*************************************************************************
*/

void
CVODEModel::applyGradientDetector(
   const Pointer< BasePatchHierarchy<NDIM> > hierarchy,
   const int level_number,
   const double time,
   const int tag_index,
   const bool initial_time,
   const bool uses_richardson_extrapolation_too)
{
   (void) uses_richardson_extrapolation_too;

   Pointer< PatchHierarchy<NDIM> > thierarchy = hierarchy;
   Pointer<PatchLevel<NDIM> > level = thierarchy->getPatchLevel(level_number);

   for (PatchLevel<NDIM>::Iterator p(level); p; p++) {
      Pointer<Patch<NDIM> > patch = level->getPatch(p());

      Pointer< CellData<NDIM,int> > tag_data = 
         patch->getPatchData(tag_index);

      // dumb implementation that tags all cells.
      tag_data->fillAll(TRUE);
   } 
}

/*
*************************************************************************
*
* Methods inherited from RefinePatchStrategy.
*
***********************************************************************
*/

void 
CVODEModel::setPhysicalBoundaryConditions(
   Patch<NDIM>& patch,
   const double time,
   const IntVector<NDIM>& ghost_width_to_fill)
{

   (void) time;

   Pointer< CellData<NDIM,double> > soln_data =
      patch.getPatchData(d_soln_scr_id);

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!soln_data.isNull());
#endif
   IntVector<NDIM> ghost_cells = soln_data->getGhostCellWidth();


#if (NDIM == 2) 

   /*
    * Set boundary conditions for cells corresponding to patch edges.
    */
   CartesianBoundaryUtilities2::
      fillEdgeBoundaryData("soln_data", soln_data,
                           patch,
                           ghost_width_to_fill,
                           d_scalar_bdry_edge_conds,
                           d_bdry_edge_val);

   /*
    *  Set boundary conditions for cells corresponding to patch nodes.
    */

   CartesianBoundaryUtilities2::
      fillNodeBoundaryData("soln_data", soln_data,
                           patch,
                           ghost_width_to_fill,
                           d_scalar_bdry_node_conds,
                           d_bdry_edge_val);

#endif // NDIM == 2

#if (NDIM == 3)

   /*
    *  Set boundary conditions for cells corresponding to patch faces.
    */
   CartesianBoundaryUtilities3::
      fillFaceBoundaryData("soln_data", soln_data,
                           patch,
                           ghost_width_to_fill,
                           d_scalar_bdry_face_conds,
                           d_bdry_face_val);

   /*
    *  Set boundary conditions for cells corresponding to patch edges.
    */
   CartesianBoundaryUtilities3::
      fillEdgeBoundaryData("soln_data", soln_data,
                           patch,
                           ghost_width_to_fill,
                           d_scalar_bdry_edge_conds,
                           d_bdry_face_val);

   /*
    *  Set boundary conditions for cells corresponding to patch nodes.
    */
   CartesianBoundaryUtilities3::
      fillNodeBoundaryData("soln_data", soln_data,
                           patch,
                           ghost_width_to_fill,
                           d_scalar_bdry_node_conds,
                           d_bdry_face_val);

#endif // NDIM == 3

//    plog << "----Boundary Conditions "  << endl;
//    soln_data->print(soln_data->getGhostBox());

}


void 
CVODEModel::preprocessRefine(
   Patch<NDIM>& fine,
   const Patch<NDIM>& coarse,
   const Box<NDIM>& fine_box,
   const IntVector<NDIM>& ratio)
{
   (void) fine;
   (void) coarse;
   (void) fine_box;
   (void) ratio;
}

void 
CVODEModel::postprocessRefine(
   Patch<NDIM>& fine,
   const Patch<NDIM>& coarse,
   const Box<NDIM>& fine_box,
   const IntVector<NDIM>& ratio)
{
   (void) fine;
   (void) coarse;
   (void) fine_box;
   (void) ratio;
}




/*************************************************************************
 *
 * Methods inherited from CoarsenPatch<NDIM>Strategy.
 *
 ************************************************************************/

void 
CVODEModel::preprocessCoarsen(
   Patch<NDIM>& coarse,
   const Patch<NDIM>& fine,
   const Box<NDIM>& coarse_box,
   const IntVector<NDIM>& ratio)
{
   (void) coarse;
   (void) fine;
   (void) coarse_box;
   (void) ratio;
}

void 
CVODEModel::postprocessCoarsen(
   Patch<NDIM>& coarse,
   const Patch<NDIM>& fine,
   const Box<NDIM>& coarse_box,
   const IntVector<NDIM>& ratio)
{
   (void) coarse;
   (void) fine;
   (void) coarse_box;
   (void) ratio;
}

/*************************************************************************
 *
 * Methods inherited from CVODEAbstractFunction
 *
 ************************************************************************/

int 
CVODEModel::evaluateRHSFunction(
   double time,
   SundialsAbstractVector* y,
   SundialsAbstractVector* y_dot)
{
   /*
    * Convert Sundials vectors to SAMRAI vectors
    */
   Pointer< SAMRAIVectorReal<NDIM,double> > y_samvect =
      Sundials_SAMRAIVector<NDIM>::getSAMRAIVector(y);
   Pointer< SAMRAIVectorReal<NDIM,double> > y_dot_samvect =
      Sundials_SAMRAIVector<NDIM>::getSAMRAIVector(y_dot);

   Pointer<PatchHierarchy<NDIM> > hierarchy = y_samvect->getPatchHierarchy();

   /*
    * Compute max norm of solution vector.
    */
   Pointer< HierarchyDataOpsReal<NDIM,double> > hierops = 
      new HierarchyCellDataOpsReal<NDIM,double>(hierarchy);
   //double max_norm = hierops->maxNorm(y_samvect->
   //                                   getComponentDescriptorIndex(0));


   if (d_print_solver_info) {
      pout << "\t\tEval RHS: " 
           << "\n   \t\t\ttime = " << time
           << "\n   \t\t\ty_maxnorm = " << y_samvect->maxNorm()
           << endl; 
   }
 
   /*
    * Allocate scratch space and fill ghost cells in the solution vector
    * 1) Create a refine algorithm
    * 2) Register with the algorithm the current & scratch space, along
    *    with a refine operator.
    * 3) Use the refine algorithm to construct a refine schedule
    * 4) Use the refine schedule to fill data on fine level. 
    */
   Pointer<RefineAlgorithm<NDIM> > bdry_fill_alg = new RefineAlgorithm<NDIM>();
   Pointer<RefineOperator<NDIM> > refine_op = d_grid_geometry->
      lookupRefineOperator(d_soln_var, "CONSERVATIVE_LINEAR_REFINE");
   bdry_fill_alg->registerRefine(d_soln_scr_id,  // dest
                                 y_samvect->
                                    getComponentDescriptorIndex(0), // src
                                 d_soln_scr_id, // scratch
                                 refine_op);

   for (int ln = hierarchy->getFinestLevelNumber(); ln >= 0; ln--) {
      Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
      if (!level->checkAllocated(d_soln_scr_id)) {
         level->allocatePatchData(d_soln_scr_id);
      }
      

      // Note:  a pointer to "this" tells the refine schedule to invoke
      // the setPhysicalBCs defined in this class.
      bdry_fill_alg->
         createSchedule(level, ln-1, hierarchy, this)->fillData(time);
   }

   /*
    * Step through the levels and compute rhs 
    */
   for (int ln = hierarchy->getFinestLevelNumber(); ln >= 0; ln--) {
      Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
   
      for (PatchLevel<NDIM>::Iterator ip(level); ip; ip++) {
         Pointer<Patch<NDIM> > patch = level->getPatch(ip());
         Pointer<CartesianPatchGeometry<NDIM> > p_geom = patch->getPatchGeometry();

         Pointer< CellData<NDIM,double> > y =
            patch->getPatchData(d_soln_scr_id);
         Pointer< SideData<NDIM,double> > diff = 
            patch->getPatchData(d_diff_id);
         Pointer< CellData<NDIM,double> > rhs =
            patch->getPatchData(y_dot_samvect->getComponentDescriptorIndex(0));

         const Index<NDIM> ifirst =patch->getBox().lower();
         const Index<NDIM> ilast  =patch->getBox().upper();

         const Pointer<CartesianPatchGeometry<NDIM> > patch_geom =
            patch->getPatchGeometry();
         const double* dx = patch_geom->getDx();

         IntVector<NDIM> ghost_cells = y->getGhostCellWidth();

         /*
          * 1 eqn radiation diffusion
          */ 
         comprhs_(ifirst(0),ilast(0),
                  ifirst(1),ilast(1),
#if (NDIM == 3)
                  ifirst(2),ilast(2),
#endif
                  ghost_cells(0), ghost_cells(1),
#if (NDIM == 3)
                  ghost_cells(2),
#endif
                  dx,
                  y->getPointer(),
                  diff->getPointer(0),
                  diff->getPointer(1),
#if (NDIM == 3)
                  diff->getPointer(2),
#endif
                  rhs->getPointer());

//         plog << "y data" << endl;
//         y->print(y->getGhostBox<NDIM>());
         
//         plog << "rhs data" << endl;
//         rhs->print(rhs->getGhostBox<NDIM>());

      } // loop over patches
   } // loop over levels


   /*
    * Deallocate scratch space.
    */
   for (int ln = hierarchy->getFinestLevelNumber(); ln >= 0; ln--) {
      Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
      level->deallocatePatchData(d_soln_scr_id);
   }

   /*
    * record current time and increment counter for number of RHS 
    * evaluations.
    */
   d_current_soln_time = time;
   d_number_rhs_eval++;

   return 0;
}

/*
 *****************************************************************
 *
 * Set up FAC preconditioner for Jacobian system.  Here we
 * use the FAC hierarchy solver in SAMRAI which automatically sets
 * up the composite grid system and uses hypre as a solver on each
 * level.  
 *
 *****************************************************************
*/

int CVODEModel::CVSpgmrPrecondSet(double t,
				  SundialsAbstractVector* y,
				  SundialsAbstractVector* fy,
				  int jok,
				  int *jcurPtr,
				  double gamma,
				  SundialsAbstractVector* vtemp1,
				  SundialsAbstractVector* vtemp2,
				  SundialsAbstractVector* vtemp3)
{
   (void) jok;
   (void) jcurPtr;
   (void) vtemp1;
   (void) vtemp2;
   (void) vtemp3;
   
#ifdef USE_FAC_PRECONDITIONER
   /*
    * Convert passed-in CVODE vectors into SAMRAI vectors
    */
   Pointer< SAMRAIVectorReal<NDIM,double> > y_samvect =
      Sundials_SAMRAIVector<NDIM>::getSAMRAIVector(y);
   Pointer< SAMRAIVectorReal<NDIM,double> > fy_samvect =
      Sundials_SAMRAIVector<NDIM>::getSAMRAIVector(fy);
   
   Pointer<PatchHierarchy<NDIM> > hierarchy = y_samvect->getPatchHierarchy();
 
   int y_indx = y_samvect->getComponentDescriptorIndex(0);
   
   /*
    * Construct refine algorithm to fill boundaries of solution vector
    */
   RefineAlgorithm<NDIM> fill_soln_vector_bounds;
   Pointer<RefineOperator<NDIM> > refine_op = d_grid_geometry->
      lookupRefineOperator(d_soln_var, "CONSERVATIVE_LINEAR_REFINE");
   fill_soln_vector_bounds.registerRefine(d_soln_scr_id,
                            y_samvect->getComponentDescriptorIndex(0),
                            d_soln_scr_id,
                            refine_op);
   
   /*
    * Construct coarsen algorithm to fill interiors on coarser levels  
    * with solution on finer level.
    */
   CoarsenAlgorithm<NDIM> fill_soln_interior_on_coarser;
   Pointer<CoarsenOperator<NDIM> > coarsen_op = d_grid_geometry->
      lookupCoarsenOperator(d_soln_var, "CONSERVATIVE_COARSEN");
   fill_soln_interior_on_coarser.registerCoarsen(y_indx,
                                                 y_indx,
                                                 coarsen_op);

   /*
    * Step through levels - largest to smallest
    */
   for (int amr_level = hierarchy->getFinestLevelNumber();
        amr_level >= 0;
        amr_level-- ) {
      Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(amr_level);
 
      /*
       * Construct a refine schedule for each level, from the algorithm 
       * constructed above, and fill ghost boundaries of solution vector.
       */
      Pointer<RefineSchedule<NDIM> > fill_soln_vector_bounds_sched =
         fill_soln_vector_bounds.createSchedule(level,
                                  amr_level-1,
                                  hierarchy,
                                  this);
      if (!level->checkAllocated(d_soln_scr_id)) {
         level->allocatePatchData(d_soln_scr_id);
      }
      fill_soln_vector_bounds_sched->fillData(t);

      /*
       * Construct a coarsen schedule for all levels larger than coarsest,  
       * and fill interiors of solution vector on coarser levels using fine 
       * data.
       */
      if (amr_level > 0) {
         Pointer<PatchLevel<NDIM> > coarser_level = 
            hierarchy->getPatchLevel(amr_level-1);
         Pointer<CoarsenSchedule<NDIM> > fill_soln_interior_on_coarser_sched = 
            fill_soln_interior_on_coarser.createSchedule(coarser_level,
                                                         level);
         fill_soln_interior_on_coarser_sched->coarsenData();
      }



      for (PatchLevel<NDIM>::Iterator p(level); p; p++) {
         Pointer<Patch<NDIM> > patch = level->getPatch(p());
 
         const Pointer<CartesianPatchGeometry<NDIM> > patch_geom =
            patch->getPatchGeometry();
 
         const Index<NDIM> ifirst =patch->getBox().lower();
         const Index<NDIM> ilast  =patch->getBox().upper();
 
         Pointer< CellData<NDIM,double> > u = patch->getPatchData(y_indx);
         Pointer< SideData<NDIM, double> > diffusion = 
            patch->getPatchData(d_diff_id);
 
         diffusion->fillAll(1.0);

#ifdef DEBUG_CHECK_ASSERTIONS
         double current_dt = t - d_current_soln_time;
         TBOX_ASSERT(current_dt >= 0.);
#endif


         /*
          * Set Neumann fluxes and flag array (if desired)
          */
         if (d_use_neumann_bcs) {

            Pointer< OuterfaceData<NDIM, int> > flag_data = 
               patch->getPatchData(d_flag_id);
            Pointer< OuterfaceData<NDIM, double> > neuf_data = 
               patch->getPatchData(d_neuf_id);

            /*
             * Outerface data access:
             *    neuf_data->getPointer(axis,face);
             * where axis specifies X, Y, or Z (0,1,2 respectively)
             * and face specifies lower or upper (0,1 respectively)
             */

#if (NDIM == 2)
            setneufluxvalues_(ifirst(0),ilast(0),
                              ifirst(1),ilast(1),
                              d_bdry_types,
                              d_bdry_edge_val.getPointer(),
                              flag_data->getPointer(0,0), // x lower
                              flag_data->getPointer(0,1), // x upper
                              flag_data->getPointer(1,0), // y lower
                              flag_data->getPointer(1,1), // y upper
                              neuf_data->getPointer(0,0), // x lower
                              neuf_data->getPointer(0,1), // x upper
                              neuf_data->getPointer(1,0), // y lower
                              neuf_data->getPointer(1,1)); // y upper
#endif
#if (NDIM == 3)
            setneufluxvalues_(ifirst(0),ilast(0),
                              ifirst(1),ilast(1),
                              ifirst(2),ilast(2),
                              d_bdry_types,
                              d_bdry_face_val.getPointer(),
                              flag_data->getPointer(0,0), // x lower
                              flag_data->getPointer(0,1), // x upper
                              flag_data->getPointer(1,0), // y lower
                              flag_data->getPointer(1,1), // y upper
                              flag_data->getPointer(2,0), // z lower
                              flag_data->getPointer(2,1), // z lower
                              neuf_data->getPointer(0,0), // x lower
                              neuf_data->getPointer(0,1), // x upper
                              neuf_data->getPointer(1,0), // y lower
                              neuf_data->getPointer(1,1), // y upper
                              neuf_data->getPointer(2,0), // z lower
                              neuf_data->getPointer(2,1)); // z upper
#endif
         }
                             
      } // patch loop

      level->deallocatePatchData(d_soln_scr_id);

   } // level loop

   
   /*
    * Setup FAC preconditioner
    */

   /* Construct a Poisson solver for a level.  The use_smg flag
    * is used to select which of HYPRE's linear solver algorithms to use
    * if true, the semicoarsening multigrid algorithm is used, and if
    * false, the ``PF'' multigrid algorithm is used.
    */
   bool use_smg = false;
   
   d_FAC_solver.setUseSMG( use_smg );

   /*
    * Set boundaries.  The "bdry_types" array holds a set of integers
    * where 0 = dirichlet and 1 = neumann boundary conditions.
    */
   if (d_use_neumann_bcs) {
      d_FAC_solver.setBoundaries("Mixed", d_neuf_id, d_flag_id, d_bdry_types);
   } else {
      d_FAC_solver.setBoundaries("Dirichlet");
   }
   
   d_FAC_solver.setCConstant(1.0/gamma);
   d_FAC_solver.setDPatchDataId(d_diff_id);

   /*
    * increment counter for number of precond setup calls
    */
   d_number_precond_setup++;

#endif
   /* 
    * We return 0 or 1 here - 0 if it passes, 1 if it fails.  For now,
    * just be optimistic and return 0. Eventually we should add some 
    * exception handling above to set what this value should be.
    */
   return 0;
}

/*
*************************************************************************
*                                                                       *
* Apply preconditioner where right-hand-side is "r" and "z" is the      *
* solution.   This routine assumes that the preconditioner setup call   *
* has already been invoked.  Return 0 if preconditioner fails;          *
* return 1 otherwise.                                                   *
*                                                                       *
*************************************************************************
*/

int CVODEModel::CVSpgmrPrecondSolve(double t,
                                   SundialsAbstractVector* y,
                                   SundialsAbstractVector* fy,
                                   SundialsAbstractVector* r,
                                   SundialsAbstractVector* z,
                                   double gamma,
                                   double delta,
                                   int lr,
				    SundialsAbstractVector* vtemp)
{
   (void) y;
   (void) fy;
   (void) vtemp;
   (void) gamma;
   (void) delta;
   (void) lr;

#ifdef USE_FAC_PRECONDITIONER

   /*
    * Convert passed-in CVODE vectors into SAMRAI vectors
    */
   Pointer< SAMRAIVectorReal<NDIM,double> > r_samvect =
      Sundials_SAMRAIVector<NDIM>::getSAMRAIVector(r);
   Pointer< SAMRAIVectorReal<NDIM,double> > z_samvect =
      Sundials_SAMRAIVector<NDIM>::getSAMRAIVector(z);

   int ret_val = 0;
 
   Pointer<PatchHierarchy<NDIM> > hierarchy = r_samvect->getPatchHierarchy();
 
   int r_indx = r_samvect->getComponentDescriptorIndex(0);
   int z_indx = z_samvect->getComponentDescriptorIndex(0);
   /******************************************************************
    *
    * We need to supply to the FAC solver a "version" of the z vector
    * that contains ghost cells.  The operations below allocate  
    * on the patches a scratch context of the solution vector z and
    * fill it with z vector data
    *
    *****************************************************************/

   /*
    * Construct a communication schedule which will fill ghosts of
    * soln_scratch with z vector data (z -> soln_scratch).
    */
   RefineAlgorithm<NDIM> fill_z_vector_bounds;
   Pointer<RefineOperator<NDIM> > refine_op = d_grid_geometry->
      lookupRefineOperator(d_soln_var, "CONSERVATIVE_LINEAR_REFINE");
   fill_z_vector_bounds.registerRefine(d_soln_scr_id,
                                       z_indx,
                                       d_soln_scr_id,
                                       refine_op);
   
   /*
    * Set initial guess for z (if applicable) and copy z data into the 
    * solution scratch context.   
    */
   int ln; 
   for (ln = hierarchy->getFinestLevelNumber(); ln >= 0; ln-- ) {
      Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);

      if (!level->checkAllocated(d_soln_scr_id)) {
         level->allocatePatchData(d_soln_scr_id);
      }

      for (PatchLevel<NDIM>::Iterator p(level); p; p++) {
         
         Pointer<Patch<NDIM> > patch = level->getPatch(p());

         const Pointer<CartesianPatchGeometry<NDIM> > patch_geom =
            patch->getPatchGeometry(); 

         Pointer< CellData<NDIM,double> > z_data = 
            patch->getPatchData(z_indx);
         
         /*
          * Set initial guess for z here.
          */
         z_data->fillAll(0.);

         /*
          * Scale RHS by 1/gamma
          */
         PatchCellDataOpsReal<NDIM,double> math_ops;
         Pointer< CellData<NDIM,double> > r_data = patch->getPatchData(r_indx);
         math_ops.scale( r_data, 1.0/gamma, r_data, r_data->getBox() );

         /*
          * Copy interior data from z vector to soln_scratch
          */
         Pointer< CellData<NDIM,double> > z_scr_data =
            patch->getPatchData(d_soln_scr_id);
         z_scr_data->copy(*z_data);
      }
      
      /*
       * Fill ghost boundaries of soln_scratch.  
       * Construct a schedule for each level, from the algorithm
       * constructed above.
       */
      Pointer<RefineSchedule<NDIM> > fill_z_vector_bounds_sched =
         fill_z_vector_bounds.createSchedule(level,
                                             ln-1,
                                             hierarchy,
                                             this);
      
      
      fill_z_vector_bounds_sched->fillData(t);
      
   }
   
   /******************************************************************
    *
    * Apply the FAC solver.  It solves the system Az=r with the
    * format "solveSystem(z, r)". A was constructed in the precondSetup()
    * method.  
    *
    ******************************************************************/

   if (d_print_solver_info) {
      pout << "\t\tBefore FAC Solve (Az=r): " 
           << "\n   \t\t\tz_l2norm = " << z_samvect->L2Norm() 
           << "\n   \t\t\tz_maxnorm = " << z_samvect->maxNorm() 
           << "\n   \t\t\tr_l2norm = " << r_samvect->L2Norm()
           << "\n   \t\t\tr_maxnorm = " << r_samvect->maxNorm() 
           << endl; 
   }
   /*
    * Set paramemters in the FAC solver.  It solves the system Az=r. 
    * Here we supply the max norm of r in order to scale the
    * residual (i.e. residual = Az - r) to properly scale the convergence
    * error.
    */ 

   d_FAC_solver.setMaxCycles(d_max_fac_its);
   d_FAC_solver.setResidualTolerance(d_fac_tol);
   const int coarsest_solve_ln = 0;
   const int finest_solve_ln = 0;
   /*
    * Note: I don't know why we are only solving on level 0 here.
    * When upgrading to the new FAC solver from the old, I noticed
    * that the old solver only solved on level 0.  BTNG.
    */
   bool converge = d_FAC_solver.solveSystem(d_soln_scr_id,
                                            r_indx,
                                            hierarchy,
                                            coarsest_solve_ln,
                                            finest_solve_ln);

   if (d_print_solver_info) {
      double avg_convergence, final_convergence;
      d_FAC_solver.getConvergenceFactors( avg_convergence, final_convergence );
      pout << "   \t\t\tFinal Residual Norm: " 
           << d_FAC_solver.getResidualNorm() << endl;
      pout << "   \t\t\tFinal Convergence Error: " 
           << final_convergence << endl;
      pout << "   \t\t\tFinal Convergence Rate: " 
           << avg_convergence << endl;
   }

   /******************************************************************
    *
    * The FAC solver has computed a solution to z but it is stored
    * in the soln_scratch data space.  Copy it from soln_scratch back
    * into the z vector.
    *          
    ******************************************************************/
   for (ln = hierarchy->getFinestLevelNumber(); ln >= 0; ln-- ) {
      Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);

      for (PatchLevel<NDIM>::Iterator p(level); p; p++) {
         Pointer<Patch<NDIM> > patch = level->getPatch(p());

         Pointer< CellData<NDIM,double> > soln_scratch =
            patch->getPatchData(d_soln_scr_id);
         Pointer< CellData<NDIM,double> > z = 
            patch->getPatchData(z_indx);

         z->copy(*soln_scratch);
      }

   }

   if (d_print_solver_info) {
      double avg_convergence, final_convergence;
      d_FAC_solver.getConvergenceFactors( avg_convergence, final_convergence );
      pout << "\t\tAfter FAC Solve (Az=r): " 
           << "\n   \t\t\tz_l2norm = " << z_samvect->L2Norm() 
           << "\n   \t\t\tz_maxnorm = " << z_samvect->maxNorm()
           << "\n   \t\t\tResidual Norm: " << d_FAC_solver.getResidualNorm()
           << "\n   \t\t\tConvergence Error: " << final_convergence
           << endl; 
   }
 
   if (converge != true) {
      ret_val = 1;
   }

   /*
    * Increment counter for number of precond solves
    */
   d_number_precond_solve++;
   
   return(ret_val);

#else
  
   return(0);
   
#endif

}

/*************************************************************************
 *
 * Methods specific to CVODEModel class.
 *
 ************************************************************************/

void
CVODEModel::setupSolutionVector(
   Pointer<PatchHierarchy<NDIM> > hierarchy)
{
   /* create SAMRAIVector */
   Pointer< SAMRAIVectorReal<NDIM,double> > soln_samvect =
      new SAMRAIVectorReal<NDIM,double>("solution", hierarchy,
          0, hierarchy->getFinestLevelNumber());
   soln_samvect->addComponent(d_soln_var,d_soln_cur_id);

   /* allocate memory for vectors. */
   soln_samvect->allocateVectorData();

   /* create SundialsAbstractVector */
   d_solution_vector = 
      Sundials_SAMRAIVector<NDIM>::createSundialsVector(soln_samvect);

#ifdef USE_FAC_PRECONDITIONER
   /* 
    * Allocate memory for preconditioner variables. 
    */

   const int nlevels = hierarchy->getNumberOfLevels();

   for (int ln = 0; ln < nlevels; ln++) {
      Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(!(level.isNull()));
#endif
      level->allocatePatchData(d_diff_id);
      if (d_use_neumann_bcs) {
         level->allocatePatchData(d_flag_id);
         level->allocatePatchData(d_neuf_id);
      }
      
   }
#endif
   
}

SundialsAbstractVector* 
CVODEModel::getSolutionVector(void)
{
   return (d_solution_vector);
}

/*
*************************************************************************
*
* Set initial conditions for CVODE solver                               *
*                                                                       *
*************************************************************************
*/
void 
CVODEModel::setInitialConditions(SundialsAbstractVector* soln_init)
{
   Pointer< SAMRAIVectorReal<NDIM,double> > soln_init_samvect =
      Sundials_SAMRAIVector<NDIM>::getSAMRAIVector(soln_init);

   Pointer<PatchHierarchy<NDIM> > hierarchy = soln_init_samvect->getPatchHierarchy();

   for (int ln = 0; ln < hierarchy->getNumberOfLevels(); ln++) {
      Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);

      for (int cn = 0; cn < soln_init_samvect->getNumberOfComponents(); cn++) {
         for (PatchLevel<NDIM>::Iterator p(level); p; p++) {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM> > p_geom = patch->getPatchGeometry();
    
            /*
             * Set initial conditions for y
             */
            Pointer< CellData<NDIM,double> > y_init =
               soln_init_samvect->getComponentPatchData(cn,*patch);
            y_init->fillAll(d_initial_value);

            /*
             * Set initial diffusion coeff values.
             * NOTE: in a "real" application, the diffusion coefficient is 
             * some function of y.  Here, we just do a simple minded 
             * approach and set it to 1.
             */
            Pointer< SideData<NDIM, double> > diffusion = 
               patch->getPatchData(d_diff_id);
 
            diffusion->fillAll(1.0);
         }
      }
   }
}
 
/*
*************************************************************************
*                                                                       *
* Return array of program counters.  Currently, the array holds the     *
* following entries:                                                    *
*    1) number of RHS evaluations                                       *
*    2) number of precond setup calls                                   *
*    3) number of precond solve calls                                   *
* More counters may be added, as desired.                               *
*                                                                       *
*************************************************************************
*/
void 
CVODEModel::getCounters(Array<int>& counters)
{
   counters.resizeArray(3);
   counters[0] = d_number_rhs_eval;
   counters[1] = d_number_precond_setup;
   counters[2] = d_number_precond_solve;
}

/*
*************************************************************************
*
* Get data from input database.                                         *
*                                                                       *
*************************************************************************
*/ 
void 
CVODEModel::getFromInput(Pointer<Database> input_db,
                         bool is_from_restart) 
{

   d_initial_value = input_db->getDoubleWithDefault("initial_value",0.0);

   IntVector<NDIM> periodic = d_grid_geometry->getPeriodicShift();
   int num_per_dirs = 0;
   for (int id = 0; id < NDIM; id++) {
      if (periodic(id)) num_per_dirs++;
   }

   if (input_db->keyExists("Boundary_data")) {
      Pointer<Database> boundary_db = input_db->getDatabase("Boundary_data");

#if (NDIM == 2)
      CartesianBoundaryUtilities2::readBoundaryInput(this,
                                                     boundary_db,
                                                     d_scalar_bdry_edge_conds,
                                                     d_scalar_bdry_node_conds,
                                                     periodic);
#endif
#if (NDIM == 3)
      CartesianBoundaryUtilities3::readBoundaryInput(this,
                                                     boundary_db,
                                                     d_scalar_bdry_face_conds,
                                                     d_scalar_bdry_edge_conds,
                                                     d_scalar_bdry_node_conds,
                                                     periodic);
#endif
         
   } else {
      TBOX_WARNING(d_object_name << ": "
                   << "Key data `Boundary_data' not found in input. " << endl);
   }

#ifdef USE_FAC_PRECONDITIONER
   d_max_fac_its = 
      input_db->getIntegerWithDefault("max_fac_its",d_max_fac_its);
   d_fac_tol = 
      input_db->getDoubleWithDefault("fac_tol",d_fac_tol);
   d_max_hypre_its = 
      input_db->getIntegerWithDefault("max_hypre_its",d_max_hypre_its);
   d_hypre_tol = 
      input_db->getDoubleWithDefault("hypre_tol",d_hypre_tol);
   d_use_neumann_bcs = 
      input_db->getBoolWithDefault("use_neumann_bcs", d_use_neumann_bcs);
   d_print_solver_info = 
      input_db->getBoolWithDefault("print_solver_info", d_print_solver_info);
#endif

}

/*
*************************************************************************
*
* Write data to  restart database.                                      *
*                                                                       *
*************************************************************************
*/
void CVODEModel::putToDatabase( Pointer<Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!db.isNull());
#endif

   db->putInteger("CVODE_MODEL_VERSION",CVODE_MODEL_VERSION);

   db->putDouble("d_initial_value",d_initial_value);

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
*
* Read data from restart database.                                      *
*                                                                       *
*************************************************************************
*/
void CVODEModel::getFromRestart()
{

   Pointer<Database> root_db = 
      RestartManager::getManager()->getRootDatabase();

   Pointer<Database> db;
   if ( root_db->isDatabase(d_object_name) ) {
      db = root_db->getDatabase(d_object_name);
   } else {
      TBOX_ERROR("Restart database corresponding to "
              << d_object_name << " not found in the restart file.");
   } 
   
   int ver  = db->getInteger("CVODE_MODEL_VERSION");
   if (ver != CVODE_MODEL_VERSION) {
      TBOX_ERROR(d_object_name << ":  "
              << "Restart file version different than class version.");
   }

   d_initial_value = db->getDouble("d_initial_value");
   
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

void CVODEModel::readDirichletBoundaryDataEntry(Pointer<Database> db,
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
                      d_bdry_edge_val);
#endif
#if (NDIM == 3)
   readStateDataEntry(db,
                      db_name,
                      bdry_location_index,
                      d_bdry_face_val);
#endif
}

void CVODEModel::readNeumannBoundaryDataEntry(Pointer<Database> db,
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
                      d_bdry_edge_val);
#endif
#if (NDIM == 3)
   readStateDataEntry(db,
                      db_name,
                      bdry_location_index,
                      d_bdry_face_val);
#endif
}

void CVODEModel::readStateDataEntry(Pointer<Database> db,
                               const string& db_name,
                               int array_indx,
                               Array<double>& val)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!db.isNull());
   TBOX_ASSERT(!db_name.empty());
   TBOX_ASSERT(array_indx >= 0);
   TBOX_ASSERT(val.getSize() > array_indx);
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
* Prints class data - writes out info in class if exception is thrown   *
*                                                                       *
*************************************************************************
*/

void CVODEModel::printClassData(ostream &os) const
{
   fflush(stdout);
   int j;

   os << "ptr CVODEModel = " << (CVODEModel*) this << endl;

   os << "d_object_name = " << d_object_name << endl;

   os << "d_soln_cur_id = " << d_soln_cur_id << endl;
   os << "d_soln_scr_id = " << d_soln_scr_id << endl;

   os << "d_initial_value = " << d_initial_value << endl;
   
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

}

void CVODEModel::setPrintSolverInfo(const bool info)
{
   d_print_solver_info = info;
}
#endif // HAVE_SUNDIALS
