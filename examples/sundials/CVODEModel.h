//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/examples/sundials/CVODEModel.h $
// Package:     SAMRAI mesh
// Copyright:   (c) 1997-2002 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1749 $
// Modified:    $LastChangedDate: 2007-12-10 15:14:02 -0800 (Mon, 10 Dec 2007) $
// Description: Example demonstrating use of CVODE vectors. 
//
 
#include "SAMRAI_config.h"

#ifndef included_iostream
#define included_iostream
#include <iostream>
using namespace std;
#endif

#if !defined(HAVE_SUNDIALS) || !defined(HAVE_HYPRE)

/*
*************************************************************************
* If the library is not compiled with CVODE, print an error.
* If we're running autotests, skip the error and compile an empty
* class.
*************************************************************************
*/
#if (TESTING != 1)
#error "This example requires SAMRAI be compiled with CVODE -and- HYPRE."
#endif

#else

/*
 * Header file for base classes.
 */
#include "BoundaryUtilityStrategy.h"
#include "StandardTagAndInitStrategy.h"
#include "RefinePatchStrategy.h"
#include "CoarsenPatchStrategy.h"
#include "CVODEAbstractFunctions.h"

/*
 * Header file for SAMRAI classes referenced in this class.
 */
#include "tbox/Array.h"
#include "Box.h"
#include "CellVariable.h"
#include "CartesianGridGeometry.h"
#include "tbox/Database.h"
#include "FaceVariable.h"
#include "Geometry.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "tbox/Pointer.h"
#include "OuterfaceVariable.h"
#include "SideVariable.h"
#include "VariableContext.h"
#include "RefineAlgorithm.h"
#include "RefineSchedule.h"

#define USE_FAC_PRECONDITIONER
// comment out line below to invoke preconditioning 
// #undef USE_FAC_PRECONDITIONER 

#ifdef USE_FAC_PRECONDITIONER
#include "CellPoissonFACSolver.h"
#endif
/*
 * Header files for CVODE wrapper classes
 */ 
#include "SundialsAbstractVector.h"
#include "Sundials_SAMRAIVector.h"

using namespace SAMRAI;
using namespace tbox;
using namespace hier;
using namespace xfer;
using namespace pdat;
using namespace math;
using namespace mesh;
using namespace geom;
using namespace solv;
using namespace appu;

#ifndef NULL
#define NULL (0)
#endif

/**
 * The cvode_Model class tests the CVODE-SAMRAI interface using
 * two problems: (1) y' = k * d/dx (dy/dx) and (2) y' = y.  
 * 
 * The choice of which problem to solve and other input parameters
 * are specified through the input database.
 * 
 * Input Parameters:
 * 


 *
 *    - \b Problem_type                
 *       1 for diffusion equation, 2 for y' = y.  By default, the heat
 *       equation is solved.
 *
 *    - \b Diffusion_coefficient       
 *       specifies the diffusion coefficient to use when the
 *       has been specified that the diffusion equation will be 
 *       solved.
 *
 *    - \b Initial_condition_type     
 *       0 for constant initial conditions, 1 for sinusoidal initial
 *       conditions
 *
 *    - \b Initial_value               
 *       specifies the initial value to be used for all grid points
 *       when constant initial conditions is specified 
 *
 *    - \b Boundary_value              
 *       specifies what value should be used for the dirichlet
 *       boundary conditions applied to this problem.
 *
 * 


 */

class CVODEModel : 
   public StandardTagAndInitStrategy<NDIM>,
   public RefinePatchStrategy<NDIM>,
   public CoarsenPatchStrategy<NDIM>,
   public BoundaryUtilityStrategy,
   public CVODEAbstractFunctions
{
public:
   /**
    * Default constructor for CVODEModel.
    */     
   CVODEModel(const string& object_name,
              tbox::Pointer<tbox::Database> input_db,
              tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geom);

   /**
    * Empty destructor for CVODEModel.
    */
   virtual ~CVODEModel();

/*************************************************************************
 *
 * Methods inherited from StandardTagAndInitStrategy<NDIM>.
 *
 ************************************************************************/
 
   /**
    * Initialize data on a new level after it is inserted into an AMR patch
    * hierarchy by the gridding algorithm.  The level number indicates
    * that of the new level.
    *
    * Generally, when data is set, it is interpolated from coarser levels
    * in the hierarchy.  If the old level pointer in the argument list is
    * non-null, then data is copied from the old level to the new level
    * on regions of intersection between those levels before interpolation
    * occurs.   In this case, the level number must match that of the old 
    * level.  The specific operations that occur when initializing level 
    * data are determined by the particular solution methods in use; i.e.,
    * in the subclass of this abstract base class.
    *
    * The boolean argument initial_time indicates whether the level is
    * being introduced for the first time (i.e., at initialization time),
    * or after some regrid process during the calculation beyond the initial
    * hierarchy construction.  This information is provided since the
    * initialization of the data may be different in each of those
    * circumstances.  The can_be_refined boolean argument indicates whether
    * the level is the finest allowable level in the hierarchy.
    */
   virtual void
   initializeLevelData(const tbox::Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                       const int level_number,
                       const double time,
                       const bool can_be_refined,
                       const bool initial_time,
                       const tbox::Pointer<BasePatchLevel<NDIM> > old_level = 
                             tbox::Pointer<BasePatchLevel<NDIM> >(NULL),
                       const bool allocate_data = true);

   /**
    * After hierarchy levels have changed and data has been initialized on 
    * the new levels, this routine can be used to reset any information 
    * needed by the solution method that is particular to the hierarchy 
    * configuration.  For example, the solution procedure may cache 
    * communication schedules to amortize the cost of data movement on the 
    * AMR patch hierarchy.  This function will be called by the gridding 
    * algorithm after the initialization occurs so that the algorithm-specific
    * subclass can reset such things.  Also, if the solution method must 
    * make the solution consistent across multiple levels after the hierarchy 
    * is changed, this process may be invoked by this routine.  Of course the 
    * details of these processes are determined by the particular solution 
    * methods in use.
    *
    * The level number arguments indicate the coarsest and finest levels
    * in the current hierarchy configuration that have changed.  It should
    * be assumed that all intermediate levels have changed as well.
    */
   virtual void resetHierarchyConfiguration(
      const tbox::Pointer< BasePatchHierarchy<NDIM> > hierarchy,
      const int coarsest_level,
      const int finest_level);

   /**
    * Set tags to the specified tag value where refinement of the given
    * level should occur using the user-supplied gradient detector.  The 
    * value "tag_index" is the index of the cell-centered integer tag
    * array on each patch in the hierarchy.  The boolean argument indicates 
    * whether cells are being tagged on the level for the first time; 
    * i.e., when the hierarchy is initially constructed.  If it is false, 
    * it should be assumed that cells are being tagged at some later time 
    * after the patch hierarchy was initially constructed.  This information
    * is provided since the application of the error estimator may be 
    * different in each of those circumstances.
    */
   virtual void 
   applyGradientDetector(const tbox::Pointer< BasePatchHierarchy<NDIM> > hierarchy,
                         const int level_number,
                         const double time,
                         const int tag_index,
                         const bool initial_time,
                         const bool uses_richardson_extrapolation_too);

   /**
    * Option to output solver info.  Set to true to turn on, false to 
    * turn off.
    */
   void setPrintSolverInfo(const bool info);

/*************************************************************************
 *
 * Methods inherited from RefinePatchStrategy.
 *
 ************************************************************************/
 
   /**
    * Set the data at patch boundaries corresponding to the physical domain
    * boundary.  The specific boundary conditions are determined by the user.
    */
   virtual void setPhysicalBoundaryConditions(
      Patch<NDIM> & patch,
      const double time,
      const IntVector<NDIM> & ghost_width_to_fill);

   /**
    * Perform user-defined refining operations.  This member function
    * is called before the other refining operators.  The preprocess
    * function should refine data from the scratch components of the
    * coarse patch into the scratch components of the fine patch on the
    * specified fine box region.  This version of the preprocess function
    * operates on a a single box at a time.  The user must define this
    * routine in the subclass.
    */
   virtual void preprocessRefine(
      Patch<NDIM> & fine,
      const Patch<NDIM> & coarse,
      const Box<NDIM> & fine_box,
      const IntVector<NDIM> & ratio);

   /**
    * Perform user-defined refining operations.  This member function
    * is called after the other refining operators.  The postprocess
    * function should refine data from the scratch components of the
    * coarse patch into the scratch components of the fine patch on the
    * specified fine box region.  This version of the postprocess function
    * operates on a a single box at a time.  The user must define this
    * routine in the subclass.
    */
   virtual void postprocessRefine(
      Patch<NDIM> & fine,
      const Patch<NDIM> & coarse,
      const Box<NDIM> & fine_box,
      const IntVector<NDIM> & ratio);

   /**
    * Return maximum stencil width needed for user-defined
    * data interpolation operations.  Default is to return
    * zero, assuming no user-defined operations provided.
    */
   virtual IntVector<NDIM>  getRefineOpStencilWidth() const
   {
      return(IntVector<NDIM>(0));
   }


/*************************************************************************
 *
 * Methods inherited from CoarsenPatchStrategy.
 *
 ************************************************************************/
 
   /**
    * Perform user-defined coarsening operations.  This member function
    * is called before the other coarsening operators.  The preprocess
    * function should copy data from the source components of the fine
    * patch into the source components of the destination patch on the
    * specified coarse box region.
    */
   virtual void preprocessCoarsen(
      Patch<NDIM> & coarse,
      const Patch<NDIM> & fine,
      const Box<NDIM> & coarse_box,
      const IntVector<NDIM> & ratio);

   /**
    * Perform user-defined coarsening operations.  This member function
    * is called after the other coarsening operators.  The postprocess
    * function should copy data from the source components of the fine
    * patch into the source components of the destination patch on the
    * specified coarse box region.
    */
   virtual void postprocessCoarsen(
      Patch<NDIM>& coarse,
      const Patch<NDIM>& fine,
      const Box<NDIM>& coarse_box,
      const IntVector<NDIM>& ratio);

   /**
    * Return maximum stencil width needed for user-defined
    * data interpolation operations.  Default is to return
    * zero, assuming no user-defined operations provided.
    */
   virtual IntVector<NDIM> getCoarsenOpStencilWidth() const
   {
      return(IntVector<NDIM>(0));
   }


/*************************************************************************
 *
 * Methods inherited from CVODEAbstractFunctions
 *
 ************************************************************************/

   /**
    * User-supplied right-hand side function evaluation.
    *
    * The function arguments are:
    * 


    * - \b t        (INPUT) {current value of the independent variable}
    * - \b y        (INPUT) {current value of dependent variable vector}
    * - \b y_dot   (OUTPUT){current value of the derivative of y}
    * 


    *
    * IMPORTANT: This function must not modify the vector y. (KTC??)
    */
   virtual int evaluateRHSFunction(double time,
                                    SundialsAbstractVector* y,
                                    SundialsAbstractVector* y_dot);

   virtual int CVSpgmrPrecondSet(double t,
				 SundialsAbstractVector* y,
				 SundialsAbstractVector* fy,
				 int jok,
				 int *jcurPtr,
				 double gamma,
				 SundialsAbstractVector* vtemp1,
				 SundialsAbstractVector* vtemp2,
				 SundialsAbstractVector* vtemp3);
   
   virtual int CVSpgmrPrecondSolve(double t,
				   SundialsAbstractVector* y,
				   SundialsAbstractVector* fy,
				   SundialsAbstractVector* r,
				   SundialsAbstractVector* z,
				   double gamma,
				   double delta,
				   int lr,
				   SundialsAbstractVector* vtemp);

/*************************************************************************
 *
 * Methods particular to CVODEModel class.
 *
 ************************************************************************/

   /** 
    * Set up solution vector.
    */
   void setupSolutionVector(tbox::Pointer<PatchHierarchy<NDIM> > hierarchy);
  
   /**
    * Get pointer to the solution vector.
    */
   SundialsAbstractVector* getSolutionVector(void);

   /**
    * Set initial conditions for problem.
    */
   void setInitialConditions(SundialsAbstractVector* y_init);

   /**
    * Return array of program counters.
    */
   void getCounters(tbox::Array<int>& counters);

   /**
    * Writes state of CVODEModel object to the specified database.
    *
    * This routine is a concrete implementation of the function
    * declared in the tbox::Serializable abstract base class.
    */
   void putToDatabase( tbox::Pointer<tbox::Database> db);

   /**
    * This routine is a concrete implementation of the virtual function
    * in the base class BoundaryUtilityStrategy.  It reads DIRICHLET
    * and NEUMANN boundary state values from the given database with the
    * given name string idenifier.  The integer location index
    * indicates the face (in 3D) or edge (in 2D) to which the boundary 
    * condition applies.
    */
   void readDirichletBoundaryDataEntry(tbox::Pointer<tbox::Database> db,
                                       string& db_name,
                                       int bdry_location_index);

   void readNeumannBoundaryDataEntry(tbox::Pointer<tbox::Database> db,
                                     string& db_name,
                                     int bdry_location_index);

   /**
    * Prints all class data members, if exception is thrown.
    */
   void printClassData(ostream &os) const;

private:
   /*
    * These private member functions read data from input and restart.
    * When beginning a run from a restart file, all data members are read
    * from the restart file.  If the boolean flag is true when reading
    * from input, some restart values may be overridden by those in the
    * input file.
    *
    * An assertion results if the database pointer is null.
    */
   virtual void getFromInput(tbox::Pointer<tbox::Database> db,
                             bool is_from_restart);

   virtual void getFromRestart();

   void readStateDataEntry(tbox::Pointer<tbox::Database> db,
                           const string& db_name,
                           int array_indx,
                           tbox::Array<double>& uval); 

   /*
    * Object name used for error/warning reporting and as a label
    * for restart database entries.
    */
   string d_object_name;

   /*
    * tbox::Pointer to solution vector 
    */
   SundialsAbstractVector* d_solution_vector;

   /*
    * Variables 
    */ 
   tbox::Pointer< CellVariable<NDIM,double> > d_soln_var;


   /*
    * Variable Contexts 
    */ 
   tbox::Pointer<VariableContext>        d_cur_cxt;
   tbox::Pointer<VariableContext>        d_scr_cxt;

   /*
    * Patch Data ids
    */ 
   int d_soln_cur_id;
   int d_soln_scr_id;

#ifdef USE_FAC_PRECONDITIONER
   tbox::Pointer< SideVariable<NDIM,double> > d_diff_var;
   tbox::Pointer< OuterfaceVariable<NDIM,int> > d_flag_var;
   tbox::Pointer< OuterfaceVariable<NDIM,double> > d_neuf_var;
   
   int d_diff_id;
   int d_flag_id;
   int d_neuf_id;
   int d_bdry_types[2*NDIM];

   solv::CellPoissonFACSolver<NDIM> d_FAC_solver;
   bool d_FAC_solver_allocated;
   bool d_level_solver_allocated;
   bool d_use_neumann_bcs;
   

   int d_max_fac_its;
   double d_fac_tol;
   int d_max_hypre_its;
   double d_hypre_tol;
   double d_current_soln_time;
#endif

   /*
    * Print CVODE solver information
    */
   bool d_print_solver_info;
   
  
   /*
    * Grid geometry 
    */
   tbox::Pointer<geom::CartesianGridGeometry<NDIM> > d_grid_geometry;

   /*
    * Initial value
    */
   double d_initial_value;

   /*
    * Program counters  
    *   1 - number of RHS evaluations
    *   2 - number of precond setups
    *   3 - number of precond solves
    */
   int d_number_rhs_eval;
   int d_number_precond_setup;
   int d_number_precond_solve;

   /*
    * Boundary condition cases and boundary values.
    * Options are: FLOW, REFLECT, DIRICHLET, NEUMANN
    * and variants for nodes and edges.
    *
    * Input file values are read into these arrays.
    */
   tbox::Array<int> d_scalar_bdry_edge_conds;
   tbox::Array<int> d_scalar_bdry_node_conds;
#if (NDIM == 3)
   tbox::Array<int> d_scalar_bdry_face_conds;
#endif

   /*
    * Boundary condition cases for scalar and vector (i.e., depth > 1)
    * variables.  These are post-processed input values and are passed
    * to the boundary routines.
    */
#if (NDIM ==2)
   tbox::Array<int> d_node_bdry_edge;
#endif
#if (NDIM == 3)
   tbox::Array<int> d_edge_bdry_face;
   tbox::Array<int> d_node_bdry_face;
#endif

   /*
    * Arrays of face (3d) or edge (2d) boundary values for DIRICHLET case.
    */
#if (NDIM ==2)
   tbox::Array<double> d_bdry_edge_val;
#endif
#if (NDIM ==3)
   tbox::Array<double> d_bdry_face_val;
#endif

};
#endif // HAVE_SUNDIALS
