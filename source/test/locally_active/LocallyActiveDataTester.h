//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/test/locally_active/LocallyActiveDataTester.h $
// Package:     SAMRAI test
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2132 $
// Modified:    $LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description: Class to test locally-active data communication 
//
 
#include "SAMRAI_config.h"

#include "LocallyActiveDataCoarsenAlgorithm.h"
#include "LocallyActiveDataCoarsenPatchStrategy.h"
#include "LocallyActiveDataCoarsenSchedule.h"
#include "LocallyActiveDataRefineAlgorithm.h"
#include "LocallyActiveDataRefinePatchStrategy.h"
#include "LocallyActiveDataRefineSchedule.h"

/*
 * Header files for SAMRAI classes referenced in this class.
 */
#include "Box.h"
#include "BoxList.h"
#include "tbox/Database.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "tbox/Pointer.h"
#include "tbox/Timer.h"
#include "Variable.h"
#include "VisItDataWriter.h"
#include "tbox/IOStream.h"

using namespace std;
using namespace SAMRAI;

/**
 * The LocallyActiveDataTester class coordinates construction of a patch hierarchy,
 * and initialization, mathematical operations for a set of locally-active
 * variables defined on that hierarchy.  LocallyActiveDataTester input data must 
 * follow this format:
 *
 * LocallyActiveDataTester {
 *   
 * Required GriddingParameters input describes AMR patch hierarchy configuration:
 *
 *    GriddingParameters { 
 *       // number of hierarchy levels [optional; default = 2]
 *       num_levels = 2
 *
 *       // integer vectors (length = NDIM) specifying ratio between index 
 *       // space of patch level to next coarser level 
 *       // [required if num_levels > 1] 
 *       ratio_to_coarser {
 *          level_1 = 2, 2, 2
 *          // all finer levels use level_1 values unless specified...
 *       }
 *    }
 *
 *    // Required test specifier; valid options are: 
 *    //    "REFINE_TEST" 
 *    //    "COARSEN_TEST" 
 *    //    "LAPLACIAN_TEST"
 *    test_to_run = string 
 *
 *    // If "LAPLACIAN_TEST" is given. optionally specify number of iterations 
 *    // to execute Laplacian computations.  
 *    laplacian_interations = integer  [optional; default is 1]
 * 
 *    // Optional flag indicating whether to check results of test.
 *    // Typically, set to false for measuring performance.
 *    check_results = bool [optional; default is TRUE]
 * 
 *    // Required test function distribution string identifier.
 *    // Options are: "EXPLICIT_FUNCTIONS", "UNIFORM_FUNCTIONS", "RANDOM_FUNCTIONS",
 *    // "FINE_PATCH_ALL", "FINE_PATCH_MOD2", "FINE_PATCH_EVEN_ODD".
 *    //
 *    // Additional input information is required as an additional input database
 *    // entry for some of the function_distribution choices.  The requirements are:
 *    // "EXPLICIT_FUNCTIONS" requires ExplicitFunctionData { ... }
 *    // "UNIFORM_FUNCTIONS" requires UniformFunctionData { ... }
 *    // "RANDOM_FUNCTIONS" requires RandomFunctionData { ... }
 *    // "FINE_PATCH_ALL", "FINE_PATCH_MOD2", or "FINE_PATCH_EVEN_ODD" require
 *    //  requires FinePatchFunctionData { ... }
 *    // The contents of these additional input databases is described below.
 *    function_distribution = string 
 *
 *    // Optional boolean flag indicating whether to dump plot file with
 *    // functions defined as constants to check which patches they are 
 *    // defined on.
 *    check_function_definitions = bool [optional; default is TRUE]
 *
 *    // Optional ghost cell width for test functions. All refine operations
 *    // use the "function" as source and destination data; thus, test results 
 *    // are checked on "function" data . All coarsen operations use the "solution" 
 *    // as source and destination data; thus, test results are checked on "solution" 
 *    // data. The "LAPLACIAN_TEST" performs both coarsen and refine operations; 
 *    // results are checked on "solution" data.
 *    nghosts_function = integer array size NDIM 
 *                       [optional gcw for function data; default is 1 ... 1]
 *    nghosts_solution = integer array size NDIM 
 *                       [optional gcw for solution data; default is 0 ... 0]
 *
 *    // Additional input information is required as an additional input database
 *    // entry for some of the function_distribution choices. This information is
 *    // used to define the values each test function may take and is summarized below.
 *    // All functions are active on all patches on level zero in the patch hierarchy.
 *    // A centroid and radius for each function determines the active patches for
 *    // finer hierarchy levels.  These parameters are illustrated below.
 *
 *    // Also, each case requires a function specification string to be given.
 *    // Valid choices are: "CONSTANT_TEST", "CONSTANT_PERF", or "GAUSSIAN".  
 *    // Both "CONSTANT_TEST" and "CONSTANT_PERF" options define functions that are
 *    // contant on the patches on which they are defined; "CONSTANT_TEST" sets 
 *    // function value based on function number, level number, and patch number, while
 *    // "CONSTANT-PERF" sets function value based on function number only.  
 *    // The "GAUSSIAN" choice defines each function as f(x) = exp( -(x-centroid)^2*alpha ) 
 *    // inside sphere of given radius about centroid and zero elsewhere. 
 *
 * ExplicitFunctionData input explicitly defines locally-active functions.
 * This data is required when function_distribution = "EXPLICIT_FUNCTIONS" is 
 * specified.
 *
 *    ExplicitFunctionData {
 *       function_spec = [required; either "CONSTANT_TEST", "CONSTANT_PERF", or "GAUSSIAN"]
 *  
 *       // explicitly-specified functions
 *       Function0 {
 *          centroid = double array size NDIM [required]
 *          radius   = double [required]
 *          alpha    = double [required; only used when function_spec = "GAUSSIAN"]
 *       }
 *       Function1 { ... }
 *       // etc.
 *    }
 * 
 * UniformFunctionData input describes locally-defined functions uniformly
 * distributed within domain.  This data is required when 
 * function_distribution = "UNIFORM_FUNCTIONS" is specified.
 *
 *    UniformFunctionData {
 *       function_spec = [required; either "CONSTANT_TEST", "CONSTANT_PERF", or "GAUSSIAN"]
 *       radius        = double   [required]
 *       alpha         = double    [required; only used when function_spec = "GAUSSIAN"]
 *       nfunctions    = int array of size NDIM  [required]
 *
 *       // optional lower and upper bounds for domain of uniform functions
 *       lower_bound   = double array size NDIM
 *       upper bound   = double array size NDIM
 *    }
 *
 * RandomFunctionData input gives number of random functions to generate. This
 * data is required when function_distribution = "RANDOM_FUNCTIONS".
 *
 *    RandomFunctionData {
 *       function_spec      = [required; either "CONSTANT_TEST", "CONSTANT_PERF", or "GAUSSIAN"]
 *       num_test_functions = int [required]
 *    } 
 * 
 * FinePatchFunctionData input describes locally-defined functions distributed
 * within domain in approximate uniform distribution based on finest level patches.  
 * This data is required when function_distribution = "FINE_PATCH_ALL", "FINE_PATCH_MOD2", or 
 * "FINE_PATCH_EVEN_ODD" is specified.
 *
 *    FinePatchFunctionData {
 *       function_spec = [required; either "CONSTANT_TEST", "CONSTANT_PERF", or "GAUSSIAN"]
 *       radius        = double   [required]
 *       alpha         = double [required; only used when function_spec = "GAUSSIAN"]
 *    }
 * 
 */

class LocallyActiveDataTester : 
   public xfer::LocallyActiveDataRefinePatchStrategy<NDIM>,
   public xfer::LocallyActiveDataCoarsenPatchStrategy<NDIM>
{
public:

   /**
    * The LocallyActiveDataTester::Function data structure contains information
    * about each function.  The functions are defined as 
    * f(x) = exp(-alpha*|x-centroid|) so we store the centroid, radius,
    * and alpha for each function. In addition, we store its name and 
    * pointers to the multi-function variables associated with it (function
    * and solution data).
    */
   struct Function {
      string                  d_name;
      tbox::Pointer< hier::Variable<NDIM> >       d_func;
      int                     d_func_data_index;
      tbox::Pointer< hier::Variable<NDIM> >       d_soln;
      int                     d_soln_data_index;
      double                  d_radius;
      double                  d_centroid[NDIM];
      double                  d_alpha;
   };
      
   /**
    * Constructor for LocallyActiveDataTester.
    */     
   LocallyActiveDataTester(const string& object_name,
                           tbox::Pointer<tbox::Database> input_db,
                           tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
                           tbox::Pointer<appu::VisItDataWriter<NDIM> > visit_writer);

   /**
    * Boring destructor for LocallyActiveDataTester.
    */
   virtual ~LocallyActiveDataTester();

   /**
    * Construct patch hierarchy based on input file information.
    */
   void buildPatchHierarchy();

   /**
    * Set active patches for each variable.
    */
   void setActivePatchesOnHierarchy();

   /**
    * Initialize data for each function
    */
   void initializeFunctionData();

   /**
    * Setup communication schedules based on test.
    */
   void setupCommunication();

   /**
    * Perform specified test.
    */
   void performTest();

   /**
    * Check results of test for pass or failure and return integer number
    * of failures.
    */
   int checkTestResult(ostream &os) const;

   /**
    * Dump the various function data on the hierarchy to the specified
    * output stream.
    */
   void printHierarchyData(ostream &os,
                           bool dump_function_data = false) const;

   /*
    * Methods inherited from xfer::LocallyActiveDataRefinePatchStrategy.
    */

   hier::IntVector<NDIM> getRefineOpStencilWidth() const {return(hier::IntVector<NDIM>(0));}
 
   void setPhysicalBoundaryConditions(
      hier::Patch<NDIM>& patch,
      const tbox::List<int>& scratch_data_ids,
      const double fill_time,
      const hier::IntVector<NDIM>& ghost_width_to_fill);

   void preprocessRefine(
      hier::Patch<NDIM>& fine,
      const hier::Patch<NDIM>& coarse,
      const tbox::List<int>& scratch_data_ids,
      const hier::Box<NDIM>& fine_box,
      const hier::IntVector<NDIM>& ratio)
   {
      (void) fine;
      (void) coarse;
      (void) scratch_data_ids;
      (void) fine_box;
      (void) ratio;
   }
  
   void postprocessRefine(
      hier::Patch<NDIM>& fine,
      const hier::Patch<NDIM>& coarse,
      const tbox::List<int>& scratch_data_ids,
      const hier::Box<NDIM>& fine_box,
      const hier::IntVector<NDIM>& ratio)
   {
      (void) fine;
      (void) coarse;
      (void) scratch_data_ids;
      (void) fine_box;
      (void) ratio;
   }

   /*
    * Methods inherited from xfer::LocallyActiveDataCoarsenPatchStrategy.
    */

   hier::IntVector<NDIM> getCoarsenOpStencilWidth() const {return(hier::IntVector<NDIM>(0));}

   void preprocessCoarsen(hier::Patch<NDIM>& coarse,
                          const hier::Patch<NDIM>& fine,
                          const tbox::List<int>& src_data_ids,
                          const hier::Box<NDIM>& coarse_box,
                          const hier::IntVector<NDIM>& ratio)
   {
      (void) coarse;
      (void) fine;
      (void) src_data_ids;
      (void) coarse_box;
      (void) ratio;
   }
   
   void postprocessCoarsen(hier::Patch<NDIM>& coarse,
                           const hier::Patch<NDIM>& fine,
                           const tbox::List<int>& src_data_ids,
                           const hier::Box<NDIM>& coarse_box,
                           const hier::IntVector<NDIM>& ratio)
   {
      (void) coarse;
      (void) fine;
      (void) src_data_ids;
      (void) coarse_box;
      (void) ratio;
   } 

   int getLevelWorkUnits(int ln) const 
   {
      int ret_val = -1;
      if (d_total_work_units.size() > ln) {
         ret_val = d_total_work_units[ln];
      }
      return(ret_val);
   }

private:
   /*
    * Read input parameters for patch hierarchy and locally-active data.
    */
   void getFromInput();
   void getGriddingParametersFromInput(tbox::Pointer<tbox::Database> gridding_db);
   void getExplicitFunctionDataFromInput(tbox::Pointer<tbox::Database> explicit_data_db);
   void getUniformFunctionDataFromInput(tbox::Pointer<tbox::Database> uniform_data_db);
   void getRandomFunctionDataFromInput(tbox::Pointer<tbox::Database> random_data_db);
   void getFinePatchFunctionDataFromInput(tbox::Pointer<tbox::Database> fine_patch_data_db);
   void getFunctionSpecification(tbox::Pointer<tbox::Database> function_db);

   /*
    * Create variables and register them with variable database. 
    */
   void registerVariables(tbox::Pointer<appu::VisItDataWriter<NDIM> > visit_writer);

   /*
    * Determine whether patch is active for given function.
    */
   bool patchIsActiveForFunction(int fn, tbox::Pointer<hier::PatchLevel<NDIM> > level, int ip);

   /*
    * Allocate active patch data and set initial values based on test case.
    */
   void setGaussianDataOnPatch(tbox::Pointer<hier::Patch<NDIM> > patch, int fn);
   void setRefineTestDataOnPatch(tbox::Pointer<hier::Patch<NDIM> > patch, int fn);
   void setCoarsenTestDataOnPatch(tbox::Pointer<hier::Patch<NDIM> > patch, int fn);

   void initializeFunctionDataForVisitCheck();

   /*
    * Return value used in either refine or coarsen test case, which is based on
    * level number, patch number and patch data index.
    */
   double getRefineTestValue(int level_number,
                             int patch_number,
                             int data_id) const;
   double getCoarsenTestValue(int level_number,
                              int patch_number,
                              int data_id) const;

   /*
    * Set boundary values for given patch, boundary fill box, and function.
    */
   void computeLaplacians();

   /*
    * Compute Laplcaian of each function, if needed for test.
    */
   void setPatchBoundaryValues(hier::Patch<NDIM>& patch,
                               const hier::Box<NDIM>& fill_box,
                               int function_id,
                               int data_index);

   /*
    * Check results for specified test.
    */
   int checkRefineTestResult(ostream &os) const;

   bool checkRefineTestOnLevelInterior(
      tbox::Pointer<hier::PatchLevel<NDIM> > check_level,
      tbox::Pointer<hier::PatchLevel<NDIM> > src_level,
      tbox::Array< tbox::Array< hier::BoxList<NDIM> > >& unchecked_boxes,
      tbox::Array<bool>& patch_is_done,
      ostream &os) const;

   bool checkRefineTestPhysicalBoundaries(tbox::Pointer<hier::Patch<NDIM> > patch,
                                          int function_id,
                                          double correct_value,
                                          hier::BoxList<NDIM>& check_boxes_to_modify,
                                          ostream &os) const; 

   bool checkRefineTestOnBox(tbox::Pointer<hier::Patch<NDIM> > patch,
                             int function_id,
                             const hier::Box<NDIM>& region,
                             double correct_value,
                             ostream &os) const;

   bool checkLevelIsDone(
      tbox::Pointer<hier::PatchLevel<NDIM> > level,
      tbox::Array< tbox::Array<hier::BoxList<NDIM> > >& unchecked_boxes,
      tbox::Array<bool>& patch_is_done) const;

   int checkCoarsenTestResult(ostream &os) const;

   /*
    * Object string name identifier for error reporting, etc.
    */
   string d_object_name;
  
   /*
    * Cached pointer to input database for LocallyActiveDataTester and input
    * data members for controlling test functions.
    */
   tbox::Pointer<tbox::Database> d_input_db;
   string d_test_to_run;
   int d_num_test_functions; 
   int d_laplacian_iterations; 
   bool d_check_results;
   bool d_check_function_definitions;

   /* 
    * Enumerated types for different test cases to be used in switch statements
    * and to avoid excessive string comparisions.
    */
   
   enum FCN_DIST_CASE {UNDEFINED_FCN_DIST = 0,
                       EXPLICIT_FUNCTIONS = 1,
                       UNIFORM_FUNCTIONS = 2,
                       RANDOM_FUNCTIONS = 3,
                       FINE_PATCH_ALL = 4,
                       FINE_PATCH_MOD2 = 5,
                       FINE_PATCH_EVEN_ODD = 6};

   enum FCN_SPEC_CASE {UNDEFINED_FCN_SPEC = 0,
                       CONSTANT_TEST = 1,
                       CONSTANT_PERF = 2,
                       GAUSSIAN = 3};

   FCN_DIST_CASE d_function_distribution;
   FCN_SPEC_CASE d_function_specification;

   /*
    * tbox::Array of functions.
    */
   tbox::Array<LocallyActiveDataTester::Function> d_functions;

   /*
    * Number of ghosts cells used for all function/solution variables.
    */
   hier::IntVector<NDIM> d_nghosts_function;
   hier::IntVector<NDIM> d_nghosts_solution;
  
   /*
    * Objects storing AMR patch hierarchy information. 
    */ 
   tbox::Pointer< hier::PatchHierarchy<NDIM> > d_hierarchy;
   int d_num_levels;
   tbox::Array< hier::IntVector<NDIM> > d_ratio_to_coarser;
   hier::IntVector<NDIM> d_npatches_on_coarsest;

   /*
    * Refine and Coarsen algorithms and schedules, used to communicate
    * the different function data over the hierarchy.  
    */
   xfer::LocallyActiveDataRefineAlgorithm<NDIM> d_bdry_fill_alg;
   xfer::LocallyActiveDataCoarsenAlgorithm<NDIM> d_coarsen_alg;
   
   tbox::Array< tbox::Pointer< xfer::LocallyActiveDataRefineSchedule<NDIM> > > d_bdry_fill_sched;
   tbox::Array< tbox::Pointer< xfer::LocallyActiveDataCoarsenSchedule<NDIM> > > d_coarsen_sched;
   
   /*
    * Arrays of Timers, length set to number of levels.
    */
  tbox::Array< tbox::Pointer<tbox::Timer> > d_coarsen_comm_setup_timers;
  tbox::Array< tbox::Pointer<tbox::Timer> > d_refine_comm_setup_timers;
  tbox::Array< tbox::Pointer<tbox::Timer> > d_coarsen_comm_execute_timers;
  tbox::Array< tbox::Pointer<tbox::Timer> > d_refine_comm_execute_timers; 
  tbox::Array< tbox::Pointer<tbox::Timer> > d_laplacian_execute_timers;
  tbox::Pointer<tbox::Timer> d_main_laplacian_timer;
  tbox::Pointer<tbox::Timer> d_sum_refine_comm_setup_timer;
  tbox::Pointer<tbox::Timer> d_sum_refine_comm_execute_timer;

  tbox::Array<int> d_total_work_units;
   
};
