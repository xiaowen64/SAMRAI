//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/examples/ConvDiff/ConvDiff.h $
// Package:     SAMRAI application
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Numerical routines for single patch in Heat equation ex.
//
 
#ifndef included_ConvDiffXD
#define included_ConvDiffXD

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#include "tbox/AbstractStream.h"
#include "Box.h"
#include "BoundaryUtilityStrategy.h"
#include "CartesianGridGeometry.h"
#include "CartesianVizamraiDataWriter.h"
#include "CellVariable.h"
#include "tbox/Database.h"
#include "IntVector.h"
#include "Index.h"
#include "MethodOfLinesIntegrator.h"
#include "MethodOfLinesPatchStrategy.h"
#include "Patch.h"
#include "tbox/Pointer.h"
#include "tbox/Serializable.h"
#include <string>
using namespace std;
#define included_String
#include "VariableContext.h"
#include "VisItDataWriter.h"

/**
 * The ConvDiff class provides numerical routines for a sample problem
 * which illustrates use of AMR for solution of a system of ODEs.
 * This class is derived from the algs::MethodOfLinesPatchStrategy<NDIM>
 * and provides implementations of the virtual functions declared in that 
 * class.  Other member functions are specific to this application.  Most 
 * member functions in ConvDiff provide numerical routines that apply to 
 * individual patches in an AMR hierarchy. 
 * 
 * The convection-diffusion equation is  du/dt + div(a*u) = mu * div^2(u) 
 * + gamma, where "u" is a scalar-valued function and "a", "mu" and "gamma"
 * are constant vectors.  Time integration of this equation is performed 
 * using the algs::MethodOfLinesIntegrator<NDIM>.  The PDE is cast as a set of ODEs 
 * (i.e. du/dt = F(u) where F(u) = -div(a*u) + mu*div^2(u) + gamma).    
 *
 * The primary numerical quantities are "u" and "F(u)", defined in the
 * code as "primitive_vars" and "function_eval", respectively.  All
 * other variables are temporary quantities used in the numerical routines.
 */

#define NEQU  (1) // depth of u

using namespace SAMRAI;

class ConvDiff : 
  public tbox::Serializable,
  public algs::MethodOfLinesPatchStrategy<NDIM>,
  public appu::BoundaryUtilityStrategy
{
public:
   /**
    * The constructor for ConvDiff sets default model parameters to
    * initialize the object. It creates variables that represent
    * the state of the solution, and initializes pertinent private
    * data members.  It also registers the object with the 
    * tbox::RestartManager.
    *
    * After setting the default values, the routine calls 
    * calls getFromRestart() if this is a restart case.  It next
    * calls getfromInput() to read values from the input database,
    * potentially overriding those in the restart file.
    */
   ConvDiff(const string& object_name,
           tbox::Pointer<tbox::Database> input_db,
           tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geom);
  
    /**
     * The destructor for ConvDiff.
     */
   ~ConvDiff();
 
   ///
   ///  The following routines:
   ///
   ///      registerModelVariables(),
   ///      initializeDataOnPatch(),
   ///      computeStableDtOnPatch(),
   ///      singleStep,
   ///      tagGradientDetectorCells(),
   ///
   ///  are concrete implementations of functions declared in the
   ///  algs::MethodOfLinesPatchStrategy<NDIM> abstract base class.
   ///

   /**
    * Register the variables with algs::MethodOfLinesIntegrator<NDIM>.  This
    * registration defines the ways in which data will be manipulated on 
    * patches.  Two variable types are available; SOLN and RHS.  For
    * instance, in the solution of du/dt = F(u), u is of SOLN type and
    * F(u) is RHS type.
    *
    * @see algs::MethodOfLinesIntegrator.
    */
   void registerModelVariables(algs::MethodOfLinesIntegrator<NDIM>* integrator);

   /**
    * Set the data on the patch interior to some initial values,
    * depending on the input parameters and numerical routines.
    * If the "initial_time" flag is false, indicating that the
    * routine is called after a regridding stepa the routine does nothing.
    */
   void initializeDataOnPatch(hier::Patch<NDIM>& patch,
                              const double time,
                              const bool initial_time) const;

   /**
    * Compute the stable time increment for a patch using a CFL-based
    * criteria.  Return computed dt.
    */
   double computeStableDtOnPatch(hier::Patch<NDIM>& patch,
                                 const double time) const;

   /**
    * Perform a single step of Runge-Kutta routine.  That is, an nth-order
    * RK scheme will perform n sub-iterations at each timestep to integrate 
    * over time dt.  The singleStep routine performs one of these 
    * sub-iterations.  
    */
   void singleStep(hier::Patch<NDIM>& patch, 
                   const double dt,
                   const double alpha_1,
                   const double alpha_2,
                   const double beta) const;

   /**
    * Tag cells which need refinement. 
    */
   void tagGradientDetectorCells(hier::Patch<NDIM>& patch, 
                                 const double regrid_time,
                                 const bool initial_error,
                                 const int tag_index,
                                 const bool uses_richardson_extrapolation_too);

   ///
   ///  The following routines:
   ///
   ///      setPhysicalBoundaryConditions()
   ///
   ///  are concrete implementations of functions declared in the
   ///  RefinePatchStrategy abstract base class.
   ///


   /**
    * Set the data in ghost cells corresponding to the physical domain 
    * boundary.  Specific boundary conditions are determined by information 
    * specified in input file and numerical routines.  
    */
   void setPhysicalBoundaryConditions(hier::Patch<NDIM>& patch,
                                      const double fill_time,
                                      const hier::IntVector<NDIM>& 
                                      ghost_width_to_fill);

   /**
    * Writes state of ConvDiff object to the specified database.
    *
    * This routine is a concrete implementation of the function
    * declared in the tbox::Serializable abstract base class.
    */
   void putToDatabase( tbox::Pointer<tbox::Database> db);

   /**
    * This routine is a concrete implementation of the virtual function
    * in the base class BoundaryUtilityStrategy.  It reads DIRICHLET
    * boundary state values from the given database with the
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
    * Register a Vizamrai data writer so this class will write
    * plot files that may be postprocessed with the Vizamrai 
    * visualization tool.
    */
   void registerVizamraiDataWriter( 
      tbox::Pointer<appu::CartesianVizamraiDataWriter<NDIM> > viz_writer);

   /**
    * Register a VisIt data writer so this class will write
    * plot files that may be postprocessed with the VisIt 
    * visualization tool.
    */
#ifdef HAVE_HDF5
   void registerVisItDataWriter( 
      tbox::Pointer<appu::VisItDataWriter<NDIM> > viz_writer);
#endif

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
    * Private member function to check correctness of boundary data.
    */
   void checkBoundaryData(int btype,
                          const hier::Patch<NDIM>& patch,
                          const hier::IntVector<NDIM>& ghost_width_to_fill,
                          const tbox::Array<int>& scalar_bconds) const;
  
   /*
    * Object name used for error/warning reporting and as a label
    * for restart database entries.
    */
   string d_object_name;

   /*
    * tbox::Pointer to the grid geometry object used (Cartesian) to setup
    * initial data and to set physical boundary conditions.  
    */
   tbox::Pointer<geom::CartesianGridGeometry<NDIM> > d_grid_geometry;

   /*
    * Vizamrai data writer writes out variables used by this patch
    * strategy.
    */
   tbox::Pointer<appu::CartesianVizamraiDataWriter<NDIM> > d_vizamrai_writer;

#ifdef HAVE_HDF5
   tbox::Pointer<appu::VisItDataWriter<NDIM> > d_visit_writer;
#endif


   /*
    * Pointers to variables.  d_primitive_vars - [u]
    *                         d_function_eval  - [F(u)]
    */
   tbox::Pointer< pdat::CellVariable<NDIM,double> > d_primitive_vars;
   tbox::Pointer< pdat::CellVariable<NDIM,double> > d_function_eval;

   /*
    * Convection-diffusion equation constant vectors
    */
   double d_convection_coeff[NDIM];
   double d_diffusion_coeff;
   double d_source_coeff;

   /*
    *  Parameters for numerical method:
    *
    *    d_cfl ................ CFL condition for timestepping.
    *
    *    d_tolerance .......... Tolerance used for tagging cells - if 
    *                           value[N] > d_tolerance[n] (where N is 
    *                           between 0 and NEQU-1) cell is tagged.
    *
    *    d_nghosts ............ number of ghost cells for cell-centered
    *                           variables
    *
    */
   double d_cfl;
   double d_tolerance[NEQU];
   hier::IntVector<NDIM> d_nghosts;
   hier::IntVector<NDIM> d_zero_ghosts;

   /*
    * Indicator for problem type and initial conditions
    */
   string d_data_problem;
   int    d_data_problem_int;

   /*
    * Input for SPHERE problem
    */
   double d_radius;
   double d_center[NDIM];
   double d_val_inside[NEQU];
   double d_val_outside[NEQU];

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

#endif
