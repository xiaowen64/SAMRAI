//
// File:        MblkEuler.h
// Package:     SAMRAI application
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2132 $
// Modified:    $LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description: Numerical routines for single patch in linear advection ex.
//
 
#ifndef included_MblkEulerXD
#define included_MblkEulerXD

#include "SAMRAI_config.h"

#include "tbox/Array.h"
#include "tbox/Pointer.h"
#include "tbox/Serializable.h"
#include "tbox/Database.h"

#include "IntVector.h"
#include "Box.h"
#include "BoundaryBox.h"
#include "Patch.h"
#include "CellVariable.h"
#include "SideVariable.h"
#include "NodeVariable.h"
#include "BlockGridGeometry.h"

#include "BoundaryUtilityStrategy.h"
#include "GriddingAlgorithm.h"

#include <string>
using namespace std;
#define included_String

#include "TimeInterpolateOperator.h"
#include "VariableContext.h"
#include "VisItDataWriter.h"

#include "MblkGeometry.h"
#include "MultiblockRefineSchedule.h"
#include "MblkHyperbolicLevelIntegrator.h"
#include "MblkHyperbolicPatchStrategy.h"

// ----------------------------------------------------------------------

using namespace SAMRAI;

class MblkEuler : 
  public tbox::Serializable,
  public MblkHyperbolicPatchStrategy,
  public appu::BoundaryUtilityStrategy
{
public:

   //
   // the constructor and destructor
   //
   MblkEuler(const string& object_name,
          tbox::Pointer<tbox::Database> input_db,
          tbox::Array< tbox::Pointer<hier::GridGeometry<NDIM> > >& grid_geoms);
 
   ~MblkEuler();

   //
   // register with the framework
   //
   void registerModelVariables(MblkHyperbolicLevelIntegrator* integrator);


   //
   // set the patch initial conditions
   //
   void initializeDataOnPatch(hier::Patch<NDIM>& patch,
                              const double data_time,
                              const bool initial_time);
 
   //
   // Compute the stable time increment for patch using a CFL
   // condition and return the computed dt.
   //
   double computeStableDtOnPatch( hier::Patch<NDIM>& patch,
				  const bool initial_time,
				  const double dt_time );
 
   //
   // compute the state extrema for debugging
   //
   void testPatchExtrema( hier::Patch<NDIM>& patch, const char *pos );
 
   //
   // compute the fluxes and the initial update in a timestep
   //
   void computeFluxesOnPatch(hier::Patch<NDIM>& patch,
                             const double time,
                             const double dt);
 
   //
   // update the state (currently only for refluxing)
   //
   void conservativeDifferenceOnPatch(hier::Patch<NDIM>& patch,
                                      const double time,
                                      const double dt,
                                      bool at_syncronization);

   //
   // Tag cells for refinement using gradient detector.
   //
   void tagGradientDetectorCells(
      hier::Patch<NDIM>& patch,
      const double regrid_time,
      const bool initial_error,
      const int tag_indexindx,
      const bool uses_richardson_extrapolation_too);

   //
   //  The following routines:
   //
   //      postprocessRefine()
   //      setPhysicalBoundaryConditions()
   //
   //  are concrete implementations of functions declared in the
   //  RefinePatchStrategy abstract base class.
   //
   
   //
   // mark the zones to track what zones are being filled
   //
   void markPhysicalBoundaryConditions( hier::Patch<NDIM>& patch, 
					const double fill_time,
					const hier::IntVector<NDIM>& ghost_width_to_fill );

   //
   // set the data in the physical ghost zones
   //
   void setPhysicalBoundaryConditions(hier::Patch<NDIM>& patch,
                                      const double fill_time,
                                      const hier::IntVector<NDIM>&
                                      ghost_width_to_fill);

   //
   // Refine operations for cell data.
   //
   void preprocessRefine( hier::Patch<NDIM>& fine,
                          const hier::Patch<NDIM>& coarse,
                          const hier::Box<NDIM>& fine_box,
                          const hier::IntVector<NDIM>& ratio);

   void postprocessRefine( hier::Patch<NDIM>& fine,
			   const hier::Patch<NDIM>& coarse,
			   const hier::Box<NDIM>& fine_box,
			   const hier::IntVector<NDIM>& ratio);

   //
   // Coarsen operations for cell data.
   //
   void preprocessCoarsen( hier::Patch<NDIM>& coarse,
			   const hier::Patch<NDIM>& fine,
			   const hier::Box<NDIM>& coarse_box,
			   const hier::IntVector<NDIM>& ratio);

   void postprocessCoarsen( hier::Patch<NDIM>& coarse,
			    const hier::Patch<NDIM>& fine,
			    const hier::Box<NDIM>& coarse_box,
			    const hier::IntVector<NDIM>& ratio);

   /**
    * Fill the singularity conditions for the multi-block case
    */
   void fillSingularityBoundaryConditions(
      hier::Patch<NDIM>& patch,
      tbox::List<xfer::MultiblockRefineSchedule<NDIM>::SingularityPatch>&
      sing_patches,
      const double fill_time,
      const hier::Box<NDIM>& fill_box,
      const hier::BoundaryBox<NDIM>& bbox);

   /**
    * Build mapped grid on patch
    */
   void setMappedGridOnPatch(const hier::Patch<NDIM>& patch, 
                             const int level_number,
                             const int block_number);   

   //
   // build the volume on a mapped grid
   //
   void setVolumeOnPatch( const hier::Patch<NDIM>& patch );   


   /**
    * Write state of MblkEuler object to the given database for restart.
    *
    * This routine is a concrete implementation of the function
    * declared in the tbox::Serializable abstract base class.
    */
   void putToDatabase(tbox::Pointer<tbox::Database> db);


   hier::IntVector<NDIM> getMultiblockRefineOpStencilWidth() const;
   hier::IntVector<NDIM> getMultiblockCoarsenOpStencilWidth();


#ifdef HAVE_HDF5
   /**
    * Register a VisIt data writer so this class will write
    * plot files that may be postprocessed with the VisIt 
    * visualization tool.
    */
   void registerVisItDataWriter( 
      tbox::Pointer<appu::VisItDataWriter<NDIM> > viz_writer);
#endif

   /**
    * Print all data members for MblkEuler class.
    */
   void printClassData(ostream& os) const;

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
   void getFromInput(tbox::Pointer<tbox::Database> db,
                     bool is_from_restart);

   void getFromRestart();

   /*
    * Private member function to check correctness of boundary data.
    */
   void checkBoundaryData(int btype,
                          const hier::Patch<NDIM>& patch,
                          const hier::IntVector<NDIM>& ghost_width_to_fill,
                          const tbox::Array<int>& scalar_bconds) const;


   /*
    * The object name is used for error/warning reporting and also as a 
    * string label for restart database entries. 
    */
   string d_object_name;

   /*
    * We cache pointers to the grid geometry and Vizamrai data writer
    * object to set up initial data, set physical boundary conditions,
    * and register plot variables.
    */
   tbox::Array< tbox::Pointer<hier::GridGeometry<NDIM> > >
   d_grid_geometries;
#ifdef HAVE_HDF5
   tbox::Pointer<appu::VisItDataWriter<NDIM> > d_visit_writer;
#endif

   //
   // Data items used for nonuniform load balance, if used.
   //
   tbox::Pointer< pdat::CellVariable<NDIM,double> > d_workload_variable;
   int d_workload_data_id;
   bool d_use_nonuniform_workload;

   //
   // =========================== State and Variable<NDIM> definitions (private) ============================
   //

   //
   // tbox::Pointer to state variable vector - [state]
   //
   int d_nState;  // depth of the state vector
   tbox::Pointer< pdat::CellVariable<NDIM,double> > d_state;

   //
   // tbox::Pointer to cell volume - [v]
   //
   tbox::Pointer< pdat::CellVariable<NDIM,double> > d_vol;

   //
   // tbox::Pointer to flux variable vector  - [F]
   //
   tbox::Pointer< pdat::SideVariable<NDIM,double> > d_flux;

   //
   // tbox::Pointer to grid - [xyz]
   //
   tbox::Pointer< pdat::NodeVariable<NDIM,double> > d_xyz;
   int d_xyz_id;

   //
   // ======================================= Initial Conditions (private) ============================
   //

   /// center of the sphere or revolution origin
   double d_center[NDIM];

   /// revolution axis
   double d_axis[NDIM];

   /// revolution radius and pos on axis of radius
   tbox::Array< tbox::Array<double> > d_rev_rad;
   tbox::Array< tbox::Array<double> > d_rev_axis;

   ///
   /// Rayleigh Taylor Shock tube experiments
   ///
   double d_dt_ampl;
   tbox::Array<double> d_amn;
   tbox::Array<double> d_m_mode;
   tbox::Array<double> d_n_mode;
   tbox::Array<double> d_phiy;
   tbox::Array<double> d_phiz;

   ///
   /// input for all the geometries
   ///

   // 
   // linear advection velocity vector for unit test
   //
   int d_advection_test;      // run the linear advection unit test 
   int d_advection_vel_type;  // type of velocity to use
   double d_advection_velocity[NDIM];

   //
   // sizing of zonal, flux, and nodal ghosts
   //
   hier::IntVector<NDIM> d_nghosts;
   hier::IntVector<NDIM> d_fluxghosts;
   hier::IntVector<NDIM> d_nodeghosts;

   //
   // Indicator for problem type and initial conditions
   //
   string d_data_problem;

   //
   // region initialization inputs
   //
   int d_number_of_regions;
   tbox::Array<double> d_front_position;

   //
   // array of initial conditions and their names [region][state]
   //
   double **d_state_ic;
   tbox::Array<string> d_state_names;

   //
   // This class stores geometry information used for constructing the
   // mapped multiblock hierarchy
   //
   MblkGeometry* d_mblk_geometry;

   /// the bound on the index space for the current block
   int d_dom_current_bounds[6];

   /// the number of boxes needed to describe the index space for the current block
   int d_dom_current_nboxes;

   /// the blocks bounding the current patch
   int d_dom_local_blocks[6];


   //
   // ======================================= Refinement Data (private) ============================
   //

   tbox::Array<string> d_refinement_criteria;

   /// history variable gradient tagging tolerance
   tbox::Array< tbox::Array<double> > d_state_grad_tol;
   tbox::Array<string> d_state_grad_names;
   tbox::Array<int>    d_state_grad_id;

   //
   // ======================================= Boundary Conditions (private) ============================
   //

   /// factors for the boundary conditions
   tbox::Array<int> d_wall_factors;

   //
   // Operators to be used with BlockGridGeometry
   //
   tbox::Pointer<xfer::TimeInterpolateOperator<NDIM> > d_cell_time_interp_op;
};

#endif
