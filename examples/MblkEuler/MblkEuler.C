//
// File:        MblkEuler.C
// Description: Template for a multiblock AMR Euler code
//
#include "MblkEuler.h" 

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

#include "tbox/ArenaManager.h"
#include "tbox/InputManager.h"
#include "tbox/MathUtilities.h"
#include "tbox/PIO.h"
#include "tbox/RestartManager.h"
#include "tbox/Utilities.h"

#include "BoxArray.h"
#include "CellData.h"
#include "CellIterator.h"
#include "CellDoubleLinearTimeInterpolateOp.h"
#include "CoarsenOperator.h"
#include "SideData.h"
#include "BlockGridGeometry.h"
#include "BlockPatchGeometry.h"
#include "Index.h"
#include "NodeData.h"
#include "NodeDoubleLinearTimeInterpolateOp.h"
#include "RefineOperator.h"
#include "TimeInterpolateOperator.h"
#include "VariableDatabase.h"

#include "CartesianBoundaryDefines.h"  // defines flags used by the boundingBox data

// Number of ghosts cells used for each variable quantity
#define CELLG           (1)
#define FLUXG           (0)
#define NODEG           (1)

// defines for cell tagging routines
#ifndef TRUE
#define TRUE (1)
#endif
#ifndef FALSE
#define FALSE (0)
#endif

// Version of MblkEuler restart file data
#define MBLKEULER_VERSION (3)

//
// some extra defines for C code
//
#define real8 double
#define POLY3(i,j,k,imin,jmin,kmin,nx,nxny) ( (i-imin) + (j-jmin)*(nx) + (k-kmin)*(nxny) )
#define MAX( a, b ) ( a > b ? a : b ) 
#define MIN( a, b ) ( a < b ? a : b ) 


//
// inline geometry functions
//
#include "GeomUtilsAMR.h"

// ================================= MblkEuler::Initialization =============================


/*
*************************************************************************
*                                                                       *
* The constructor for MblkEuler class sets data members to defualt values, *
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

MblkEuler::MblkEuler(
   const string& object_name,
   tbox::Pointer<tbox::Database> input_db,
   tbox::Array< tbox::Pointer<hier::GridGeometry<NDIM> > >& grid_geoms)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(!input_db.isNull());
   TBOX_ASSERT(grid_geoms.getSize() > 0);
#endif

   d_object_name = object_name;
   tbox::RestartManager::getManager()->registerRestartItem(d_object_name, this);

   d_grid_geometries = grid_geoms;

   //
   // Setup MblkGeometry object to manage construction of mapped grids
   //
   d_mblk_geometry = new MblkGeometry("MblkGeometry",
                                      input_db,
                                      grid_geoms.getSize());
   
   d_use_nonuniform_workload = false;

   tbox::Pointer<tbox::Database> mbe_db = input_db->getDatabase("MblkEuler");

   //
   // Default parameters for the numerical method.
   //
   d_nghosts    = hier::IntVector<NDIM>(CELLG);
   d_fluxghosts = hier::IntVector<NDIM>(FLUXG);
   d_nodeghosts = hier::IntVector<NDIM>(NODEG);
   
   //
   // zero out initial condition variables
   //
   tbox::MathUtilities<double>::setArrayToSignalingNaN( d_center, NDIM );
   tbox::MathUtilities<double>::setArrayToSignalingNaN( d_axis,   NDIM );

   d_front_position.resizeArray(0);

   d_rev_rad.resizeArray(0);
   d_rev_axis.resizeArray(0);

   d_amn.resizeArray(0);
   d_m_mode.resizeArray(0);
   d_n_mode.resizeArray(0);
   d_phiy.resizeArray(0);
   d_phiz.resizeArray(0);

   //
   // Initialize object with data read from given input/restart databases.
   //
   bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();
   if (is_from_restart){
      getFromRestart();
   }
   getFromInput(input_db, is_from_restart);

   //
   // quantities that define state of linear advection problem.
   //
   d_state = new pdat::CellVariable<NDIM,double>( "state", d_nState );
   d_vol   = new pdat::CellVariable<NDIM,double>( "vol",  1 );
   d_flux  = new pdat::SideVariable<NDIM,double>( "flux", d_nState );
   d_xyz   = new pdat::NodeVariable<NDIM,double>( "xyz",  NDIM );

   //
   // drop the region layout as a table
   //
   tbox::plog << "region layout follows:" << endl;

   tbox::plog << "field";
   for ( int ir = 0 ; ir < d_number_of_regions ; ir++ )
      tbox::plog << "\t" << ir;
   tbox::plog << endl;

   for ( int ii = 0 ;  ii < d_number_of_regions*10; ii++ ) 
      tbox::plog << "-";
   tbox::plog << endl;
      
   for ( int istate = 0 ;  istate < d_nState; istate++ ) {
      tbox::plog << d_state_names[istate];
      for ( int ir = 0 ; ir < d_number_of_regions ; ir++ )
	 tbox::plog << "\t" << d_state_ic[ir][istate];
      tbox::plog << endl;
   }
}

/*
*************************************************************************
*                                                                       *
* Empty destructor for MblkEuler class.                                    *
*                                                                       *
*************************************************************************
*/

MblkEuler::~MblkEuler() {
   if (d_mblk_geometry) delete d_mblk_geometry;
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

void MblkEuler::registerModelVariables(
   MblkHyperbolicLevelIntegrator* integrator)
{

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(integrator != (MblkHyperbolicLevelIntegrator*)NULL);
#endif

   //
   // zonal data and its fluxes
   //
   d_cell_time_interp_op = 
      new pdat::CellDoubleLinearTimeInterpolateOp<NDIM>();

   integrator->registerVariable(d_state, d_nghosts,
                                MblkHyperbolicLevelIntegrator::TIME_DEP,
                                NULL,
                                NULL,
                                d_cell_time_interp_op);

   integrator->registerVariable(d_vol, d_nghosts,
                                MblkHyperbolicLevelIntegrator::TIME_DEP,
                                NULL,
                                NULL,
                                d_cell_time_interp_op);  

   integrator->registerVariable(d_flux, d_fluxghosts,
				MblkHyperbolicLevelIntegrator::FLUX ,
                                NULL );

   //
   // The nodal position data
   //
   tbox::Pointer<xfer::TimeInterpolateOperator<NDIM> > node_time_interp_op = 
      new pdat::NodeDoubleLinearTimeInterpolateOp<NDIM>();

   integrator->registerVariable(d_xyz, d_nodeghosts,
                                MblkHyperbolicLevelIntegrator::TIME_DEP,
                                NULL,
                                NULL,
                                node_time_interp_op);
   
   hier::VariableDatabase<NDIM>* vardb = hier::VariableDatabase<NDIM>::getDatabase();



   if (d_visit_writer.isNull()) {
      TBOX_WARNING(d_object_name << ": registerModelVariables()"
                   << "\nVisit data writer was"
                   << "\nregistered.  Consequently, no plot data will"
                   << "\nbe written." << endl);
   }

   for ( int n = 0 ; n < d_nState; n++ ) {
      string vname = d_state_names[n];
      d_visit_writer->registerPlotQuantity( vname, "SCALAR",
					    vardb->mapVariableAndContextToIndex( d_state, integrator->getPlotContext() ),
					    n );
   }

   d_visit_writer->registerPlotQuantity( "vol", "SCALAR",
					 vardb->mapVariableAndContextToIndex( d_vol, integrator->getPlotContext() ) );
   
   d_visit_writer->registerNodeCoordinates( vardb->mapVariableAndContextToIndex( d_xyz, integrator->getPlotContext() ) );
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
void MblkEuler::initializeDataOnPatch(hier::Patch<NDIM>& patch,
                                   const double data_time,
                                   const bool initial_time)
{
   //(void) data_time;

   //
   // Build the mapped grid on the patch.
   //
   int block_number = getBlockNumber();
   int level_number = patch.getPatchLevelNumber();
   setMappedGridOnPatch(patch, level_number, block_number);

   if (initial_time) {
   
      tbox::Pointer< pdat::CellData<NDIM,double> > state = patch.getPatchData( d_state, getDataContext() );
      tbox::Pointer< pdat::CellData<NDIM,double> > vol   = patch.getPatchData( d_vol,   getDataContext() );
      tbox::Pointer< pdat::NodeData<NDIM,double> > xyz   = patch.getPatchData( d_xyz,   getDataContext() );

#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT( !state.isNull() );
      TBOX_ASSERT( !vol.isNull()   );
      TBOX_ASSERT( !xyz.isNull()   );
      TBOX_ASSERT( state->getGhostCellWidth() == vol->getGhostCellWidth() );
#endif 

      hier::IntVector<NDIM> state_ghosts = state->getGhostCellWidth();
      hier::IntVector<NDIM> xyz_ghosts   = xyz->getGhostCellWidth();

      const hier::Index<NDIM> ifirst = patch.getBox().lower();
      const hier::Index<NDIM> ilast  = patch.getBox().upper();

      int imin = ifirst(0) - state_ghosts(0);
      int imax = ilast(0)  + state_ghosts(0);
      int jmin = ifirst(1) - state_ghosts(1);
      int jmax = ilast(1)  + state_ghosts(1);
      int kmin = ifirst(2) - state_ghosts(2);
      int kmax = ilast(2)  + state_ghosts(2);
      int nx   = imax - imin + 1;
      int ny   = jmax - jmin + 1;
      int nxny = nx*ny;
      
      int nd_imin = ifirst(0)     - xyz_ghosts(0);
      int nd_imax = ilast(0)  + 1 + xyz_ghosts(0);
      int nd_jmin = ifirst(1)     - xyz_ghosts(1);
      int nd_jmax = ilast(1)  + 1 + xyz_ghosts(1);
      int nd_kmin = ifirst(2)     - xyz_ghosts(2);
      int nd_nx   = nd_imax - nd_imin + 1;
      int nd_ny   = nd_jmax - nd_jmin + 1;
      int nd_nxny = nd_nx*nd_ny;
      
      //
      // get the pointers
      // 
      double *cvol = vol->getPointer();
      
      double *x = xyz->getPointer(0);
      double *y = xyz->getPointer(1);
      double *z = xyz->getPointer(2);
      
      hier::IntVector<NDIM> ghost_cells  = state->getGhostCellWidth();
      pdat::CellData<NDIM,double> elemCoords(      patch.getBox(), 3, ghost_cells );  // storage for the average of the element coordinates
      pdat::CellData<NDIM,int>    region_ids_data( patch.getBox(), 1, ghost_cells );  // storage for the slopes

      int *region_ids = region_ids_data.getPointer();
      double *xc = elemCoords.getPointer(0);
      double *yc = elemCoords.getPointer(1);
      double *zc = elemCoords.getPointer(2);

      //
      //  ---------------- compute the element coordinates
      //
      for ( int k = kmin; k <= kmax; k++ ) {
	 for ( int j = jmin; j <= jmax; j++ ) {
	    for ( int i = imin; i <= imax; i++ ) {
	       int ind  = POLY3(i,j,k,imin,jmin,kmin,nx,nxny);

	       int n1 = POLY3(i,j,k,nd_imin,nd_jmin,nd_kmin,nd_nx,nd_nxny);
	       int n2 = n1 + 1;
	       int n3 = n1 + 1 + nd_nx;
	       int n4 = n1 + nd_nx;
	       
	       int n5 = n1 + nd_nxny;
	       int n6 = n1 + nd_nxny + 1;
	       int n7 = n1 + nd_nxny + 1 + nd_nx;
	       int n8 = n1 + nd_nxny + nd_nx;
	       
	       double xavg = 0.125*( x[n1] + x[n2] + x[n3] + x[n4] +
				    x[n5] + x[n6] + x[n7] + x[n8] );
	       
	       double yavg = 0.125*( y[n1] + y[n2] + y[n3] + y[n4] +
				    y[n5] + y[n6] + y[n7] + y[n8] );
	       
	       double zavg = 0.125*( z[n1] + z[n2] + z[n3] + z[n4] +
				    z[n5] + z[n6] + z[n7] + z[n8] );
	       xc[ind] = xavg;
	       yc[ind] = yavg;
	       zc[ind] = zavg;
	    }
	 }
      }


      //
      //  ---------------- compute the element volume
      //
      for ( int k = ifirst(2); k <= ilast(2); k++ ) {
	 for ( int j = ifirst(1); j <= ilast(1); j++ ) {
	    for ( int i = ifirst(0); i <= ilast(0); i++ ) {
	       
	       int cind = POLY3(i,j,k,imin,jmin,kmin,nx,nxny);
	       
	       int n1 = POLY3(i,j,k,nd_imin,nd_jmin,nd_kmin,nd_nx,nd_nxny);
	       int n2 = n1 + 1;
	       int n3 = n1 + 1 + nd_nx;
	       int n4 = n1 + nd_nx;
	       
	       int n5 = n1 + nd_nxny;
	       int n6 = n1 + nd_nxny + 1;
	       int n7 = n1 + nd_nxny + 1 + nd_nx;
	       int n8 = n1 + nd_nxny + nd_nx;
	       
	       double lvol = UpwindVolume( x[n1], x[n2], x[n3], x[n4],
					   x[n5], x[n6], x[n7], x[n8],
					  
					   y[n1], y[n2], y[n3], y[n4],
					   y[n5], y[n6], y[n7], y[n8],
					   
					   z[n1], z[n2], z[n3], z[n4],
					   z[n5], z[n6], z[n7], z[n8] );
	       
	       cvol[cind] = lvol;
	    }
	 }
      }

      //
      //  ---------------- process the different initialization regions
      //
      if ( d_data_problem == "REVOLUTION" ) {

	for ( int m = 0 ; m < d_number_of_regions; m++ ) {  // loop over the regions and shape in data

	  tbox::Array<double> &lrad  = d_rev_rad[m];
	  tbox::Array<double> &laxis = d_rev_axis[m];
	  int naxis = laxis.getSize();

	  for ( int k = kmin ; k <= kmax ; k++ ) {
	     for ( int j = jmin ; j <= jmax ; j++ ) {
		for ( int i = imin ; i <= imax ; i++ ) {
		   int ind  = POLY3(i,j,k,imin,jmin,kmin,nx,nxny);
		   
		   double x2 = zc[ind] - d_center[2];
		   double x1 = yc[ind] - d_center[1];
		   double x0 = xc[ind] - d_center[0];
		   
		   double Dz = x0*d_axis[0] + x1*d_axis[1] + x2*d_axis[2];
		   double r2 = x0*x0 + x1*x1 + x2*x2;
		   double Dr = sqrt(r2 - Dz*Dz);
		   
		   if ( laxis[0] <= Dz && Dz <= laxis[naxis-1] ) {  // test to see if we are contained
		      
		      int lpos = 0;
		      while ( Dz > laxis[lpos] )
			 lpos++;
		      
		      double a  = (Dz - laxis[lpos-1]) / (laxis[lpos] - laxis[lpos-1] );
		      double lr = (a)*lrad[lpos] + (1.0-a)*lrad[lpos-1];
		      
		      if ( Dr <= lr ) { // if we are within the radius set the region id
			 region_ids[ind] = m;
		      }
		      
		   } // laxis
		   
		}
	     }
	  } // k
	  
	} // end of region loop
      }

      //
      //  the spherical initialization
      //
      else if ( d_data_problem == "SPHERE" ) {
	 double zero = 0.0;
	 double angle = 0.0;
	 double *front     = d_front_position.getPointer();
	 for ( int k = kmin ; k <= kmax ; k++ ) {
	    for ( int j = jmin ; j <= jmax ; j++ ) {
	       for ( int i = imin ; i <= imax ; i++ ) {
		  int ind  = POLY3(i,j,k,imin,jmin,kmin,nx,nxny);
		  
		  double x2 = zc[ind] - d_center[2];
		  double x1 = yc[ind] - d_center[1];
		  double x0 = xc[ind] - d_center[0];
		  
		  if (x1 == zero && x0 == zero) {
		     angle = zero;
		  } else {
		     angle = atan2(x1,x0);
		  }
		  
		  double rad2 = sqrt( x0*x0 + x1*x1);
		  // double phi  = atan2(rad2,x2);  
		  double rad3 = sqrt(rad2*rad2 +x2*x2);
		  
		  int ifr = 0;  // find the region we draw from (ifr=0 is always origin)
		  while ( rad3 > front[ifr+1] ) {
		     ifr++;
		  }
		  region_ids[ind] = ifr;
	       }
	    }
	 }
      }
      

      //
      //  the planar initialization
      //
      else if ( ( d_data_problem == "PIECEWISE_CONSTANT_X" )
		|| ( d_data_problem == "PIECEWISE_CONSTANT_Y" )
		|| ( d_data_problem == "PIECEWISE_CONSTANT_Z" )) {
	 
	 double *front     = d_front_position.getPointer();
	 double *xx = xc;
	 if ( d_data_problem == "PIECEWISE_CONSTANT_Y" ) {
	    xx = yc;
	 }
	 if ( d_data_problem == "PIECEWISE_CONSTANT_Z" ) {
	    xx = zc;
	 }
	 
	 for ( int k = kmin ; k <= kmax ; k++ ) {
	    for ( int j = jmin ; j <= jmax ; j++ ) {
	       for ( int i = imin ; i <= imax ; i++ ) {
		  int ind  = POLY3(i,j,k,imin,jmin,kmin,nx,nxny);
		  
		  int ifr  = 0;  // the geion we are in
		  while (xx[ind] > front[ifr+1]) {
		     ifr++;
		  }
		  region_ids[ind] = ifr;
	       }
	    }
	 }
      }

      //
      //  the planar initialization
      //
      else if ( ( d_data_problem == "RT_SHOCK_TUBE" ) ) {
	
	 double *front    = d_front_position.getPointer();
	 double shock_pos = front[1]; // the shock front
	 double front_pos = front[2]; // the sinusoidal perturbation between the two fluids

	 double dt_ampl  = d_dt_ampl;
	 int nmodes      = d_amn.getSize();
	 double *amn     = d_amn.getPointer();
	 double *n_mode  = d_n_mode.getPointer();
	 double *m_mode  = d_m_mode.getPointer();
	 double *phiy    = d_phiy.getPointer();
	 double *phiz    = d_phiz.getPointer();

	 // ... this is a cartesian problem by definition
	 const double *xdlo = 0;  // d_cart_xlo[0][0];
	 const double *xdhi = 0;  // d_cart_xhi[0][0];
	 TBOX_ASSERT(0);  // xdlo, xdhi wrong

	 double l_y      = xdhi[1] - xdlo[1]; // width of the domain in y and z
	 double l_z      = xdhi[2] - xdlo[2];
	 double lpi      = 3.14159265358979310862446895044;

	 for ( int k = kmin ; k <= kmax ; k++ ) {
	    for ( int j = jmin ; j <= jmax ; j++ ) {
	       for ( int i = imin ; i <= imax ; i++ ) {
		  int ind  = POLY3(i,j,k,imin,jmin,kmin,nx,nxny);
		  
		  double lpert = 0.0;
		  double ly   = yc[ind];
		  double lz   = zc[ind];
		  for ( int m = 0; m < nmodes; m++ ) {
		     double lphiy = 2.0*lpi*n_mode[m]*ly/l_y + phiy[m];
		     double lphiz = 2.0*lpi*m_mode[m]*lz/l_z + phiz[m];
		     double cy   = cos(lphiy);
		     double cz   = cos(lphiz);
		     lpert += amn[m]*cy*cz;
		  }
		  lpert *= dt_ampl;
		  
		  int ifr  = 0;
		  if ( xc[ind] > shock_pos )
		     ifr = 1;
		  if ( xc[ind] > front_pos + lpert ) 
		     ifr = 2;

		  region_ids[ind] = ifr;
	       }
	    }
	 }
      }

      else if ( d_data_problem == "BLIP" ) {
	TBOX_ASSERT(0);
      }

      //
      // ---------------- state vector
      //
      int depth = state->getDepth();

      for ( int idepth = 0 ; idepth < depth ; idepth++ ) {
	 double *psi = state->getPointer(idepth);

	 for ( int k = ifirst(2); k <= ilast(2); k++ ) {
	    for ( int j = ifirst(1); j <= ilast(1); j++ ) {
	       for ( int i = ifirst(0); i <= ilast(0); i++ ) {
		  
		  int cind = POLY3(i,j,k,imin,jmin,kmin,nx,nxny);
#if 0
		  // psi[cind] = (double) i;  // set values to coordinate ids as debug tool
		  // psi[cind] = (double) j;
		  // psi[cind] = (double) k;
#endif 
#if 1
		  int ireg  = region_ids[cind];
		  psi[cind] = d_state_ic[ireg][idepth];
#endif
	       }
	    }
	 }

      } // end of depth loop

   } // end of initial time if test

   if (d_use_nonuniform_workload) {
      if (!patch.checkAllocated(d_workload_data_id)) {
         patch.allocatePatchData(d_workload_data_id);
      }
      tbox::Pointer< pdat::CellData<NDIM,double> > workload_data =
         patch.getPatchData(d_workload_data_id);
      workload_data->fillAll(1.0);
   }

}


// ================================= MblkEuler::Integration =============================


/*
*************************************************************************
*                                                                       *
* Compute stable time increment for patch.  Return this value.          *
*                                                                       *
*************************************************************************
*/

double MblkEuler::computeStableDtOnPatch( hier::Patch<NDIM>& patch,
					  const bool initial_time,
					  const double dt_time)
{
   (void) initial_time;
   (void) dt_time;

   //
   // Build the mapped grid on the patch.
   //
   int block_number = getBlockNumber();
   int level_number =  patch.getPatchLevelNumber();
   setMappedGridOnPatch(patch, level_number, block_number);

   const hier::Index<NDIM> ifirst=patch.getBox().lower();
   const hier::Index<NDIM> ilast =patch.getBox().upper();

   tbox::plog << "--------------------- start stableDtOnPatch on patch (" << level_number << ")" << endl;
   tbox::plog << "level = " << level_number   << endl;
   tbox::plog << "box   = " << patch.getBox() << endl;

   //
   // get the cell data and their bounds
   // 
   tbox::Pointer< pdat::CellData<NDIM,double> > state = patch.getPatchData( d_state, getDataContext() );
   tbox::Pointer< pdat::NodeData<NDIM,double> > xyz   = patch.getPatchData( d_xyz,   getDataContext() );
   hier::IntVector<NDIM> state_ghosts = state->getGhostCellWidth();
   hier::IntVector<NDIM> xyz_ghosts   = xyz->getGhostCellWidth();

   pdat::CellData<NDIM,double> Aii( patch.getBox(), 9, 0 );

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!state.isNull());
   TBOX_ASSERT(!xyz.isNull());
#endif 

   int imin = ifirst(0);
   int imax = ilast(0);
   int jmin = ifirst(1);
   int jmax = ilast(1);
   int kmin = ifirst(2);
   int nx   = imax - imin + 1;
   int ny   = jmax - jmin + 1;
   int nxny = nx*ny;

   int nd_imin = ifirst(0)     - xyz_ghosts(0);
   int nd_imax = ilast(0)  + 1 + xyz_ghosts(0);
   int nd_jmin = ifirst(1)     - xyz_ghosts(1);
   int nd_jmax = ilast(1)  + 1 + xyz_ghosts(1);
   int nd_kmin = ifirst(2)     - xyz_ghosts(2);
   int nd_nx   = nd_imax - nd_imin + 1;
   int nd_ny   = nd_jmax - nd_jmin + 1;
   int nd_nxny = nd_nx*nd_ny;

   double *x = xyz->getPointer(0);
   double *y = xyz->getPointer(1);
   double *z = xyz->getPointer(2);
   
   int pind = 0;
   double *a11 = Aii.getPointer(pind); pind++;
   double *a12 = Aii.getPointer(pind); pind++;
   double *a13 = Aii.getPointer(pind); pind++;
   double *a21 = Aii.getPointer(pind); pind++;
   double *a22 = Aii.getPointer(pind); pind++;
   double *a23 = Aii.getPointer(pind); pind++;
   double *a31 = Aii.getPointer(pind); pind++;
   double *a32 = Aii.getPointer(pind); pind++;
   double *a33 = Aii.getPointer(pind); pind++;

   //
   // compute direction cosines 
   //
   for ( int k = ifirst(2); k <= ilast(2); k++ ) {
      for ( int j = ifirst(1); j <= ilast(1); j++ ) {
         for ( int i = ifirst(0); i <= ilast(0); i++ ) {

            int ind  = POLY3(i,j,k,imin,jmin,kmin,nx,nxny);

            int n1 = POLY3(i,j,k,nd_imin,nd_jmin,nd_kmin,nd_nx,nd_nxny);
            int n2 = n1 + 1;
            int n3 = n1 + 1 + nd_nx;
            int n4 = n1 + nd_nx;

            int n5 = n1 + nd_nxny;
            int n6 = n1 + nd_nxny + 1;
            int n7 = n1 + nd_nxny + 1 + nd_nx;
            int n8 = n1 + nd_nxny + nd_nx;

	    // ------------------------------------------------ x
	    
	    double x1 = 0.25*( x[n1] + x[n4] + x[n5] + x[n8] );  // xi
	    double x2 = 0.25*( x[n2] + x[n3] + x[n6] + x[n7] );
	    
	    double x3 = 0.25*( x[n1] + x[n2] + x[n5] + x[n6] );  // eta
	    double x4 = 0.25*( x[n3] + x[n4] + x[n7] + x[n8] );
	    
	    double x5 = 0.25*( x[n1] + x[n2] + x[n3] + x[n4] );  // zeta
	    double x6 = 0.25*( x[n5] + x[n6] + x[n7] + x[n8] );
	    
	    // ------------------------------------------------ y
	    
	    double y1 = 0.25*( y[n1] + y[n4] + y[n5] + y[n8] );
	    double y2 = 0.25*( y[n2] + y[n3] + y[n6] + y[n7] );
	    
	    double y3 = 0.25*( y[n1] + y[n2] + y[n5] + y[n6] );
	    double y4 = 0.25*( y[n3] + y[n4] + y[n7] + y[n8] );
	    
	    double y5 = 0.25*( y[n1] + y[n2] + y[n3] + y[n4] );
	    double y6 = 0.25*( y[n5] + y[n6] + y[n7] + y[n8] );
	    
	    // ------------------------------------------------ z
	    
	    double z1 = 0.25*( z[n1] + z[n4] + z[n5] + z[n8] );
	    double z2 = 0.25*( z[n2] + z[n3] + z[n6] + z[n7] );
	    
	    double z3 = 0.25*( z[n1] + z[n2] + z[n5] + z[n6] );
	    double z4 = 0.25*( z[n3] + z[n4] + z[n7] + z[n8] );
	    
	    double z5 = 0.25*( z[n1] + z[n2] + z[n3] + z[n4] );
	    double z6 = 0.25*( z[n5] + z[n6] + z[n7] + z[n8] );
	    
	    //
	    // the components of the matrices that we want to invert
	    //
	    double dx_xi   = 0.5*(x2 - x1);
	    double dy_xi   = 0.5*(y2 - y1);
	    double dz_xi   = 0.5*(z2 - z1);
	    
	    double dx_eta  = 0.5*(x4 - x3);
	    double dy_eta  = 0.5*(y4 - y3);
	    double dz_eta  = 0.5*(z4 - z3);
	    
	    double dx_zeta = 0.5*(x6 - x5);
	    double dy_zeta = 0.5*(y6 - y5);
	    double dz_zeta = 0.5*(z6 - z5);
	    
	    //
	    // invert M = dx/dxi to find the matrix needed to convert 
	    // displacements in x into displacements in xi (dxi/dx)
	    // via kramer's rule
	    //
	    double detM = ( dx_xi   * dy_eta  * dz_zeta +
			   dx_eta  * dy_zeta * dz_xi   +
			   dx_zeta * dy_xi   * dz_eta  -
			   dx_zeta * dy_eta  * dz_xi   -
			   dx_eta  * dy_xi   * dz_zeta -
			   dx_xi   * dy_zeta * dz_eta  );
	    
	    double detB11 = dy_eta  * dz_zeta - dy_zeta * dz_eta;
	    double detB21 = dy_zeta * dz_xi   - dy_xi   * dz_zeta;
	    double detB31 = dy_xi   * dz_eta  - dy_eta  * dz_xi;
	    
	    double detB12 = dx_zeta * dz_eta  - dx_eta  * dz_zeta;
	    double detB22 = dx_xi   * dz_zeta - dx_zeta * dz_xi;
	    double detB32 = dx_eta  * dz_xi   - dx_xi   * dz_eta;
	    
	    double detB13 = dx_eta  * dy_zeta - dx_zeta * dy_eta;
	    double detB23 = dx_zeta * dy_xi   - dx_xi   * dy_zeta;
	    double detB33 = dx_xi   * dy_eta  - dx_eta  * dy_xi;
	    
	    a11[ind] = detB11/detM;
	    a21[ind] = detB21/detM;
	    a31[ind] = detB31/detM;
	    
	    a12[ind] = detB12/detM;
	    a22[ind] = detB22/detM;
	    a32[ind] = detB32/detM;
	    
	    a13[ind] = detB13/detM;
	    a23[ind] = detB23/detM;
	    a33[ind] = detB33/detM;
	 }
      }
   }

   //
   // print out patch extrema
   //
   testPatchExtrema( patch, "before timestep evaluation" );

   //
   // the timestep evaluation
   //
   double stabdt = 1.e20;;
   if ( d_advection_test ) {
      // ------------------------------------- for linear advection
      double u = d_advection_velocity[0];
      double v = d_advection_velocity[1];
      double w = d_advection_velocity[2];
      
      for ( int k = ifirst(2); k <= ilast(2); k++ ) {
	 for ( int j = ifirst(1); j <= ilast(1); j++ ) {
	    for ( int i = ifirst(0); i <= ilast(0); i++ ) {
	       int ind  = POLY3(i,j,k,imin,jmin,kmin,nx,nxny);
	       
	       double uxi   = a11[ind]*u + a12[ind]*v + a13[ind]*w;  // max parametric signal speed
	       double ueta  = a21[ind]*u + a22[ind]*v + a23[ind]*w;
	       double uzeta = a31[ind]*u + a32[ind]*v + a33[ind]*w;
	       
	       uxi   = fabs( uxi   );
	       ueta  = fabs( ueta  );
	       uzeta = fabs( uzeta );
	       
	       double dtxi   = 2.0 / MAX( 1.e-80, uxi   );
	       double dteta  = 2.0 / MAX( 1.e-80, ueta  );
	       double dtzeta = 2.0 / MAX( 1.e-80, uzeta );
	       
	       double ldelt = MIN( dtxi, MIN( dteta, dtzeta ) );
	       stabdt       = MIN( stabdt, ldelt );
	       
	    }
	 }
      }
   } 
   else {
      TBOX_ASSERT( d_advection_test );
   }

   //
   // process the timestep constraints
   // 
   double dt_fixed = -1;
   double returned_dt = stabdt;
   if ( dt_fixed > 0.0 )
      returned_dt = dt_fixed;

   tbox::plog << "stabdt      = " << stabdt      << endl;
   tbox::plog << "dt_fixed    = " << dt_fixed    << endl;
   tbox::plog << "returned_dt = " << returned_dt << endl;
   
   tbox::plog << "--------------------- end stableDtOnPatch on patch" << endl;
   
   return returned_dt;
}


// ---------------------------------------------------------------------------

//
// Test the extrema on the patch
//
void MblkEuler::testPatchExtrema( hier::Patch<NDIM>& patch, const char *pos )
{
   tbox::Pointer< pdat::CellData<NDIM,double> > state = patch.getPatchData( d_state, getDataContext() );
   hier::IntVector<NDIM> state_ghosts = state->getGhostCellWidth();


   const hier::Index<NDIM> ifirst = patch.getBox().lower();
   const hier::Index<NDIM> ilast  = patch.getBox().upper();

   int imgn = ifirst(0) - state_ghosts(0);
   int imgx = ilast(0)  + state_ghosts(0);
   int jmgn = ifirst(1) - state_ghosts(1);
   int jmgx = ilast(1)  + state_ghosts(1);
   int kmgn = ifirst(2) - state_ghosts(2);
   int nxg   = imgx - imgn + 1;
   int nyg   = jmgx - jmgn + 1;
   int nxnyg = nxg*nyg;

   //
   //  compute field extrema to tag along with this and print out
   //
   double *psi_min = new double[d_nState];
   double *psi_max = new double[d_nState];
   for ( int ii = 0; ii < d_nState ; ii++ ) {
      psi_max[ii] = -1.e80;
      psi_min[ii] =  1.e80;
   }

   for ( int ii = 0; ii < d_nState ; ii++ ) {
      double *lstate = state->getPointer(ii);

      for ( int k = ifirst(2); k <= ilast(2) ; k++ ) {  // just loop over interior elements
	 for ( int j = ifirst(1); j <= ilast(1) ; j++ ) {
	    for ( int i = ifirst(0); i <= ilast(0) ; i++ ) {
	       
	       int gind = POLY3(i,j,k,imgn,jmgn,kmgn,nxg,nxnyg);

	       //
	       // some extra bounds checks for sanity
	       // 
	       psi_max[ii] = MAX( psi_max[ii], lstate[gind] );
	       psi_min[ii] = MIN( psi_min[ii], lstate[gind] );
	    }
	 }
      }
   }
   
   tbox::plog << endl << "extrema for the state follow " << pos << " (min,max) = " << endl;
   for ( int ii = 0; ii < d_nState ; ii++ ) {
      tbox::plog << d_state_names[ii] << " (min,max) = ";
      tbox::plog << psi_min[ii]       << " " << psi_max[ii] << endl;
   }

   delete [] psi_min;
   delete [] psi_max;
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

void MblkEuler::computeFluxesOnPatch( hier::Patch<NDIM>& patch,
				      const double time,
				      const double dt)
{
   (void) time;

   //
   // process the SAMRAI data
   //
   int block_number = getBlockNumber();
   int level_number =  patch.getPatchLevelNumber();
   setMappedGridOnPatch( patch, level_number, block_number );

   hier::Box<NDIM> pbox = patch.getBox();
   const hier::Index<NDIM> ifirst = patch.getBox().lower();
   const hier::Index<NDIM> ilast  = patch.getBox().upper();

   tbox::Pointer< pdat::CellData<NDIM,double> > state = patch.getPatchData( d_state, getDataContext() );
   tbox::Pointer< pdat::CellData<NDIM,double> > vol   = patch.getPatchData( d_vol,   getDataContext() );
   tbox::Pointer< pdat::SideData<NDIM,double> > flux  = patch.getPatchData( d_flux,  getDataContext() );
   tbox::Pointer< pdat::NodeData<NDIM,double> > xyz   = patch.getPatchData( d_xyz,   getDataContext() );

   TBOX_ASSERT( !state.isNull() );
   TBOX_ASSERT( !flux.isNull()  );
   TBOX_ASSERT( !vol.isNull()   );
   TBOX_ASSERT( state->getGhostCellWidth() == d_nghosts );
   TBOX_ASSERT( vol->getGhostCellWidth()   == d_nghosts );
   TBOX_ASSERT( flux->getGhostCellWidth()  == d_fluxghosts );
   TBOX_ASSERT( xyz->getGhostCellWidth()   == d_nodeghosts );

   tbox::plog << "--------------------- start computeFluxesOnPatch on patch (";
   tbox::plog << level_number << ")" << endl;
   tbox::plog << "TIMESTEP for level = " << level_number << ", ";
   tbox::plog << "dt    = " << dt << ", stime  = " << time << ", ftime = " << time+dt << endl;
   tbox::plog << "box   = " << pbox << endl;

   //
   // ------------------------------- the upwind bounds ----------------------------
   //

   int fx_imn  = ifirst(0)     - d_fluxghosts(0);
   int fx_imx  = ilast(0)  + 1 + d_fluxghosts(0);
   int fx_jmn  = ifirst(1)     - d_fluxghosts(1);
   int fx_jmx  = ilast(1)      + d_fluxghosts(1);
   int fx_kmn  = ifirst(2)     - d_fluxghosts(2);
   int fx_nx   = fx_imx - fx_imn + 1;
   int fx_ny   = fx_jmx - fx_jmn + 1;
   int fx_nxny = fx_nx*fx_ny;

   int fy_imn  = ifirst(0)     - d_fluxghosts(0);
   int fy_imx  = ilast(0)      + d_fluxghosts(0);
   int fy_jmn  = ifirst(1)     - d_fluxghosts(1);
   int fy_jmx  = ilast(1)  + 1 + d_fluxghosts(1);
   int fy_kmn  = ifirst(2)     - d_fluxghosts(2);
   int fy_nx   = fy_imx - fy_imn + 1;
   int fy_ny   = fy_jmx - fy_jmn + 1;
   int fy_nxny = fy_nx*fy_ny;

   int fz_imn  = ifirst(0)      - d_fluxghosts(0);
   int fz_imx  = ilast(0)       + d_fluxghosts(0);
   int fz_jmn  = ifirst(1)      - d_fluxghosts(1);
   int fz_jmx  = ilast(1)       + d_fluxghosts(1);
   int fz_kmn  = ifirst(2)      - d_fluxghosts(2);
   int fz_nx   = fz_imx - fz_imn + 1;
   int fz_ny   = fz_jmx - fz_jmn + 1;
   int fz_nxny = fz_nx*fz_ny;

   int imin = ifirst(0) - d_nghosts(0);
   int imax = ilast(0)  + d_nghosts(0);
   int jmin = ifirst(1) - d_nghosts(1);
   int jmax = ilast(1)  + d_nghosts(1);
   int kmin = ifirst(2) - d_nghosts(2);
   int nx   = imax - imin + 1;
   int ny   = jmax - jmin + 1;
   int nxny = nx*ny;

   int nd_imin = ifirst(0)     - d_nodeghosts(0);
   int nd_imax = ilast(0)  + 1 + d_nodeghosts(0);
   int nd_jmin = ifirst(1)     - d_nodeghosts(1);
   int nd_jmax = ilast(1)  + 1 + d_nodeghosts(1);
   int nd_kmin = ifirst(2)     - d_nodeghosts(2);
   int nd_nx   = nd_imax - nd_imin + 1;
   int nd_ny   = nd_jmax - nd_jmin + 1;
   int nd_nxny = nd_nx*nd_ny;

   //
   // get the pointers
   // 
   double *cvol = vol->getPointer();

   double *x = xyz->getPointer(0);
   double *y = xyz->getPointer(1);
   double *z = xyz->getPointer(2);

   double u = d_advection_velocity[0];
   double v = d_advection_velocity[1];
   double w = d_advection_velocity[2];
   double u0 = sqrt( u*u + v*v + w*w );

   // 0, cartesian, 1 R translation, 2 rigid body theta rotation, 4, rigid body phi rotation
   int VEL_TYPE = d_advection_vel_type;

   //
   // note on areas, we set the area vector Avector,
   //
   //     Avector = -0.5*( x3 - x1 ) cross product ( x4 - x2 ) 
   //
   // where the x1 .. x4 are the right hand rule circulation for the face
   //

   int depth = state->getDepth();
   
   for ( int idepth = 0 ; idepth < depth ; idepth++ ) {

      double *psi = state->getPointer(idepth); // assumed single depth here !!!!
      
      double *fx = flux->getPointer(0,idepth);
      double *fy = flux->getPointer(1,idepth);
      double *fz = flux->getPointer(2,idepth);

      //
      // compute the fluxes for the upwind method
      //
      for ( int k = ifirst(2); k <= ilast(2); k++ ) {
	 for ( int j = ifirst(1); j <= ilast(1); j++ ) {
	    for ( int i = ifirst(0); i <= ilast(0)+1; i++ ) {

	       // --------- get the indices
	       int ifx  = POLY3( i, j, k, fx_imn, fx_jmn, fx_kmn, fx_nx, fx_nxny );

	       int ind  = POLY3( i, j, k, imin, jmin, kmin, nx, nxny );
	       int ib   = ind-1;

	       int n1 = POLY3( i, j, k, nd_imin, nd_jmin, nd_kmin, nd_nx, nd_nxny );
	       int n4 = n1 + nd_nx;
	       int n5 = n1 + nd_nxny;
	       int n8 = n1 + nd_nxny + nd_nx;

	       // --------- get the positions // 1 - 4 - 8 - 5
	       double x1 = x[n1];  double y1 = y[n1];  double z1 = z[n1];
	       double x2 = x[n4];  double y2 = y[n4];  double z2 = z[n4];
	       double x3 = x[n8];  double y3 = y[n8];  double z3 = z[n8];
	       double x4 = x[n5];  double y4 = y[n5];  double z4 = z[n5];

	       double xm = 0.25*(x1 + x2 + x3 + x4);
	       double ym = 0.25*(y1 + y2 + y3 + y4);
	       double zm = 0.25*(z1 + z2 + z3 + z4);

	       double xn1 = sqrt( xm*xm + ym*ym );
	       double R   = sqrt( xm*xm + ym*ym + zm*zm );
	    
	       // --------- compute the flux
	       double u1 = u0 * xm / R;
	       double v1 = u0 * ym / R;
	       double w1 = u0 * zm / R;
  
	       double dx31 = x3 - x1;
	       double dx42 = x4 - x2;
	    
	       double dy31 = y3 - y1;
	       double dy42 = y4 - y2;
	    
	       double dz31 = z3 - z1;
	       double dz42 = z4 - z2;
	    
	       double Ax = -0.5 * (dy42 * dz31 - dz42 * dy31);
	       double Ay = -0.5 * (dz42 * dx31 - dx42 * dz31);
	       double Az = -0.5 * (dx42 * dy31 - dy42 * dx31);

	       double Audotn;
	       if (  VEL_TYPE == 0 ) {
		  Audotn = Ax*u + Ay*v + Az*w;
	       }
	       if (  VEL_TYPE == 1 ) {
		  Audotn = Ax*u1 + Ay*v1 + Az*w1;
	       }
	       if (  VEL_TYPE == 2 ) {
		  double u2  = - u0 * ym;
		  double v2  =   u0 * xm;
		  double w2  =   0.0;
		  Audotn = Ax*u2 + Ay*v2 + Az*w2;
	       }
	       if (  VEL_TYPE == 3 ) {
		  double cosphi = xn1/R;
		  double sinphi = sqrt( 1 - xn1*xn1/(R*R)) *( zm > 0.0 ? 1.0 : -1.0 );
		  double u2  = - u0 * xm * sinphi / cosphi;
		  double v2  = - u0 * ym * sinphi / cosphi;
		  double w2  =   u0 * cosphi / sinphi;
		  Audotn = Ax*u2 + Ay*v2 + Az*w2;
	       }

	       double lflux = ( Audotn > 0.0 ? psi[ib] : psi[ind] )*Audotn;

	       fx[ifx] = lflux;  // (Ax,Ay,Az) vector points towards ind
	    }
	 }
      }


      for ( int k = ifirst(2); k <= ilast(2); k++ ) {
	 for ( int j = ifirst(1); j <= ilast(1)+1; j++ ) {
	    for ( int i = ifirst(0); i <= ilast(0); i++ ) {

	       // --------- get the indices
	       int ify  = POLY3( i, j, k, fy_imn, fy_jmn, fy_kmn, fy_nx, fy_nxny );

	       int ind  = POLY3( i, j, k, imin, jmin, kmin, nx, nxny );
	       int jb   = ind-nx;

	       int n1 = POLY3( i, j, k, nd_imin, nd_jmin, nd_kmin, nd_nx, nd_nxny );
	       int n2 = n1 + 1;
	       int n5 = n1 + nd_nxny;
	       int n6 = n1 + nd_nxny + 1;

	       // --------- get the positions // 1 - 5 - 6 - 2
	       double x1 = x[n1];  double y1 = y[n1];  double z1 = z[n1];
	       double x2 = x[n5];  double y2 = y[n5];  double z2 = z[n5];
	       double x3 = x[n6];  double y3 = y[n6];  double z3 = z[n6];
	       double x4 = x[n2];  double y4 = y[n2];  double z4 = z[n2];

	       double xm = 0.25*(x1 + x2 + x3 + x4);
	       double ym = 0.25*(y1 + y2 + y3 + y4);
	       double zm = 0.25*(z1 + z2 + z3 + z4);
	    
	       double xn1 = sqrt( xm*xm + ym*ym );
	       double R   = sqrt( xm*xm + ym*ym + zm*zm );

	       // --------- compute the flux
	       double u1 = u0 * xm / R;
	       double v1 = u0 * ym / R;
	       double w1 = u0 * zm / R;
  
	       double dx31 = x3 - x1;
	       double dx42 = x4 - x2;
	    
	       double dy31 = y3 - y1;
	       double dy42 = y4 - y2;
	    
	       double dz31 = z3 - z1;
	       double dz42 = z4 - z2;
	    
	       double Ax = -0.5 * (dy42 * dz31 - dz42 * dy31);
	       double Ay = -0.5 * (dz42 * dx31 - dx42 * dz31);
	       double Az = -0.5 * (dx42 * dy31 - dy42 * dx31);

	       double Audotn;
	       if (  VEL_TYPE == 0 ) {
		  Audotn = Ax*u + Ay*v + Az*w;
	       }
	       if (  VEL_TYPE == 1 ) {
		  Audotn = Ax*u1 + Ay*v1 + Az*w1;
	       }
	       if (  VEL_TYPE == 2 ) {
		  double u2  = - u0 * ym;
		  double v2  =   u0 * xm;
		  double w2  =   0.0;
		  Audotn = Ax*u2 + Ay*v2 + Az*w2;
	       }
	       if (  VEL_TYPE == 3 ) {
		  double cosphi = xn1/R;
		  double sinphi = sqrt( 1 - xn1*xn1/(R*R)) *( zm > 0.0 ? 1.0 : -1.0 );
		  double u2  = - u0 * xm * sinphi / cosphi;
		  double v2  = - u0 * ym * sinphi / cosphi;
		  double w2  =   u0 * cosphi / sinphi;
		  Audotn = Ax*u2 + Ay*v2 + Az*w2;
	       }

	       double lflux = ( Audotn > 0.0 ? psi[jb] : psi[ind] )*Audotn;
 
	       fy[ify] = lflux;  // (Ax,Ay,Az) vector points towards ind
	    }
	 }
      }

      for ( int k = ifirst(2); k <= ilast(2)+1; k++ ) {
	 for ( int j = ifirst(1); j <= ilast(1); j++ ) {
	    for ( int i = ifirst(0); i <= ilast(0); i++ ) {

	       // --------- get the indices
	       int ifz  = POLY3( i,j,k,fz_imn,fz_jmn,fz_kmn,fz_nx,fz_nxny );

	       int ind  = POLY3( i,j,k,imin,jmin,kmin,nx,nxny );
	       int kb   = ind-nxny;

	       int n1 = POLY3( i,j,k,nd_imin,nd_jmin,nd_kmin,nd_nx,nd_nxny );
	       int n2 = n1 + 1;
	       int n3 = n1 + 1 + nd_nx;
	       int n4 = n1 + nd_nx;

	       // --------- get the positions // 1 - 2 - 3 - 4
	       double x1 = x[n1];  double y1 = y[n1];  double z1 = z[n1];
	       double x2 = x[n2];  double y2 = y[n2];  double z2 = z[n2];
	       double x3 = x[n3];  double y3 = y[n3];  double z3 = z[n3];
	       double x4 = x[n4];  double y4 = y[n4];  double z4 = z[n4];

	       double xm = 0.25*(x1 + x2 + x3 + x4);
	       double ym = 0.25*(y1 + y2 + y3 + y4);
	       double zm = 0.25*(z1 + z2 + z3 + z4);

	       double xn1 = sqrt( xm*xm + ym*ym );
	       double R   = sqrt( xm*xm + ym*ym + zm*zm );
	    
	       // --------- compute the flux
	       double u1 = u0 * xm / R;
	       double v1 = u0 * ym / R;
	       double w1 = u0 * zm / R;
  
	       double dx31 = x3 - x1;
	       double dx42 = x4 - x2;
	    
	       double dy31 = y3 - y1;
	       double dy42 = y4 - y2;
	    
	       double dz31 = z3 - z1;
	       double dz42 = z4 - z2;
	    
	       double Ax = -0.5 * (dy42 * dz31 - dz42 * dy31);
	       double Ay = -0.5 * (dz42 * dx31 - dx42 * dz31);
	       double Az = -0.5 * (dx42 * dy31 - dy42 * dx31);

	       double Audotn;
	       if (  VEL_TYPE == 0 ) {
		  Audotn = Ax*u + Ay*v + Az*w;
	       }
	       if (  VEL_TYPE == 1 ) {
		  Audotn = Ax*u1 + Ay*v1 + Az*w1;
	       }
	       if (  VEL_TYPE == 2 ) {
		  double u2  = - u0 * ym;
		  double v2  =   u0 * xm;
		  double w2  =   0.0;
		  Audotn = Ax*u2 + Ay*v2 + Az*w2;
	       }
	       if (  VEL_TYPE == 3 ) {
		  double cosphi = xn1/R;
		  double sinphi = sqrt( 1 - xn1*xn1/(R*R)) *( zm > 0.0 ? 1.0 : -1.0 );
		  double u2  = - u0 * xm * sinphi / cosphi;
		  double v2  = - u0 * ym * sinphi / cosphi;
		  double w2  =   u0 * cosphi / sinphi;
		  Audotn = Ax*u2 + Ay*v2 + Az*w2;
	       }

	       double lflux = ( Audotn > 0.0 ? psi[kb] : psi[ind] )*Audotn;

	       fz[ifz] = lflux;  // (Ax,Ay,Az) vector points towards ind
	    }
	 }
      }



      //
      // compute the source due to the upwind method
      //
      for ( int k = ifirst(2); k <= ilast(2); k++ ) {
	 for ( int j = ifirst(1); j <= ilast(1); j++ ) {
	    for ( int i = ifirst(0); i <= ilast(0); i++ ) {

	       int ind  = POLY3(i,j,k,imin,jmin,kmin,nx,nxny);

	       int ib   = POLY3(i,  j,k,fx_imn,fx_jmn,fx_kmn,fx_nx,fx_nxny);
	       int ie   = POLY3(i+1,j,k,fx_imn,fx_jmn,fx_kmn,fx_nx,fx_nxny);

	       int jb   = POLY3(i,j,  k,fy_imn,fy_jmn,fy_kmn,fy_nx,fy_nxny);
	       int je   = POLY3(i,j+1,k,fy_imn,fy_jmn,fy_kmn,fy_nx,fy_nxny);

	       int kb   = POLY3(i,j,k,  fz_imn,fz_jmn,fz_kmn,fz_nx,fz_nxny);
	       int ke   = POLY3(i,j,k+1,fz_imn,fz_jmn,fz_kmn,fz_nx,fz_nxny);


	       //   have set up the normal so that it always points in the positive, i, j, and k directions
	       int n1 = POLY3(i,j,k,nd_imin,nd_jmin,nd_kmin,nd_nx,nd_nxny);
	       int n2 = n1 + 1;
	       int n3 = n1 + 1 + nd_nx;
	       int n4 = n1 + nd_nx;

	       int n5 = n1 + nd_nxny;
	       int n6 = n1 + nd_nxny + 1;
	       int n7 = n1 + nd_nxny + 1 + nd_nx;
	       int n8 = n1 + nd_nxny + nd_nx;

	       double lvol = UpwindVolume( x[n1], x[n2], x[n3], x[n4],
					   x[n5], x[n6], x[n7], x[n8],
				      
					   y[n1], y[n2], y[n3], y[n4],
					   y[n5], y[n6], y[n7], y[n8],
					   
					   z[n1], z[n2], z[n3], z[n4],
					   z[n5], z[n6], z[n7], z[n8] );
	    
	       cvol[ind] = lvol;

	       psi[ind] -= dt*( ( fx[ie] - fx[ib] ) +
				( fy[je] - fy[jb] ) +
				( fz[ke] - fz[kb] ) )/lvol;


	    }
	 }
      }

   } // end of field loop

   tbox::plog << "--------------------- end computeFluxesOnPatch on patch" << endl;
}

/*
*************************************************************************
*                                                                       *
* Update solution variables by performing a conservative                *
* difference with the fluxes calculated in computeFluxesOnPatch().      *
*                                                                       *
*************************************************************************
*/

void MblkEuler::conservativeDifferenceOnPatch(hier::Patch<NDIM>& patch,
                                           const double time,
                                           const double dt,
                                           bool at_syncronization)
{
   (void) time;
   (void) dt;
   (void) at_syncronization;

   return;
}


 
/*
*************************************************************************
*                                                                       *
* this routine initializes all the zonal data that one needs to fill to dummy
* values as a check
*                                                                       *
*************************************************************************
*/


void MblkEuler::markPhysicalBoundaryConditions( hier::Patch<NDIM>& patch, 
						const double fill_time,
						const hier::IntVector<NDIM>& ghost_width_to_fill )
{
   tbox::Pointer< pdat::CellData<NDIM,double> > state = patch.getPatchData(d_state, getDataContext());

   //
   // the domain and its ghost box
   //
   const hier::Box<NDIM>& interior = patch.getBox();

   const hier::Box<NDIM>&  ghost_box      = state->getGhostBox();
   const hier::IntVector<3>& ghost_cells  = state->getGhostCellWidth();

   hier::IntVector<3> gcw_to_fill = hier::IntVector<3>::min( ghost_cells,
							     ghost_width_to_fill );


   int bc_types[3] = { FACE3D_BDRY_TYPE, EDGE3D_BDRY_TYPE, NODE3D_BDRY_TYPE };

   int imin = ghost_box.lower(0);
   int imax = ghost_box.upper(0);
   int jmin = ghost_box.lower(1);
   int jmax = ghost_box.upper(1);
   int kmin = ghost_box.lower(2);
   int nx   = imax - imin + 1;
   int ny   = jmax - jmin + 1;
   int nxny = nx*ny;

   const tbox::Pointer<geom::BlockPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();

   for ( int ii = 0 ; ii < 3 ; ii++ ) {
      
      const tbox::Array<hier::BoundaryBox<3> >& bc_bdry =
	 pgeom->getCodimensionBoundary(bc_types[ii]);
   
      for (int jj = 0; jj < bc_bdry.getSize(); jj++) {

	 hier::Box<3> fill_box = pgeom->getBoundaryFillBox( bc_bdry[jj],
							    interior,
							    gcw_to_fill );
	 int l_imin = fill_box.lower(0);
	 int l_jmin = fill_box.lower(1);
	 int l_kmin = fill_box.lower(2);

	 int l_imax = fill_box.upper(0);
	 int l_jmax = fill_box.upper(1);
	 int l_kmax = fill_box.upper(2);
	 
	 int nd = state->getDepth();

	 for ( int n = 0 ; n < nd ; n++ ) {
	    double *sptr = state->getPointer(n);

	    for ( int k = l_kmin ; k <= l_kmax ; k++ ) {
	       for ( int j = l_jmin ; j <= l_jmax ; j++ ) {
		  for ( int i = l_imin ; i <= l_imax ; i++ ) {
		     
		     int ind  = POLY3(i,j,k,imin,jmin,kmin,nx,nxny);
		     
		     sptr[ind] = ii;

		  }
	       }
	    }
	 } // end of marking loops

      } // loop over boxes in each type

   } // loop over box types

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

#if 1
#define BCMULTIBLOCK bcmultiblock_
#else
#define BCMULTIBLOCK BCMULTIBLOCK
#endif

extern "C" { 
   void BCMULTIBLOCK( const int *, const int *, const int *, 
		       const int *, const int *, const int *, 
		       double *,
		       const int *, const int *, const int *, 
		       const int *, const int *, const int *, 
		       double *,
		       const int* ncomp,
		       const int* dlo, const int* dhi,
		       const int* glo, const int* ghi,
		       const int* lo,  const int* hi,
		       const int* dir, const int* side );
}

void MblkEuler::setPhysicalBoundaryConditions( hier::Patch<NDIM>& patch, 
					       const double fill_time,
					       const hier::IntVector<NDIM>& ghost_width_to_fill )
{
   markPhysicalBoundaryConditions( patch, fill_time, ghost_width_to_fill );

   (void) fill_time;
   
   tbox::Pointer< pdat::NodeData<NDIM,double> > position = patch.getPatchData(d_xyz, getDataContext());
   tbox::Pointer< pdat::CellData<NDIM,double> > state    = patch.getPatchData(d_state, getDataContext());

   //
   // the patch and its ghost box
   //
   const hier::Box<NDIM>& patch_box  = patch.getBox();
   const hier::Box<NDIM>&  ghost_box = state->getGhostBox();
   const hier::Box<NDIM>& nghost_box = position->getGhostBox();

   double *state_ptr    = state->getPointer();
   int state_depth      = state->getDepth();
   double *position_ptr = position->getPointer();

   //
   // the index space of this block and its neighbors
   //
   int block_number = getBlockNumber();
   const tbox::Pointer<geom::BlockPatchGeometry<NDIM> > patch_geom = patch.getPatchGeometry();
   const hier::IntVector<NDIM> ratio = patch_geom->getRatio();
   hier::BoxArray<NDIM> domain_boxes;
   d_grid_geometries[block_number]->computePhysicalDomain(domain_boxes, ratio); 

   const hier::IntVector<NDIM>& periodic = d_grid_geometries[block_number]->getPeriodicShift();

   const hier::Box<NDIM>& domain_box = domain_boxes(0);
   // domain_box.refine(patch_geom->getRatio());

   d_mblk_geometry->buildLocalBlocks( patch_box,
				      domain_box,
				      block_number,
				      d_dom_local_blocks );
   //
   // loop over the dimensions, filling in boundary conditions where needed
   // note that here we have a check for is this really a physical boundary or is it
   // just a periodic boundary condition, or is it just an internal block boundary 
   //
   int imin = ghost_box.lower(0);
   int imax = ghost_box.upper(0);
   int jmin = ghost_box.lower(1);
   int jmax = ghost_box.upper(1);
   int kmin = ghost_box.lower(2);
   int kmax = ghost_box.upper(2);
   
   int nd_imin = nghost_box.lower(0);
   int nd_imax = nghost_box.upper(0);
   int nd_jmin = nghost_box.lower(1);
   int nd_jmax = nghost_box.upper(1);
   int nd_kmin = nghost_box.lower(2);
   int nd_kmax = nghost_box.upper(2);

   for (int dir = 0; dir < NDIM; dir++) {
      if ( !periodic(dir) ) {
	 
	 if ( ( ghost_box.lower(dir) < domain_box.lower(dir) ) && 
	      ( d_dom_local_blocks[dir] == block_number      ) ) {
	    int iside = 0;
	    BCMULTIBLOCK( &nd_imin, &nd_imax, &nd_jmin, &nd_jmax, &nd_kmin, &nd_kmax, position_ptr,
			   &imin, &imax, &jmin, &jmax, &kmin, &kmax, state_ptr, &state_depth,
			   domain_box.lower(),
			   domain_box.upper(),
			   ghost_box.lower(),
			   ghost_box.upper(),
			   patch_box.lower(),
			   patch_box.upper(),
			   &dir, &iside );
	 }
	 
	 if ( ( ghost_box.upper(dir) > domain_box.upper(dir) ) && 
	      ( d_dom_local_blocks[NDIM + dir] == block_number      ) ) {
	    int iside = 1;
	    BCMULTIBLOCK( &nd_imin, &nd_imax, &nd_jmin, &nd_jmax, &nd_kmin, &nd_kmax, position_ptr,
			   &imin, &imax, &jmin, &jmax, &kmin, &kmax, state_ptr, &state_depth,
			   domain_box.lower(),
			   domain_box.upper(),
			   ghost_box.lower(),
			   ghost_box.upper(),
			   patch_box.lower(),
			   patch_box.upper(),
			   &dir, &iside );
	 }
	 
      } // end of periodic check
      
   } // end of dimension loop
   
}


/*
*************************************************************************
*                                                                       *
* Refine operations 
*                                                                       *
*************************************************************************
*/

void MblkEuler::preprocessRefine( 
   hier::Patch<NDIM>& fine,
   const hier::Patch<NDIM>& coarse,
   const hier::Box<NDIM>& fine_box,
   const hier::IntVector<NDIM>& ratio )
{

   int xyz_id =  hier::VariableDatabase<NDIM>::getDatabase()->
      mapVariableAndContextToIndex(d_xyz, getDataContext());

   int fln = fine.getPatchLevelNumber();
   int cln = coarse.getPatchLevelNumber();   
   if (fln < 0) {
      fln = cln + 1;
      if (!fine.checkAllocated(xyz_id)) {
         TBOX_ERROR(d_object_name << ":preprocessRefine()"
                    << "\nfine xyz data not allocated" << endl);
      }
   }
   if (cln < 0) {
      cln = fln - 1;
      if (!coarse.checkAllocated(xyz_id)) {
         TBOX_ERROR(d_object_name << ":preprocessRefine()"
                    << "\ncoarse xyz data not allocated" << endl);
      }
   }
   int block_number = getBlockNumber();
   setMappedGridOnPatch(coarse, cln, block_number);
   setMappedGridOnPatch(fine, fln, block_number);
   setVolumeOnPatch( coarse );
   setVolumeOnPatch( fine );
}



// ---------------------------------------------------------------------------

//
// the refinement operator
//
void MblkEuler::postprocessRefine( hier::Patch<NDIM>& fine,
				   const hier::Patch<NDIM>& coarse,
				   const hier::Box<NDIM>& fine_box, // where the fine data is needed
				   const hier::IntVector<NDIM>& ratio )
{
   tbox::Pointer< pdat::CellData<NDIM,double> > cstate = coarse.getPatchData( d_state, getDataContext() );
   tbox::Pointer< pdat::CellData<NDIM,double> > cvol   = coarse.getPatchData( d_vol,   getDataContext() );

   tbox::Pointer< pdat::CellData<NDIM,double> > fstate = fine.getPatchData( d_state, getDataContext() );
   tbox::Pointer< pdat::CellData<NDIM,double> > fvol   = fine.getPatchData( d_vol,   getDataContext() );

   int depth = cstate->getDepth();

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( !cstate.isNull() );
   TBOX_ASSERT( !fstate.isNull() );
   TBOX_ASSERT( !cvol.isNull()   );
   TBOX_ASSERT( !fvol.isNull()   );
   TBOX_ASSERT( cstate->getDepth() == fstate->getDepth() );
#endif

   //
   // get the boxes and bounds
   //

   // ... the bounds of the data
   const hier::Box<NDIM> cgbox( cstate->getGhostBox() ); 
   const hier::Box<NDIM> fgbox( fstate->getGhostBox() ); 

   const hier::Index<NDIM> cilo = cgbox.lower();
   const hier::Index<NDIM> cihi = cgbox.upper();
   const hier::Index<NDIM> filo = fgbox.lower(); 
   const hier::Index<NDIM> fihi = fgbox.upper();

   // ... the bounds we actually work on
   const hier::Box<NDIM> coarse_box = hier::Box<NDIM>::coarsen(fine_box, ratio ); 
   const hier::Index<NDIM> ifirstc = coarse_box.lower();
   const hier::Index<NDIM> ilastc  = coarse_box.upper();
   const hier::Index<NDIM> ifirstf = fine_box.lower();
   const hier::Index<NDIM> ilastf  = fine_box.upper();

   int flev = fine.getPatchLevelNumber();
   int clev = coarse.getPatchLevelNumber();

   tbox::plog << "--------------------- start postprocessRefineData (" << flev << "," << clev << ")" << endl;
   tbox::plog << "flevel     = " << flev << endl;
   tbox::plog << "clevel     = " << clev << endl << endl;
   tbox::plog << "fine_box         = " << fine_box        << endl;
   tbox::plog << "fine_patch_box   = " << fine.getBox()   << endl;
   tbox::plog << "fine_ghost_box   = " << fgbox           << endl << endl;

   tbox::plog << "coarse_box       = " << coarse_box      << endl;
   tbox::plog << "coarse_patch_box = " << coarse.getBox() << endl;
   tbox::plog << "coarse_ghost_box = " << coarse_box      << endl;

   //
   // setup work variables 
   //
   const hier::IntVector<NDIM> tmp_ghosts(0 );
   const int numInterp = 1;
   pdat::CellData<NDIM,double> val_vals(    coarse_box, numInterp, tmp_ghosts );
   pdat::CellData<NDIM,double> slope0_vals( coarse_box, numInterp, tmp_ghosts );
   pdat::CellData<NDIM,double> slope1_vals( coarse_box, numInterp, tmp_ghosts );
   pdat::CellData<NDIM,double> slope2_vals( coarse_box, numInterp, tmp_ghosts );
   double *val      	  = val_vals.getPointer();
   double *slope0 	  = slope0_vals.getPointer();
   double *slope1 	  = slope1_vals.getPointer();
   double *slope2 	  = slope2_vals.getPointer();

   //
   // setup coarse strides
   //
   int imin = ifirstc(0);  // the box of coarse elements being refined
   int imax = ilastc(0);   // the work data is sized to this box
   int jmin = ifirstc(1);
   int jmax = ilastc(1);
   int kmin = ifirstc(2);
   int kmax = ilastc(2);
   int nx   = imax - imin + 1;
   int ny   = jmax - jmin + 1;
   int nz   = kmax - kmin + 1;
   int nxny = nx*ny;
   int nel  = nx*ny*nz;

   int cimin = cilo(0);  // the coarse data bounds
   int cimax = cihi(0);
   int cjmin = cilo(1);
   int cjmax = cihi(1);
   int ckmin = cilo(2);
   int cnx   = cimax - cimin + 1;
   int cny   = cjmax - cjmin + 1;
   int cnxny = cnx*cny;

   int fimin = filo(0);  // the fine data bounds
   int fimax = fihi(0);
   int fjmin = filo(1);
   int fjmax = fihi(1);
   int fkmin = filo(2);
   int fnx   = fimax - fimin + 1;
   int fny   = fjmax - fjmin + 1;
   int fnxny = fnx*fny;
   
   double rat0 = ratio[0];
   double rat1 = ratio[1];
   double rat2 = ratio[2];
   double fact = 2.0;  // xi varies from -1 to 1

   //
   // ================================= state variable refinement ====================
   //
   double *cdata = NULL; // keeps pointers around till end of loop
   double *fdata = NULL; 

   for ( int n = 0 ; n < depth; n++ ) {   

      cdata = cstate->getPointer(n);
      fdata = fstate->getPointer(n);

      for ( int l = 0; l < nel ; l++ ) {    // default slopes are zero
	 slope0[l] = 0.0;  // this yields piecewise constant interpolation
	 slope1[l] = 0.0;  // and makes a handy initializer
	 slope2[l] = 0.0;
      }

      for ( int k = ifirstc(2); k <= ilastc(2); k++ ) {  
	 for ( int j = ifirstc(1); j <= ilastc(1); j++ ) {
	    for ( int i = ifirstc(0); i <= ilastc(0); i++ ) {
	       
	       int ind  = POLY3(i,j,k,imin,jmin,kmin,nx,nxny);
	       int cind = POLY3(i,j,k,cimin,cjmin,ckmin,cnx,cnxny);

	       double aii = cdata[ cind ];
	       val[ind] = aii;
	       
#ifdef DEBUG_CHECK_ASSERTIONS
	       TBOX_ASSERT( ind >= 0 ); // debug assertions
	       TBOX_ASSERT( ind < nel );
#endif
#if 0 // turn to zero for simple interp
	       int im1 = cind-1;
	       int ip1 = cind+1;
	       int jm1 = cind-cnx;
	       int jp1 = cind+cnx;
	       int km1 = cind-cnxny;
	       int kp1 = cind+cnxny;
	       
               double w_i  = cvolume[ cind ];
               double w_ip = cvolume[ ip1 ];
               double w_im = cvolume[ im1 ];
	       double w_jp = cvolume[ jp1 ];
               double w_jm = cvolume[ jm1 ];
               double w_kp = cvolume[ kp1 ];
               double w_km = cvolume[ km1 ];
	       
	       double aip = cdata[ ip1 ];
	       double aim = cdata[ im1 ];
	       double ajp = cdata[ jp1 ];
	       double ajm = cdata[ jm1 ];
	       double akp = cdata[ kp1 ];
	       double akm = cdata[ km1 ];
	       
	       my_slopes( aii, aip, aim, ajp, ajm, akp, akm,
                          w_i, w_ip, w_im, w_jp, w_jm, w_kp, w_km,
                          slope0[ind], 
                          slope1[ind], 
                          slope2[ind] );
#endif
	    }
	 }
      }

      //
      // compute the interpolated data from the cached slopes, looping
      // over the fine zones
      //
      for ( int k = ifirstf(2); k <= ilastf(2); k++ ) {  
	 for ( int j = ifirstf(1); j <= ilastf(1); j++ ) {
	    for ( int i = ifirstf(0); i <= ilastf(0); i++ ) {	       

	       int find = POLY3(i,j,k,fimin,fjmin,fkmin,fnx,fnxny);
	       
	       double ric   = (double(i) + 0.5)/rat0 - 0.5;
	       double rjc   = (double(j) + 0.5)/rat1 - 0.5;
	       double rkc   = (double(k) + 0.5)/rat2 - 0.5;
	       
	       int ic    = (int) (ric + ( ric >= 0 ? 0.5 : -0.5 ) ); 
	       int jc    = (int) (rjc + ( rjc >= 0 ? 0.5 : -0.5 ) );
	       int kc    = (int) (rkc + ( rkc >= 0 ? 0.5 : -0.5 ) ); 
	       
	       double ldx   = fact*( ric - ic );
	       double ldy   = fact*( rjc - jc );
	       double ldz   = fact*( rkc - kc );
	       
	       int ind  = POLY3(ic,jc,kc,imin,jmin,kmin,nx,nxny);  // work pos

#ifdef DEBUG_CHECK_ASSERTIONS
	       TBOX_ASSERT( 0   <= ind );
	       TBOX_ASSERT( ind <  nel );
#endif 

	       fdata[find] =
		  val[ind] + 
		  ldx*slope0[ind] +
		  ldy*slope1[ind] +
		  ldz*slope2[ind];
	    }
	 }
      } // end of i,j,k loops for finding the fine state variables

   } // end of state loop

   tbox::plog << "--------------------- end postprocessRefine" << endl;
}


/*
*************************************************************************
*                                                                       *
* Coarsen operations 
*                                                                       *
*************************************************************************
*/
void MblkEuler::preprocessCoarsen( hier::Patch<NDIM>& coarse,
                                   const hier::Patch<NDIM>& fine,
                                   const hier::Box<NDIM>& coarse_box,
                                   const hier::IntVector<NDIM>& ratio)
{
   int xyz_id =  hier::VariableDatabase<NDIM>::getDatabase()->
      mapVariableAndContextToIndex(d_xyz, getDataContext());

   int fln = fine.getPatchLevelNumber();
   int cln = coarse.getPatchLevelNumber();   
   if (fln < 0) {
      fln = cln + 1;
      if (!fine.checkAllocated(xyz_id)) {
         TBOX_ERROR(d_object_name << ":preprocessCoarsen()"
                    << "\nfine xyz data not allocated" << endl);
      }
   }
   if (cln < 0) {
      cln = fln - 1;
      if (!coarse.checkAllocated(xyz_id)) {
         TBOX_ERROR(d_object_name << ":preprocessCoarsen()"
                    << "\ncoarse xyz data not allocated" << endl);
      }
   }
   int block_number = getCoarsenBlockNumber();
   setMappedGridOnPatch(coarse, cln, block_number);
   setMappedGridOnPatch(fine, fln, block_number);
   setVolumeOnPatch( coarse );
   setVolumeOnPatch( fine );
}

// ---------------------------------------------------------------------------

//
// the coarsening function
//
void MblkEuler::postprocessCoarsen( hier::Patch<NDIM>& coarse,
				    const hier::Patch<NDIM>& fine,
				    const hier::Box<NDIM>& coarse_box,
				    const hier::IntVector<NDIM>& ratio )
{
   tbox::Pointer< pdat::CellData<NDIM,double> > cstate     = 
      coarse.getPatchData(d_state, getDataContext());
   tbox::Pointer< pdat::CellData<NDIM,double> > cvol      = 
      coarse.getPatchData(d_vol, getDataContext());

   tbox::Pointer< pdat::CellData<NDIM,double> > fstate     = 
      fine.getPatchData(d_state, getDataContext());
   tbox::Pointer< pdat::CellData<NDIM,double> > fvol      = 
      fine.getPatchData(d_vol, getDataContext());

   int depth = cstate->getDepth();
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!cstate.isNull());
   TBOX_ASSERT(!cvol.isNull());
   TBOX_ASSERT(!fstate.isNull());
   TBOX_ASSERT(!fvol.isNull());
   TBOX_ASSERT(cstate->getDepth() == fstate->getDepth());
#endif

   //
   // box and geometry information
   //
   const hier::Index<NDIM> filo = fstate->getGhostBox().lower();
   const hier::Index<NDIM> fihi = fstate->getGhostBox().upper();  
   const hier::Index<NDIM> cilo = cstate->getGhostBox().lower();
   const hier::Index<NDIM> cihi = cstate->getGhostBox().upper();

   const hier::Box<NDIM> fine_box = hier::Box<NDIM>::refine(coarse_box, ratio );
   const hier::Index<NDIM> ifirstc = coarse_box.lower();  // coarse basis
   const hier::Index<NDIM> ilastc  = coarse_box.upper();
   const hier::Index<NDIM> ifirstf = fine_box.lower();    // fine data
   const hier::Index<NDIM> ilastf  = fine_box.upper();

   int flev = fine.getPatchLevelNumber();
   int clev = coarse.getPatchLevelNumber();

   tbox::plog << "--------------------- start postprocessCoarsenData (" << flev << "," << clev << ")" << endl;
   tbox::plog << "flevel     = " << flev << endl;
   tbox::plog << "clevel     = " << clev << endl;
   tbox::plog << "fine       = " << fine.getBox()   << endl;
   tbox::plog << "coarse     = " << coarse.getBox() << endl;
   tbox::plog << "fine_box   = " << fine_box        << endl;
   tbox::plog << "coarse_box = " << coarse_box      << endl;

   tbox::plog << "filo = " << filo << ", fihi = " << fihi << endl;
   tbox::plog << "cilo = " << cilo << ", cihi = " << cihi << endl;

   //
   // work variables
   //

   int cimin = cilo(0);
   int cimax = cihi(0);
   int cjmin = cilo(1);
   int cjmax = cihi(1);
   int ckmin = cilo(2);
   int ckmax = cihi(2);
   int cnx   = cimax - cimin + 1;
   int cny   = cjmax - cjmin + 1;
   int cnz   = ckmax - ckmin + 1;
   int cnxny = cnx*cny;
   int cnel  = cnx*cny*cnz;

   int fimin = filo(0);
   int fimax = fihi(0);
   int fjmin = filo(1);
   int fjmax = fihi(1);
   int fkmin = filo(2);
   int fkmax = fihi(2);
   int fnx   = fimax - fimin + 1;
   int fny   = fjmax - fjmin + 1;
   int fnz   = fkmax - fkmin + 1;
   int fnxny = fnx*fny;
   int fnel  = fnx*fny*fnz;
   
   double rat0 = (double)(ratio(0));
   double rat1 = (double)(ratio(1));
   double rat2 = (double)(ratio(2));

   double *cvolume = cvol->getPointer();
   double *fvolume = fvol->getPointer();
   double *cdata   = cstate->getPointer();
   double *fdata   = fstate->getPointer();
      
   //
   // average the data
   //
   for ( int n = 0 ; n < depth; n++ ) {
      
      //
      // zero out the underlying coarse data to serve as a counter
      //
      for ( int k = ifirstc(2); k <= ilastc(2); k++ ) {  // loop over the coarse zones
	 for ( int j = ifirstc(1); j <= ilastc(1); j++ ) {
	    for ( int i = ifirstc(0); i <= ilastc(0); i++ ) {
	       
	       int chind     = POLY3(i,j,k,cimin,cjmin,ckmin,cnx,cnxny) + n*cnel;

	       cdata[chind] = 0.0; 
	    }
	 }
      }
     
      //
      // compute the interpolated data from the cached slopes
      //
      
      for ( int k = ifirstf(2); k <= ilastf(2); k++ ) {  // loop over the coarse zones
	 for ( int j = ifirstf(1); j <= ilastf(1); j++ ) {
	    for ( int i = ifirstf(0); i <= ilastf(0); i++ ) {
	       
	       int vol_ind = POLY3(i,j,k,fimin,fjmin,fkmin,fnx,fnxny);
	       int fhind   = vol_ind + n*fnel;
	       
	       double ric   = (double(i) + 0.5)/rat0 - 0.5;
	       double rjc   = (double(j) + 0.5)/rat1 - 0.5;
	       double rkc   = (double(k) + 0.5)/rat2 - 0.5;
	   
	       int ic    = (int) (ric + ( ric >= 0 ? 0.5 : -0.5 ) ); // a round operation
	       int jc    = (int) (rjc + ( rjc >= 0 ? 0.5 : -0.5 ) ); // shift up and truncate if ic > 0
	       int kc    = (int) (rkc + ( rkc >= 0 ? 0.5 : -0.5 ) ); // shift down and truncate if ic < 0
	       
	       int chind  = POLY3(ic,jc,kc,cimin,cjmin,ckmin,cnx,cnxny) + n*cnel; // post + state offset

#ifdef DEBUG_CHECK_ASSERTIONS
	       TBOX_ASSERT( cimin <= ic && ic <= cimax );
	       TBOX_ASSERT( cjmin <= jc && jc <= cjmax );
	       TBOX_ASSERT( ckmin <= kc && kc <= ckmax );
#endif
	   
	       double fmass = fvolume[vol_ind];

	       cdata[chind] += fmass*fdata[fhind];  // sum of fine extensives now
	    }
	 }
      }

      //
      // normalize the completed sum by converting back from extensive to intensive
      //
      for ( int k = ifirstc(2); k <= ilastc(2); k++ ) {  // loop over the coarse zones
	 for ( int j = ifirstc(1); j <= ilastc(1); j++ ) {
	    for ( int i = ifirstc(0); i <= ilastc(0); i++ ) {
	       
	       int vol_ind = POLY3(i,j,k,cimin,cjmin,ckmin,cnx,cnxny);
	       int chind   = vol_ind + n*cnel;

	       double cmass = cvolume[vol_ind];

	       cdata[chind] /= cmass + 1.e-80;
	    }
	 }
      }

   }
   tbox::plog << "--------------------- end postprocessCoarsen" << endl;
}


// -------------------------------------------------------------------


/*
*************************************************************************
*                                                                       *
* Tag cells for refinement using gradient detector.  Tagging criteria   *
* defined in input.                                                     *
*                                                                       *
*************************************************************************
*/

// 
// tag cells for refinement
//
void MblkEuler::tagGradientDetectorCells( hier::Patch<NDIM>& patch,
					   const double regrid_time,
					   const bool initial_error,
					   const int tag_indx,
					   const bool uses_richardson_extrapolation_too)
{
   int block_number = getBlockNumber();
   int level_number = patch.getPatchLevelNumber();
   setMappedGridOnPatch(patch, level_number, block_number);

   //
   // get geometry data
   //
   const int error_level_number = patch.getPatchLevelNumber();

   tbox::Pointer< pdat::NodeData<NDIM,double> > xyz = patch.getPatchData( d_xyz, getDataContext() );
   double *x = xyz->getPointer(0);
   double *y = xyz->getPointer(1);
   double *z = xyz->getPointer(2);

   hier::Box<NDIM> pbox           = patch.getBox(); 
   const hier::Index<NDIM> ifirst = patch.getBox().lower();
   const hier::Index<NDIM> ilast  = patch.getBox().upper();
   int level                      = patch.getPatchLevelNumber();
   
   tbox::plog << "--------------------- start tagGradientCells (" << level << ")" << endl;
   tbox::plog << "level  = " << level << endl;
   tbox::plog << "box    = " << patch.getBox() << endl;
   
   tbox::Pointer< pdat::CellData<NDIM,int>    > tags = patch.getPatchData( tag_indx );
   tbox::Pointer< pdat::CellData<NDIM,double> > var  = patch.getPatchData( d_state, getDataContext() );
   
   //
   // Create a set of temporary tags and set to untagged value.
   //
   tbox::Pointer< pdat::CellData<NDIM,int> > temp_tags  = new pdat::CellData<NDIM,int>(pbox, 1, d_nghosts );
   temp_tags->fillAll(FALSE );
	 
   hier::IntVector<NDIM> tag_ghost = tags->getGhostCellWidth();
      
   hier::IntVector<NDIM> nghost_cells = xyz->getGhostCellWidth();
   int nd_imin = ifirst(0)     - nghost_cells(0);
   int nd_imax = ilast(0)  + 1 + nghost_cells(0);
   int nd_jmin = ifirst(1)     - nghost_cells(1);
   int nd_jmax = ilast(1)  + 1 + nghost_cells(1);
   int nd_kmin = ifirst(2)     - nghost_cells(2);
   int nd_nx   = nd_imax - nd_imin + 1;
   int nd_ny   = nd_jmax - nd_jmin + 1;
   int nd_nxny = nd_nx*nd_ny;
      
   hier::IntVector<NDIM> v_ghost = var->getGhostCellWidth();     // has ghost zones
   int imin = ifirst(0) - v_ghost(0); // the polynomial for the field
   int imax = ilast(0)  + v_ghost(0);
   int jmin = ifirst(1) - v_ghost(1);
   int jmax = ilast(1)  + v_ghost(1);
   int kmin = ifirst(2) - v_ghost(2);
   int nx   = imax - imin + 1;
   int ny   = jmax - jmin + 1;
   int nxny = nx*ny;
      
   hier::IntVector<NDIM> temp_tag_ghost = temp_tags->getGhostCellWidth();
   int imn = ifirst(0) - temp_tag_ghost(0);  // the polynomial for temp_tags
   int imx = ilast(0)  + temp_tag_ghost(0);
   int jmn = ifirst(1) - temp_tag_ghost(1);
   int jmx = ilast(1)  + temp_tag_ghost(1);
   int kmn = ifirst(2) - temp_tag_ghost(2);
   int tnx   = imx - imn + 1;
   int tny   = jmx - jmn + 1;
   int tnxny = tnx*tny;

   int *ltags   = temp_tags->getPointer();
   double  dv_x[3];
   double dv_xi[3];
   double xfact = 0.5;


   //
   // Possible tagging criteria includes 
   //    DENSITY_DEVIATION, DENSITY_GRADIENT, DENSITY_SHOCK
   //    PRESSURE_DEVIATION, PRESSURE_GRADIENT, PRESSURE_SHOCK
   // The criteria is specified over a time interval.
   //
   // Loop over criteria provided and check to make sure we are in the
   // specified time interval.  If so, apply appropriate tagging for
   // the level.
   //
   for (int ncrit = 0; ncrit < d_refinement_criteria.getSize(); ncrit++) {

#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT( !var.isNull() );
#endif
      string ref = d_refinement_criteria[ncrit];

      if (ref == "GRADIENT" ) {
	 int nStateLocal = d_state_grad_names.getSize();
	 for ( int id = 0 ; id < nStateLocal ; id++ ) {
	    
	    double *lvar = var->getPointer( d_state_grad_id[id] );

	    int size = d_state_grad_tol[id].getSize();  // max depth of gradient tolerance
	    double tol = ( ( error_level_number < size)  // find the tolerance
			   ? d_state_grad_tol[id][error_level_number] 
			   : d_state_grad_tol[id][size-1]  );
	    
	    for ( int k = ifirst(2) ; k <= ilast(2) ; k++ ) {
	       for ( int j = ifirst(1) ; j <= ilast(1) ; j++ ) {
		  for ( int i = ifirst(0) ; i <= ilast(0) ; i++ ) {
		     
		     int  ind = POLY3(i,j,k,imin,jmin,kmin,nx,nxny);
		     int tind = POLY3(i,j,k,imn,jmn,kmn,tnx,tnxny);
		     
		     int ib = ind-1;
		     int ie = ind+1;
		     int jb = ind-nx;
		     int je = ind+nx;
		     int kb = ind-nxny;
		     int ke = ind+nxny;
		     
		     //
		     // vector is now a max gradient in xi, eta, zeta 
		     //
		     dv_xi[0] = xfact * MAX( fabs( lvar[ind]-lvar[ib] ), fabs( lvar[ind]-lvar[ie] ) );
		     dv_xi[1] = xfact * MAX( fabs( lvar[ind]-lvar[jb] ), fabs( lvar[ind]-lvar[je] ) );
		     dv_xi[2] = xfact * MAX( fabs( lvar[ind]-lvar[kb] ), fabs( lvar[ind]-lvar[ke] ) );
		     
		     //
		     // convert vector to max gradine in x, y, z
		     //
		     int n1 = POLY3(i,j,k,nd_imin,nd_jmin,nd_kmin,nd_nx,nd_nxny);  // -1, -1, -1
		     int n2 = n1 + 1;                                              //  1, -1, -1
		     int n3 = n1 + 1 + nd_nx;                                      //  1,  1, -1
		     int n4 = n1 + nd_nx;                                          // -1,  1, -1
		     
		     int n5 = n1 + nd_nxny;                                        // -1, -1,  1
		     int n6 = n1 + nd_nxny + 1;		   	             //  1, -1,  1
		     int n7 = n1 + nd_nxny + 1 + nd_nx;			     //  1,  1,  1
		     int n8 = n1 + nd_nxny + nd_nx;	         		     // -1,  1,  1
		     
		     // ------------------------------------------------ x
		     
		     double x1 = 0.25*( x[n1] + x[n4] + x[n5] + x[n8] );  // xi
		     double x2 = 0.25*( x[n2] + x[n3] + x[n6] + x[n7] );
		     
		     double x3 = 0.25*( x[n1] + x[n2] + x[n5] + x[n6] );  // eta
		     double x4 = 0.25*( x[n3] + x[n4] + x[n7] + x[n8] );
		     
		     double x5 = 0.25*( x[n1] + x[n2] + x[n3] + x[n4] );  // zeta
		     double x6 = 0.25*( x[n5] + x[n6] + x[n7] + x[n8] );
		     
		     // ------------------------------------------------ y
		     
		     double y1 = 0.25*( y[n1] + y[n4] + y[n5] + y[n8] );
		     double y2 = 0.25*( y[n2] + y[n3] + y[n6] + y[n7] );
		     
		     double y3 = 0.25*( y[n1] + y[n2] + y[n5] + y[n6] );
		     double y4 = 0.25*( y[n3] + y[n4] + y[n7] + y[n8] );
		     
		     double y5 = 0.25*( y[n1] + y[n2] + y[n3] + y[n4] );
		     double y6 = 0.25*( y[n5] + y[n6] + y[n7] + y[n8] );
		     
		     // ------------------------------------------------ z
		     
		     double z1 = 0.25*( z[n1] + z[n4] + z[n5] + z[n8] );
		     double z2 = 0.25*( z[n2] + z[n3] + z[n6] + z[n7] );
		     
		     double z3 = 0.25*( z[n1] + z[n2] + z[n5] + z[n6] );
		     double z4 = 0.25*( z[n3] + z[n4] + z[n7] + z[n8] );
		     
		     double z5 = 0.25*( z[n1] + z[n2] + z[n3] + z[n4] );
		     double z6 = 0.25*( z[n5] + z[n6] + z[n7] + z[n8] );
		     
		     //
		     // the components of the matrices that we want to invert
		     //
		     double dx_xi   = 0.5*(x2 - x1);
		     double dy_xi   = 0.5*(y2 - y1);
		     double dz_xi   = 0.5*(z2 - z1);
		     
		     double dx_eta  = 0.5*(x4 - x3);
		     double dy_eta  = 0.5*(y4 - y3);
		     double dz_eta  = 0.5*(z4 - z3);
		     
		     double dx_zeta = 0.5*(x6 - x5);
		     double dy_zeta = 0.5*(y6 - y5);
		     double dz_zeta = 0.5*(z6 - z5);
		     
		     //
		     // invert dx/dxi as in dx/dxi d/dx = d/dxi, note this
		     // is the transpose of the above matrix M, also via
		     // Kramer's rule
		     //
		     double detMt = ( dx_xi   * dy_eta  * dz_zeta +
				     dx_eta  * dy_zeta * dz_xi   +
				     dx_zeta * dy_xi   * dz_eta  -
				     dx_zeta * dy_eta  * dz_xi   -
				     dx_eta  * dy_xi   * dz_zeta -
				     dx_xi   * dy_zeta * dz_eta  );
		     
		     double detC11 = dy_eta  * dz_zeta - dz_eta  * dy_zeta;
		     double detC21 = dx_eta  * dz_zeta - dz_eta  * dx_zeta;
		     double detC31 = dx_eta  * dy_zeta - dy_eta  * dx_zeta;
		     
		     double detC12 = dy_xi   * dz_zeta - dz_xi   * dy_zeta;
		     double detC22 = dx_xi   * dz_zeta - dz_xi   * dx_zeta;
		     double detC32 = dx_xi   * dy_zeta - dy_xi   * dx_zeta;
		     
		     double detC13 = dy_xi   * dz_eta  - dz_xi   * dy_eta;
		     double detC23 = dx_xi   * dz_eta  - dz_xi   * dx_eta;
		     double detC33 = dx_xi   * dy_eta  - dy_xi   + dx_eta;
		     
		     // -------------------
		     
		     double b11 = detC11/detMt;
		     double b21 = detC21/detMt;
		     double b31 = detC31/detMt;
		     
		     double b12 = detC12/detMt;
		     double b22 = detC22/detMt;
		     double b32 = detC32/detMt;
		     
		     double b13 = detC13/detMt;
		     double b23 = detC23/detMt;
		     double b33 = detC33/detMt;
		     
		     //
		     // determine the maximum gradient in x, y and z (nice orthonormal basis)
		     //
		     dv_x[0] = b11*dv_xi[0] + b12*dv_xi[1] + b13*dv_xi[2];
		     dv_x[1] = b21*dv_xi[0] + b22*dv_xi[1] + b23*dv_xi[2];
		     dv_x[2] = b31*dv_xi[0] + b32*dv_xi[1] + b33*dv_xi[2];
		     
		     double vmax = MAX( dv_x[0], MAX( dv_x[1], dv_x[2] ) );
		     
		     if ( vmax > tol )
			ltags[tind] = TRUE;

		  } // i, j, k loops
	       }
	    }

	 } // end of state position loop
      } // criteria = STATE_GRADIENT

      //
      // For user-defined fixed refinement, access refine box data from the MblkGeometry
      // class.
      //
      if (ref == "USER_DEFINED") {
            
         hier::BoxArray<NDIM> refine_boxes(1); 
         if (d_mblk_geometry->getRefineBoxes(refine_boxes,
                                             block_number,
                                             level_number) ) {
            for (int b = 0; b < refine_boxes.getNumberOfBoxes(); b++) {
               hier::Box<NDIM> intersect = pbox * refine_boxes(b);
               if (!intersect.empty()) {
                  temp_tags->fill(TRUE,intersect);
               }
            }
         }

      } // criteria = USER_DEFINED      
      
   }  // loop over criteria
   
   //
   // Update tags
   //
   for (pdat::CellIterator<NDIM> ic(pbox ); ic; ic++) {
      (*tags)(ic(),0) = (*temp_tags)(ic(),0 );
   }
   
   tbox::plog << "--------------------- end tagGradientCells" << endl;

}


/*
*************************************************************************
*                                                                       *
* Fill the singularity conditions for the multi-block case
*                                                                       *
*************************************************************************
*/
void MblkEuler::fillSingularityBoundaryConditions(
   hier::Patch<NDIM>& patch,
   tbox::List<xfer::MultiblockRefineSchedule<NDIM>::SingularityPatch>&
      sing_patches,
   const double fill_time,
   const hier::Box<NDIM>& fill_box,
   const hier::BoundaryBox<NDIM>& bbox)
{

   NULL_USE(patch);
   NULL_USE(sing_patches);
   NULL_USE(fill_time);
   NULL_USE(fill_box);
   NULL_USE(bbox);
   
   /* 
    * Here we need to specify how to fill the singularities...
    */
#if 0
   for (int i = 0; i < d_variables.getSize(); i++) {

      tbox::Pointer< pdat::CellData<NDIM,double> > cell_data =
         patch.getPatchData(d_variables[i], getDataContext());

      hier::Box<NDIM> sing_fill_box(cell_data->getGhostBox() * fill_box);
      cell_data->fillAll(0.0, sing_fill_box);

      int depth = cell_data->getDepth();

      /*
       * If sing_patches is not empty, that means there is enhanced
       * connectivity, and we get data from other blocks
       */

      if (sing_patches.size()) {

         for (tbox::List
              <xfer::MultiblockRefineSchedule<NDIM>::SingularityPatch>::
              Iterator sp(sing_patches); sp; sp++) {
            tbox::Pointer< pdat::CellData<NDIM,double> > sing_data =
               sp().d_patch->getPatchData(d_variables[i], getDataContext());
            int sing_neighbor_id = sp().d_id;
            for (pdat::CellIterator<NDIM> ci(sing_fill_box); ci; ci++) {
               for (int d = 0; d < depth; d++) {
                  (*cell_data)(ci(),d) += sing_neighbor_id;
               }
            }
         }

         for (pdat::CellIterator<NDIM> ci(sing_fill_box); ci; ci++) {
            for (int d = 0; d < depth; d++) {
               (*cell_data)(ci(),d) /= sing_patches.size();
            }
         }

      /*
       * In cases of reduced connectivity, there are no other blocks
       * from which to acquire data.
       */

      } else {

         cell_data->fillAll(
            (double)bbox.getLocationIndex()+200.0, fill_box);  
 
      }
   }
#endif

}


/*
*************************************************************************
*                                                                       *
* Private method to build XYZ coordinates on a patch                    *
*                                                                       *
*************************************************************************
*/
void MblkEuler::setMappedGridOnPatch(const hier::Patch<NDIM>& patch,
                                      const int level_number,
                                      const int block_number)
{
   //
   // compute level domain
   //
   const tbox::Pointer<geom::BlockPatchGeometry<NDIM> > 
      patch_geom = patch.getPatchGeometry();
   hier::IntVector<NDIM> ratio = patch_geom->getRatio();
   hier::BoxArray<NDIM> domain_boxes;
   d_grid_geometries[block_number]->computePhysicalDomain(domain_boxes, ratio); 

   //
   // statistics on the level domain
   //
   d_dom_current_nboxes = domain_boxes.getNumberOfBoxes();
   
   d_dom_current_bounds[0] = domain_boxes(0).lower(0);
   d_dom_current_bounds[1] = domain_boxes(0).lower(1);
   d_dom_current_bounds[2] = domain_boxes(0).lower(2);
   d_dom_current_bounds[3] = domain_boxes(0).upper(0);
   d_dom_current_bounds[4] = domain_boxes(0).upper(1);
   d_dom_current_bounds[5] = domain_boxes(0).upper(2);
   
   for ( int i = 1 ; i < d_dom_current_nboxes ; i++ ) {
      d_dom_current_bounds[0] = MIN( d_dom_current_bounds[0], domain_boxes(i).lower(0) );
      d_dom_current_bounds[1] = MIN( d_dom_current_bounds[1], domain_boxes(i).lower(1) );
      d_dom_current_bounds[2] = MIN( d_dom_current_bounds[2], domain_boxes(i).lower(2) );
      
      d_dom_current_bounds[3] = MAX( d_dom_current_bounds[3], domain_boxes(i).upper(0) );
      d_dom_current_bounds[4] = MAX( d_dom_current_bounds[4], domain_boxes(i).upper(1) );
      d_dom_current_bounds[5] = MAX( d_dom_current_bounds[5], domain_boxes(i).upper(2) );
   }

   //
   // now build the mesh
   // 
   int xyz_id =  hier::VariableDatabase<NDIM>::getDatabase()->
      mapVariableAndContextToIndex(d_xyz, getDataContext());

   d_mblk_geometry->buildGridOnPatch(patch,
                                     domain_boxes(0),
                                     xyz_id,
                                     level_number,
                                     block_number,
				     d_dom_local_blocks );
}



/*
*************************************************************************
*                                                                       *
* Private method to build the volume on a patch                         *
*                                                                       *
*************************************************************************
*/
void MblkEuler::setVolumeOnPatch( const hier::Patch<NDIM>& patch )
{
   tbox::Pointer< pdat::CellData<NDIM,double> > vol = 
      patch.getPatchData(d_vol, getDataContext());

   tbox::Pointer< pdat::NodeData<NDIM,double> > xyz = 
      patch.getPatchData(d_xyz, getDataContext());

   hier::IntVector<NDIM> vol_ghosts = vol->getGhostCellWidth();
   hier::IntVector<NDIM> xyz_ghosts = xyz->getGhostCellWidth();

   const hier::Index<NDIM> ifirst=patch.getBox().lower();
   const hier::Index<NDIM> ilast =patch.getBox().upper();
   
   int imin = ifirst(0) - vol_ghosts(0);
   int imax = ilast(0)  + vol_ghosts(0);
   int jmin = ifirst(1) - vol_ghosts(1);
   int jmax = ilast(1)  + vol_ghosts(1);
   int kmin = ifirst(2) - vol_ghosts(2);
   int nx   = imax - imin + 1;
   int ny   = jmax - jmin + 1;
   int nxny = nx*ny;
   
   int nd_imin = ifirst(0)     - xyz_ghosts(0);
   int nd_imax = ilast(0)  + 1 + xyz_ghosts(0);
   int nd_jmin = ifirst(1)     - xyz_ghosts(1);
   int nd_jmax = ilast(1)  + 1 + xyz_ghosts(1);
   int nd_kmin = ifirst(2)     - xyz_ghosts(2);
   int nd_nx   = nd_imax - nd_imin + 1;
   int nd_ny   = nd_jmax - nd_jmin + 1;
   int nd_nxny = nd_nx*nd_ny;

   double *cvol = vol->getPointer();
   
   double *x = xyz->getPointer(0);
   double *y = xyz->getPointer(1);
   double *z = xyz->getPointer(2);
   
   for ( int k = ifirst(2); k <= ilast(2); k++ ) {
      for ( int j = ifirst(1); j <= ilast(1); j++ ) {
	 for ( int i = ifirst(0); i <= ilast(0); i++ ) {
	    
	    int cind = POLY3(i,j,k,imin,jmin,kmin,nx,nxny);
	    
	    int n1 = POLY3(i,j,k,nd_imin,nd_jmin,nd_kmin,nd_nx,nd_nxny);
	    int n2 = n1 + 1;
	    int n3 = n1 + 1 + nd_nx;
	    int n4 = n1 + nd_nx;
	    
	    int n5 = n1 + nd_nxny;
	    int n6 = n1 + nd_nxny + 1;
	    int n7 = n1 + nd_nxny + 1 + nd_nx;
	    int n8 = n1 + nd_nxny + nd_nx;
	    
	    double lvol = UpwindVolume( x[n1], x[n2], x[n3], x[n4],
					x[n5], x[n6], x[n7], x[n8],
					
					y[n1], y[n2], y[n3], y[n4],
					y[n5], y[n6], y[n7], y[n8],
					
					z[n1], z[n2], z[n3], z[n4],
					z[n5], z[n6], z[n7], z[n8] );
	    
	    cvol[cind] = lvol;
	 }
      }
   }
}


// ================================= MblkEuler::Visualization and IO =============================


/*
*************************************************************************
*                                                                       *
* Register VisIt data writer to write data to plot files that may       *
* be postprocessed by the VisIt tool.                                   *
*                                                                       *
*************************************************************************
*/

#ifdef HAVE_HDF5
void MblkEuler::registerVisItDataWriter(
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
* Write MblkEuler object state to specified stream.                        *
*                                                                       *
*************************************************************************
*/

void MblkEuler::printClassData(ostream &os) const 
{
   int j;

   os << "\nMblkEuler::printClassData..." << endl;
   os << "MblkEuler: this = " << (MblkEuler*)this << endl;
   os << "d_object_name = " << d_object_name << endl;
   os << "d_grid_geometries = " << endl;
   for (j=0; j < d_grid_geometries.getSize(); j++) {
//      os << (geom::BlockGridGeometry<NDIM>*)d_grid_geometries[j] << endl;
   }

   // ----------------------------------------------

   os << "Parameters for numerical method ..." << endl;
   os << "   d_advection_velocity = ";
   for (j=0;j<NDIM;j++) os << d_advection_velocity[j] << " ";
   os << endl;

   os << "   d_nghosts    = " << d_nghosts << endl;
   os << "   d_fluxghosts = " << d_fluxghosts << endl;

   os << "Problem description and initial data..." << endl;
   os << "   d_data_problem = " << d_data_problem << endl;
}


/*
*************************************************************************
*                                                                       *
* Read data members from input.  All values set from restart can be	*
* overridden by values in the input database.
*                                                                       *
*************************************************************************
*/
void MblkEuler::getFromInput( 
   tbox::Pointer<tbox::Database> input_db,
   bool is_from_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!input_db.isNull());
#endif

   tbox::Pointer<tbox::Database> db = input_db->getDatabase("MblkEuler");   

   //
   // --------------- load balancing inputs
   //
   // Note: if we are restarting, then we only allow nonuniform
   // workload to be used if nonuniform workload was used originally.
   //
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

   //
   // --------------- initialize boundary condition factors
   //
   if ( db->keyExists( "wall_factors" )) {
      d_wall_factors = db->getIntegerArray( "wall_factors" );
   } else {
      d_wall_factors.resizeArray(6);
      for ( int i = 0 ; i < 6 ; i++ ) d_wall_factors[i] = 0;
   }

   //
   // --------------- process the linear advection test ---------------------
   //
   d_advection_test = 0;
   d_advection_velocity[0] = d_advection_velocity[1] = d_advection_velocity[2] = FLT_MAX; 
   d_advection_vel_type = 0;
   if ( db->keyExists( "advection_test" ) ) {
      d_advection_test = db->getInteger( "advection_test" );
      if (db->keyExists("advection_velocity")){
	 db->getDoubleArray("advection_velocity",
			    d_advection_velocity, NDIM);
	 d_advection_vel_type = db->getInteger("advection_vel_type");
      } else {
	 TBOX_ERROR(d_object_name << ":  "
		    << "Key data `advection_velocity' not found in input.");
      } 
   }

   //
   // --------------- The state names inputs ---------------------
   //
   if ( d_advection_test ) {
      if ( db->keyExists("state_names") ) {
	 d_state_names = db->getStringArray( "state_names" );
	 d_nState      = d_state_names.getSize();
      } else {
	 TBOX_ERROR( "missing 'state_names' input for sizing the state" << endl );
      }
   }
   else {
      TBOX_ASSERT( d_advection_test );
   }

   //
   // --------------- region Initialization inputs ---------------------
   //
   if (!is_from_restart) { 
 
      if (db->keyExists("data_problem")) {
         d_data_problem = db->getString("data_problem");
      } else {
         TBOX_ERROR(d_object_name << ": "
                  << "`data_problem' value not found in input." << endl);
      }

      int problem_1d = 1;

      //
      //  axis of revolution inputs
      //
      if ( d_data_problem == "REVOLUTION" ) {
	 
	 if (db->keyExists( "center" )) {
	    db->getDoubleArray( "center", d_center, NDIM );
	 } else {
	    TBOX_ERROR( "`center' input required for REVOLUTION problem." << endl );
	 }
	 
	 if (db->keyExists( "axis" )) {
	    db->getDoubleArray( "axis", d_axis, NDIM );
	 } else {
	    TBOX_ERROR( "`axis' input required for REVOLUTION problem." << endl );
	 }
	 
	 // normalize the axis to a unit vector
	 double anorm = sqrt( d_axis[0]*d_axis[0] + d_axis[1]*d_axis[1] + d_axis[2]*d_axis[2] );
	 d_axis[0] /= anorm;
	 d_axis[1] /= anorm;
	 d_axis[2] /= anorm;

	 d_rev_rad.resizeArray(  d_number_of_regions );
	 d_rev_axis.resizeArray( d_number_of_regions );

	 for ( int i = 0 ; i < d_number_of_regions; i++ ) {
	    
	    char tmp[20];
	    sprintf( tmp, "region_%d", i+1 );  // 
 	    string lkey = tmp;
	    tbox::Pointer<tbox::Database> region_db = db->getDatabase(lkey);
	    
	    d_rev_rad[i]  = region_db->getDoubleArray( "radius" );
	    d_rev_axis[i] = region_db->getDoubleArray( "axis"   );
	    
	 }
	 
	 problem_1d = 0;
      }
      
      //
      //  Spherical inputs
      //
      if ( d_data_problem == "SPHERE" ) {
	 
	 if (db->keyExists( "center" )) {
	    db->getDoubleArray( "center", d_center, NDIM );
	 } else {
	    TBOX_ERROR( "`center' input required for SPHERE problem." << endl );
	 }
      }

      //
      //  Rayleigh tayler shock tube inputs
      //
      if ( d_data_problem == "RT_SHOCK_TUBE" ) {
	 
	 if (db->keyExists( "amn" )) {
	    d_dt_ampl = db->getDouble( "ampl" );
	    d_amn     = db->getDoubleArray( "amn"    );
	    d_m_mode  = db->getDoubleArray( "m_mode" );
	    d_n_mode  = db->getDoubleArray( "n_mode" );
	    d_phiy    = db->getDoubleArray( "phiy"   );
	    d_phiz    = db->getDoubleArray( "phiz"   );
	 } else {
	    TBOX_ERROR( "missing input for RT_SHOCK_TUBE problem." << endl );
	 }
      }
      
      //
      //  shared inputs for the 1d style problems
      //
      if ( problem_1d ) {
	 if (db->keyExists( "front_position" )) {
	    d_front_position = db->getDoubleArray( "front_position" );
	    d_number_of_regions = d_front_position.getSize()-1;
	    TBOX_ASSERT( d_number_of_regions > 0 );
	 } else {
	    TBOX_ERROR( "Missing`front_position' input required" << endl );
	 }
      }

      //
      // the state data entries for the initial conditions
      //
      int llen = d_number_of_regions*d_nState;
      double *tmp = new double[llen];
      for ( int ii = 0 ; ii < llen ; ii++ ) {
	 tmp[ii] = FLT_MAX;
      }

      d_state_ic  = new double *[d_number_of_regions];
      for ( int iReg = 0 ; iReg < d_number_of_regions; iReg++ ) {
	 d_state_ic[iReg] = &tmp[d_nState*iReg];
      } 

      //
      // pull in the data for each region
      //
      if ( d_advection_test ) {
	 if (db->keyExists( "state_data" )) {
	    tbox::Pointer<tbox::Database> state_db = db->getDatabase("state_data");
	    tbox::Array<double> lpsi;
	    for ( int iState = 0 ; iState < d_nState; iState++ ) {
	       lpsi = state_db->getDoubleArray( d_state_names[iState] );
	       for ( int iReg = 0 ; iReg < d_number_of_regions; iReg++ ) {
		  d_state_ic[iReg][iState] = lpsi[iReg];
	       }
	    }
	 } else {
	    TBOX_ERROR( "missing 'state_data' input for initial conditions" << endl );
	 }
      }
      else {
	 TBOX_ASSERT( d_advection_test );
      }

   } // if !is_from_restart read in problem data

   //
   //  --------------- refinement criteria inputs
   //
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
            
	    //
	    // allow only valid refinement criteria as the remaining keys
	    //
            if ( !(error_key == "GRADIENT") ) { 

               TBOX_ERROR( "Unknown refinement criteria: " << error_key
			   << " in input." << endl );

            } else {
               error_db = refine_db->getDatabase(error_key );
               ref_keys_defined[def_key_cnt] = error_key;
               def_key_cnt++;
            }
               
	    //
	    // process the specific keys
	    //
            if (!error_db.isNull() && error_key == "GRADIENT" ) {

	       d_state_grad_names = error_db->getStringArray( "names" );
	       int nStateLocal    = d_state_grad_names.getSize();

	       d_state_grad_tol.resizeArray(nStateLocal);
	       d_state_grad_id.resizeArray(nStateLocal);

	       for ( int id = 0 ; id < nStateLocal; id++ ) {
		  string grad_name =  d_state_grad_names[id];

		  // ... the index needed
		  d_state_grad_id[id] = -1;
		  bool found = false;
		  for ( int idj = 0 ; idj < d_nState && !found ; idj++ ) {
		     if ( grad_name == d_state_names[idj] ) {
			found = true;
			d_state_grad_id[id] = idj;
		     }
		  }

		  // ... the tolerance array needed
		  if (error_db->keyExists( grad_name ) ) {
		     d_state_grad_tol[id] = error_db->getDoubleArray( grad_name );
		  } else {
		     TBOX_ERROR( "No tolerance array " << grad_name << "found for gradient detector" << endl );
		  }
	       
	       }
	    }
	    

         } // refine criteria if test
         
      } // loop over refine criteria
   
   } // refine db entry exists 


   //
   // --------------- boundary condition inputs ---------------------
   //
   hier::IntVector<NDIM> periodic = d_grid_geometries[0]->getPeriodicShift();
   int num_per_dirs = 0;
   for (int id = 0; id < NDIM; id++) {
      if (periodic(id)) num_per_dirs++;
   }

   /*
    * If there are multiple blocks, periodicity is not currently supported.
    */
   if ((d_grid_geometries.getSize() > 1) && (num_per_dirs > 0)) {
      TBOX_ERROR(d_object_name << ": cannot have periodic BCs when there"
                 << "\nare multiple blocks." << endl);
   }
}


/*
*************************************************************************
*                                                                       *
* Routines to put/get data members to/from restart database.            *
*                                                                       *
*************************************************************************
*/

void MblkEuler::putToDatabase( 
   tbox::Pointer<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!db.isNull());
#endif

   db->putInteger("MBLKEULER_VERSION",MBLKEULER_VERSION);

   db->putDoubleArray("d_advection_velocity",d_advection_velocity,NDIM);

   db->putIntegerArray("d_nghosts", (int*)d_nghosts, NDIM);
   db->putIntegerArray("d_fluxghosts", (int*)d_fluxghosts, NDIM);

   db->putString("d_data_problem", d_data_problem);
}


/*
*************************************************************************
*                                                                       *
*    Access class information from restart database.                    *
*                                                                       *
*************************************************************************
*/
void MblkEuler::getFromRestart() 
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

   int ver = db->getInteger("MBLKEULER_VERSION");
   if (ver != MBLKEULER_VERSION){
      TBOX_ERROR(d_object_name << ":  "
              << "Restart file version different than class version.");
   }

   d_data_problem = db->getString("d_data_problem");
}


/*
*************************************************************************
*                                                                       *
* Routine to check boundary data when debugging.                        *
*                                                                       *
*************************************************************************
*/

void MblkEuler::checkBoundaryData(int btype,
                              const hier::Patch<NDIM>& patch,
                              const hier::IntVector<NDIM>& ghost_width_to_check,
                              const tbox::Array<int>& scalar_bconds) const
{
}

hier::IntVector<NDIM> MblkEuler::getMultiblockRefineOpStencilWidth() const
{
   return (hier::IntVector<NDIM>(1));
}
                                                                                
hier::IntVector<NDIM> MblkEuler::getMultiblockCoarsenOpStencilWidth()
{
   return (hier::IntVector<NDIM>(0));
}



