#ifndef included_solv_CartesianRobinBcHelper_C
#define included_solv_CartesianRobinBcHelper_C

/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/solvers/poisson/CartesianRobinBcHelper.C $
 * Package:     SAMRAI application utilities
 * Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 2147 $
 * Modified:    $LastChangedDate: 2008-04-23 16:48:12 -0700 (Wed, 23 Apr 2008) $
 * Description: Robin boundary condition support on cartesian grids.
 */


#include "CartesianPatchGeometry.h"
#include "Index.h"
#include "VariableDatabase.h"
#include "PatchCellDataOpsReal.h"
#include "CellVariable.h"
#include "tbox/Array.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"

#include IOMANIP_HEADER_FILE

#include "CartesianRobinBcHelper.h"


extern "C" {
void settype1cells2d_(
   double *data,
   const int &difirst, const int &dilast,
   const int &djfirst, const int &djlast,
   const double *a, const double *b, const double *g,
   const int &ifirst, const int &ilast,
   const int &jfirst, const int &jlast,
   const int &ibeg, const int &iend,
   const int &jbeg, const int &jend,
   const int &face, const int &ghos, const int &inte, const int &location,
   const double &h, const int &zerog );
void settype2cells2d_(
   double *data,
   const int &difirst, const int &dilast,
   const int &djfirst, const int &djlast,
   const int *lower, const int *upper, const int &location );
void settype1cells3d_(
   double *data,
   const int &difirst, const int &dilast,
   const int &djfirst, const int &djlast,
   const int &dkfirst, const int &dklast,
   const double *a, const double *b, const double *g,
   const int &ifirst, const int &ilast,
   const int &jfirst, const int &jlast,
   const int &kfirst, const int &klast,
   const int &ibeg, const int &iend,
   const int &jbeg, const int &jend,
   const int &kbeg, const int &kend,
   const int &face, const int &ghos, const int &inte, const int &location,
   const double &h, const int &zerog );
void settype2cells3d_(
   double *data,
   const int &difirst, const int &dilast,
   const int &djfirst, const int &djlast,
   const int &dkfirst, const int &dklast,
   const int *lower, const int *upper, const int &location );
void settype3cells3d_(
   double *data,
   const int &difirst, const int &dilast,
   const int &djfirst, const int &djlast,
   const int &dkfirst, const int &dklast,
   const int *lower, const int *upper, const int &location );
}

namespace SAMRAI {
    namespace solv {



/*
************************************************************************
* Constructor                                                          *
************************************************************************
*/

template<int DIM>  CartesianRobinBcHelper<DIM>::CartesianRobinBcHelper(
   std::string object_name,
   RobinBcCoefStrategy<DIM> *coef_strategy)
  : d_object_name(object_name),
    d_coef_strategy(NULL),
    d_target_data_id(-1),
    d_homogeneous_bc(false)
{

   NULL_USE(coef_strategy);

   t_set_boundary_values_in_cells = tbox::TimerManager::getManager()->
      getTimer("solv::CartesianRobinBcHelper::setBoundaryValuesInCells()");
   t_use_set_bc_coefs = tbox::TimerManager::getManager()->
      getTimer("solv::CartesianRobinBcHelper::setBoundaryValuesInCells()_setBcCoefs");

   return;
}



/*
************************************************************************
* Destructor                                                           *
************************************************************************
*/

template<int DIM>  CartesianRobinBcHelper<DIM>::~CartesianRobinBcHelper(void) {}





/*
************************************************************************
*     Set physical boundary conditions in cells.                       *
************************************************************************
*/

template<int DIM> void CartesianRobinBcHelper<DIM>::setBoundaryValuesInCells (
   hier::Patch<DIM> &patch ,
   const double fill_time ,
   const hier::IntVector<DIM> &ghost_width_to_fill ,
   int target_data_id ,
   bool homogeneous_bc ) const
{
   NULL_USE(fill_time);

   t_set_boundary_values_in_cells->start();

#ifdef DEBUG_CHECK_ASSERTIONS
   if ( !d_coef_strategy ) {
      TBOX_ERROR(d_object_name << ": coefficient strategy is not set.\n"
                 << "Use setCoefImplementation() to set it.\n");
   }
#endif

   if (DIM == 1) {
      TBOX_ERROR(d_object_name << ": DIM = 1 not supported");
   }
   math::PatchCellDataOpsReal<DIM,double> cops;

   /*
    * Get info on the data.
    */
   hier::VariableDatabase<DIM> *vdb =
      hier::VariableDatabase<DIM>::getDatabase();
   tbox::Pointer< hier::Variable<DIM> > variable_ptr;
   vdb->mapIndexToVariable( target_data_id, variable_ptr );
   tbox::Pointer<pdat::CellVariable<DIM,double> > cell_variable_ptr
      = variable_ptr;
   if ( !variable_ptr ) {
      TBOX_ERROR(d_object_name << ": No variable for index "
                 << target_data_id );
   }
   if ( !cell_variable_ptr ) {
      TBOX_ERROR(d_object_name << ": hier::Patch data index " << target_data_id
                 << " does not correspond to\n"
                 << "a cell-centered double variable.\n");
   }

   /*
    * Get the data.
    */
   tbox::Pointer< hier::PatchData<DIM> >
      data_ptr = patch.getPatchData(target_data_id);
   tbox::Pointer<pdat::CellData<DIM,double> >
      cell_data_ptr = data_ptr;
   if ( !data_ptr ) {
      TBOX_ERROR(d_object_name << ": No data for index " << target_data_id );
   }
   if ( !cell_data_ptr ) {
      TBOX_ERROR(d_object_name << ": hier::Patch data index " << target_data_id
                 << " does not correspond to\n"
                 << "cell-centered double data.\n");
   }
   pdat::CellData<DIM,double> &data = *cell_data_ptr;

   const hier::IntVector<DIM> &ghost_cells =
      cell_data_ptr->getGhostCellWidth();
   hier::IntVector<DIM> gcw_to_fill = hier::IntVector<DIM>::min(ghost_cells,
                                                      ghost_width_to_fill);
   if ( gcw_to_fill == hier::IntVector<DIM>(0) ) {
      return;
   }


   /*
    * Given a and g in a*u + (1-a)*un = g,
    * where un is the derivative in the outward normal direction,
    * and ui (the value of u in the first interior cell),
    * we compute the value on the outer face
    * uf = ...
    * and the normal derivative on the outer face
    * un = ...
    * and the uo (the value in the first ghost cell)
    * uo = ...
    */
   const hier::Box<DIM> &patch_box(patch.getBox());

   /*
    * These definitions can go in the next block.
    * They are kept her for debugging.
    */
   tbox::Pointer< geom::CartesianPatchGeometry<DIM> > pg =
      patch.getPatchGeometry();

   const tbox::Array< hier::BoundaryBox<DIM> > &codim1_boxes =
      pg->getCodimensionBoundaries(1);

   const int n_codim1_boxes=codim1_boxes.getSize();


   const hier::Box<DIM> &ghost_box = data.getGhostBox();
   const double *h = pg->getDx();
   const int num_coefs( homogeneous_bc ? 1 : 2 );
   const int zerog = num_coefs == 1;

   for (int n=0; n<n_codim1_boxes; ++n ) {

      const int location_index = codim1_boxes[n].getLocationIndex();
      const int normal_dir = location_index/2;
      if ( ! gcw_to_fill(normal_dir) ) {
         // Zero ghost width to fill for this boundary box.
         continue;
      }
      hier::IntVector<DIM> extension_amount(1);
      extension_amount(normal_dir) = 0;
      const hier::BoundaryBox<DIM> boundary_box = 
         d_coef_strategy->numberOfExtensionsFillable() >= extension_amount ?
         trimBoundaryBox(codim1_boxes[n],ghost_box) :
         trimBoundaryBox(codim1_boxes[n],patch_box);
      const hier::Index<DIM> &lower = boundary_box.getBox().lower();
      const hier::Index<DIM> &upper = boundary_box.getBox().upper();
      const hier::Box<DIM> coefbox = makeFaceBoundaryBox(boundary_box);
      tbox::Pointer<pdat::ArrayData<DIM,double> >
        acoef_data = new pdat::ArrayData<DIM,double>( coefbox, 1 ),
        bcoef_data = new pdat::ArrayData<DIM,double>( coefbox, 1 ),
        gcoef_data = homogeneous_bc ? NULL :
                     new pdat::ArrayData<DIM,double>( coefbox, 1 );
      t_use_set_bc_coefs->start();
      d_coef_strategy->setBcCoefs( acoef_data,
                                   bcoef_data,
                                   gcoef_data,
                                   variable_ptr,
                                   patch,
                                   boundary_box ,
                                   fill_time );
      t_use_set_bc_coefs->stop();

      int igho, ifac, iint, ibeg, iend;
      double dx;
      int jgho, jfac, jint, jbeg, jend;
      double dy;
      int kgho, kfac, kint, kbeg, kend;
      double dz;


      if (DIM == 2) {
	 switch (location_index) {
	    case 0:
	       // min i edge
	       dx = h[0];
	       igho = lower[0]; // Lower and upper are the same.
	       ifac = igho + 1;
	       iint = igho+1;
	       jbeg = lower[1]; jend = upper[1];
	       settype1cells2d_( data.getPointer(0),
				 ghost_box.lower()[0], ghost_box.upper()[0],
				 ghost_box.lower()[1], ghost_box.upper()[1],
				 acoef_data->getPointer(),
				 bcoef_data->getPointer(),
				 gcoef_data?gcoef_data->getPointer():NULL,
				 coefbox.lower()[0], coefbox.upper()[0],
				 coefbox.lower()[1], coefbox.upper()[1],
				 igho, igho, jbeg, jend,
				 ifac, igho, iint, location_index, dx, zerog
		  );
	       break;
	    case 1:
	       // max i edge
	       dx = h[0];
	       igho = lower[0]; // Lower and upper are the same.
	       ifac = igho;
	       iint = igho-1;
	       jbeg = lower[1]; jend = upper[1];
	       settype1cells2d_( data.getPointer(0),
				 ghost_box.lower()[0], ghost_box.upper()[0],
				 ghost_box.lower()[1], ghost_box.upper()[1],
				 acoef_data->getPointer(),
				 bcoef_data->getPointer(),
				 gcoef_data?gcoef_data->getPointer():NULL,
				 coefbox.lower()[0], coefbox.upper()[0],
				 coefbox.lower()[1], coefbox.upper()[1],
				 igho, igho, jbeg, jend,
				 ifac, igho, iint, location_index, dx, zerog
		  );
	       break;
	    case 2:
	       // min j edge
	       dy = h[1];
	       jgho = lower[1]; // Lower and upper are the same.
	       jfac = jgho + 1;
	       jint = jgho+1;
	       ibeg = lower[0]; iend = upper[0];
	       settype1cells2d_( data.getPointer(0),
				 ghost_box.lower()[0], ghost_box.upper()[0],
				 ghost_box.lower()[1], ghost_box.upper()[1],
				 acoef_data->getPointer(),
				 bcoef_data->getPointer(),
				 gcoef_data?gcoef_data->getPointer():NULL,
				 coefbox.lower()[0], coefbox.upper()[0],
				 coefbox.lower()[1], coefbox.upper()[1],
				 ibeg, iend, jgho, jgho,
				 jfac, jgho, jint, location_index, dy, zerog
		  );
	       break;
	    case 3:
	       // max j edge
	       dy = h[1];
	       jgho = lower[1]; // Lower and upper are the same.
	       jfac = jgho;
	       jint = jgho-1;
	       ibeg = lower[0]; iend = upper[0];
	       settype1cells2d_( data.getPointer(0),
				 ghost_box.lower()[0], ghost_box.upper()[0],
				 ghost_box.lower()[1], ghost_box.upper()[1],
				 acoef_data->getPointer(),
				 bcoef_data->getPointer(),
				 gcoef_data?gcoef_data->getPointer():NULL,
				 coefbox.lower()[0], coefbox.upper()[0],
				 coefbox.lower()[1], coefbox.upper()[1],
				 ibeg, iend, jgho, jgho,
				 jfac, jgho, jint, location_index, dy, zerog
		  );
	       break;
	    default:
	       TBOX_ERROR(d_object_name << ": Invalid location index ("
			  << location_index << ") in\n"
			  << "setBoundaryValuesInCells" );
	 }
      } else if (DIM == 3) {
	 switch (location_index) {
	    case 0:
	       // min i face
	       dx = h[0];
	       igho = lower[0]; // Lower and upper are the same.
	       ifac = igho + 1;
	       iint = igho+1;
	       jbeg = lower[1]; jend = upper[1];
	       kbeg = lower[2]; kend = upper[2];
	       settype1cells3d_( data.getPointer(0),
				 ghost_box.lower()[0], ghost_box.upper()[0],
				 ghost_box.lower()[1], ghost_box.upper()[1],
				 ghost_box.lower()[2], ghost_box.upper()[2],
				 acoef_data->getPointer(),
				 bcoef_data->getPointer(),
				 gcoef_data?gcoef_data->getPointer():NULL,
				 coefbox.lower()[0], coefbox.upper()[0],
				 coefbox.lower()[1], coefbox.upper()[1],
				 coefbox.lower()[2], coefbox.upper()[2],
				 igho, igho, jbeg, jend, kbeg, kend,
				 ifac, igho, iint, location_index, dx, zerog
		  );
	       break;
	    case 1:
	       // max i face
	       dx = h[0];
	       igho = lower[0]; // Lower and upper are the same.
	       ifac = igho;
	       iint = igho-1;
	       jbeg = lower[1]; jend = upper[1];
	       kbeg = lower[2]; kend = upper[2];
	       settype1cells3d_( data.getPointer(0),
				 ghost_box.lower()[0], ghost_box.upper()[0],
				 ghost_box.lower()[1], ghost_box.upper()[1],
				 ghost_box.lower()[2], ghost_box.upper()[2],
				 acoef_data->getPointer(),
				 bcoef_data->getPointer(),
				 gcoef_data?gcoef_data->getPointer():NULL,
				 coefbox.lower()[0], coefbox.upper()[0],
				 coefbox.lower()[1], coefbox.upper()[1],
				 coefbox.lower()[2], coefbox.upper()[2],
				 igho, igho, jbeg, jend, kbeg, kend,
                           ifac, igho, iint, location_index, dx, zerog
		  );
	       break;
	    case 2:
	       // min j face
	       dy = h[1];
	       jgho = lower[1]; // Lower and upper are the same.
	       jfac = jgho + 1;
	       jint = jgho+1;
	       ibeg = lower[0]; iend = upper[0];
	       kbeg = lower[2]; kend = upper[2];
	       settype1cells3d_( data.getPointer(0),
				 ghost_box.lower()[0], ghost_box.upper()[0],
				 ghost_box.lower()[1], ghost_box.upper()[1],
				 ghost_box.lower()[2], ghost_box.upper()[2],
				 acoef_data->getPointer(),
				 bcoef_data->getPointer(),
				 gcoef_data?gcoef_data->getPointer():NULL,
				 coefbox.lower()[0], coefbox.upper()[0],
				 coefbox.lower()[1], coefbox.upper()[1],
				 coefbox.lower()[2], coefbox.upper()[2],
				 ibeg, iend, jgho, jgho, kbeg, kend,
				 jfac, jgho, jint, location_index, dy, zerog
		  );
	       break;
	    case 3:
	       // max j face
	       dy = h[1];
	       jgho = lower[1]; // Lower and upper are the same.
	       jfac = jgho;
	       jint = jgho-1;
	       ibeg = lower[0]; iend = upper[0];
	       kbeg = lower[2]; kend = upper[2];
	       settype1cells3d_( data.getPointer(0),
				 ghost_box.lower()[0], ghost_box.upper()[0],
				 ghost_box.lower()[1], ghost_box.upper()[1],
				 ghost_box.lower()[2], ghost_box.upper()[2],
				 acoef_data->getPointer(),
				 bcoef_data->getPointer(),
				 gcoef_data?gcoef_data->getPointer():NULL,
				 coefbox.lower()[0], coefbox.upper()[0],
				 coefbox.lower()[1], coefbox.upper()[1],
				 coefbox.lower()[2], coefbox.upper()[2],
				 ibeg, iend, jgho, jgho, kbeg, kend,
				 jfac, jgho, jint, location_index, dy, zerog
		  );
	       break;
	    case 4:
	       // min k face
	       dz = h[2];
	       kgho = lower[2]; // Lower and upper are the same.
	       kfac = kgho + 1;
	       kint = kgho+1;
	       ibeg = lower[0]; iend = upper[0];
	       jbeg = lower[1]; jend = upper[1];
	       settype1cells3d_( data.getPointer(0),
				 ghost_box.lower()[0], ghost_box.upper()[0],
				 ghost_box.lower()[1], ghost_box.upper()[1],
				 ghost_box.lower()[2], ghost_box.upper()[2],
				 acoef_data->getPointer(),
				 bcoef_data->getPointer(),
				 gcoef_data?gcoef_data->getPointer():NULL,
				 coefbox.lower()[0], coefbox.upper()[0],
				 coefbox.lower()[1], coefbox.upper()[1],
				 coefbox.lower()[2], coefbox.upper()[2],
				 ibeg, iend, jbeg, jend, kgho, kgho,
				 kfac, kgho, kint, location_index, dz, zerog
		  );
	       break;
	    case 5:
	       // max k face
	       dz = h[2];
	       kgho = lower[2]; // Lower and upper are the same.
	       kfac = kgho;
	       kint = kgho-1;
	       ibeg = lower[0]; iend = upper[0];
	       jbeg = lower[1]; jend = upper[1];
	       settype1cells3d_( data.getPointer(0),
				 ghost_box.lower()[0], ghost_box.upper()[0],
				 ghost_box.lower()[1], ghost_box.upper()[1],
				 ghost_box.lower()[2], ghost_box.upper()[2],
				 acoef_data->getPointer(),
				 bcoef_data->getPointer(),
				 gcoef_data?gcoef_data->getPointer():NULL,
				 coefbox.lower()[0], coefbox.upper()[0],
				 coefbox.lower()[1], coefbox.upper()[1],
				 coefbox.lower()[2], coefbox.upper()[2],
				 ibeg, iend, jbeg, jend, kgho, kgho,
				 kfac, kgho, kint, location_index, dz, zerog
		  );
	       break;
	    default:
	       TBOX_ERROR(d_object_name << ": Invalid location index ("
			  << location_index << ") in\n"
			  << "setBoundaryValuesInCells" );
	 }
      }
   }


   /*
    * Now that the surface boundary boxes have been set,
    * the rest of this function set the lower-dimensional
    * boundary boxes.  Users may not need to have these
    * set, but refiners may.
    */


   if (DIM == 2) {
      /*
       * The node boundary conditions are set from a linear interpolation
       * through the nearest interior cell and the two nearest edge values.
       * This data may be used by refinement operators to do interpolation.
       */
      
      const tbox::Array< hier::BoundaryBox<DIM> > &node_boxes =
         pg->getNodeBoundaries();
      const int n_node_boxes=node_boxes.getSize();
      for ( int n=0; n<n_node_boxes; ++n ) {
	 const hier::BoundaryBox<DIM> &bb = node_boxes[n];
#ifdef DEBUG_CHECK_ASSERTIONS
	 TBOX_ASSERT( bb.getBoundaryType() == 2 );	// Must be a node boundary.
#endif
	 const hier::Box<DIM> &bb_box = bb.getBox();
	 const hier::Index<DIM> &lower = bb_box.lower();
	 const hier::Index<DIM> &upper = bb_box.upper();
	 const int location_index = bb.getLocationIndex();
      settype2cells2d_( data.getPointer(0),
                        ghost_box.lower()[0], ghost_box.upper()[0],
                        ghost_box.lower()[1], ghost_box.upper()[1],
                        lower, upper, location_index );
      }
   } else if (DIM == 3) {
      /*
       * The edge boundary conditions are set from a linear interpolation
       * through the nearest interior cell and the two nearest side values.
       * This data may be used by refinement operators to do interpolation.
       */
      const tbox::Array< hier::BoundaryBox<DIM> > &edge_boxes =
         pg->getEdgeBoundaries();
      const int n_edge_boxes=edge_boxes.getSize();
      for ( int n=0; n<n_edge_boxes; ++n ) {
	 const int location_index = edge_boxes[n].getLocationIndex();
	 const int edge_dir = location_index/4;
	 hier::IntVector<DIM> extension_amount(0);
	 extension_amount(edge_dir) = 1;
	 const hier::BoundaryBox<DIM> boundary_box = 
	    d_coef_strategy->numberOfExtensionsFillable() >= extension_amount ?
	    trimBoundaryBox(edge_boxes[n],ghost_box) :
            trimBoundaryBox(edge_boxes[n],patch_box);
#ifdef DEBUG_CHECK_ASSERTIONS
	 TBOX_ASSERT( boundary_box.getBoundaryType() == 2 );
#endif
	 const hier::Index<DIM> &lower = boundary_box.getBox().lower();
	 const hier::Index<DIM> &upper = boundary_box.getBox().upper();
	 settype2cells3d_( data.getPointer(0),
			   ghost_box.lower()[0], ghost_box.upper()[0],
			   ghost_box.lower()[1], ghost_box.upper()[1],
			   ghost_box.lower()[2], ghost_box.upper()[2],
			   lower, upper, location_index );
      }

      /*
       * The node boundary conditions are set from a linear interpolation
       * through the nearest interior cell and the three nearest edge values.
       * This data may be used by refinement operators to do interpolation.
       */
      const tbox::Array< hier::BoundaryBox<DIM> > &node_boxes = pg->getNodeBoundaries();
      const int n_node_boxes=node_boxes.getSize();
      for ( int n=0; n<n_node_boxes; ++n ) {
	 const hier::BoundaryBox<DIM> &bb = node_boxes[n];
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT( bb.getBoundaryType() == 3 );	// Must be an node boundary.
#endif
      const hier::Box<DIM> &bb_box = bb.getBox();
      const hier::Index<DIM> &lower = bb_box.lower();
      const hier::Index<DIM> &upper = bb_box.upper();
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT( lower == upper );
#endif
      const int location_index = bb.getLocationIndex();
      settype3cells3d_( data.getPointer(0),
                        ghost_box.lower()[0], ghost_box.upper()[0],
                        ghost_box.lower()[1], ghost_box.upper()[1],
                        ghost_box.lower()[2], ghost_box.upper()[2],
                        lower, upper, location_index );
      }
   } else {
      TBOX_ERROR("CartesianRobinBcHelper<DIM>::setBoundaryValuesInCells error ..."
		 << "\n not implemented for DIM>3" << std::endl);
   }

   
   t_set_boundary_values_in_cells->stop();

   return;
}





/*
************************************************************************
* Set physical boundary conditions in cells, for all patches in a      *
* given level.                                                         *
************************************************************************
*/

template<int DIM> void CartesianRobinBcHelper<DIM>::setBoundaryValuesInCells (
   hier::PatchLevel<DIM> &level ,
   const double fill_time ,
   const hier::IntVector<DIM> &ghost_width_to_fill ,
   int target_data_id ,
   bool homogeneous_bc ) const
{

   for ( typename hier::PatchLevel<DIM>::Iterator p(level); p; p++ ) {
      tbox::Pointer< hier::Patch<DIM> > patch = level.getPatch(p());
      setBoundaryValuesInCells( *patch,
                                fill_time,
                                ghost_width_to_fill ,
                                target_data_id,
                                homogeneous_bc );
   }
   return;
}





/*
************************************************************************
* Set physical boundary conditions at nodes.                           *
************************************************************************
*/

template<int DIM> void CartesianRobinBcHelper<DIM>::setBoundaryValuesAtNodes (
   hier::Patch<DIM> &patch ,
   const double fill_time ,
   int target_data_id ,
   bool homogeneous_bc ) const
{
   NULL_USE(patch);
   NULL_USE(fill_time);
   NULL_USE(target_data_id);
   NULL_USE(homogeneous_bc);

   TBOX_ERROR(d_object_name << ": Using incomplete implementation"
              << "CartesianRobinBcHelper<DIM>::setBoundaryValuesAtNodes"
              << "is not implemented because there is not a need for it (yet)"
              << std::endl);
   return;
}





/*
************************************************************************
* Set the coefficient strategy pointer that will be used to get        *
* Robin bc coefficients.  It should be some external implementation.   *
* This function implies that the simple mappings for                   *
* parallelpiped domains are not used and resets those arrays to null.  *
* is a parallelpiped (not checked) and that the boundary condition     *
* coefficients are functions only of the location index of the         *
* boundary.                                                            *
************************************************************************
*/

template<int DIM> void CartesianRobinBcHelper<DIM>::setCoefImplementation (
   const RobinBcCoefStrategy<DIM> *coef_strategy )
{
   if ( ! coef_strategy ) {
      TBOX_ERROR(d_object_name << ": Invalid pointer value"
                 << std::endl);
   }
   d_coef_strategy = coef_strategy;
   return;
}




template<int DIM> void CartesianRobinBcHelper<DIM>::setTargetDataId(
   int target_data_id )
{
   d_target_data_id = target_data_id;
   return;
}




template<int DIM> void CartesianRobinBcHelper<DIM>::setHomogeneousBc(
   bool is_homogeneous )
{
   d_homogeneous_bc = is_homogeneous;
}





/*
***********************************************************************
*                                                                     *
*  Virtual functions or xfer::RefinePatchStrategy<DIM>.                    *
*                                                                     *
***********************************************************************
*/

template<int DIM> void CartesianRobinBcHelper<DIM>::setPhysicalBoundaryConditions (
   hier::Patch<DIM> &patch ,
   const double fill_time ,
   const hier::IntVector<DIM> &ghost_width_to_fill ) {

   setBoundaryValuesInCells( patch ,
                             fill_time ,
                             ghost_width_to_fill ,
                             d_target_data_id ,
                             d_homogeneous_bc );
   return;
}


template<int DIM> hier::IntVector<DIM> CartesianRobinBcHelper<DIM>::getRefineOpStencilWidth() const
{
   return hier::IntVector<DIM>(0);
}

template<int DIM> void CartesianRobinBcHelper<DIM>::preprocessRefineBoxes (
      hier::Patch<DIM> &fine ,
      const hier::Patch<DIM> &coarse ,
      const hier::BoxList<DIM> &fine_boxes ,
      const hier::IntVector<DIM> &ratio )
{
   NULL_USE(fine);
   NULL_USE(coarse);
   NULL_USE(fine_boxes);
   NULL_USE(ratio);
   return;
}
template<int DIM> void CartesianRobinBcHelper<DIM>::preprocessRefine (
      hier::Patch<DIM> &fine ,
      const hier::Patch<DIM> &coarse ,
      const hier::Box<DIM> &fine_box ,
      const hier::IntVector<DIM> &ratio )
{
   NULL_USE(fine);
   NULL_USE(coarse);
   NULL_USE(fine_box);
   NULL_USE(ratio);
   return;
}
template<int DIM> void CartesianRobinBcHelper<DIM>::postprocessRefineBoxes (
      hier::Patch<DIM> &fine ,
      const hier::Patch<DIM> &coarse ,
      const hier::BoxList<DIM> &fine_box ,
      const hier::IntVector<DIM> &ratio )
{
   NULL_USE(fine);
   NULL_USE(coarse);
   NULL_USE(fine_box);
   NULL_USE(ratio);
  return;
}
template<int DIM> void CartesianRobinBcHelper<DIM>::postprocessRefine (
      hier::Patch<DIM> &fine ,
      const hier::Patch<DIM> &coarse ,
      const hier::Box<DIM> &fine_boxes ,
      const hier::IntVector<DIM> &ratio )
{
   NULL_USE(fine);
   NULL_USE(coarse);
   NULL_USE(fine_boxes);
   NULL_USE(ratio);
   return;
}




/*
************************************************************************
* Trim a boundary box so it does not stick out past the corners of a   *
* given box.  This removes the extension parallel to the boundary,     *
* past the corner of the limit box.                                    *
************************************************************************
*/

template<int DIM> hier::BoundaryBox<DIM> CartesianRobinBcHelper<DIM>::trimBoundaryBox(
   const hier::BoundaryBox<DIM> &boundary_box ,
   const hier::Box<DIM> &limit_box ) const
{
   if ( boundary_box.getBoundaryType() == DIM ) {
      // This is a node boundary box and cannot be trimmed anymore.
      return boundary_box;
   }
   const hier::Box<DIM> &bbox = boundary_box.getBox();
   const hier::Index<DIM> &plo = limit_box.lower();
   const hier::Index<DIM> &pup = limit_box.upper();
   const hier::Index<DIM> &blo = bbox.lower();
   const hier::Index<DIM> &bup = bbox.upper();
   hier::Index<DIM> newlo, newup;
   int key_direction;
   int d;
   switch ( boundary_box.getBoundaryType() ) {
   case 2:
      key_direction = boundary_box.getLocationIndex()/4;
      for ( d=0; d<DIM; ++d ) {
         if ( d == key_direction ) {
           newlo(d) = tbox::MathUtilities<int>::Max( blo(d) , plo(d) );
           newup(d) = tbox::MathUtilities<int>::Min( bup(d) , pup(d) );
         }
         else {
            newlo(d) = blo(d);
            newup(d) = bup(d);
         }
      }
      break;
   case 1:
      key_direction = boundary_box.getLocationIndex()/2;
      /*
       * Loop through directions.
       * Preserve box size in direction normal to boundary.
       * Trim box size in direction transverse to boundary.
       */
      for ( d=0; d<DIM; ++d ) {
         if ( d == key_direction ) {
            newlo(d) = blo(d);
            newup(d) = bup(d);
         }
         else {
            // Min side.  Use max between boundary and patch boxes.
            newlo(d) = tbox::MathUtilities<int>::Max( blo(d) , plo(d) );
            // Max side.  Use min between boundary and patch boxes.
            newup(d) = tbox::MathUtilities<int>::Min( bup(d) , pup(d) );
         }
      }
      break;
   }
   const hier::Box<DIM> newbox(newlo,newup);
   const hier::BoundaryBox<DIM> newbbox( newbox,
                                    boundary_box.getBoundaryType(),
                                    boundary_box.getLocationIndex() );
   return newbbox;
}



/*
************************************************************************
* Make surface box on boundary using standard boundary box             *
************************************************************************
*/

template<int DIM> hier::Box<DIM> CartesianRobinBcHelper<DIM>::makeFaceBoundaryBox(
   const hier::BoundaryBox<DIM> &boundary_box ) const
{
   if ( boundary_box.getBoundaryType() != 1 ) {
      TBOX_ERROR(d_object_name << ": makeFaceBoundaryBox called with\n"
                 << "improper boundary box\n"
                 << "for " << d_object_name );
   }
   hier::Box<DIM> face_indices = boundary_box.getBox();
   int location_index = boundary_box.getLocationIndex();
   if ( location_index%2 == 0 ) {
      /*
       * On the min index side, the face indices are one higher
       * than the boundary cell indices, in the direction normal
       * to the boundary.
       */
      face_indices.shift(location_index/2,1);
   }
   return face_indices;
}



/*
************************************************************************
* Make node box on boundary using standard boundary box                *
************************************************************************
*/

template<int DIM> hier::Box<DIM> CartesianRobinBcHelper<DIM>::makeNodeBoundaryBox(
   const hier::BoundaryBox<DIM> &boundary_box ) const
{
   if ( boundary_box.getBoundaryType() != 1 ) {
      TBOX_ERROR(d_object_name << ": makeNodeBoundaryBox called with\n"
                 << "improper boundary box\n"
                 << "for " << d_object_name );
   }
   hier::Box<DIM> node_indices = boundary_box.getBox();
   int location_index = boundary_box.getLocationIndex();
   if ( location_index%2 == 0 ) {
      /*
       * On the min index side, the node indices are one higher
       * than the boundary cell indices, in the direction normal
       * to the boundary.
       */
      node_indices.shift(location_index/2,1);
   }
   /*
    * The node indices range one higher than the cell indices,
    * in the directions parallel to the boundary.
    */
   hier::IntVector<DIM> parallel_growth(1);
   parallel_growth(location_index/2) = 0;
   node_indices.growUpper( parallel_growth );
   return node_indices;
}



}
}
#endif
