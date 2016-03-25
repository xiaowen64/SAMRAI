#ifndef included_solv_GhostCellRobinBcCoefs_C
#define included_solv_GhostCellRobinBcCoefs_C

/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/solvers/poisson/GhostCellRobinBcCoefs.C $
 * Package:     SAMRAI application utilities
 * Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 1988 $
 * Modified:    $LastChangedDate: 2008-02-14 10:04:49 -0800 (Thu, 14 Feb 2008) $
 * Description: Robin boundary condition support on cartesian grids.
 */


#include "VariableDatabase.h"
#include "Variable.h"
#include "CartesianPatchGeometry.h"
#include "ArrayDataBasicOps.h"
#include "CellData.h"
#include "CellVariable.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"
#include IOMANIP_HEADER_FILE

#include "GhostCellRobinBcCoefs.h"



namespace SAMRAI {
    namespace solv {



/*
************************************************************************
* Constructor                                                          *
************************************************************************
*/

template<int DIM>  GhostCellRobinBcCoefs<DIM>::GhostCellRobinBcCoefs(
   std::string object_name )
   : d_object_name(object_name),
     d_ghost_data_id(-1)
{

   t_set_bc_coefs = tbox::TimerManager::getManager()->
      getTimer("solv::GhostCellRobinBcCoefs::setBcCoefs()");

   return;
}



/*
************************************************************************
* Destructor                                                           *
************************************************************************
*/

template<int DIM>  GhostCellRobinBcCoefs<DIM>::~GhostCellRobinBcCoefs(void) {}





/*
************************************************************************
* Set the index of the data providing ghost cell values                *
************************************************************************
*/

template<int DIM> void GhostCellRobinBcCoefs<DIM>::setGhostDataId(
   int ghost_data_id,
   hier::IntVector<DIM> extensions_fillable)
{
   d_ghost_data_id = ghost_data_id;
   d_extensions_fillable = extensions_fillable;
   /*
    * Check for correctness of data index.
    * Unfortunately, the ghost width is not provided by the
    * variable database, so we cannot check that also.
    */
   if ( d_ghost_data_id != -1 ) {
      hier::VariableDatabase<DIM> *vdb = hier::VariableDatabase<DIM>::getDatabase();
      tbox::Pointer< hier::Variable<DIM> > variable_ptr;
      vdb->mapIndexToVariable( ghost_data_id, variable_ptr );
      if ( !variable_ptr ) {
         TBOX_ERROR(d_object_name << ": hier::Index " << ghost_data_id
                    << " does not correspond to any variable.");
      }
      tbox::Pointer<pdat::CellVariable<DIM,double> >
         cell_variable_ptr = variable_ptr;
      if ( !cell_variable_ptr ) {
         TBOX_ERROR(d_object_name << ": hier::Index " << ghost_data_id
                    << " does not correspond to a cell-double variable.");
      }
   }
   return;
}





/*
************************************************************************
* Set the bc coefficients reflect the value at the ghost cell centers. *
* The factor 1.0/(1+0.5*h) appears in a and g.  This factor comes      *
* from a linear approximation of the data through the patch boundary,  *
* going through the centers of the first interior and ghost cells      *
* and having the specified values there.                               *
************************************************************************
*/

template<int DIM> void GhostCellRobinBcCoefs<DIM>::setBcCoefs (
   tbox::Pointer<pdat::ArrayData<DIM,double> > &acoef_data ,
   tbox::Pointer<pdat::ArrayData<DIM,double> > &bcoef_data ,
   tbox::Pointer<pdat::ArrayData<DIM,double> > &gcoef_data ,
   const tbox::Pointer< hier::Variable<DIM> > &variable ,
   const hier::Patch<DIM> &patch ,
   const hier::BoundaryBox<DIM> &bdry_box ,
   double fill_time ) const
{
   NULL_USE(variable);
   NULL_USE(fill_time);

   t_set_bc_coefs->start();

   tbox::Pointer< geom::CartesianPatchGeometry<DIM> >
      patch_geom = patch.getPatchGeometry();
   const int norm_dir = bdry_box.getLocationIndex()/2;
   const double *dx = patch_geom->getDx();
   const double h = dx[ norm_dir ];

   /*
    * Set acoef_data to 1.0/(1+0.5*h) uniformly.  This value
    * corresponds to the fact that the solution is fixed at
    * the ghost cell centers.  bcoef_data is 1-acoef_data.
    */
   if ( !acoef_data.isNull() ) {
      acoef_data->fill( 1.0/(1+0.5*h) );
   }
   if ( !bcoef_data.isNull() ) {
      bcoef_data->fill( 0.5*h/(1+0.5*h) );
   }

   if ( !gcoef_data.isNull() ) {

     if ( d_ghost_data_id == -1 ) {
       TBOX_ERROR(d_object_name << ": Coefficient g requested without\n"
                  << "having valid ghost data id.\n");
     }

     /*
      * Fill in gcoef_data with data from d_ghost_data_id.
      * The data is first looked for in a pdat::OutersideData<DIM> object
      * and a pdat::CellData<DIM> object in that order.  Data from the
      * first place with allocated storage is used.
      */
     tbox::Pointer< hier::PatchData<DIM> > patch_data
       = patch.getPatchData( d_ghost_data_id );
     if ( patch_data.isNull() ) {
        TBOX_ERROR(d_object_name << ": hier::Patch data for index "
                   << d_ghost_data_id << " does not exist.");
     }
     tbox::Pointer<pdat::CellData<DIM,double> > cell_data = patch_data;
     if ( cell_data.isNull() ) {
        TBOX_ERROR(d_object_name << ": hier::Patch data for index "
                   << d_ghost_data_id << " is not cell double data.");
     }
     const int location_index = bdry_box.getLocationIndex();
     const hier::IntVector<DIM> &gw = cell_data->getGhostCellWidth();
     if ( gw[norm_dir] < 1 ) {
        TBOX_ERROR(d_object_name << ": hier::Patch data for index "
                   << d_ghost_data_id << " has zero ghost width.");
     }
     const pdat::ArrayData<DIM,double> &cell_array_data
       = cell_data->getArrayData();
     hier::IntVector<DIM> shift_amount(0);
     if ( location_index%2 == 0 ) shift_amount[location_index/2] = 1;
     gcoef_data->copy( cell_array_data,
                       makeSideBoundaryBox(bdry_box),
                       shift_amount );
     math::ArrayDataBasicOps<DIM,double> aops;
     /*
      * To convert from the value at the ghost cell center
      * to the coefficient g, we must scale the data by
      * 1/(1+h/2), according to our linear approximation
      * of the data at the patch boundary.
      */
     aops.scale( *gcoef_data, 1.0/(1+0.5*h), *gcoef_data,
                 makeSideBoundaryBox(bdry_box) );

   }

   t_set_bc_coefs->stop();

   return;
}



/*
***********************************************************************
* This class can only set coeficients for boundary boxes that extend  *
* no more than what the data it uses provides.                        *
***********************************************************************
*/
template<int DIM> hier::IntVector<DIM> GhostCellRobinBcCoefs<DIM>::numberOfExtensionsFillable()
   const
{
   return d_extensions_fillable;
}



/*
************************************************************************
* Make surface box on boundary using standard boundary box             *
************************************************************************
*/

template<int DIM> hier::Box<DIM> GhostCellRobinBcCoefs<DIM>::makeSideBoundaryBox(
   const hier::BoundaryBox<DIM> &boundary_box ) const
{
   if ( boundary_box.getBoundaryType() != 1 ) {
      TBOX_ERROR(d_object_name << ": CartesianRobinBcHelper<DIM>::makeSideBoundaryBox called with\n"
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



}
}
#endif
