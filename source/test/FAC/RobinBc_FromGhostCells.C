/*
  File:		$RCSfile$
  Copyright:	(c) 1997-2005 The Regents of the University of California
  Revision:	$Revision: 173 $
  Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
  Description:	RobinBc_FromGhostCells class implementation
*/

#include "SAMRAI_config.h"

#include "printObject.h"

#include "RobinBc_FromGhostCells.h"

#include "tbox/Pointer.h"
#include "tbox/ConstPointer.h"
#include "tbox/Array.h"
#include "tbox/Utilities.h"

#include "CartesianPatchGeometry.h"
#include "Index.h"
#include "Patch.h"
#include "VariableDatabase.h"
#include "ArrayDataBasicOps.h"
#include "CellData.h"
#include "CellVariable.h"
#include "OutersideData.h"
#include "OutersideVariable.h"

using namespace SAMRAI;

/*
******************************************************************
Constructor.
******************************************************************
*/
template<int NDIM>
RobinBc_FromGhostCells<NDIM>::RobinBc_FromGhostCells(
  const string &object_name ,
  ostream *out_stream ,
  ostream *log_stream
)
: solv::CartesianRobinBc<NDIM>() ,
  d_ostream(out_stream) ,
  d_lstream(log_stream)
{
  return;
}

/*
******************************************************************
Set descriptor for ghost data to use.
******************************************************************
*/
template<int NDIM>
void RobinBc_FromGhostCells<NDIM>::setGhostDataId( int ghost_data_id ) {
   /*
     Check that the id is appropriate.  Do not check that the data
     is allocated because it does not have to be yet.
   */
   hier::VariableDatabase<NDIM> *vdb = hier::VariableDatabase<NDIM>::getDatabase();
   tbox::Pointer< hier::Variable<NDIM> > variable;
   vdb->mapIndexToVariable( ghost_data_id, variable );
   if ( variable.isNull() ) {
      TBOX_ERROR("Index " << ghost_data_id << " does not correspond\n"
		 << "to any database variable in " << d_object_name << ".\n");
   }
   tbox::Pointer<pdat::CellVariable<NDIM,double> > cell_variable = variable;
   tbox::Pointer<pdat::OutersideVariable<NDIM,double> > outerside_variable = variable;
   if ( cell_variable || outerside_variable ) {
      d_ghost_data_id = ghost_data_id;
   }
   else {
      TBOX_ERROR("Index " << ghost_data_id << " does not correspond\n"
		 << "to an appropriate database variable in\n"
		 << d_object_name << ".\n");
   }
   return;
}

/*
******************************************************************
Set the boundary condition coefficients on sides.
******************************************************************
*/
template<int NDIM>
void RobinBc_FromGhostCells<NDIM>::setBcCoefOnBoundaryBoxAtSides (
  tbox::Pointer<pdat::ArrayData<NDIM,double> > &acoef_data ,
  tbox::Pointer<pdat::ArrayData<NDIM,double> > &gcoef_data ,
  const hier::Patch<NDIM> &patch ,
  const hier::BoundaryBox<NDIM> &bdry_box ,
  double fill_time ) const
{
   tbox::Pointer<geom::CartesianPatchGeometry<NDIM> >
      patch_geom = patch.getPatchGeometry();
   const int norm_dir = bdry_box.getLocationIndex()/2;
   const double *dx = patch_geom->getDx();
   const double h = dx[ norm_dir ];

   /*
     Set acoef_data to 1.0/(1+0.5*h) uniformly.
   */
   setBcCoefOnBoundaryBoxAtSides_Uniform(
      acoef_data,
      gcoef_data,
      patch ,
      bdry_box,
      1.0/(1+0.5*h) );

   if ( !gcoef_data ) return;

   /*
     Fill in gcoef_data with data from d_ghost_data_id.
     The data is first looked for in a pdat::OutersideData<NDIM> object
     and a pdat::CellData<NDIM> object in that order.  Data from the
     first place with allocated storage is used.
   */
   bool finished = false;
   tbox::Pointer< hier::PatchData<NDIM> > patch_data
      = patch.getPatchData( d_ghost_data_id );
   tbox::Pointer<pdat::CellData<NDIM,double> > cell_data = patch_data;
   tbox::Pointer<pdat::OutersideData<NDIM,double> > outerside_data = patch_data;
   const int location_index = bdry_box.getLocationIndex();
   if ( ! outerside_data.isNull() ) {
      const pdat::ArrayData<NDIM,double> &outerside_array_data
	 = outerside_data->getArrayData(location_index/2,location_index%2);
      gcoef_data->copy( outerside_array_data,
			makeSideBoundaryBox(bdry_box) );
      finished = true;
   }
   else if ( ! cell_data.isNull() ) {
      const hier::IntVector<NDIM> &gw = cell_data->getGhostCellWidth();
      if ( gw[norm_dir] > 0 ) {
	 /* Ghost cells exist. */
	 const pdat::ArrayData<NDIM,double> &cell_array_data
	    = cell_data->getArrayData();
	 hier::IntVector<NDIM> shift_amount(0);
	 if ( location_index%2 == 0 ) shift_amount[location_index/2] = 1;
	 gcoef_data->copy( cell_array_data,
			   makeSideBoundaryBox(bdry_box),
			   shift_amount );
	 math::ArrayDataBasicOps<NDIM,double> aops;
	 aops.scale( *gcoef_data, 1.0/(1+0.5*h), *gcoef_data,
		     makeSideBoundaryBox(bdry_box) );
	 tbox::plog << "Cell data:\n";
	 cell_data->print( cell_data->getGhostBox(), 0 );
	 tbox::plog << "gcoef data:\n";
	 printObject<NDIM>(plog,*gcoef_data);
	 finished = true;
      }
   }
   else {
      TBOX_ERROR("No variable for index " << d_ghost_data_id
		 << " found in " << d_object_name);
   }

   if ( finished == false ) {
      TBOX_ERROR("Unable to set boundary coefficient data for "
		 << d_object_name);
   }

   return;
}
