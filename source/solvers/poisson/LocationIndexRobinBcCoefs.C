#ifndef included_solv_LocationIndexRobinBcCoefs_C
#define included_solv_LocationIndexRobinBcCoefs_C

/*
 * File:        LocationIndexRobinBcCoefs.C
 * Package:     SAMRAI application utilities
 * Copyright:   (c) 1997-2005 The Regents of the University of California
 * Revision:    $Revision: 724 $
 * Modified:    $Date: 2005-11-10 14:55:14 -0800 (Thu, 10 Nov 2005) $
 * Description: Robin boundary condition support on cartesian grids.
 */


#include "LocationIndexRobinBcCoefs.h"

#include "CartesianPatchGeometry.h"
#include "Index.h"
#include "tbox/Array.h"
#include "tbox/IEEE.h"
#include "tbox/Utilities.h"
#include IOMANIP_HEADER_FILE
#include <stdio.h>

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif


namespace SAMRAI {
    namespace solv {



/*
************************************************************************
* Constructor                                                          *
************************************************************************
*/

template<int DIM>  LocationIndexRobinBcCoefs<DIM>::LocationIndexRobinBcCoefs(
   const string &object_name,
   tbox::Pointer<tbox::Database> database
)
   : d_object_name(object_name)
{
   int i;
   for ( i=0; i<2*DIM; ++i ) {
      d_a_map[i] = d_g_map[i] = tbox::IEEE::getSignalingNaN();
   }
   if ( !database.isNull() ) {
      getFromInput(database);
   }
   return;
}



/*
************************************************************************
* Destructor                                                           *
************************************************************************
*/

template<int DIM>  LocationIndexRobinBcCoefs<DIM>::~LocationIndexRobinBcCoefs(void) {}




/*
********************************************************************
* Set state from input database                                    *
********************************************************************
*/

template<int DIM> void LocationIndexRobinBcCoefs<DIM>::getFromInput(
   tbox::Pointer<tbox::Database> database )
{
   if ( database ) {
      int i;
      char buf[20];
      for ( i=0; i<2*DIM; ++i ) {
	 sprintf( buf, "boundary_%d", i );
	 string name(buf);
	 if ( database->isString(name) ) {
	    d_a_map[i] = 1.0;
	    d_g_map[i] = 0.0;
	    tbox::Array<string> specs = database->getStringArray(name);
	    if ( specs[0] == "value" ) {
	       d_a_map[i] = 1.0;
	       if ( specs.size() > 1 ) d_g_map[i] = atof(specs[1].c_str());
	    }
	    else if ( specs[0] == "slope" ) {
	       d_a_map[i] = 0.0;
	       if ( specs.size() > 1 ) d_g_map[i] = atof(specs[1].c_str());
	    }
	    else if ( specs[0] == "coefficients" ) {
	       if ( specs.size() > 1 ) d_a_map[i] = atof(specs[1].c_str());
	       if ( specs.size() > 2 ) d_g_map[i] = atof(specs[2].c_str());
	    }
	    else {
	       TBOX_ERROR(d_object_name << ": Bad boundary specifier\n"
			  << "'" << specs[0] << "'.  Use either 'value'\n"
			  << "'slope' or 'coefficients'.\n");
	    }
	 }
      }
   }
   return;
}





/*
************************************************************************
* Set the boundary value for a Dirichlet boundary condition.           *
************************************************************************
*/

template<int DIM> void LocationIndexRobinBcCoefs<DIM>::setBoundaryValue (
   int location_index,
   double value)
{
   if ( location_index < 0 || location_index >= 2*DIM ) {
      TBOX_ERROR("Location index in " << DIM << "D must be\n"
                 <<"in [0," << 2*DIM-1 << "].\n");
   }
   d_a_map[location_index] = 1.0;
   d_g_map[location_index] = value;
   return;
}





/*
************************************************************************
* Set the slpe for a Neumann boundary condition.                       *
************************************************************************
*/

template<int DIM> void LocationIndexRobinBcCoefs<DIM>::setBoundarySlope (
   int location_index,
   double slope)
{
   if ( location_index >= 2*DIM ) {
      TBOX_ERROR("Location index in " << DIM << "D must be\n"
                 <<"in [0," << 2*DIM-1 << "].\n");
   }
   d_a_map[location_index] = 0.0;
   d_g_map[location_index] = slope;
   return;
}





/*
************************************************************************
* Set the raw bc coefficients.                                         *
************************************************************************
*/

template<int DIM> void LocationIndexRobinBcCoefs<DIM>::setRawCoefficients (
   int location_index,
   double a,
   double g)
{
   if ( location_index >= 2*DIM ) {
      TBOX_ERROR("Location index in " << DIM << "D must be\n"
                 <<"in [0," << 2*DIM-1 << "].\n");
   }
   if ( a < 0.0 || a > 1.0 ) {
      TBOX_ERROR("a coefficient must be between 0 and 1\n");
   }
   d_a_map[location_index] = a;
   d_g_map[location_index] = g;
   return;
}





/*
************************************************************************
* Set the bc coefficients to their mapped values.                      *
************************************************************************
*/

template<int DIM> void LocationIndexRobinBcCoefs<DIM>::setBcCoefs (
   tbox::Pointer<pdat::ArrayData<DIM,double> > &acoef_data ,
   tbox::Pointer<pdat::ArrayData<DIM,double> > &gcoef_data ,
   const tbox::Pointer< hier::Variable<DIM> > &variable ,
   const hier::Patch<DIM> &patch ,
   const hier::BoundaryBox<DIM> &bdry_box ,
   double fill_time ) const
{
   int location = bdry_box.getLocationIndex();
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( location >= 0 && location < 2*DIM );
#endif
   if ( acoef_data ) {
      acoef_data->fill( d_a_map[location] );
   }
   if ( gcoef_data ) {
      gcoef_data->fill( d_g_map[location] );
   }
   return;
}




template<int DIM> hier::IntVector<DIM> LocationIndexRobinBcCoefs<DIM>::numberOfExtensionsFillable()
   const
{
   /*
    * Return some really big number.  We have no limits.
    */
   return hier::IntVector<DIM>(1<<(sizeof(int)-1));
}




template<int DIM> void LocationIndexRobinBcCoefs<DIM>::getCoefficients( int i,
                                                                        double &a,
                                                                        double &g)
   const
{
   a = d_a_map[i];
   g = d_g_map[i];
   return;
}



}
}
#endif
