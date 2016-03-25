#ifndef included_solv_LocationIndexRobinBcCoefs_C
#define included_solv_LocationIndexRobinBcCoefs_C

/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/solvers/poisson/LocationIndexRobinBcCoefs.C $
 * Package:     SAMRAI application utilities
 * Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 2172 $
 * Modified:    $LastChangedDate: 2008-05-02 11:02:08 -0700 (Fri, 02 May 2008) $
 * Description: Robin boundary condition support on cartesian grids.
 */

#include <stdlib.h>

#include "LocationIndexRobinBcCoefs.h"

#include "CartesianPatchGeometry.h"
#include "Index.h"
#include "tbox/Array.h"
#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"
#include IOMANIP_HEADER_FILE

namespace SAMRAI {
    namespace solv {



/*
************************************************************************
* Default constructor                                                  *
************************************************************************
*/

template<int DIM>  LocationIndexRobinBcCoefs<DIM>::LocationIndexRobinBcCoefs(
)
   : d_object_name("")
{
   int i;
   for ( i=0; i<2*DIM; ++i ) {
      d_a_map[i] = tbox::MathUtilities<double>::getSignalingNaN();
      d_b_map[i] = tbox::MathUtilities<double>::getSignalingNaN();
      d_g_map[i] = tbox::MathUtilities<double>::getSignalingNaN();
   }
   return;
}



/*
************************************************************************
* Constructor using database 
************************************************************************
*/

template<int DIM>  LocationIndexRobinBcCoefs<DIM>::LocationIndexRobinBcCoefs(
   const std::string &object_name,
   tbox::Pointer<tbox::Database> database
)
   : d_object_name(object_name)
{
   int i;
   for ( i=0; i<2*DIM; ++i ) {
      d_a_map[i] = tbox::MathUtilities<double>::getSignalingNaN();
      d_b_map[i] = tbox::MathUtilities<double>::getSignalingNaN();
      d_g_map[i] = tbox::MathUtilities<double>::getSignalingNaN();
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
      for ( i=0; i<2*DIM; ++i ) {
	 std::string name = "boundary_" + tbox::Utilities::intToString(i);
	 if ( database->isString(name) ) {
	    d_a_map[i] = 1.0;
	    d_g_map[i] = 0.0;
	    tbox::Array<std::string> specs = database->getStringArray(name);
	    if ( specs[0] == "value" ) {
	       d_a_map[i] = 1.0;
	       d_b_map[i] = 0.0;
	       if ( specs.size() > 1 ) d_g_map[i] = atof(specs[1].c_str());
	    }
	    else if ( specs[0] == "slope" ) {
	       d_a_map[i] = 0.0;
	       d_b_map[i] = 1.0;
	       if ( specs.size() > 1 ) d_g_map[i] = atof(specs[1].c_str());
	    }
	    else if ( specs[0] == "coefficients" ) {
	       if ( specs.size() > 1 ) d_a_map[i] = atof(specs[1].c_str());
	       if ( specs.size() > 2 ) d_b_map[i] = atof(specs[2].c_str());
	       if ( specs.size() > 3 ) d_g_map[i] = atof(specs[3].c_str());
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
   d_b_map[location_index] = 0.0;
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
   d_b_map[location_index] = 1.0;
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
   double b,
   double g)
{
   if ( location_index >= 2*DIM ) {
      TBOX_ERROR("Location index in " << DIM << "D must be\n"
                 <<"in [0," << 2*DIM-1 << "].\n");
   }
   d_a_map[location_index] = a;
   d_b_map[location_index] = b;
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
   tbox::Pointer<pdat::ArrayData<DIM,double> > &bcoef_data ,
   tbox::Pointer<pdat::ArrayData<DIM,double> > &gcoef_data ,
   const tbox::Pointer< hier::Variable<DIM> > &variable ,
   const hier::Patch<DIM> &patch ,
   const hier::BoundaryBox<DIM> &bdry_box ,
   double fill_time ) const
{
   NULL_USE(variable);
   NULL_USE(patch);
   NULL_USE(fill_time);

   int location = bdry_box.getLocationIndex();
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( location >= 0 && location < 2*DIM );
#endif
   if ( acoef_data ) {
      acoef_data->fill( d_a_map[location] );
   }
   if ( bcoef_data ) {
      bcoef_data->fill( d_b_map[location] );
   }
   if ( gcoef_data ) {
      gcoef_data->fill( d_g_map[location] );
   }
   return;
}




template<int DIM> hier::IntVector<DIM>
LocationIndexRobinBcCoefs<DIM>::numberOfExtensionsFillable() const
{
   /*
    * Return some really big number.  We have no limits.
    */
   return hier::IntVector<DIM>(1<<(sizeof(int)-1));
}




template<int DIM> void LocationIndexRobinBcCoefs<DIM>::getCoefficients(
   int i,
   double &a,
   double &b,
   double &g) const
{
   a = d_a_map[i];
   b = d_b_map[i];
   g = d_g_map[i];
   return;
}

/*
************************************************************************
* Assignment operator                                                  *
************************************************************************
*/

template<int DIM>
const LocationIndexRobinBcCoefs<DIM> &LocationIndexRobinBcCoefs<DIM>::operator=(
   const LocationIndexRobinBcCoefs<DIM> &r )
{
   d_object_name = r.d_object_name;
   for ( size_t i=0; i<2*DIM; ++i ) {
      d_a_map[i] = r.d_a_map[i];
      d_b_map[i] = r.d_b_map[i];
      d_g_map[i] = r.d_g_map[i];
   }
   return *this;
}




}
}
#endif
