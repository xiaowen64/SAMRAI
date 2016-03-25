/*
 * File:        LocationIndexRobinBcCoefs.h
 * Package:     SAMRAI solver package
 * Copyright:   (c) 1997-2005 The Regents of the University of California
 * Revision:    $Revision: 724 $
 * Modified:    $Date: 2005-11-10 14:55:14 -0800 (Thu, 10 Nov 2005) $
 * Description: Robin boundary condition problem-dependent interfaces
 */

#ifndef included_solv_LocationIndexRobinBcCoefs
#define included_solv_LocationIndexRobinBcCoefs

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

/*
 * SAMRAI classes
 */

#ifndef included_hier_BoundaryBox
#include "BoundaryBox.h"
#endif

#ifndef included_hier_Patch
#include "Patch.h"
#endif

#ifndef included_tbox_ArrayData
#include "ArrayData.h"
#endif

#ifndef included_solv_RobinBcCoefStrategy
#include "RobinBcCoefStrategy.h"
#endif

#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif


namespace SAMRAI {
    namespace solv {


/*!
 * @brief A prefabricated Robin boundary condition coefficients
 * for coefficients that are entirely specified by the boundary
 * box location index.
 *
 * This implementation of the strategy
 * class RobinBcCoefStrategy<DIM> may be used when your
 * Robin boundary condition coefficients are completely determined
 * by the location index of the boundary box.
 *
 * Before this class is used in to provide the boundary condition
 * coefficients, you must specify what boundary conditions to
 * associate with what location index.  Methods for specifying
 * these are setBoundaryValue(), setBoundarySlope() and
 * setRawCoefficients().  The first two are for Dirichlet
 * and Neumann boundary conditions, respectively.  If the boundary
 * condition is the more general Robin boundary condition,
 * the third function should be used to set the coefficients
 * a and g directly (see RobinBcCoefStrategy) for the
 * meanings of a and g.
 *
 * @b Inputs:
 * You can specify the boundary conditions for any location index
 * through the input database.  One line is required for each
 * location index.  The input parameters are "boundary_N", where
 * N is the index of the location.  Each parameter must be
 * a string array so that all boundary types can be accomodated
 * the same way.  The first string must be one of "value",
 * "slope" or "coefficients".  If the string is "value" or "slope"
 * the next string is the value you want to set, defaulting to
 * zero if not specified.  If the first string is "coefficients",
 * the next two strings specifies the values of a and g.
 *
 * @b Examples inputs:
 * @verbatim
 * boundary_0 = "value", "0.0"
 * boundary_1 = "value", "1.0"
 * boundary_2 = "slope", "0.0"
 * boundary_4 = "coefficients", "1.0", "0.0"
 * @endverbatim
 */
template<int DIM> class LocationIndexRobinBcCoefs
  : public RobinBcCoefStrategy<DIM>
{

public:

   /*!
    * @brief Constructor
    */
   LocationIndexRobinBcCoefs( const string &object_name,
				    tbox::Pointer<tbox::Database> database );


   /*!
    * @brief Destructor.
    */
   virtual ~LocationIndexRobinBcCoefs<DIM>(void);


   /*!
    * @brief Function to fill arrays of Robin boundary
    * condition coefficients at a patch boundary.
    *
    * This implementation of the virtual function
    * RobinBcCoefStrategy<DIM>::setBcCoefs()
    * fills the coefficient arrays with constant values
    * set according to the location index of the boundary box.
    *
    * @param acoef_data boundary coefficient data.
    *        This is defined to include index range for
    *        the boundary faces on the boundary box @c bdry_box.
    * @param gcoef_data boundary coefficient data.
    *        This is defined to include index range for
    *        the boundary faces on the boundary box @c bdry_box.
    *        If this is a null pointer, then the calling function
    *        is not interested in g, and you can disregard it.
    * @param variable variable to set the coefficients for.
    * @param patch patch requiring bc coefficients
    * @param bdry_box boundary box showing where on the boundary
    *        the coefficient data is needed.
    * @param fill_time Solution time corresponding to filling,
    *        for use when coefficients are time-dependent.
    */
   void setBcCoefs (
      tbox::Pointer<pdat::ArrayData<DIM,double> > &acoef_data ,
      tbox::Pointer<pdat::ArrayData<DIM,double> > &gcoef_data ,
      const tbox::Pointer< hier::Variable<DIM> > &variable ,
      const hier::Patch<DIM> &patch ,
      const hier::BoundaryBox<DIM> &bdry_box ,
      double fill_time=0.0 ) const;

   /*
    * @brief Return how many cells past the edge or corner of the
    * patch the object can fill.
    */
   hier::IntVector<DIM> numberOfExtensionsFillable() const;


   /*!
    * @brief Set the boundary value at a given location index.
    *
    * @param location_index Set coefficients for this index.
    * @param value Boundary value at @c location_index.
    */
   void setBoundaryValue(int location_index,
                         double value);


   /*!
    * @brief Set the boundary slope at a given location index.
    *
    * @param location_index Set coefficients for this index.
    * @param slope Boundary slope at @c location_index.
    */
   void setBoundarySlope(int location_index,
                         double slope);


   /*!
    * @brief Set the values of coefficients a and g at a
    * given location index.
    *
    * See RobinBcCoefStrategy<DIM> for the definitions
    * of coefficients a and g.
    *
    * If the boundary condition is neither Dirichlet nor
    * Neumann (a general Robin boundary condition), use
    * this function to set the values of the bc coefficients.
    *
    * @param location_index Set coefficients for this index.
    * @param a Value of coefficient a at given location index.
    * @param g Value of coefficient g at given location index.
    */
   void setRawCoefficients(int location_index,
                           double a,
                           double g);


   /*!
    * @brief Access coefficients.
    */
   void getCoefficients( int location_index, double &a, double &g ) const;


private:


   /*
    * @brief Set state from input database.
    */
   void getFromInput( tbox::Pointer<tbox::Database> database );


   /*
    * @brief Object name.
    */
   string d_object_name;

   /*
    * @brief Mapping for a coefficient.
    */
   double d_a_map[2*DIM];
   /*
    * @brief Mapping for g coefficient.
    */
   double d_g_map[2*DIM];

};

}
}

#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "LocationIndexRobinBcCoefs.C"
#endif
