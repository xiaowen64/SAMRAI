//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/geometry/cartesian/grid_geom/CartesianGridGeometry.h $
// Package:	SAMRAI geometry package
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Simple Cartesian grid geometry for an AMR hierarchy.
//

#ifndef included_geom_CartesianGridGeometry
#define included_geom_CartesianGridGeometry

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_hier_BoundaryBox
#include "BoundaryBox.h"
#endif
#ifndef included_hier_BoxArray
#include "BoxArray.h"
#endif
#ifndef included_hier_IntVector
#include "IntVector.h"
#endif
#ifndef included_hier_Patch
#include "Patch.h"
#endif
#ifndef included_hier_PatchLevel
#include "PatchLevel.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif
#ifndef included_tbox_Database
#include "tbox/Database.h"
#endif
#ifndef included_tbox_Serializable
#include "tbox/Serializable.h"
#endif
#ifndef included_String
#include <string>
#define included_String
#endif
#ifndef included_xfer_Geometry
#include "Geometry.h"
#endif

namespace SAMRAI {
    namespace geom {

/**
 * Class CartesianGridGeometry<DIM> provides simple Cartesian mesh geometry
 * management on an AMR hierarchy.  The mesh coordinates on each hierarchy 
 * level are limited to mesh increments specified as DIM-tuple 
 * (dx[0],...,dx[DIM-1]) and spatial coordinates of the lower and upper
 * corners of the smallest parallelepiped bounding the entire computational
 * domain.  The mesh increments on each level are defined with respect to 
 * the coarsest hierarchy level and multiplying those values by the proper
 * refinement ratio.  This class sets geometry information on each patch in
 * an AMR hierarchy.  This class is derived from the xfer::Geometry<DIM> base 
 * class which is further derived from the hier::GridGeometry<DIM> base class.
 *
 * An object of this class requires numerous parameters to be read from
 * input.  Also, data must be written to and read from files for restart.
 * The input and restart data are summarized as follows:
 *
 * Required input keys and data types: 
 * 


 *    
 *    - \b    domain_boxes   
 *       tbox::Array of boxes representing the index space for the entire
 *       domain (on the coarsest refinement level).
 *
 *    - \b    x_lo   
 *       tbox::Array of double values representing the spatial coordinates of
 *       the lower corner of the physical domain.
 *
 *    - \b    x_up   
 *       tbox::Array of double values representing the spatial coordinates of
 *       the upper corner of the physical domain.
 *
 * 


 * 
 * Optional input keys, data types, and defaults:
 * 


 *
 *    - \b    periodic_dimension   
 *       tbox::Array of integer values representing the directions in which
 *       the physical domain is periodic.  A non-zero value indicates
 *       that the direction is periodic.  A zero value indicates that  
 *       the direction is not periodic.  If no values are specified, then
 *       the array is initialized to all zeros (no periodic directions).
 * 
 *    - \b    use_original_location_indices 
 *       Boolean argument to handle backward compatibility with new
 *       location index scheme for codimension 2 in 3 dimensions.
 *       Set to true to use the location index scheme that was used in
 *       SAMRAI v. 1.4 and earlier versions, and false to use the new
 *       scheme.  In 3 dimensions, this key defaults to true if not
 *       present in the input.  In all other dimensions, this key
 *       is irrelevant.
 *

 * 
 * No input values can overwrite restart values.
 *
 * A sample input file for a two-dimensional problem might look like:
 *
 * @verbatim
 *
 *    domain_boxes = [(0,0) , (49,39)]
 *    x_lo = 0.0 , 0.0
 *    x_up = 50.0 , 40.0
 *    periodic_dimension = 0, 1  // periodic in y only
 *
 * @endverbatim
 *
 * This generates a two-dimensional rectangular domain periodic in the
 * y-direction, and having 50 cells in the x-direction and 40 cells in
 * the y-direction, with the cell size 1 unit in each direction.
 *
 * @see xfer::Geometry
 * @see hier::GridGeometry
 */

template<int DIM> class CartesianGridGeometry 
: public xfer::Geometry<DIM>,
  public tbox::Serializable
{
public:
   /**
    * Constructor for CartesianGridGeometry<DIM> initializes data
    * members based on parameters read from the specified input database
    * or from the restart database corresponding to the specified 
    * object name.  The constructor also registers this object 
    * for restart using the specified object name when the boolean
    * argument is true.  Whether object will write its state to restart 
    * files during program execution is determined by this argument.  
    * Note that it has a default state of true.
    *
    * Errors: passing in a null database pointer or an empty string
    * will result in an unrecoverable assertion.
    */
   CartesianGridGeometry(const std::string& object_name,
                               tbox::Pointer<tbox::Database> input_db,
                               bool register_for_restart = true);

   /**
    * Constructor for CartesianGridGeometry<DIM> sets data members 
    * based on arguments.  The constructor also registers this object
    * for restart using the specified object name when the boolean
    * argument is true.  Whether object will write its state to restart 
    * files during program execution is determined by this argument.  
    * Note that it has a default state of true. 
    * 
    * Errors: passing in an empty string, or null data pointers will
    * result in an unrecoverable assertion.
    */
   CartesianGridGeometry(const std::string& object_name,
                               const double* x_lo,
                               const double* x_up,
                               const hier::BoxArray<DIM>& domain,
                               bool register_for_restart = true);

   /**
    * Destructor for CartesianGridGeometry<DIM> deallocates 
    * data describing grid geometry and unregisters the object with
    * the restart manager if previously registered.
    */
   virtual ~CartesianGridGeometry<DIM>();

   
   /**
    * Create and return a pointer to a refined version of this Cartesian grid 
    * geometry object. This function is pure virtual in the hier::GridGeometry<DIM> base class.
    */
   tbox::Pointer<hier::GridGeometry<DIM> > makeRefinedGridGeometry(
      const std::string& fine_geom_name,
      const hier::IntVector<DIM>& refine_ratio,
      bool register_for_restart) const;

   /**
    * Create and return a pointer to a coarsened version of this Cartesian grid
    * geometry object. This function is pure virtual in the hier::GridGeometry<DIM> base class.
    */
   tbox::Pointer<hier::GridGeometry<DIM> > makeCoarsenedGridGeometry(
      const std::string& coarse_geom_name,
      const hier::IntVector<DIM>& coarsen_ratio,
      bool register_for_restart) const;

   /*
    * Compute grid data for patch and assign new geom_CartesianPatchGeometry<DIM>
    * object to patch.  This function is pure virtual in the hier::GridGeometry<DIM> 
    * base class.
    */
   void setGeometryDataOnPatch(
      hier::Patch<DIM>& patch,
      const hier::IntVector<DIM>& ratio_to_level_zero,
      const tbox::Array< tbox::Array<bool> >& touches_regular_bdry,
      const tbox::Array< tbox::Array<bool> >& touches_periodic_bdry) const;

   /**
    * Set data members for this CartesianGridGeometry<DIM> object.
    */
   void setGeometryData(const double* x_lo,
                        const double* x_up,
                        const hier::BoxArray<DIM>& domain);

   /**
    * Return const pointer to dx array for reference level in hierarchy.
    */
   const double* getDx() const; 

   /**
    * Return const pointer to lower spatial coordinate for reference 
    * level in hierarchy.
    */
   const double* getXLower() const; 

   /**
    * Return const pointer to upper spatial coordinate for reference 
    * level in hierarchy.
    */
   const double* getXUpper() const; 

   /**
    * Print class data representation.
    */
   virtual void printClassData(std::ostream& os) const;

   /**
    * Writes the state of the CartesianGridGeometry object to the
    * database.  
    * 
    * When assertion checking is active, db cannot be a null database pointer.
    */
   virtual void putToDatabase(tbox::Pointer<tbox::Database> db);

private:

   /*
    * Reads in domain_boxes, x_lo, and x_up from the input database.
    * Data is read from input only if the simulation is not from restart.
    * Otherwise, all values specified in the input database are ignored.
    *
    * Arguments: is_from_restart is true when simulation is from restart
    * Assertions: db must not be a null pointer.
    */
   void getFromInput(tbox::Pointer<tbox::Database> db, 
                     bool is_from_restart);

   /*
    * Read object state from the restart file and initialize class data
    * members.  The database from which the restart data is read is
    * determined by the object_name specified in the constructor.
    *
    * Unrecoverable Errors:
    *
    *    -The database corresponding to object_name is not found
    *     in the restart file.
    *
    *    -The class version number and restart version number do not
    *     match.
    *
    */
   void getFromRestart();

   /*
    * Create default operator for Cartesian grid geometry.
    */
   void makeStandardOperators();

   /*
    * String name for object used in error reporting and restart operations.
    * Boolean is set in constructor and determines whether object should
    * dump its state to restart files during program execution.
    */
   std::string d_object_name;
   bool d_registered_for_restart;

   double d_dx[DIM];                // mesh increments for level 0.
   double d_x_lo[DIM];              // spatial coordinates of lower corner
                                     // (i.e., box corner) of problem domain.
   double d_x_up[DIM];              // spatial coordinates of upper corner
                                     // (i.e., box corner) of problem domain.

   hier::Box<DIM> d_domain_box;           // smallest box covering coarsest level
                                     // (i.e., reference level) index space.

   bool d_using_original_locations;
 
};

}
}
#ifndef DEBUG_NO_INLINE
#include "CartesianGridGeometry.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "CartesianGridGeometry.C"
#endif
