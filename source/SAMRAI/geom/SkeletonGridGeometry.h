/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Skeleton grid geometry for an AMR hierarchy. 
 *
 ************************************************************************/

#ifndef included_geom_SkeletonGridGeometry
#define included_geom_SkeletonGridGeometry

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/tbox/Serializable.h"
#include "SAMRAI/hier/GridGeometry.h"

namespace SAMRAI {
namespace geom {

/**
 * Class SkeletonGridGeometry is a concrete grid geometry class that
 * contains no information about the physical domain characteristics of an AMR
 * mesh apart from the index space.  The purpose of this class is to allow
 * an application that needs to use a special mesh to manage the physical
 * domain within the application code.  The skeleton grid geometry only
 * manages the index space of the adaptive grid.  This class sets geometry
 * information on each patch in an AMR hierarchy.  This class is derived
 * from the hier::GridGeometry base class.
 *
 * An object of this class requires parameters to be read from input to create
 * the hier::BoxList that stores the index space in the hier::GridGeometry
 * superclass.  Also, data must be written to and read from files for restart.
 * The input and restart data are summarized as follows:
 *
 *
 * Required input keys and data types:
 *
 *    - \b    domain_boxes
 *       Array of boxes representing the index space for the entire
 *       domain (on the coarsest refinement level).
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
 * No input values can overwrite restart values.
 *
 * A sample input file for a two-dimensional problem might look like:
 *
 * @verbatim
 *
 *    domain_boxes = [(0,0) , (49,39)]
 *    periodic_dimension = 0, 1  // periodic in y only
 *
 * @endverbatim
 *
 * This generates a two-dimensional domain periodic in the
 * y-direction, and having 50 cells in the x-direction and 40 cells in
 * the y-direction.
 *
 * @see hier::GridGeometry
 */

class SkeletonGridGeometry:
   public hier::GridGeometry
{

   typedef hier::PatchGeometry::TwoDimBool TwoDimBool;

public:
   /**
    * Constructor for SkeletonGridGeometry initializes data
    * members based on parameters read from the specified input and
    * restart databases.  The constructor also registers this object
    * for restart using the specified object name, when the boolean
    * argument is true.  Whether object will write its state to restart
    * files during program execution is determined by this argument.
    * Note that it has a default state of true.
    *
    * Errors: passing in a null database pointer or an empty std::string
    * will result in an unrecoverable assertion.
    */
   explicit SkeletonGridGeometry(
      const tbox::Dimension& dim,
      const std::string& object_name,
      tbox::Pointer<tbox::Database> input_db,
      bool register_for_restart = true);

   /**
    * Constructor for SkeletonGridGeometry sets index space domain
    * based on arguments. The constructor also registers this object
    * for restart using the specified object name, when the boolean
    * argument is true.  Whether object will write its state to restart
    * files during program execution is determined by this argument.
    * Note that it has a default state of true.
    *
    * Errors: passing in an empty std::string, or null data pointers will
    * result in an unrecoverable assertion.
    */
   explicit SkeletonGridGeometry(
      const std::string& object_name,
      const hier::BoxList& level_domain,
      bool register_for_restart = true);

   /**
    * Destructor for SkeletonGridGeometry unregisters the object
    * with the restart manager if previously registered.
    */
   virtual ~SkeletonGridGeometry();

   /**
    * Create and return a pointer to a refined version of this Cartesian grid
    * geometry object. This function is pure virtual in the hier_GridGeometry base class.
    */
   tbox::Pointer<hier::GridGeometry>
   makeRefinedGridGeometry(
      const std::string& fine_geom_name,
      const hier::IntVector& refine_ratio,
      bool register_for_restart) const;

   /**
    * Create and return a pointer to a coarsened version of this Cartesian grid
    * geometry object. This function is pure virtual in the hier_GridGeometry base class.
    */
   tbox::Pointer<hier::GridGeometry>
   makeCoarsenedGridGeometry(
      const std::string& coarse_geom_name,
      const hier::IntVector& coarsen_ratio,
      bool register_for_restart) const;

   /*
    * Compute grid data for patch and assign new geom_SkeletonPatchGeometry
    * object to patch.  This function is pure virtual in the hier_GridGeometry
    * base class.
    */
   void
   setGeometryDataOnPatch(
      hier::Patch& patch,
      const hier::IntVector& ratio_to_level_zero,
      const TwoDimBool& touches_regular_bdry,
      const TwoDimBool& touches_periodic_bdry) const;

   /**
    * Print class data representation.
    */
   virtual void
   printClassData(
      std::ostream& os) const;

   /**
    * Writes the state of the SkeletonGridGeometry object to the database.
    *
    * When assertion checking is active, db cannot be a null database pointer.
    */
   virtual void
   putToDatabase(
      tbox::Pointer<tbox::Database> db);

private:
   /*
    * Static integer constant describing class's version number.
    */
   static const int GEOM_SKELETON_GRID_GEOMETRY_VERSION;

   /*
    * Reads in d_physical_domain (from hier::GridGeometry superclass),
    * from the specified input database.  If the simulation is from restart,
    * these values are taken from restart and newly specified values in the
    * input file are ignored.
    * Arguments: restart_flag is true when simulation is from restart
    * Assertions: db must not be a NULL pointer.
    */
   void
   getFromInput(
      tbox::Pointer<tbox::Database> db,
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
   void
   getFromRestart();

   /*
    * Flag to determine whether this instance is registered for restart.
    */
   bool d_registered_for_restart;


};

}
}

#endif
