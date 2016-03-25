//
// File:	SkeletonGridGeometry.h
// Package:	SAMRAI geometry package
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Skeleton grid geometry for an AMR hierarchy.
//

#ifndef included_geom_SkeletonGridGeometry
#define included_geom_SkeletonGridGeometry

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_tbox_Serializable
#include "tbox/Serializable.h"
#endif
#ifndef included_xfer_Geometry
#include "Geometry.h"
#endif

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
 * from the xfer::Geometry<DIM> base class which is further derived from
 * the hier::GridGeometry<DIM> base class.
 *
 * An object of this class requires parameters to be read from input to
 * create the hier::BoxArray that stores the index space in the hier::GridGeometry
 * superclass.  Also, data must be written to and read from files for restart.
 * The input and restart data are summarized as follows:
 *
 * \verbatim
 * Input: 
 *    Required keyword assignment: domain_boxes  
 *    Optional keyword assignments: using_original_location_indices
 * 
 * A sample input file entry might look like:
 *
 *    domain_boxes = [(0,0) , (50,50)], [(20,25), (39, 17)]
 *    use_original_location_indices = TRUE
 *
 * \endverbatim
 *
 * The optional keyword use_original_location_indices is used to
 * tell whether to construct BoundaryBox objects with the original
 * location index scheme that has existed in SAMRAI v. 1.4.0 and earlier
 * or to use a new scheme that is compatible with the
 * hier_BoundaryLookupTable.  This only makes a difference for BoundaryBox
 * objects in 3D with codimension 2.  If this keyword is omitted from
 * input, the default value is true for  dimension 3 and is irrelevant
 * for any other dimension.
 *
 * @see xfer::Geometry
 * @see hier::GridGeometry
 */

template<int DIM> class SkeletonGridGeometry 
: public xfer::Geometry<DIM>,
  public tbox::Serializable
{
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
    * Errors: passing in a null database pointer or an empty string
    * will result in an unrecoverable assertion.
    */
   SkeletonGridGeometry(const string& object_name,
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
    * Errors: passing in an empty string, or null data pointers will
    * result in an unrecoverable assertion.
    */
   SkeletonGridGeometry(const string& object_name,
                              const hier::BoxArray<DIM>& level_domain,
                              bool register_for_restart = true);

   /**
    * Destructor for SkeletonGridGeometry unregisters the object 
    * with the restart manager if previously registered.
    */
   virtual ~SkeletonGridGeometry<DIM>();
   
   /**
    * Create and return a pointer to a refined version of this Cartesian grid
    * geometry object. This function is pure virtual in the hier_GridGeometry base class.
    */
   tbox::Pointer<hier::GridGeometry<DIM> > makeRefinedGridGeometry(
      const string& fine_geom_name,
      const hier::IntVector<DIM>& refine_ratio,
      bool register_for_restart) const;

   /**
    * Create and return a pointer to a coarsened version of this Cartesian grid
    * geometry object. This function is pure virtual in the hier_GridGeometry base class.
    */
   tbox::Pointer<hier::GridGeometry<DIM> > makeCoarsenedGridGeometry(
      const string& coarse_geom_name,
      const hier::IntVector<DIM>& coarsen_ratio,
      bool register_for_restart) const;

   /*
    * Compute grid data for patch and assign new geom_SkeletonPatchGeometry
    * object to patch.  This function is pure virtual in the hier_GridGeometry
    * base class.
    */
   void setGeometryDataOnPatch(
      hier::Patch<DIM>& patch,
      const hier::IntVector<DIM>& ratio_to_level_zero,
      const tbox::Array< tbox::Array<bool> >& touches_regular_bdry,
      const tbox::Array< tbox::Array<bool> >& touches_periodic_bdry) const;
   /**
    * Print class data representation.
    */
   void printClassData(ostream& os) const;

   /**
    * Writes the state of the SkeletonGridGeometry object to the database.
    * 
    * When assertion checking is active, db cannot be a null database pointer.
    */
   virtual void putToDatabase(tbox::Pointer<tbox::Database> db);

private:

   /*
    * Reads in d_physical_domain (from hier::GridGeometry superclass),
    * from the specified input database.  If the simulation is from restart,
    * these values are taken from restart and newly specified values in the
    * input file are ignored.
    * Arguments: restart_flag is true when simulation is from restart
    * Assertions: db must not be a NULL pointer.
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
    * Create default operator for Skeleton grid geometry.
    */
   void makeStandardOperators();

   /*
    * String name for object used in error reporting and restart operations.
    * Boolean is set in constructor and determines whether object should
    * dump its state to restart files during program execution.
    */
   string d_object_name;
   bool d_registered_for_restart;

   bool d_using_original_locations;

};

}
}

#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "SkeletonGridGeometry.C"
#endif
