//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/trunk/source/geometry/skeleton/grid_geom/BlockGridGeometry.h $
// Package:	SAMRAI multiblock package
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 878 $
// Modified:	$LastChangedDate: 2006-01-09 16:55:30 -0800 (Mon, 09 Jan 2006) $
// Description: Block grid geometry for an AMR hierarchy.
//

#ifndef included_geom_BlockGridGeometry
#define included_geom_BlockGridGeometry

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
 */

template<int DIM> class BlockGridGeometry 
: public xfer::Geometry<DIM>,
  public tbox::Serializable
{
public:
   /**
    * Constructor for BlockGridGeometry initializes data
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
   BlockGridGeometry(const std::string& object_name,
                     tbox::Pointer<tbox::Database> input_db,
                     const int block_number,
                     bool register_for_restart = true);

   /**
    * Constructor for BlockGridGeometry sets index space domain 
    * based on arguments. The constructor also registers this object
    * for restart using the specified object name, when the boolean
    * argument is true.  Whether object will write its state to restart 
    * files during program execution is determined by this argument.  
    * Note that it has a default state of true.
    *
    * Errors: passing in an empty string, or null data pointers will
    * result in an unrecoverable assertion.
    */
   BlockGridGeometry(const std::string& object_name,
                           const hier::BoxArray<DIM>& level_domain,
                           bool register_for_restart = true);

   /**
    * Destructor for BlockGridGeometry unregisters the object 
    * with the restart manager if previously registered.
    */
   virtual ~BlockGridGeometry<DIM>();
   
   /**
    * Create and return a pointer to a refined version of this Cartesian grid
    * geometry object. This function is pure virtual in the hier_GridGeometry base class.
    */
   tbox::Pointer<hier::GridGeometry<DIM> > makeRefinedGridGeometry(
      const std::string& fine_geom_name,
      const hier::IntVector<DIM>& refine_ratio,
      bool register_for_restart) const;

   /**
    * Create and return a pointer to a coarsened version of this Cartesian grid
    * geometry object. This function is pure virtual in the hier_GridGeometry base class.
    */
   tbox::Pointer<hier::GridGeometry<DIM> > makeCoarsenedGridGeometry(
      const std::string& coarse_geom_name,
      const hier::IntVector<DIM>& coarsen_ratio,
      bool register_for_restart) const;

   /*
    * Compute grid data for patch and assign new BlockPatchGeometry
    * object to patch.  This function is pure virtual in the hier_GridGeometry
    * base class.
    */
   void setGeometryDataOnPatch(
      hier::Patch<DIM>& patch,
      const hier::IntVector<DIM>& ratio_to_level_zero,
      const tbox::Array< tbox::Array<bool> >& touches_regular_bdry,
      const tbox::Array< tbox::Array<bool> >& touches_periodic_bdry) const;

   int getBlockNumber() const
   {
      return (d_block_number);
   }

   /**
    * Print class data representation.
    */
   virtual void printClassData(std::ostream& os) const;

   /**
    * Writes the state of the BlockGridGeometry object to the database.
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
    * Create default operator for Block grid geometry.
    */
   void makeStandardOperators();

   /*
    * String name for object used in error reporting and restart operations.
    * Boolean is set in constructor and determines whether object should
    * dump its state to restart files during program execution.
    */
   std::string d_object_name;
   bool d_registered_for_restart;

   bool d_using_original_locations;

   int d_block_number;

};

}
}

#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "BlockGridGeometry.C"
#endif
