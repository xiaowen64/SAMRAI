//
// File:	PatchHierarchy.h
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 647 $
// Modified:	$Date: 2005-10-03 17:37:15 -0700 (Mon, 03 Oct 2005) $
// Description:	An AMR hierarchy of patch levels
//

#ifndef included_hier_PatchHierarchy
#define included_hier_PatchHierarchy

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_tbox_Array
#include "tbox/Array.h"
#endif
#ifndef included_String
#include <string>
using namespace std;
#define included_String
#endif
#ifndef included_hier_BasePatchHierarchy
#include "BasePatchHierarchy.h"
#endif
#ifndef included_hier_ComponentSelector
#include "ComponentSelector.h"
#endif
#ifndef included_hier_GridGeometry
#include "GridGeometry.h"
#endif
#ifndef included_hier_PatchDescriptor
#include "PatchDescriptor.h"
#endif
#ifndef included_hier_PatchFactory
#include "PatchFactory.h"
#endif
#ifndef included_hier_PatchLevel
#include "PatchLevel.h"
#endif
#ifndef included_hier_PatchLevelFactory
#include "PatchLevelFactory.h"
#endif
#ifndef included_hier_VariableDatabase
#include "VariableDatabase.h"
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
#ifndef included_tbox_DescribedClass
#include "tbox/DescribedClass.h"
#endif


namespace SAMRAI {
    namespace hier {

/**
 * Class PatchHierarchy<DIM> maintains the array of patch levels that
 * describe the AMR hierarchy.  All patches on all levels in the patch
 * hierarchy use the same patch descriptor instance set by the constructor.
 *
 * @see hier::PatchLevel
 * @see hier::PatchDescriptor
 */

template<int DIM> class PatchHierarchy 
:
public BasePatchHierarchy<DIM>
{
public:

   /**
    * The constructor for the PatchHierarchy initializes the number of
    * levels to zero, sets the geometry for the PatchHierarchy, and
    * registers the PatchHierarchy for restart with the specified name
    * when the boolean argument is true.  Whether the hierarchy object 
    * will write its state to restart files during program execution 
    * is determined by this argument.  Note that it has a default 
    * state of true.
    *
    * Errors: passing in an empty string or a null grid geometry pointer 
    * will result in an unrecoverable assertion when assertion checking is
    * active.
    */ 
   PatchHierarchy(
      const string& object_name,
      tbox::Pointer< GridGeometry<DIM> > geometry,
      bool register_for_restart = true);

   /**
    * Destructor for patch hierarchy objects.
    */
   ~PatchHierarchy<DIM>();


  /**
    * Create and return pointer to a patch hierarchy that is a refined
    * version of this patch hierarchy object.  That is, the data
    * members of the returned patch hierarchy are set by refining
    * information on this hierarchy by the given ratio.  The refined
    * hierarchy will cover the same physical space as this hierarchy
    * and will have the same number of levels and same mapping of
    * patches to processors on each level.  However, the index space
    * of each level will be refined by the specified ratio.  Note that
    * this function does not allocate patch data so this must be done
    * before any data operations can be performed on the new
    * hierarchy.
    */
   tbox::Pointer<hier::PatchHierarchy<DIM> > makeRefinedPatchHierarchy(
      const string& fine_hierarchy_name,
      const hier::IntVector<DIM>& refine_ratio,
      bool register_for_restart) const;

   /**
    * Create and return pointer to a patch hierarchy that is a
    * coarsened version of this patch hierarchy object.  That is, the
    * data members of the returned patch hierarchy are set by
    * coarsening information on this hierarchy by the given ratio.
    * The coarsened hierarchy will cover the same physical space as
    * this hierarchy and will have the same number of levels and same
    * mapping of patches to processors on each level.  However, the
    * index space of each level will be coarsened by the specified
    * ratio.  Note that this function does not allocate patch data so
    * this must be done before any data operations can be performed on
    * the new hierarchy.
    */
   tbox::Pointer<hier::PatchHierarchy<DIM> > makeCoarsenedPatchHierarchy(
      const string& coarse_hierarchy_name,
      const hier::IntVector<DIM>& coarsen_ratio,
      bool register_for_restart) const;

   /**
    * Construct new patch level in hierarchy at given level number.
    * Level construction involves making patches whose boxes correspond 
    * to given box arrays and assigning those patches to the patch level.  
    * Boxes are assigned to processors based on the processor mapping.
    * The argument ratio to coarsest gives the ratio to the reference 
    * refinement level (typically the coarsest level; i.e., level zero). 
    */ 
   void makeNewPatchLevel(
      const int l,
      const IntVector<DIM>& ratio_to_coarsest,
      const BoxArray<DIM>& patch_boxes,
      const ProcessorMapping& mapping);

   /**
    * Remove PatchLevel and adjust number of levels accordingly.
    */
   void removePatchLevel(const int l);

   /**
    * Return a pointer to the specified patch level.
    */
   tbox::Pointer< BasePatchLevel<DIM> > getPatchLevel(const int l) const;

   /**
    * Return a pointer to the patch descriptor used for the patches in
    * the patch hierarchy.
    */
   tbox::Pointer< PatchDescriptor<DIM> > getPatchDescriptor() const;

   /**
    * Returns true if the array of patch levels contains the specified
    * patch level.
    */
   bool levelExists(const int l) const;

   /**
    * Returns true if the array of patch levels contains a patch level
    * finer than the specified patch level. Otherwise, false is returned.
    */
   bool finerLevelExists(const int l) const;

   /**
    * Return the number of levels that currently exist in the hierarchy.
    */ 
   int getNumberLevels() const;

   /**
    * Return the level number of the finest resolution patch level residing
    * in the hierarchy.
    */
   int getFinestLevelNumber() const;

   /**
    * Set the factory used to create patch objects.  If a factory is not
    * specified, then the default factory will create patch objects of
    * type Patch<DIM>.
    */
   void setPatchFactory(tbox::Pointer< PatchFactory<DIM> > factory);

   /**
    * Set the factory used to create patch level objects.  If a factory
    * is not specified, then the default factory will create patch level
    * objects of type PatchLevel<DIM>.
    */
   void setPatchLevelFactory(tbox::Pointer< PatchLevelFactory<DIM> > factory);

   /**
    * Return a pointer to the grid geometry object.
    */
   tbox::Pointer< GridGeometry<DIM> > getGridGeometry() const;

   /**
    * Writes the state of the PatchHierarchy object and the PatchLevels
    * it contains to the database.  It should be noted that
    * only those patch data which have been registered for restart with 
    * the VariableDatabase<DIM> will be written to the database.
    * This method implements the pure virtual method in tbox::Serializable 
    * class which is used by the tbox::RestartManager for writing the 
    * PatchHierarchy to a restart file.
    *
    * When assertion checking is active, the database pointer must be non-null.
    */
   void putToDatabase(tbox::Pointer<tbox::Database> database);

   /**
    * Writes the state of the PatchHierarchy object and the PatchLevels
    * it contains to the database.  Only those patchdata corresponding to 
    * the set bits in the ComponentSelector are written to the 
    * specified database.  
    *
    * When assertion checking is active, the database pointer must be non-null.
    */
   void putToDatabase(tbox::Pointer<tbox::Database> database,
                      const ComponentSelector& patchdata_write_table);

   /**
    * Read in the entire hierarchy from the restart file.  The database 
    * from which the restart data is read is determined by the object_name
    * specified in the constructor.
    * 
    * Notes:
    * 
    *    
    *    -
    *        This method handles the memory allocation for each PatchLevel 
    *        it reads in.  
    * 
    *    - 
    *        The number of levels read in is the minimum of the max_levels 
    *        argument and the number of levels stored in the database. 
    *
    * 
    *  
    * When assertion checking is active, the max_levels argument must be 
    * greater than zero.  An unrecoverable exception will result if the
    * database cannot be found in the restart file or the data in the 
    * restart file is erroneous.
    */
   void getFromRestart(const int max_levels);

   /**
    * Read in the entire hierarchy from the specified database.  The
    * component_selector argument specifies which patch data components 
    * are to be read in from the database.  The max_levels argument 
    * indicates the highest level in the hierarchy that is read in 
    * from the database.
    * 
    * Notes:
    * 
    *    
    *    -
    *        Warning messages will be printed to the log file if 
    *        any patch data component specified in the 
    *        component_selector cannot be found in the database.
    *
    *    -
    *        By default the max_levels argument is set to -1 which
    *        indicates that all patch levels in hierarchy stored
    *        in the database should be read in.
    * 
    *    -
    *        This method handles the memory allocation for each PatchLevel 
    *        it reads in.  
    * 
    *
    * Assertion checks:
    * 
    *
    *    -
    *        The database argument must not be null.
    * 
    *    -
    *        The number of levels (if given) must be greater than zero.
    *
    */
   void getFromDatabase(tbox::Pointer<tbox::Database> database,
                        const ComponentSelector component_selector,
                        const int max_levels = -1);



   /*!
    * @brief Print a patch hierarchy to varying details.
    *
    * If depth>0, print function will be called for each level in the hierarchy.
    */
int recursivePrint( ostream &os ,
                    const string &border=string() ,
                    unsigned short depth=0 );

private:
   string                                       d_object_name;
   bool                                         d_registered_for_restart;
   int                                          d_number_levels;
   tbox::Array< tbox::Pointer< PatchLevel<DIM> > > d_patch_levels;
   tbox::Pointer< GridGeometry<DIM> >             d_grid_geometry;
   tbox::Pointer< PatchDescriptor<DIM> >          d_patch_descriptor;
   tbox::Pointer< PatchFactory<DIM> >             d_patch_factory;
   tbox::Pointer< PatchLevelFactory<DIM> >        d_patch_level_factory;
};

}
}

#ifndef DEBUG_NO_INLINE
#include "PatchHierarchy.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "PatchHierarchy.C"
#endif
