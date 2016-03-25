//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/patches/PatchHierarchy.C $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2141 $
// Modified:	$LastChangedDate: 2008-04-23 08:36:33 -0700 (Wed, 23 Apr 2008) $
// Description:	An AMR hierarchy of patch levels
//

#ifndef included_hier_PatchHierarchy_C
#define included_hier_PatchHierarchy_C

#include "PatchHierarchy.h"

#include <stdio.h>

#include "tbox/RestartManager.h"
#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"

#define HIER_PATCH_HIERARCHY_VERSION (2)

#ifdef DEBUG_NO_INLINE
#include "PatchHierarchy.I"
#endif

namespace SAMRAI {
    namespace hier {

/*
*************************************************************************
*									*
* Instantiate the patch hierarchy and set default values.		*
* Initialize from restart if necessary.                                 *
*									*
*************************************************************************
*/

template<int DIM>  PatchHierarchy<DIM>::PatchHierarchy(
   const std::string& object_name,
   tbox::Pointer< GridGeometry<DIM> > geometry,
   bool register_for_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(!geometry.isNull());
#endif
   d_object_name         = object_name;
   d_registered_for_restart = register_for_restart;
   d_number_levels       = 0;
   d_grid_geometry       = geometry;
   d_patch_descriptor    = VariableDatabase<DIM>::getDatabase()->
                                                   getPatchDescriptor();
   d_patch_factory       = new PatchFactory<DIM>;
   d_patch_level_factory = new PatchLevelFactory<DIM>;

   if (d_registered_for_restart) {
      tbox::RestartManager::getManager()->
         registerRestartItem(d_object_name, this);
   }
}

/*
*************************************************************************
*									*
* The destructor tells the tbox::RestartManager to remove this hierarchy *
* from the list of restart items and automatically deletes all          *
* allocated resources through smart pointers and arrays.                *
*									*
*************************************************************************
*/

template<int DIM>  PatchHierarchy<DIM>::~PatchHierarchy()
{
   if (d_registered_for_restart) {
      tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
   }
}


/*
*************************************************************************
*									*
* Create a copy of this patch hierarchy with each level refined by      *
* the given ratio and return a pointer to it.                           *
*									*
*************************************************************************
*/

template<int DIM> tbox::Pointer<hier::PatchHierarchy<DIM> > 
PatchHierarchy<DIM>::makeRefinedPatchHierarchy(
   const std::string& fine_hierarchy_name,
   const hier::IntVector<DIM>& refine_ratio,
   bool register_for_restart) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!fine_hierarchy_name.empty());
   TBOX_ASSERT(fine_hierarchy_name != d_object_name);
   TBOX_ASSERT(refine_ratio > hier::IntVector<DIM>(0));
#endif

   tbox::Pointer<hier::GridGeometry<DIM> > fine_geometry =
      d_grid_geometry->makeRefinedGridGeometry(fine_hierarchy_name + "GridGeometry",
                                               refine_ratio,
                                               register_for_restart);

   hier::PatchHierarchy<DIM>* fine_hierarchy = 
      new hier::PatchHierarchy<DIM>(fine_hierarchy_name,
                               fine_geometry,
                               register_for_restart);

   for (int ln = 0; ln < d_number_levels; ln++) {

      tbox::Pointer<hier::PatchLevel<DIM> > new_level = new hier::PatchLevel<DIM>();
      new_level->setRefinedPatchLevel(d_patch_levels[ln],
                                      refine_ratio,
                                      fine_geometry);

      new_level->setLevelNumber(ln);
      new_level->setNextCoarserHierarchyLevelNumber(ln-1);
      new_level->setLevelInHierarchy(true);
      new_level->setRatioToCoarserLevel(d_patch_levels[ln]->getRatioToCoarserLevel());
      if (ln >= fine_hierarchy->d_number_levels) {
         fine_hierarchy->d_number_levels = ln+1;
         fine_hierarchy->d_patch_levels.resizeArray(d_number_levels);
      }
      fine_hierarchy->d_patch_levels[ln] = new_level;
   }

   return(tbox::Pointer<hier::PatchHierarchy<DIM> >(fine_hierarchy));

}

/*
*************************************************************************
*									*
* Create a copy of this patch hierarchy with each level coarsened by    *
* the given ratio and return a pointer to it.                           *
*									*
*************************************************************************
*/

template<int DIM> tbox::Pointer<hier::PatchHierarchy<DIM> > 
hier::PatchHierarchy<DIM>::makeCoarsenedPatchHierarchy(
   const std::string& coarse_hierarchy_name,
   const hier::IntVector<DIM>& coarsen_ratio,
   bool register_for_restart) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!coarse_hierarchy_name.empty());
   TBOX_ASSERT(coarse_hierarchy_name != d_object_name);
   TBOX_ASSERT(coarsen_ratio > hier::IntVector<DIM>(0));
#endif

   tbox::Pointer<hier::GridGeometry<DIM> > coarse_geometry =
      d_grid_geometry->makeCoarsenedGridGeometry(coarse_hierarchy_name + "GridGeometry",
                                                 coarsen_ratio,
                                                 register_for_restart);

   hier::PatchHierarchy<DIM>* coarse_hierarchy = 
      new hier::PatchHierarchy<DIM>(coarse_hierarchy_name,
                               coarse_geometry,
                               register_for_restart);

   for (int ln = 0; ln < d_number_levels; ln++) {
      tbox::Pointer<hier::PatchLevel<DIM> > new_level = new hier::PatchLevel<DIM>();
      new_level->setCoarsenedPatchLevel(d_patch_levels[ln],
                                        coarsen_ratio,
                                        coarse_geometry);
      new_level->setLevelNumber(ln);
      new_level->setNextCoarserHierarchyLevelNumber(ln-1);
      new_level->setLevelInHierarchy(true);
      new_level->setRatioToCoarserLevel(d_patch_levels[ln]->getRatioToCoarserLevel());
      if (ln >= coarse_hierarchy->d_number_levels) {
         coarse_hierarchy->d_number_levels = ln+1;
         coarse_hierarchy->d_patch_levels.resizeArray(d_number_levels);
      }
      coarse_hierarchy->d_patch_levels[ln] = new_level;
   }

   return(tbox::Pointer<hier::PatchHierarchy<DIM> >(coarse_hierarchy));

}


/*
*************************************************************************
*									*
* Create a new patch level in the hierarchy.				*
*									*
*************************************************************************
*/

template<int DIM> void PatchHierarchy<DIM>::makeNewPatchLevel(
   const int l, 
   const IntVector<DIM>& ratio_to_coarsest,
   const BoxArray<DIM>& patch_boxes,
   const ProcessorMapping& mapping,
   const bool defer_boundary_box_creation)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(l >= 0);
   for (int i = 0; i < DIM; i++) {
      TBOX_ASSERT( ratio_to_coarsest(i) > 0 );
   }
#endif

   if (l >= d_number_levels) {
      d_number_levels = l+1;
      d_patch_levels.resizeArray(d_number_levels);
   }

   d_patch_levels[l] = d_patch_level_factory->allocate(
      patch_boxes, mapping, ratio_to_coarsest,
      d_grid_geometry, d_patch_descriptor, d_patch_factory,
      defer_boundary_box_creation);

   d_patch_levels[l]->setLevelNumber(l);
   d_patch_levels[l]->setNextCoarserHierarchyLevelNumber(l-1);
   d_patch_levels[l]->setLevelInHierarchy(true);

   if ((l > 0) && (!d_patch_levels[l-1].isNull())) {
      d_patch_levels[l]->setRatioToCoarserLevel(
         ratio_to_coarsest / (d_patch_levels[l-1]->getRatio()) );
   }

}

/*
*************************************************************************
*									*
* Remove the specified patch level from the hierarchy.			*
*									*
*************************************************************************
*/

template<int DIM> void PatchHierarchy<DIM>::removePatchLevel(const int l)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((l >= 0) && (l < d_number_levels));
#endif
   d_patch_levels[l].setNull();
   if (d_number_levels == l+1) {
     d_number_levels--;
   }
}

/*
*************************************************************************
*									*
* Writes out the class version number and the number of levels in the 	*
* hierarchy and has each patch_level write itself out.			*
* The database keys for the patch levels are given by                   *
* "level#" where # is the level number for the patch_level.             *
* The patchdata that are written to the database are determined by      *
* which those bits in the VariableDatabase<DIM> restart table.          *
*									*
* Asserts that the database pointer passed in is not NULL.		*
*									*
*************************************************************************
*/

template<int DIM> void PatchHierarchy<DIM>::putToDatabase(
   tbox::Pointer<tbox::Database> database)
{
   putToDatabase(database, 
                 VariableDatabase<DIM>::getDatabase()->getPatchDataRestartTable());
}

/*
*************************************************************************
*									*
* Writes out the class version number and the number of levels in the 	*
* hierarchy and has each patch_level write itself out.			*
* The database keys for the patch levels are given by                   *
* "level#" where # is the level number for the patch_level.             *
* The patchdata that are written to the database are determined by      *
* which those bits in the specified ComponentSelector that are     *
* set.                                                                  *
*									*
* Asserts that the database pointer passed in is not NULL.		*
*									*
*************************************************************************
*/

template<int DIM> void PatchHierarchy<DIM>::putToDatabase(
   tbox::Pointer<tbox::Database> database,
   const ComponentSelector& patchdata_write_table)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!database.isNull());
#endif

   database->putInteger("HIER_PATCH_HIERARCHY_VERSION",
                         HIER_PATCH_HIERARCHY_VERSION);
   database->putInteger("d_number_levels", d_number_levels);


   for (int i = 0; i < d_number_levels; i++) {
      std::string level_name = "level_" + tbox::Utilities::levelToString(i);

      tbox::Pointer<tbox::Database> level_database = 
         database->putDatabase(level_name);

      d_patch_levels[i]->putToDatabase(level_database, patchdata_write_table);
   }
}

/*
*************************************************************************
*									*
* Gets the database in the root database that corresponds to the object *
* name.  This method then checks the class version against restart      * 
* file version.  If they match, it creates each hierarchy level and     *
* reads in the level data.   The number of levels read from restart is  *
* the minimum of the argument max levels and the number of levels in    *
* the restart file.                                                     *
*									*
*************************************************************************
*/
template<int DIM> void PatchHierarchy<DIM>::getFromRestart(
   const int max_levels)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(max_levels > 0);
#endif

   tbox::Pointer<tbox::Database> restart_db = 
      tbox::RestartManager::getManager()->getRootDatabase();

   tbox::Pointer<tbox::Database> database;

   if ( restart_db->isDatabase(d_object_name) ) {
      database = restart_db->getDatabase(d_object_name);
   } else {
      TBOX_ERROR("PatchHierarchy<DIM>::getFromRestart() error...\n"
              << "   Restart database with name "
              << d_object_name << " not found in restart file" << std::endl);
   }

   int ver = database->getInteger("HIER_PATCH_HIERARCHY_VERSION");
   if (ver != HIER_PATCH_HIERARCHY_VERSION) {
      TBOX_ERROR("PatchHierarchy<DIM>::getFromRestart error...\n" 
          << "  object name = " << d_object_name 
          << " : Restart file version different than class version" << std::endl);
   }

   getFromDatabase(
      database, 
      VariableDatabase<DIM>::getDatabase()->getPatchDataRestartTable(),
      max_levels);
}

template<int DIM> void PatchHierarchy<DIM>::getFromDatabase(
   tbox::Pointer<tbox::Database> database,
   const ComponentSelector component_selector,
   const int max_levels)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!database.isNull());
   TBOX_ASSERT( (max_levels == -1) || (max_levels > 0) );
#endif

   d_number_levels = database->getInteger("d_number_levels");
   if (d_number_levels <= 0) {
      TBOX_ERROR("PatchHierarchy<DIM>::getFromDatabase error ...\n"
         << "  object name = " << d_object_name
         << " : `d_number_levels' is <= zero in restart file");
   }

   if (max_levels != -1) {
      d_number_levels = 
         tbox::MathUtilities<int>::Min(d_number_levels, max_levels);
   }

   d_patch_levels.resizeArray(d_number_levels);

   for (int i = 0; i< d_number_levels; i++) {
      std::string level_name = "level_" + tbox::Utilities::levelToString(i);

      tbox::Pointer<tbox::Database> level_database =
         database->getDatabase(level_name);
      d_patch_levels[i] = d_patch_level_factory->allocate(
           level_database, 
           d_grid_geometry, 
           d_patch_descriptor, 
           d_patch_factory,
           component_selector,
           false);
   }
}


template<int DIM> int PatchHierarchy<DIM>::recursivePrint( std::ostream &os ,
                                                           const std::string &border ,
                                                           unsigned short depth )
{
  int totl_npatches = 0;
  int totl_ncells = 0;
  int nlevels = getNumberOfLevels();
  os << border << "Number of levels = " << nlevels << "\n";
  if ( depth > 0 ) {
    int ln;
    for ( ln=0; ln<nlevels; ++ln ) {
      os << border << "Level " << ln << '/' << nlevels << "\n";
      tbox::Pointer<PatchLevel<DIM> > level = getPatchLevel(ln);
      level->recursivePrint( os, border+"\t", depth-1 );
      totl_npatches += level->getNumberOfPatches();
      const BoxArray<DIM> &level_boxes = level->getBoxes();
      for ( int pn=0; pn<level_boxes.size(); ++pn ) {
         totl_ncells += level_boxes[pn].size();
      } 
    }
  }
  os << border << "Total number of patches = " << totl_npatches << "\n";
  os << border << "Total number of cells = " << totl_ncells << "\n";
  return 0;
}

}
}

#endif
