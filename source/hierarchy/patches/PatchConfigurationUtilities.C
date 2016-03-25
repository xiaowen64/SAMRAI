//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/patches/PatchConfigurationUtilities.C $
// Package:     SAMRAI hierarchy
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2141 $
// Modified:    $LastChangedDate: 2008-04-23 08:36:33 -0700 (Wed, 23 Apr 2008) $
// Description: Utilities class with routines for understanding spatial
//              relationships among patches
//

#ifndef included_hier_PatchConfigurationUtilities_C
#define included_hier_PatchConfigurationUtilities_C

#include "PatchConfigurationUtilities.h"
#include "BoundaryLookupTable.h"


#include "tbox/Utilities.h"

#include "BoxTree.h"
#include "BoundaryLookupTable.h"

namespace SAMRAI {
    namespace hier {

/*
*************************************************************************
*                                                                       *
* Patch operation utilities constructor and destructor.                 *
*                                                                       *
*************************************************************************
*/

template<int DIM>
PatchConfigurationUtilities<DIM>::PatchConfigurationUtilities(
   const std::string& object_name,
   tbox::Pointer< PatchHierarchy<DIM> > hierarchy)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!object_name.empty());
#endif
   d_object_name     = object_name; 
   d_patch_hierarchy = hierarchy;
}

template<int DIM>
PatchConfigurationUtilities<DIM>::~PatchConfigurationUtilities()
{
   for (int ln = 0; ln < d_patch_level_info.size(); ++ln) {
      clearPatchLevelInfo(ln);
   }
   d_patch_level_info.resizeArray(0);
   d_patch_hierarchy.setNull(); 
}

/*
*************************************************************************
*                                                                       *
* Determine whether information is set for level or patch.              *
*                                                                       *
*************************************************************************
*/

template<int DIM>
bool PatchConfigurationUtilities<DIM>::levelIsSet(
   const LevelNumber& level_num) const
{
   int ln = level_num.ln;
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( (0 <= ln) && (ln < d_patch_level_info.size()) );
#endif
   return( d_patch_level_info[ln].d_is_all_set );
}

template<int DIM>
bool PatchConfigurationUtilities<DIM>::patchIsSet(
   const PatchNumber& patch_num,
   const LevelNumber& level_num) const
{
   bool ret_val = false;
   if ( levelIsSet(level_num) ) {
      int ln = level_num.ln;  
      int pn = patch_num.pn;  
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT( (0 <= pn) &&
              (pn < d_patch_level_info[ln].d_patch_level->getNumberOfPatches()) );
#endif
      ret_val = ( d_patch_level_info[ln].d_patch_info[pn] ? true : false );
   }
   return( ret_val );
}

/*
*************************************************************************
*                                                                       *
* Initialize patch data info for single level or entire hierarchy.      *
*                                                                       *
*************************************************************************
*/

template<int DIM>
void PatchConfigurationUtilities<DIM>::initialize(
   tbox::Pointer< PatchLevel<DIM> > level)
{
   if ( d_patch_hierarchy.isNull() ) {

      /*
       *  No patch hierarchy was passed to ctor. Attempt to initialize
       *  object with single patch level (if level is non-null).
       * 
       *  Check to see if object is already initialized with a different level;
       *  if it is, print warning message and do nothing.
       */

      if ( !level.isNull() ) {

         if ( d_patch_level_info.size() > 0 ) {

            if ( !levelMatch(level, 0) ) {

               TBOX_WARNING("PatchConfigurationUtilities<DIM>::initialize warning..."
                            << "\n  object named = " << d_object_name
                            << "\n  Object is set up to handle only a single level"
                            << " and it is already initialized with a level."
                            << "\n  But, that level does not match argument level."
                            << "\n  Thus, initialize() routine does nothing."
                            << std::endl);

            } 

         } else {

            allocatePatchLevelInfoArray(1);
            setPatchLevelInfo( level, 0 );

         }

      }  // else level is null and there is nothing to do...

   } else {

      /*
       *  A patch hierarchy was passed to ctor. Attempt to initialize
       *  information for given patch level if non-null.  Otherwise,
       *  set all uninitialized patch level information for hierarchy.
       * 
       *  If a non-null patch level argument is passed, check to see if
       *  the level lives in the patch hierarchy or if the object is already 
       *  initialized with a different level at the corresponding PatchLevelInfo 
       *  array slot.  If the level is not in the hierarchy or differs from 
       *  the one with which the object is already initialized,
       *  print warning message and do nothing.
       */

      if ( d_patch_level_info.size() == 0 ) {
         allocatePatchLevelInfoArray( d_patch_hierarchy->getNumberOfLevels() );
      }

      if ( level.isNull() ) {

         for (int ln = d_patch_level_info.size() - 1; ln >= 0; --ln) {
            tbox::Pointer< PatchLevel<DIM> > t_level = 
               d_patch_hierarchy->getPatchLevel(ln);
            setPatchLevelInfo( t_level );
         }

      } else {

         if ( !levelMatch(level) ) {

            TBOX_WARNING("PatchConfigurationUtilities<DIM>::initialize warning..."
                         << "\n  object named = " << d_object_name
                         << "\n  Object is set up to handle a hierarchy of levels."
                         << "\n  But, argument level does not live in the patch" 
                         << " hierarchy or it does not match the level with" 
                         << "\n  level number = " << level->getLevelNumber() 
                         << " that was used to initialize object earlier."
                         << "\n  Thus, initialize() routine does nothing."
                         << std::endl);

         } else {

            setPatchLevelInfo( level );

         }

      }

   }  // else, hierarchy is non-null
}

/*
*************************************************************************
*                                                                       *
* Clear patch data info for single level or entire hierarchy.           *
*                                                                       *
*************************************************************************
*/

template<int DIM>
void PatchConfigurationUtilities<DIM>::clear(
   tbox::Pointer< PatchLevel<DIM> > level)
{
   if ( level.isNull() ) {

      for (int ln = 0; ln < d_patch_level_info.size(); ++ln) {
         clearPatchLevelInfo(ln);
      }

   } else {

      if ( d_patch_level_info.size() > 0 ) {

         if ( d_patch_hierarchy.isNull() ) {

            /*
             *  No patch hierarchy was passed to ctor. Attempt to clear
             *  single patch level from object (if level is non-null).
             *
             *  Check to see if level matches that which object is initialized;
             *  if not, print warning message and do nothing.
             */

            if ( !levelMatch(level, 0) ) {

               TBOX_WARNING("PatchConfigurationUtilities<DIM>::clear warning..."
                            << "\n  object named = " << d_object_name
                            << "\n  Object is set up to handle only a single level"
                            << " and is initialized with a level."
                            << "\n  But, that level does not match argument level."
                            << "\n  Thus, clear() routine does nothing."
                            << std::endl);

            } else {
      
               clearPatchLevelInfo(0);
               d_patch_level_info.resizeArray(0);

            }

         } else {

            /*
             *  A patch hierarchy was passed to ctor. Attempt to clear
             *  information for given patch level.
             *
             *  Check to see if level matches that which object is initialized;
             *  if not, print warning message and do nothing.
             * 
             */

            if ( !levelMatch(level) ) {

               TBOX_WARNING("PatchConfigurationUtilities<DIM>::clear warning..."
                            << "\n  object named = " << d_object_name
                            << "\n  Object is set up to handle a hierarchy of levels."
                            << "\n  But, argument level does not live in the patch" 
                            << " hierarchy or it does not match the level with" 
                            << "\n  level number = " << level->getLevelNumber() 
                            << " that was used to initialize object earlier."
                            << "\n  Thus, clear() routine does nothing."
                            << std::endl);

            } else {

               clearPatchLevelInfo( level->getLevelNumber() );

            }

         }  // else hierarchy is non-null

      }  // if some patch level info initialized

   } // level is non-null 

   /*
    *  Clear patch level info array if necessary.
    */

   bool some_level_set = false;
   for (int ln = 0; ln < d_patch_level_info.size(); ++ln) {
      some_level_set |= 
         ( d_patch_level_info[ln].d_neighbor_patch_info_is_set ||
           d_patch_level_info[ln].d_finer_level_info_is_set ||
           d_patch_level_info[ln].d_coarser_level_info_is_set );
   }

   if ( !some_level_set ) {
      d_patch_level_info.resizeArray(0);
   }

}

/*
*************************************************************************
*                                                                       *
* Functions to retrieve patch neighbor information.                     *
*                                                                       *
*************************************************************************
*/

template<int DIM>
const tbox::Array< typename PatchConfigurationUtilities<DIM>::NeighborPatchInfo >&
PatchConfigurationUtilities<DIM>::getNodeNeighborInfo(
   const PatchNumber& patch_num,
   const LevelNumber& level_num) const
{
   return( getCodimensionNeighborPatchInfo(DIM, patch_num, level_num) );
}

template<int DIM>
const tbox::Array< typename PatchConfigurationUtilities<DIM>::NeighborPatchInfo >&
PatchConfigurationUtilities<DIM>::getEdgeNeighborInfo(
   const PatchNumber& patch_num,
   const LevelNumber& level_num) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(DIM > 1);
#endif

   return( getCodimensionNeighborPatchInfo(DIM-1, patch_num, level_num) );
}

template<int DIM>
const tbox::Array< typename PatchConfigurationUtilities<DIM>::NeighborPatchInfo >&
PatchConfigurationUtilities<DIM>::getFaceNeighborInfo(
   const PatchNumber& patch_num,
   const LevelNumber& level_num) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(DIM > 2);
#endif

   return( getCodimensionNeighborPatchInfo(DIM-2, patch_num, level_num) );
}

template<int DIM>
const tbox::Array< typename PatchConfigurationUtilities<DIM>::NeighborPatchInfo >&
PatchConfigurationUtilities<DIM>::getCodimensionNeighborPatchInfo(
   int codim, 
   const PatchNumber& patch_num,
   const LevelNumber& level_num) const
{
   int ln = ( d_patch_hierarchy.isNull() ? 0 : level_num.ln );
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( (codim > 0) && (codim <= DIM) );
   TBOX_ASSERT( (0 <= ln) && (ln < d_patch_level_info.size()) );
   TBOX_ASSERT( d_patch_level_info[ln].d_neighbor_patch_info_is_set );
   TBOX_ASSERT( (0 <= patch_num.pn) && 
           (patch_num.pn < d_patch_level_info[ln].d_patch_level->getNumberOfPatches()) );
#endif
   return( d_patch_level_info[ln].d_patch_info[patch_num.pn]->
              d_neighbors[codim-1] );
}

/*
*************************************************************************
*                                                                       *
* Functions to retrieve information about finer level patches that      *
* overlay a given patch.                                                *
*                                                                       *
*************************************************************************
*/

template<int DIM>
const tbox::Array<int>&
PatchConfigurationUtilities<DIM>::getFinerLevelOverlapPatchIndices(
   const PatchNumber& patch_num,
   const LevelNumber& level_num) const
{
   if ( d_patch_hierarchy.isNull() ||
        (level_num.ln == d_patch_hierarchy->getFinestLevelNumber()) ) {
      return( d_empty_array );
   } else {
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( (0 <= level_num.ln) && (level_num.ln < d_patch_level_info.size()) );
   TBOX_ASSERT( d_patch_level_info[level_num.ln].d_finer_level_info_is_set );
   TBOX_ASSERT( (0 <= patch_num.pn) && 
           (patch_num.pn < 
            d_patch_level_info[level_num.ln].d_patch_level->getNumberOfPatches()) );
#endif
      return( d_patch_level_info[level_num.ln].d_patch_info[patch_num.pn]->
                 d_finer_level_overlap_patch_indices );
   }
}

template<int DIM>
const Box<DIM>&
PatchConfigurationUtilities<DIM>::getCoarsenedFinerLevelPatchBox(
   const PatchNumber& finer_level_patch_num,
   const LevelNumber& level_num) const
{
   if ( d_patch_hierarchy.isNull() ||
        (level_num.ln == d_patch_hierarchy->getFinestLevelNumber()) ) {
      return( d_empty_box );
   } else {
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( (0 <= level_num.ln) && (level_num.ln < d_patch_level_info.size() - 1) );
   TBOX_ASSERT( d_patch_level_info[level_num.ln].d_finer_level_info_is_set );
   TBOX_ASSERT( (0 <= finer_level_patch_num.pn) &&
           (finer_level_patch_num.pn <
            d_patch_level_info[level_num.ln + 1].d_patch_level->getNumberOfPatches()) );
#endif
      return( d_patch_level_info[level_num.ln].
                 d_coarsened_fine_patch_boxes[finer_level_patch_num.pn] );
   }
}

/*
*************************************************************************
*                                                                       *
* Functions to retrieve information about coarser level patches that    *
* overlay a given patch.                                                *
*                                                                       *
*************************************************************************
*/

template<int DIM>
const tbox::Array<int>&
PatchConfigurationUtilities<DIM>::getCoarserLevelOverlapPatchIndices(
   const PatchNumber& patch_num,
   const LevelNumber& level_num) const
{
   if ( d_patch_hierarchy.isNull() || (level_num.ln == 0) ) {
      return( d_empty_array );
   } else {
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( (0 <= level_num.ln) && (level_num.ln < d_patch_level_info.size()) );
   TBOX_ASSERT( d_patch_level_info[level_num.ln].d_coarser_level_info_is_set );
   TBOX_ASSERT( (0 <= patch_num.pn) && 
           (patch_num.pn < 
            d_patch_level_info[level_num.ln].d_patch_level->getNumberOfPatches()) );
#endif
      return( d_patch_level_info[level_num.ln].d_patch_info[patch_num.pn]->
                 d_coarser_level_overlap_patch_indices );
   }
}

template<int DIM>
const Box<DIM>&
PatchConfigurationUtilities<DIM>::getRefinedCoarserLevelPatchBox(
   const PatchNumber& coarser_level_patch_num,
   const LevelNumber& level_num) const
{
   if ( d_patch_hierarchy.isNull() || (level_num.ln == 0) ) {
      return( d_empty_box );
   } else {
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( (0 < level_num.ln) && (level_num.ln < d_patch_level_info.size()) );
   TBOX_ASSERT( d_patch_level_info[level_num.ln].d_coarser_level_info_is_set );
   TBOX_ASSERT( (0 <= coarser_level_patch_num.pn) &&
           (coarser_level_patch_num.pn <
            d_patch_level_info[level_num.ln - 1].d_patch_level->getNumberOfPatches()) );
#endif
      return( d_patch_level_info[level_num.ln].
                 d_refined_coarse_patch_boxes[coarser_level_patch_num.pn] );
   }
}

/*
*************************************************************************
*                                                                       *
* Function to print internal data structures for a single level or      *
* entire patch hierarchy.                                               *
*                                                                       *
*************************************************************************
*/

template<int DIM>
void PatchConfigurationUtilities<DIM>::printPatchLevelData(
   const LevelNumber& level_num,
   std::ostream& os) const
{
   os << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
   os << "PatchConfigurationUtilities<DIM> data for object named ";
   os << d_object_name;
   os << "\n----------------------------------------------------------" << std::endl;

   int level_number_in = level_num.ln; 

   int start_ln = 0;
   int end_ln = -1;
   if ( !d_patch_hierarchy.isNull() ) {
      if ( (0 <= level_number_in) && 
           (level_number_in < d_patch_hierarchy->getNumberOfLevels()) &&
           (level_number_in < d_patch_level_info.size()) ) {
         start_ln = level_number_in; 
         end_ln = level_number_in; 
      } else {
         if (d_patch_level_info.size() > 0) {
            end_ln = d_patch_level_info.size() - 1;
         }
      }
   } else {
      if (d_patch_level_info.size() > 0) {
         end_ln = 0;
      }
   }

   os << "d_patch_hierarchy = " 
      << (PatchHierarchy<DIM>*)d_patch_hierarchy << std::endl;
   if ( !d_patch_hierarchy.isNull() ) {
      os << "Number of levels in patch hierarchy = "
         << d_patch_hierarchy->getNumberOfLevels() << std::endl;
   }
   os << "d_patch_level_info size = " << d_patch_level_info.size() << std::endl;

   for (int ln = start_ln; ln <= end_ln; ++ln) {
      const PatchLevelInfo& level_info = d_patch_level_info[ln];
      os << "\n  Patch level info for level number = " << ln << std::endl;
      os << "    d_neighbor_patch_info_is_set = " 
         << level_info.d_neighbor_patch_info_is_set << std::endl;
      os << "    d_finer_level_info_is_set = " 
         << level_info.d_finer_level_info_is_set << std::endl;
      os << "    d_coarser_level_info_is_set = " 
         << level_info.d_coarser_level_info_is_set << std::endl;
      os << "    d_is_all_set = " 
         << level_info.d_is_all_set << std::endl;
      if ( level_info.d_is_all_set || 
           level_info.d_neighbor_patch_info_is_set ) {
         printPatchLevelInfo(os, ln); 
      }
      if ( level_info.d_finer_level_info_is_set ) {
         os << "    d_coarsened_fine_patch_boxes number = " << std::endl;
         level_info.d_coarsened_fine_patch_boxes.print(os);
         os << std::endl; 
      }
      if ( level_info.d_coarser_level_info_is_set ) {
         os << "    d_refined_coarse_patch_boxes number = " << std::endl;
         level_info.d_refined_coarse_patch_boxes.print(os);
         os << std::endl;
      }
   }
   os << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
}

/*
*************************************************************************
*                                                                       *
* Function to print internal data structures for a single patch on      *
* a patch level.                                                        *
*                                                                       *
*************************************************************************
*/

template<int DIM>
void PatchConfigurationUtilities<DIM>::printPatchData(
   const PatchNumber& patch_num,
   const LevelNumber& level_num,
   std::ostream& os) const
{
   int ln = level_num.ln;
   int pn = patch_num.pn;
   if ( d_patch_hierarchy.isNull() ) {
      ln = 0;
   } else {
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT( (0 <= ln) && (ln < d_patch_level_info.size()) );
#endif
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( (0 <= pn) &&
           (pn < d_patch_level_info[ln].d_patch_info.size()) );
#endif

   const ProcessorMapping& mapping = 
      d_patch_level_info[ln].d_patch_level->getProcessorMapping();

   if ( mapping.isMappingLocal(pn) ) {

      if ( d_patch_level_info[ln].d_neighbor_patch_info_is_set ) {
         os << "\n Neighbors for patch " << pn << ":" << std::endl;
         for (int codim = 1; codim <= DIM; ++codim) {
            const
            tbox::Array< typename PatchConfigurationUtilities<DIM>::NeighborPatchInfo >&
               neighbors = getCodimensionNeighborPatchInfo(codim,
                                                           patch_num,
                                                           level_num);
            os << "   number of neighbors of codim " << codim
               << " is " << neighbors.size() << std::endl;
            for (int ni = 0; ni < neighbors.size(); ++ni) {
               os << "neighbor[" << ni << "] = " << std::endl;
               os << "      d_neighbor_patch_number = "
                  << neighbors[ni].d_neighbor_patch_number << std::endl;
               os << "      d_neighbor_type = "
                  << neighbors[ni].d_neighbor_type << std::endl;
               os << "      d_location_index = "
                  << neighbors[ni].d_location_index << std::endl;
               os << "      d_neighbor_shift = "
                  << neighbors[ni].d_neighbor_shift << std::endl;
            }
            os << std::endl;
         }
      } else {
         os << "\n Neighbor patch info not set for level "  
            << ln << std::endl << std::endl;
      }

      if ( d_patch_level_info[ln].d_finer_level_info_is_set ) {

         const tbox::Array<int>& fine_overlap =
            d_patch_level_info[ln].d_patch_info[pn]->
               d_finer_level_overlap_patch_indices;
         os << "   number of fine overlap patches is " << fine_overlap.size()
            << " : fine patch ids = ...." << std::endl;
         int ipf = 0;
         for ( ; ipf < fine_overlap.size() - 1; ++ipf) {
            os << fine_overlap[ipf] << " , ";
         }
         if ( ipf < fine_overlap.size() ) {
            os << fine_overlap[ipf] << std::endl;
         }

      } else {
         os << "\n Finer level patch overlap info not set for level "  
            << ln << std::endl << std::endl;
      }

      if ( d_patch_level_info[ln].d_coarser_level_info_is_set ) {

         const tbox::Array<int>& coarse_overlap =
            d_patch_level_info[ln].d_patch_info[pn]->
               d_coarser_level_overlap_patch_indices;
         os << "   number of coarse overlap patches is " << coarse_overlap.size()
            << " : coarse patch ids = ...." << std::endl;
         int ipc = 0;
         for ( ; ipc < coarse_overlap.size() - 1; ++ipc) {
            os << coarse_overlap[ipc] << " , ";
         }
         if ( ipc < coarse_overlap.size() ) {
            os << coarse_overlap[ipc] << std::endl;
         }

      } else {
         os << "\n Coarser level patch overlap info not set for level "
            << ln << std::endl << std::endl;
      }

   } else {
      os << "\n patch " << pn << " is on processor "
         << mapping.getProcessorAssignment(pn) << std::endl << std::endl;
   }

}

/*
*************************************************************************
*                                                                       *
* Private function to print PatchLevelInfo data structures for a        *
* single patch level.                                                   *
*                                                                       *
*************************************************************************
*/

template<int DIM>
void PatchConfigurationUtilities<DIM>::printPatchLevelInfo(std::ostream& os,
                                                           int ln) const
{
   os << "\n==========================================================" << std::endl;
   const PatchLevelInfo& level_info = d_patch_level_info[ln];
   tbox::Pointer< PatchLevel<DIM> > level = level_info.d_patch_level;
   os << "\n  Patch info for level number = " << ln << std::endl;
   os << "    d_patch_level number = " 
      << level->getLevelNumber() << std::endl;
   os << "    d_patch_level number of patches = " 
      << level->getNumberOfPatches() << std::endl;
   os << "    d_patch_info size = " << level_info.d_patch_info.size() << std::endl;

   os << "\n  Patch Neighbor Info (rank = " << tbox::SAMRAI_MPI::getRank() << ")...." << std::endl;
   for (int ip = 0; ip < level_info.d_patch_info.size(); ++ip) {
      printPatchData(PatchNumber(ip), LevelNumber(ln), os); 
   }
   os << "==========================================================" << std::endl;
}

/*
*************************************************************************
*                                                                       *
* Private method to set info for single patch level.                    *
*                                                                       *
*************************************************************************
*/

template<int DIM>
void PatchConfigurationUtilities<DIM>::setPatchLevelInfo( 
   tbox::Pointer< PatchLevel<DIM> > level,
   int array_slot) 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!level.isNull());
#endif

   const int ln_slot = 
      ( (array_slot < 0) ? level->getLevelNumber() : array_slot );

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( (0 <= ln_slot) && (ln_slot < d_patch_level_info.size()) );
#endif

   PatchLevelInfo& level_info = d_patch_level_info[ln_slot]; 

   if ( !level_info.d_is_all_set ) {

      if ( level_info.d_patch_level.isNull() ) { 
         level_info.d_patch_level = level;
      } 

      setNeighborPatchInfo( ln_slot ); 

      if ( d_patch_hierarchy.isNull() ) {

         level_info.d_finer_level_info_is_set = false;
         level_info.d_coarser_level_info_is_set = false;
         level_info.d_is_all_set = true;

      } else {

         setFinerPatchLevelInfo( ln_slot );
         setCoarserPatchLevelInfo( ln_slot );
         if ( ( level_info.d_finer_level_info_is_set ||
                (d_patch_hierarchy->getFinestLevelNumber() == ln_slot) ) &&
              ( level_info.d_coarser_level_info_is_set || (ln_slot == 0) ) ) {
            level_info.d_is_all_set = true;
         }

         int coarser_ln_slot = ln_slot - 1;
         if ( coarser_ln_slot >= 0 ) {
            PatchLevelInfo& coarser_level_info = 
               d_patch_level_info[ coarser_ln_slot ];
            if ( coarser_level_info.d_neighbor_patch_info_is_set &&
                 !coarser_level_info.d_finer_level_info_is_set ) {
               setFinerPatchLevelInfo( coarser_ln_slot );
            }
         }

         int finer_ln_slot = ln_slot + 1;
         if ( finer_ln_slot <= d_patch_hierarchy->getFinestLevelNumber() ) {
            PatchLevelInfo& finer_level_info = 
               d_patch_level_info[ finer_ln_slot ];
            if ( finer_level_info.d_neighbor_patch_info_is_set &&
                 !finer_level_info.d_coarser_level_info_is_set ) {
               setCoarserPatchLevelInfo( finer_ln_slot );
            }
         } 

      }

   }

}

/*
*************************************************************************
*                                                                       *
* Private method to set patch neighbor info for single patch level.     *
*                                                                       *
*************************************************************************
*/

template<int DIM>
void PatchConfigurationUtilities<DIM>::setNeighborPatchInfo(
   int array_slot) 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( (0 <= array_slot) && (array_slot < d_patch_level_info.size()) );
   TBOX_ASSERT( !d_patch_level_info[array_slot].d_patch_level.isNull() );
#endif

   PatchLevelInfo& level_info = d_patch_level_info[array_slot]; 

   if ( !level_info.d_neighbor_patch_info_is_set ) { 

      tbox::Pointer< PatchLevel<DIM> > level = level_info.d_patch_level;

      const int npatches = level->getNumberOfPatches();

      level_info.d_patch_info.resizeArray(0);
      level_info.d_patch_info.resizeArray(npatches);

      tbox::Pointer< BoxTree<DIM> > box_tree = level->getBoxTree();
      const IntVector<DIM> pboxgrow(1);

      const ProcessorMapping& mapping = level->getProcessorMapping();
      for (int ip = 0; ip < npatches; ++ip) {

         if ( mapping.isMappingLocal(ip) ) {

            level_info.d_patch_info[ip] = new PatchInfo;

            tbox::Pointer< Patch<DIM> > patch = level->getPatch(ip);
            Box<DIM> pghostbox( patch->getBox() );
            pghostbox.grow(pboxgrow);

            tbox::Array<int> patch_neighbor_indices;
            box_tree->findOverlapIndices(
               patch_neighbor_indices, pghostbox);

            tbox::Array<bool> done_with_neighbor( patch_neighbor_indices.size() );
            for (int ni = 0; ni < patch_neighbor_indices.size(); ++ni) {
               done_with_neighbor[ni] = false;
            }  

            for (int codim = DIM; codim > 0; --codim) {

               findPatchNeighbors(*level_info.d_patch_info[ip],
                                  patch,
                                  level,
                                  patch_neighbor_indices,
                                  done_with_neighbor,
                                  codim);
            }

         } else {  // else patch is not local to processor set patch info to null

            level_info.d_patch_info[ip] = (PatchInfo*)0;

         }

      }  // iterate over patch indices on patch level

      level_info.d_neighbor_patch_info_is_set = true;

   }  // if level neighbor patch info not already set

}

/*
*************************************************************************
*                                                                       *
* Private method to set info about finer level patches that overlay     *
* patches on level associated with given patch level info array slot.   *
*                                                                       *
*************************************************************************
*/

template<int DIM>
void PatchConfigurationUtilities<DIM>::setFinerPatchLevelInfo(
   int array_slot)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( (0 <= array_slot) && (array_slot < d_patch_level_info.size()) );
   TBOX_ASSERT( !d_patch_level_info[array_slot].d_patch_level.isNull() );
#endif

   PatchLevelInfo& level_info = d_patch_level_info[array_slot];

   if ( !d_patch_hierarchy.isNull() &&
        !level_info.d_finer_level_info_is_set ) {

      const int finer_ln_slot = array_slot + 1;

      if ( (array_slot < d_patch_hierarchy->getFinestLevelNumber()) &&
           (!d_patch_hierarchy->getPatchLevel( array_slot ).isNull()) &&
           (!d_patch_hierarchy->getPatchLevel( finer_ln_slot ).isNull()) ) {

         tbox::Pointer< PatchLevel<DIM> > level = 
            d_patch_hierarchy->getPatchLevel( array_slot );

         const BoxArray<DIM>& boxes = level->getBoxes();

#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT( level_info.d_patch_info.size() == boxes.size() );
#endif

         tbox::Pointer< PatchLevel<DIM> > finer_level =
            d_patch_hierarchy->getPatchLevel( finer_ln_slot );

         level_info.d_coarsened_fine_patch_boxes = 
            finer_level->getBoxes();
         level_info.d_coarsened_fine_patch_boxes.
            coarsen( finer_level->getRatioToCoarserLevel() ); 

         BoxTree<DIM> box_tree( level_info.d_coarsened_fine_patch_boxes );

         const ProcessorMapping& mapping = level->getProcessorMapping();
         for (int ip = 0; ip < boxes.size(); ++ip) {

            if ( mapping.isMappingLocal(ip) ) {
#ifdef DEBUG_CHECK_ASSERTIONS
               TBOX_ASSERT( level_info.d_patch_info[ip] );
#endif
               tbox::Array<int> fine_overlap_indices;
               box_tree.findOverlapIndices(
                   level_info.d_patch_info[ip]->
                      d_finer_level_overlap_patch_indices,
                   boxes[ip] );

            }

         }  // iterate over patch level boxes

         level_info.d_finer_level_info_is_set = true;

      } else {  // there is a no finer hierarchy level

         level_info.d_finer_level_info_is_set = false;

      }

   }  // if hierarchy is non-null, and finer level info not already set

}

/*
*************************************************************************
*                                                                       *
* Private method to set info about coarser level patches that overlay   *
* patches on level associated with given patch level info array slot.   *
*                                                                       *
*************************************************************************
*/

template<int DIM>
void PatchConfigurationUtilities<DIM>::setCoarserPatchLevelInfo(
   int array_slot)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( (0 <= array_slot) && (array_slot < d_patch_level_info.size()) );
   TBOX_ASSERT( !d_patch_level_info[array_slot].d_patch_level.isNull() );
#endif

   PatchLevelInfo& level_info = d_patch_level_info[array_slot];

   if ( !d_patch_hierarchy.isNull() &&
        !level_info.d_coarser_level_info_is_set ) {

      const int coarser_ln_slot = array_slot - 1;

      if ( (array_slot > 0) &&
           (!d_patch_hierarchy->getPatchLevel( array_slot ).isNull()) &&
           (!d_patch_hierarchy->getPatchLevel( coarser_ln_slot ).isNull()) ) {

         tbox::Pointer< PatchLevel<DIM> > level = 
            d_patch_hierarchy->getPatchLevel( array_slot );

         const BoxArray<DIM>& boxes = level->getBoxes();

#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT( level_info.d_patch_info.size() == boxes.size() );
#endif

         tbox::Pointer< PatchLevel<DIM> > coarser_level =
            d_patch_hierarchy->getPatchLevel( coarser_ln_slot );

         level_info.d_refined_coarse_patch_boxes = 
            coarser_level->getBoxes();
         level_info.d_refined_coarse_patch_boxes.
            refine( level->getRatioToCoarserLevel() ); 

         BoxTree<DIM> box_tree( level_info.d_refined_coarse_patch_boxes );

         const ProcessorMapping& mapping = level->getProcessorMapping();
         for (int ip = 0; ip < boxes.size(); ++ip) {

            if ( mapping.isMappingLocal(ip) ) {
#ifdef DEBUG_CHECK_ASSERTIONS
               TBOX_ASSERT( level_info.d_patch_info[ip] );
#endif
               tbox::Array<int> coarse_overlap_indices;
               box_tree.findOverlapIndices(
                   level_info.d_patch_info[ip]->
                      d_coarser_level_overlap_patch_indices,
                   boxes[ip] );

            }

         }  // iterate over patch level boxes

         level_info.d_coarser_level_info_is_set = true;

      } else {  // there is a no coarser hierarchy level

         level_info.d_coarser_level_info_is_set = false;

      }

   }  // if hierarchy is non-null, and coarser level info not already set

}

/*
*************************************************************************
*                                                                       *
* Private method to clear info for single patch level.                  *
*                                                                       *
*************************************************************************
*/

template<int DIM>
void PatchConfigurationUtilities<DIM>::clearPatchLevelInfo(int array_slot) 
{
   d_patch_level_info[array_slot].d_patch_level.setNull();

   for (int ip = 0; 
        ip < d_patch_level_info[array_slot].d_patch_info.size();
        ++ip) {
      if ( d_patch_level_info[array_slot].d_patch_info[ip] ) {
         delete d_patch_level_info[array_slot].d_patch_info[ip];
      }
   }

   d_patch_level_info[array_slot].d_patch_info.resizeArray(0);
   
   d_patch_level_info[array_slot].d_coarsened_fine_patch_boxes.resizeBoxArray(0);

   d_patch_level_info[array_slot].d_neighbor_patch_info_is_set = false; 
   d_patch_level_info[array_slot].d_finer_level_info_is_set = false; 
   d_patch_level_info[array_slot].d_coarser_level_info_is_set = false; 
   d_patch_level_info[array_slot].d_is_all_set = false; 
}

/*
*************************************************************************
*                                                                       *
* Private method to find patch neighbors for given codim.               *
*                                                                       *
*************************************************************************
*/

template<int DIM> 
void PatchConfigurationUtilities<DIM>::findPatchNeighbors(
   PatchInfo& patch_info,
   tbox::Pointer< Patch<DIM> > patch,
   tbox::Pointer< PatchLevel<DIM> > level,
   const tbox::Array<int>& neighbor_patch_ids,
   tbox::Array<bool>& done_with_neighbor,
   int codim)
{

   const int blut_loc_idx = codim - 1;
   BoundaryLookupTable<DIM>* blut =
      BoundaryLookupTable<DIM>::getLookupTable();
   const int num_locations = 
      blut->getMaxLocationIndices()[blut_loc_idx];

   patch_info.d_neighbors[blut_loc_idx].resizeArray(num_locations);

   const int pnum = patch->getPatchNumber();
   const Box<DIM>& pbox = patch->getBox();
   const Box<DIM> pgbox( Box<DIM>::grow(patch->getBox(), IntVector<DIM>(1)) );

   const BoxArray<DIM>& level_boxes = level->getBoxes();

   /*
    * If patch has non-zero shifts, it abuts a periodic domain boundary.
    * Only the shifts with codim non-zero entries need to be checked.
    */
 
   tbox::List< IntVector<DIM> > check_shifts;
   typename tbox::List< IntVector<DIM> >::Iterator 
     is( level->getShiftsForPatch(pnum) );
   for ( ; is; is++) {
      const IntVector<DIM>& shift = is();
      int ndirs = 0;
      for (int id = 0; id < DIM; ++id) {
         if ( shift[id] != 0 ) {
            ndirs++;
         }
      }
      if (ndirs <= codim) {
         check_shifts.appendItem(shift);
      }
   }

   tbox::Array<bool> done_with_location(num_locations);
   for (int loc = 0; loc < num_locations; ++loc) {
      done_with_location[loc] = false;
   }

   int num_neighbors = 0;

   typename tbox::List< IntVector<DIM> >::Iterator 
      check_shift( check_shifts );
   bool zero_shift = true;

   while (check_shift || zero_shift) {

      hier::IntVector<DIM> nshift(0);
      if (!zero_shift) {
         nshift = -check_shift();
      }

      for (int loc = 0; loc < num_locations; ++loc) {

         if ( !done_with_location[loc] ) {

            const tbox::Array<int>& dirs = blut->getDirections(loc, codim);

            Box<DIM> bregion( pbox );
            for (int i = 0; i < dirs.size(); ++i) {
               if ( blut->isUpper(loc, codim, i) ) {
                  bregion.lower(dirs[i]) = pbox.upper(dirs[i]) + 1;
                  bregion.upper(dirs[i]) = pbox.upper(dirs[i]) + 1;
               } else {
                  bregion.lower(dirs[i]) = pbox.lower(dirs[i]) - 1;
                  bregion.upper(dirs[i]) = pbox.lower(dirs[i]) - 1;
               }
            }

            for (int ni = 0; ni < neighbor_patch_ids.size(); ++ni) {

               if ( !done_with_neighbor[ni] ) {

                  const int npnum = neighbor_patch_ids[ni];

                  if ( !zero_shift || (npnum != pnum) ) {

                     const Box<DIM> shifted_nbox( 
                        Box<DIM>::shift(level_boxes[npnum], nshift) );

                     const Box<DIM> intersection = shifted_nbox * bregion;

                     if ( !intersection.empty() ) {

                        NeighborPatchInfo& info = 
                           patch_info.d_neighbors[blut_loc_idx][num_neighbors];
                        info.d_neighbor_patch_number = npnum;
                        info.d_neighbor_type = codim;
                        info.d_location_index = loc;
                        info.d_neighbor_shift = nshift;

                        num_neighbors++;

                        if ( intersection == ( pgbox * bregion) ) {
                           done_with_location[loc] = true; 
                           if ( check_shifts.size() == 0 ) {
                              done_with_neighbor[ni] = true;
                           }
                        }

                     }  // add neighbor info to codim neighbor array

                  }  // patch may only be a neighbor of itself when shift is non-zero

               }  // skip neighbor if done with it
            
            } // iterate over potential patch neighbors

         }  // skip location if done with it

      }  // iterate over locations for codimension

      if (!zero_shift) {
         check_shift++;
      } else {
         zero_shift = false;
      } 

   }  // iterate over valid shifts for codimension

   patch_info.d_neighbors[blut_loc_idx].resizeArray(num_neighbors);

}

/*
*************************************************************************
*                                                                       *
* Private method to check whether patch level matches level managed     *
* by patch op utilities object.                                         *
*                                                                       *
*************************************************************************
*/

template<int DIM> 
bool PatchConfigurationUtilities<DIM>::levelMatch(
   tbox::Pointer< PatchLevel<DIM> > level,
   int array_slot) 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!level.isNull());
#endif

   if ( d_patch_level_info.size() == 0 ) {
      return false;
   } else {

      const int check_ln =
         ( (array_slot < 0) ? level->getLevelNumber() : array_slot );

      if ( (check_ln >= 0) && (check_ln < d_patch_level_info.size()) ) {

         if ( !d_patch_hierarchy.isNull() ) {

            tbox::Pointer< PatchLevel<DIM> > check_level =
               d_patch_hierarchy->getPatchLevel(check_ln);

            if ( level.getPointer() != check_level.getPointer() ) {
               return false;
            }

         }

         const bool check_level_null = 
            d_patch_level_info[check_ln].d_patch_level.isNull();

         return ( check_level_null || 
                  (level.getPointer() == 
                   d_patch_level_info[check_ln].
                      d_patch_level.getPointer()) );

      } else {
         return false;
      }

   }
}

/*
*************************************************************************
*                                                                       *
* Private method to allocate PatchLevelInfo array.                      *
*                                                                       *
*************************************************************************
*/

template<int DIM> 
void PatchConfigurationUtilities<DIM>::allocatePatchLevelInfoArray(
   int nlevels) 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(d_patch_level_info.size() == 0);
   if ( d_patch_hierarchy.isNull() ) {
      TBOX_ASSERT( nlevels == 1 );
   } else {
      TBOX_ASSERT( nlevels == d_patch_hierarchy->getNumberOfLevels() ); 
   }    
#endif
   d_patch_level_info.resizeArray(nlevels);
   for (int ln = 0; ln < nlevels; ++ln) {
      d_patch_level_info[ln].d_patch_level.setNull();    
      d_patch_level_info[ln].d_patch_info.resizeArray(0);    
      d_patch_level_info[ln].d_coarsened_fine_patch_boxes.resizeBoxArray(0);    
      d_patch_level_info[ln].d_neighbor_patch_info_is_set = false;
      d_patch_level_info[ln].d_finer_level_info_is_set = false;
      d_patch_level_info[ln].d_coarser_level_info_is_set = false;
      d_patch_level_info[ln].d_is_all_set = false;
   }
}

}
}
#endif
