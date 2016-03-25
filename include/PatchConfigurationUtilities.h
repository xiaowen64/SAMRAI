//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/patches/PatchConfigurationUtilities.h $
// Package:     SAMRAI hierarchy
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Utilities class with routines for understanding spatial
//              relationships among patches
//

#ifndef included_hier_PatchConfigurationUtilities
#define included_hier_PatchConfigurationUtilities

#include "SAMRAI_config.h"

#include "PatchLevel.h"
#include "PatchHierarchy.h"

#include "ErrorCheckIntTypes.h"

#include "tbox/Pointer.h"
#include "tbox/Array.h"


#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif

namespace SAMRAI {
    namespace hier {

/*!
 * @brief Class PatchConfigurationUtilities contains routines that 
 * provide information about spatial relationships among patches 
 * on an AMR patch hierarchy.  Such information is useful for building 
 * operations on patch data that must account for these relationships.  
 * Examples include knowing which patches are neighbors of a given patch 
 * on the same patch level and the positions of those neighbors with 
 * respect to the patch, and knowing which patches on a finer patch 
 * level overlap a given patch, etc.
 *
 * @see hier::PatchLevel
 * @see hier::PatchHierarchy
 */

template<int DIM>
class PatchConfigurationUtilities 
{
public:
   /*!
    * The <TT>PatchConfigurationUtilities<DIM>::NeighborPatchInfo</TT> data 
    * structure contains items that describe the relationship between
    * a patch and one of its neigbors on the same patch level.  
    *
    * The patch number indicates the number of the neighbor patch in
    * the array of patches on the patch level.  The neighbor type and 
    * location index describe the orientation of the neighbor with 
    * respect to the patch.  These values use the conventions of the 
    * boundary type and location index designations, respectively, in 
    * the hier::BoundaryBox<DIM> class.  See that class header file for 
    * details.  The shift vector indicates how the neighbor patch box 
    * must be shifted to to be a neighbor of the patch.  This vector
    * contains non-zero entries only when the patch abuts a periodic
    * boundary.
    */
   struct NeighborPatchInfo {
      int d_neighbor_patch_number;
      int d_neighbor_type;
      int d_location_index;
      IntVector<DIM> d_neighbor_shift;
   };

   /*!
    * The constructor for PatchConfigurationUtilities initializes the 
    * utilities object to manage data operations for the given patch 
    * hierarchy.  
    *
    * Note that the constructed utilities object cannot do 
    * anything useful until hierarchy patch level information is
    * initialized via the initialize() or reset() methods.
    * 
    * @param object_name  String name identifier for utilities object. 
    *                     Cannot be empty when assertion checking is turned on.
    *                     Should be unique among names used for other utilities
    *                     objects for proper error checking and restart functionality.
    * @param hierarchy    Optional pointer to patch hierarchy containing patch levels 
    *                     for which requests for patch information will be made.  
    *                     If null (default), then this object can only be used to 
    *                     obtain information about relationships between patches on a 
    *                     single level, which is set in the call to initialize() or
    *                     reset. 
    */
   PatchConfigurationUtilities(
      const std::string& object_name,
      tbox::Pointer< PatchHierarchy<DIM> > hierarchy = 0);

   /*!
    * Virtual destructor for utilities objects.
    */
   virtual ~PatchConfigurationUtilities<DIM>();

   /*!
    * Check whether information is set for given level number. 
    * 
    * @return Boolean true if level is set; false otherwise.
    * 
    * @param level_num  Const reference to hier::LevelNumber type indicating 
    *                   level number of interest.  When assertion checking is 
    *                   active, an assertion is thrown if number does not 
    *                   correspond to level in patch hierarchy.
    */
   bool levelIsSet(const LevelNumber& level_num) const;
 
   /*!
    * Check whether information is set for given patch and level number. 
    * 
    * Note that this routine will always return false if the patch is
    * not mapped to this processor.
    * 
    * @return Boolean true if patch is is set; false otherwise.
    * 
    * @param patch_num  Const reference to hier::PatchNumber type indicating 
    *                   patch number of interest.  When assertion checking
    *                   is active, an assertion is thrown if number does not
    *                   correspond to a patch on the patch level.
    * @param level_num  Const reference to hier::LevelNumber type indicating 
    *                   level number of interest.  When assertion checking is 
    *                   active, an assertion is thrown if number does not 
    *                   correspond to level in patch hierarchy.
    */
   bool patchIsSet(const PatchNumber& patch_num, 
                   const LevelNumber& level_num) const; 

   /*!
    * Initialize utilities object so that it can service requests for 
    * patch information for the given patch level or for all patch 
    * levels in the hierarchy.
    *
    * If this method is called more than once, then patch information is built 
    * only for patch levels for which information has not been built already.
    * In particular, if the method is called with the same patch level more than
    * once, each call after the first is essentially a no-op.
    *
    * This method cannot be used to replace patch info associated with a given 
    * patch level number with that for a new level with the same level number.  
    * For that to occur, the clear() method must be called first. 
    *
    * For consistency and to prevent unexpected behavior, this routine issues 
    * warnings and will not perform operations in certain situations.  These
    * situations include: the given patch level is not in the patch hierarchy passed 
    * to the constructor when this object was created, and the given patch level 
    * does not match the level with the same level number used to build patch 
    * information in the most recent call to this initialization routine. Either 
    * of these situations indicates potentially erroneous usage such as mixing
    * levels from different patch hierarchied in a single utilities object.
    * 
    * @param level Optional pointer to patch level. If a non-null level pointer is
    *              supplied, patch information is built only for the level. 
    *              Otherwise, information is built for all hierarchy levels 
    *              for which it has not been built already.
    */
   void initialize(tbox::Pointer< PatchLevel<DIM> > level = 0);

   /*!
    * Clear patch information for the given patch level or for all patch levels 
    * in the hierarchy.  
    *
    * For consistency and to prevent unexpected behavior, this routine issues 
    * warnings and will not perform operations in certain situations.  These
    * situations include: the given patch level is not in the patch hierarchy passed 
    * to the constructor when this object was created, and the given patch level 
    * does not match the level with the same level number used to build patch 
    * information in the most recent call to the initialization routine. Either 
    * of these situations indicates potentially erroneous usage such as mixing
    * levels from different patch hierarchied in a single utilities object.
    *
    * If no patch level argument is supplied. this method returns state of object 
    * to default state at construction. In particular, only the object name and 
    * hierarchy (if passed to constructor originally) remain intact.
    * 
    * @param level Optional pointer to patch level. If a non-null level pointer is
    *              supplied, patch information is cleared only for the level.
    *              Otherwise, information is cleared for all hierarchy levels.
    */
   void clear(tbox::Pointer< PatchLevel<DIM> > level = 0);

   /*!
    * Get an array of NeighborPatchInfo structs describing the patches on
    * the same patch level that are node neighbors of the given patch.  
    * A node neighbor intersects the patch at a single corner point.
    *
    * @return Const reference to array of NeighborPatchInfo for node neighbors.
    *         If the patch has no node neighbors, the array will have size zero.
    *  
    * @param patch_num  Const reference to hier::PatchNumber type indicating patch
    *                   number of patch of interest.  When assertion checking
    *                   is active, an assertion is thrown if number does not
    *                   correspond to a patch on the patch level.
    * @param level_num  Optional const reference to hier::LevelNumber type indicating 
    *                   level number of patch of interest.  The level number must be
    *                   given when this object is used to manage information for a
    *                   patch hierarchy (as opposed to a single level).  When 
    *                   assertion checking is active, an assertion is thrown if 
    *                   number does not correspond to level in patch hierarchy.
    *                   If no patch hierarchy was passed to the constructor when 
    *                   this object was created, then the level number is ignored. 
    */
   const tbox::Array< typename PatchConfigurationUtilities<DIM>::NeighborPatchInfo >& 
   getNodeNeighborInfo(const PatchNumber& patch_num,
                       const LevelNumber& level_num = LevelNumber(-1)) const;

   /*!
    * Get an array of NeighborPatchInfo structs describing the patches on
    * the same patch level that are edge neighbors of the given patch.
    * An edge neighbor intersects the patch along a 1-dimensional edge.
    *
    * @return Const reference to array of NeighborPatchInfo for node neighbors.
    *         If the patch has no edge neighbors, or when DIM < 2, this routine 
    *         returns an array of size zero.
    *
    * @param patch_num  Const reference to hier::PatchNumber type indicating patch
    *                   number of patch of interest.  When assertion checking
    *                   is active, an assertion is thrown if number does not
    *                   correspond to a patch on the patch level.
    * @param level_num  Optional const reference to hier::LevelNumber type indicating
    *                   level number of patch of interest.  The level number must be
    *                   given when this object is used to manage information for a
    *                   patch hierarchy (as opposed to a single level).  When
    *                   assertion checking is active, an assertion is thrown if
    *                   number does not correspond to level in patch hierarchy.
    *                   If no patch hierarchy was passed to the constructor when
    *                   this object was created, then the level number is ignored.
    */
   const tbox::Array< typename PatchConfigurationUtilities<DIM>::NeighborPatchInfo >&
   getEdgeNeighborInfo(const PatchNumber& patch_num,
                       const LevelNumber& level_num = LevelNumber(-1)) const;

   /*!
    * Get an array of NeighborPatchInfo structs describing the patches on
    * the same patch level that are face neighbors of the given patch.
    * A face neighbor intersects the patch along a 2-dimensional face.
    *
    * @return Const reference to array of NeighborPatchInfo for node neighbors.
    *         If the patch has no face neighbors, or when DIM < 3, this routine 
    *         returns an array of size zero.
    *
    * @param patch_num  Const reference to hier::PatchNumber type indicating patch
    *                   number of patch of interest. When assertion checking
    *                   is active, an assertion is thrown if number does not
    *                   correspond to a patch on the patch level.
    * @param level_num  Optional const reference to hier::LevelNumber type indicating
    *                   level number of patch of interest.  The level number must be
    *                   given when this object is used to manage information for a
    *                   patch hierarchy (as opposed to a single level).  When
    *                   assertion checking is active, an assertion is thrown if
    *                   number does not correspond to level in patch hierarchy.
    *                   If no patch hierarchy was passed to the constructor when
    *                   this object was created, then the level number is ignored.
    */
   const tbox::Array< typename PatchConfigurationUtilities<DIM>::NeighborPatchInfo >&
   getFaceNeighborInfo(const PatchNumber& patch_num,
                       const LevelNumber& level_num = LevelNumber(-1)) const;

   /*!
    * Get an array of NeighborPatchInfo structs describing the patches on
    * the same patch level that are neighbors along an (DIM - codim)-dimensional 
    * region of the given patch.
    * 
    * For spatial dimensions less than 4, this routine operates similarly to the
    * functions above.  Specifically,
    *
    * if DIM == 1: (codim == 1) => same components as getNodeBoundary.
    *
    * if DIM == 2, (codim == 1) => same components as getEdgeBoundary.
    *              (codim == 2) => same components as getNodeBoundary.
    *
    * if DIM == 3, (codim == 1) => same components as getFaceBoundary.
    *              (codim == 2) => same components as getEdgeBoundary.
    *              (codim == 3) => same components as getNodeBoundary.
    *
    * @return Const reference to array of NeighborPatchInfo for node neighbors.
    *         If the patch has no neighbors corresponding to the codimension, 
    *         or when codim <= 0 or codim > DIM, the array will have size zero.
    *
    * @param codim      Integer codimension of neighbor boundary.  The value
    *                   should satisfy 0 < codim <= DIM.  If not, an array
    *                   of size zero will be returned.
    * @param patch_num  Const reference to hier::PatchNumber type indicating patch
    *                   number of patch of interest. When assertion checking
    *                   is active, an assertion is thrown if number does not
    *                   correspond to a patch on the patch level.
    * @param level_num  Optional const reference to hier::LevelNumber type indicating
    *                   level number of patch of interest.  The level number must be
    *                   given when this object is used to manage information for a
    *                   patch hierarchy (as opposed to a single level).  When
    *                   assertion checking is active, an assertion is thrown if
    *                   number does not correspond to level in patch hierarchy.
    *                   If no patch hierarchy was passed to the constructor when
    *                   this object was created, then the level number is ignored.
    */
   const tbox::Array< typename PatchConfigurationUtilities<DIM>::NeighborPatchInfo >&
   getCodimensionNeighborPatchInfo(int codim,
                                   const PatchNumber& patch_num,
                                   const LevelNumber& level_num = LevelNumber(-1)) const;

   /*!
    * Get array of integer patch numbers of patches on the next finer hierarchy 
    * level that overlap the specified patch.  
    *
    * @return Const reference to int array of patch numbers of finer level
    *         patches that overlap given patch.  When the level number 
    *         corresponds to the finest hierarchy level, or when a null 
    *         patch hierarchy pointer was passed to the constructor, the
    *         array will be size zero.
    * 
    * @param patch_num  Const reference to hier::PatchNumber type indicating patch
    *                   number of patch of interest. When assertion checking
    *                   is active, an assertion is thrown if number does not
    *                   correspond to a patch on the patch level.
    * @param level_num  Const reference to hier::LevelNumber type indicating level
    *                   number of patch of interest. When assertion checking
    *                   is active, an assertion is thrown if number does not
    *                   correspond to a level in patch hierarchy if a non-null
    *                   patch hierarchy pointer was passed to the constructor.
    */
   const tbox::Array<int>& 
   getFinerLevelOverlapPatchIndices(const PatchNumber& patch_num,
                                    const LevelNumber& level_num) const;

   /*!
    * Get box representing the coarsening of a patch box on the next finer 
    * hierarchy level to the index space of the specified level. 
    *
    * @return Const reference to finer level patch box coarsened to index
    *         space of specified next coarser level.  When the level number 
    *         corresponds to the finest hierarchy level, or when a null 
    *         patch hierarchy pointer was passed to the constructor, the
    *         box will be empty.
    *
    * @param finer_level_patch_num  Const reference to hier::PatchNumber type 
    *                   indicating patch number of patch on next finer level. 
    *                   When assertion checking is active, an assertion is thrown 
    *                   if number does not correspond to a patch on the patch level.
    * @param level_num  Const reference to hier::LevelNumber type indicating level
    *                   number of interest. When assertion checking
    *                   is active, an assertion is thrown if number does not
    *                   correspond to level in patch hierarchy.
    */
   const Box<DIM>& 
   getCoarsenedFinerLevelPatchBox(const PatchNumber& finer_level_patch_num,
                                  const LevelNumber& level_num) const;

   /*!
    * Get array of integer patch numbers of patches on the next coarser hierarchy 
    * level that overlap the specified patch.
    *
    * @return Const reference to int array of patch numbers of coarser level
    *         patches that overlap given patch.  When the level number 
    *         corresponds to the coarsest hierarchy level, or when a null 
    *         patch hierarchy pointer was passed to the constructor, the
    *         array will be size zero.
    * 
    * @param patch_num  Const reference to hier::PatchNumber type indicating patch
    *                   number of patch of interest. When assertion checking
    *                   is active, an assertion is thrown if number does not
    *                   correspond to a patch on the patch level.
    * @param level_num  Const reference to hier::LevelNumber type indicating level
    *                   number of patch of interest. When assertion checking
    *                   is active, an assertion is thrown if number does not
    *                   correspond to a level in patch hierarchy if a non-null
    *                   patch hierarchy pointer was passed to the constructor.
    */
   const tbox::Array<int>& 
   getCoarserLevelOverlapPatchIndices(const PatchNumber& patch_num,
                                      const LevelNumber& level_num) const;

   /*!
    * Get box representing the refinement of a patch box on the next coarser 
    * hierarchy level to the index space of the specified level. 
    *
    * @return Const reference to coarser level patch box refined to index
    *         space of specified next finer level.  When the level number 
    *         corresponds to the coarsest hierarchy level (i.e., level 0), 
    *         or when a null patch hierarchy pointer was passed to the 
    *         constructor, the box will be empty.
    *
    * @param coarser_level_patch_num  Const reference to hier::PatchNumber type 
    *                   indicating patch number of patch on next coarser level. 
    *                   When assertion checking is active, an assertion is thrown 
    *                   if number does not correspond to a patch on the patch level.
    * @param level_num  Const reference to hier::LevelNumber type indicating level
    *                   number of interest. When assertion checking
    *                   is active, an assertion is thrown if number does not
    *                   correspond to level in patch hierarchy.
    */
   const Box<DIM>& 
   getRefinedCoarserLevelPatchBox(const PatchNumber& coarser_level_patch_num,
                                  const LevelNumber& level_num) const;

   /*!
    * Print internal data structures for patch level to given output stream.
    *
    * @param level_num  Optional const reference to hier::LevelNumber type 
    *                   indicating number of patch level of interest. 
    *                   If not given and this utilities object is set up for 
    *                   a patch hierarchy, data for all patch levels is output.  
    *                   If a valid level number is given and this utilities 
    *                   object is set up for a patch hierarchy, only data for 
    *                   the specified level is output.  If the value is not a 
    *                   valid hierarchy level number, the routine does nothing.  
    *                   In the case that this utilities object is set up for a 
    *                   single patch level, the value is ignored and data for 
    *                   the level is output.
    * @param os  Optional output stream.  If not given, tbox::plog is used.
    */
   void printPatchLevelData(const LevelNumber& level_num = LevelNumber(-1),
                            std::ostream& os = tbox::plog) const;

   /*!
    * Print internal data structures for patch to given output stream.
    *
    * @param patch_num  Const reference to hier::PatchNumber type
    *                   indicating number of patch of interest.  When
    *                   assertion checking is active, an assertion results 
    *                   when this value is not a valid patch number for the
    *                   level of interest.
    * @param level_num  Optional const reference to hier::LevelNumber type
    *                   indicating number of patch level of interest.  When
    *                   this utilities object is set up for a single patch 
    *                   level, the value is ignored and level passed in to 
    *                   the constructor of this object is assumed.
    *                   If this object is set up for a patch hierarchy, this
    *                   value must be given and must be a valid level number for
    *                   the patch hierarchy.  If not, an assertion results, when
    *                   assertion checking is active.
    * @param os  Optional output stream.  If not given, tbox::plog is used.
    */
   void printPatchData(const PatchNumber& patch_num,
                       const LevelNumber& level_num = LevelNumber(-1),
                       std::ostream& os = tbox::plog) const;

private:
   // The following three methods are not implemented; they are declared
   // private here to prevent the compiler from creating them implicitly.
   PatchConfigurationUtilities();
   PatchConfigurationUtilities(const PatchConfigurationUtilities<DIM>&);
   void operator=(const PatchConfigurationUtilities<DIM>&);

   /*
    * Struct used to maintain information about neighbor patches
    * on the same level and indices overlapping patches on the next
    * finer level for an individual patch.
    */
   struct PatchInfo {
#ifdef LACKS_NAMESPACE_IN_DECLARE
      tbox::Array< NeighborPatchInfo > d_neighbors[DIM];
#else
      tbox::Array< typename PatchConfigurationUtilities<DIM>::NeighborPatchInfo > 
         d_neighbors[DIM];
#endif
      tbox::Array<int> d_finer_level_overlap_patch_indices;
      tbox::Array<int> d_coarser_level_overlap_patch_indices;
   };

   /*
    * Struct used to maintain information about an individual patch level.
    */
   struct PatchLevelInfo {
      tbox::Pointer< PatchLevel<DIM> > d_patch_level;
#ifdef LACKS_NAMESPACE_IN_DECLARE
      tbox::Array< PatchInfo* > d_patch_info;
#else
      tbox::Array< typename PatchConfigurationUtilities<DIM>::PatchInfo* > 
         d_patch_info;
#endif
      BoxArray<DIM> d_coarsened_fine_patch_boxes;
      BoxArray<DIM> d_refined_coarse_patch_boxes;
      bool d_neighbor_patch_info_is_set;
      bool d_finer_level_info_is_set;
      bool d_coarser_level_info_is_set;
      bool d_is_all_set;
   };

   std::string d_object_name;

#ifdef LACKS_NAMESPACE_IN_DECLARE
   tbox::Array< PatchLevelInfo > d_patch_level_info;
#else
   tbox::Array< typename PatchConfigurationUtilities<DIM>::PatchLevelInfo > 
      d_patch_level_info;
#endif

   tbox::Pointer< PatchHierarchy<DIM> > d_patch_hierarchy;

   /*
    * Dummy objects for returning empty info from member functions.
    */
   tbox::Array<int> d_empty_array;
   Box<DIM> d_empty_box;

   /*
    * Private method to set PatchLevelInfo array entry for given
    * patch level.  Method assumes patch level pointer is non-null
    * and that PatchLevelInfo array is large enough for the level entry.
    */
   void setPatchLevelInfo(tbox::Pointer< PatchLevel<DIM> > level, 
                          int array_slot = -1);

   /*
    * Private method to set PatchInfo array for given PatchLevelInfo
    * array slot (ln_slot).  Method assumes patch level pointer 
    * is set in corresponding PatchLevelInfo object.
    */
   void setNeighborPatchInfo(int array_slot);

   /*
    * Private method to set information about configuration of finer level 
    * patches with respect to patches on level with given PatchLevelInfo array 
    * slot.  Method assumes patch level pointer is set in corresponding 
    * PatchLevelInfo object.
    */
   void setFinerPatchLevelInfo(int array_slot);

   /*
    * Private method to set information about configuration of coarser level 
    * patches with respect to patches on level with given PatchLevelInfo array 
    * slot.  Method assumes patch level pointer is set in corresponding 
    * PatchLevelInfo object.
    */
   void setCoarserPatchLevelInfo(int array_slot);

   /*
    * Private method to clear PatchLevelInfo array entry for given
    * array slot.  Method assumes that array slot value is valid.
    */
   void clearPatchLevelInfo(int array_slot);

   /* 
    * Private method to set neighbor patch info for given codimension.
    */
   void findPatchNeighbors(PatchInfo& patch_info,
                           tbox::Pointer< Patch<DIM> > patch,
                           tbox::Pointer< PatchLevel<DIM> > level,
                           const tbox::Array<int>& neighbor_patch_ids,
                           tbox::Array<bool>& done_with_neighbor,
                           int codim);
   
   /*
    * Private method to check whether a different level is assigned
    * to a given slot in the PatchLevelInfo array.  Method assumes 
    * patch level pointer is non-null.
    */
   bool levelMatch(tbox::Pointer< PatchLevel<DIM> > level, 
                   int array_slot = -1);

   /*
    * Private method to allocate array of PatchLevelInfo structures and
    * initialize each to default state.
    */
   void allocatePatchLevelInfoArray(int nlevels);

   /*
    * Private method to print patch info for level.
    */
   void printPatchLevelInfo(std::ostream& os, int ln) const;

};

}
}

#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "PatchConfigurationUtilities.C"
#endif
