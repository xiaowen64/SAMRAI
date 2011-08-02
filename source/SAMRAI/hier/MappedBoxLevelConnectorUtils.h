/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Utilities for working on DLBG edges. 
 *
 ************************************************************************/
#ifndef included_hier_MappedBoxLevelConnectorUtils
#define included_hier_MappedBoxLevelConnectorUtils

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Connector.h"
#include "SAMRAI/hier/MappedBoxLevel.h"
#include "SAMRAI/hier/MultiblockMappedBoxTree.h"

namespace SAMRAI {
namespace hier {

/*!
 * @brief Utilities for common operating on MappedBoxLevels.
 *
 * Objects of this class can be set to perform certain sanity checks
 * on the pre and post conditions of the methods.  See
 * setSanityCheckMethodPreconditions() and
 * setSanityCheckMethodPostconditions().
 */
class MappedBoxLevelConnectorUtils
{

public:
   /*!
    * @brief Default constructor.
    *
    * By default, sanity checks are disabled.  To enable them, see
    * setSanityCheckMethodPreconditions() and
    * setSanityCheckMethodPostconditions().
    */
   MappedBoxLevelConnectorUtils();

   /*!
    * @brief Set whether to run expensive sanity checks on input parameters.
    *
    * Mainly for debugging.
    *
    * @param[in] do_check
    */
   void
   setSanityCheckMethodPreconditions(
      bool do_check);

   /*!
    * @brief Set whether to run expensive sanity checks on output parameters.
    *
    * Mainly for debugging.
    *
    * @param[in] do_check
    */
   void
   setSanityCheckMethodPostconditions(
      bool do_check);


   //@{

   //! @name Comparing boxes of two MappedBoxLevels

   /*!
    * @brief Given an overlap Connector, determine the extent to which
    * the Connector's base nests in its head.
    *
    * This method returns true if the base, grown by @c base_swell,
    * nests inside the head, grown by @c head_swell, by a margin of @c
    * head_nesting_margin.  @c base_swell and @c head_swell should be
    * non-negative and specified in the base and head index spaces,
    * respectively.  @c head_nesting_margin should be in the head
    * index space.
    *
    * The Connector width must be at least the sum of the @c
    * base_swell and the (appropriately converted) @c head_swell and
    * @c head_nesting_margin.  We require and assume without verifying
    * that the Connector is complete.
    *
    * If the domain is given, non-nesting parts outside of the domain
    * are disregarded.
    *
    * @param[out] locally_nests Whether the local parts of the base
    * nests in the head.  This output may vary among the processes.
    *
    * @param[in] connector
    *
    * @param[in] base_swell the amount that the base is grown by, given in the
    * base index space and non-negative
    *
    * @param[in] head_swell the amount that the head is grown by, given in the
    * head index space and non-negative
    *
    * @param[in] head_nesting_margin given in the head index space.
    *
    * @param[in] domain Domain description, in reference index space,
    * in search tree format.
    *
    * @return True if the given base MappedBoxLevel nests in the head,
    * otherwise False.
    */
   bool
   baseNestsInHead(
      bool* locally_nests,
      const Connector& connector,
      const IntVector& base_swell,
      const IntVector& head_swell,
      const IntVector& head_nesting_margin,
      const MappedBoxTree* domain = NULL) const;

   /*!
    * @brief Multiblock version of baseNestsInHead.
    */
   bool
   baseNestsInHeadForMultiblock(
      bool* locally_nests,
      const Connector& connector,
      const IntVector& base_swell,
      const IntVector& head_swell,
      const IntVector& head_nesting_margin,
      const MultiblockMappedBoxTree* domain = NULL) const;

   /*!
    * @brief Given base and head MappedBoxLevels, determine the extent
    * to which the base nests in the head.
    *
    * This method is similar to the version taking a Connector instead
    * of the base and head MappedBoxLevels, except that it will build
    * the overlap Connector from the base to the head.  It should be
    * used when you don't have a pre-built Connector or don't have one
    * with sufficient width.
    *
    * This method simply builds a Connector from head to base with
    * sufficient width and calls baseNestsInHead().  Building of the
    * Connector is not scalable.  If you already have such a
    * Connector, you should call baseNestsInHead(bool*, const Connector&,
    * const IntVector&, const IntVector&, const IntVector&,
    * const MappedBoxTree*) directly.
    *
    * TODO: Is this method really needed now that we can automatically
    * build the base--->head Connector automatically using
    * PersistentOverlapConnectors?  Using PersistentOverlapConnector
    * does have the side-effect that the generated Connector persists
    * longer than the immediate need.
    *
    * @param domain Domain description, in reference index space, in
    * search tree format.
    *
    * @return Whether the given base MappedBoxLevel nests in the head.
    *
    * @param[out] locally_nests Whether the local parts of the base
    * nests in the head.  This output may vary among the processes.
    *
    * @param[in] base
    *
    * @param[in] head
    *
    * @param[in] base_swell the amount that the base is grown by, given in the
    * base index space and non-negative
    *
    * @param[in] head_swell the amount that the head is grown by, given in the
    * head index space and non-negative
    *
    * @param[in] head_nesting_margin given in the head index space.
    *
    * @param[in] domain Domain description, in reference index space,
    * in search tree format.
    *
    * @return Whether the given base MappedBoxLevel nests in the head.
    */
   bool
   baseNestsInHead(
      bool* locally_nests,
      const MappedBoxLevel& base,
      const MappedBoxLevel& head,
      const IntVector& base_swell,
      const IntVector& head_swell,
      const IntVector& head_margin,
      const MappedBoxTree* domain = NULL) const;
   bool

   /*!
    * @brief Multiblock version of baseNestsInHead.
    */
   baseNestsInHeadForMultiblock(
      bool* locally_nests,
      const MappedBoxLevel& base,
      const MappedBoxLevel& head,
      const IntVector& base_swell,
      const IntVector& head_swell,
      const IntVector& head_margin,
      const MultiblockMappedBoxTree* domain = NULL) const;

   /*!
    * @brief Compute the parts of one MappedBoxLevel that are external
    * to another MappedBoxLevel.
    *
    * This is the singleblock version of computeExternalPartsForMultiblock.
    *
    * @see computeExternalPartsForMultiblock
    */
   void
   computeExternalParts(
      MappedBoxLevel& external,
      Connector& input_to_external,
      const Connector& input_to_reference,
      const IntVector& nesting_width,
      const MappedBoxTree& domain) const;

   /*!
    * @brief Compute the parts of one MappedBoxLevel that are internal
    * to another MappedBoxLevel.
    *
    * This is the singleblock version of computeInternalPartsForMultiblock.
    *
    * @see computeInternalPartsForMultiblock
    */
   void
   computeInternalParts(
      MappedBoxLevel& internal,
      Connector& input_to_internal,
      const Connector& input_to_reference,
      const IntVector& nesting_width,
      const MappedBoxTree& domain) const;

   /*!
    * @brief Compute the parts of one MappedBoxLevel that are external
    * to another MappedBoxLevel.
    *
    * Compare an input MappedBoxLevel to a "reference" MappedBoxLevel.
    * Compute the parts of the input that are external to the
    * reference.  Build the "external" MappedBoxLevel representing the
    * external parts.  Build a mapping Connector with the input as its
    * base and the external as its head.
    *
    * A partially external input cell (possible when input is coarser
    * than reference) is considered to be external.
    *
    * For the purpose of defining what is external, the reference
    * level can be grown by nesting_width before comparing.  This
    * feature can be used to determing which parts of the input does
    * not nest in the reference comparison.  A negative growth
    * indicates shrinking of the reference level.
    *
    * This method does not require any communication.
    *
    * @param[out] external.  The existing state will be discarded.
    *
    * @param[out] input_to_external.  The existing state will be
    * discarded.
    *
    * @param[in] input_to_reference Overlap Connector from input to
    * reference MappedBoxLevel.
    *
    * @param[in] nesting_width Growth of the reference MappedBoxLevel
    * for the purpose of comparing to input.  Must be in resolution of
    * input MappedBoxLevel.  Must be either non-negative or
    * non-positive but not mixed.
    *
    * @param[in] domain The domain representation, without periodic
    * images, in search tree form.  These boxes should be in the
    * reference index space.  If domain is given, do not shrink the
    * reference MappedBoxLevel where it touches the domain boundary.
    */
   void
   computeExternalPartsForMultiblock(
      MappedBoxLevel& external,
      Connector& input_to_external,
      const Connector& input_to_reference,
      const IntVector& nesting_width,
      const MultiblockMappedBoxTree& domain = MultiblockMappedBoxTree() ) const;

   /*!
    * @brief Compute the parts of one MappedBoxLevel that are internal
    * to another MappedBoxLevel.
    *
    * Compare an input MappedBoxLevel to a "reference" mapped_box_level.
    * Identify parts of the input that are internal to the reference
    * MappedBoxLevel, and store the internal parts in a
    * MappedBoxLevel.  Set up a mapping Connector between the input
    * and its internal parts.
    *
    * A partially internal input cell (possible when input is coarser
    * than reference) is considered to be internal.
    *
    * For the purpose of defining what is external, the reference
    * level can be grown by nesting_width before comparing.  This
    * feature can be used to determing which parts of the input does
    * not nest in the reference comparison.  A negative growth
    * indicates shrinking of the reference level.
    *
    * This method does not require any communication.
    *
    * @param[out] internal.  The existing state will be discarded.
    *
    * @param[out] input_to_internal.  The existing state will be
    * discarded.
    *
    * @param[in] input_to_reference Overlap Connector from input to
    * reference MappedBoxLevel.
    *
    * @param[in] nesting_width Growth of the reference MappedBoxLevel
    * for the purpose of comparing to input.  Must be in resolution of
    * input MappedBoxLevel.  Must be either non-negative or
    * non-positive but not mixed.
    *
    * @param[in] domain The domain representation, without periodic
    * images, in search tree form.  These boxes should be in the
    * reference index space.  If domain is given, do not shrink the
    * reference MappedBoxLevel where it touches the domain boundary.
    */
   void
   computeInternalPartsForMultiblock(
      MappedBoxLevel& internal,
      Connector& input_to_internal,
      const Connector& input_to_reference,
      const IntVector& nesting_width,
      const MultiblockMappedBoxTree& domain = MultiblockMappedBoxTree() ) const;

   /*!
    * @brief Compute the parts of one MappedBoxLevel that is internal
    * to another MappedBoxLevel.
    *
    * This version of computeInternalParts() does not require a domain.
    * The domain is taken to be big enough that it does not affect
    * the definition of "internal".
    * @see computeInternalParts(MappedBoxLevel&,Connector&,const Connector&,const IntVector&,const MappedBoxTree&)const
    *
    * @param[out] internal Any internal part is owned by the process
    * owning the input MappedBox that generated it.
    *
    * @param[out] input_to_internal Relationships from input
    * MapppedBoxes to their internal parts.  This is a local map.
    *
    * @param[in] input_to_reference Overlap Connector from input to
    * reference MappedBoxLevel.
    *
    * @param[in] nesting_width Growth of the reference MappedBoxLevel
    * for the purpose of comparing to input.  Must be in coordinate
    * system of input.
    */
   void
   computeInternalParts(
      MappedBoxLevel& internal,
      Connector& input_to_internal,
      const Connector& input_to_reference,
      const IntVector& nesting_width) const;

   //@}


   /*!
    * @brief Given a set of MappedBoxes, compute its boundary as a set
    * of boxes located just outside it.
    *
    * @param boundary_boxes[o] Boundary boxes, sorted into BoxLists
    * according to the BlockId.
    *
    * @param mapped_boxes[i] MappedBoxes to find the boundary for.
    *
    * @param refinement_ratio[i] Refinement ratio of mapped_boxes.
    *
    * @param grid_geometry[i]
    *
    * @param simplify_boundary_boxes Whether to simplify the boundary
    * boxes after computing them.
    */
   void computeBoxesAroundBoundary(
      std::map<BlockId,BoxList> &boundary,
      const IntVector &refinement_ratio,
      const tbox::ConstPointer<GridGeometry> &grid_geometry,
      const bool simplify_boundary_boxes = true ) const;


   //@{

   //! @name Setting up common mapping Connectors

   /*
    * @brief Sort the MappedBoxes in MappedBoxLevel and make a mapping
    * Connector from the unsorted MappedBoxLevel to the sorted one.
    * The sorting can renumber the LocalIndices of the MappedBoxes
    * or put the MappedBoxes in spatial ordering, or both.
    *
    * The Connector map created is local (no MappedBox is mapped to a new
    * owner).
    *
    * If @c sort_mapped_boxes_by_corner is true, the map will reorder
    * local MappedBoxes by their box corners.  This is useful for
    * making random box ordering into something deterministic.
    *
    * If @c sequentialize_global_indices is true, determine the lowest
    * index the local processor should use for the output
    * MappedBoxLevel so that the global set of MappedBoxes have
    * sequential indices.  (This requires communication.)  If false,
    * each processor start numbering with LocalIndex zero.  If true,
    * @c initial_sequential_index can specify the first index of first
    * MappedBox of the lowest rank processor.
    *
    * For more information on mapping Connectors, see
    * MappingConnectorAlgorithm.
    *
    * @param[out] sorted_mapped_box_level Sorted version of the input
    * unsorted_mapped_box_level.
    *
    * @param[out] output_map Mapping from @c unsorted_mapped_box_level
    * to @c sorted_mapped_box_level.
    *
    * @param[in] unsorted_mapped_box_level
    *
    * @param[in] sort_mapped_boxes_by_corner Whether to sort local
    * MappedBoxes by their indices to make their ordering solely a
    * function of box positions.
    *
    * @param[in] sequentialize_global_indices Whether to renumber the
    * LocalIndices into a globally sequential numbering.
    *
    * @param[in] initial_sequential_index The first index of first
    * MappedBox of the lowest rank process.  This parameter is
    * disregarded when not globally sequentializing the indices.
    */
   void
   makeSortingMap(
      MappedBoxLevel& sorted_mapped_box_level,
      Connector& output_map,
      const MappedBoxLevel& unsorted_mapped_box_level,
      bool sort_mapped_boxes_by_corner = true,
      bool sequentialize_global_indices = true,
      LocalId initial_sequential_index = LocalId::getZero()) const;

   /*
    * @brief Given a mapping from an original MappedBoxLevel to parts
    * to be removed (rejected), construct the remainder MappedBoxLevel
    * and the mapping from the original to a remainder.
    *
    * @see MappingConnectorAlgorithm.
    *
    * @param[out] remainder The new MappedBoxLevel resulting from
    * removing the rejected parts from the original MappedBoxLevel.
    *
    * @param[out] orig_to_remainder The output mapping.  This is a
    * local map.
    *
    * @param[in] orig_to_rejections Mapping from original
    * MappedBoxLevel to its parts that should be be removed.  This
    * must be a local map.
    */
   void
   makeRemainderMap(
      MappedBoxLevel& remainder,
      Connector& orig_to_remainder,
      const Connector& orig_to_rejections) const;

   //@}


   //@{

   //! @name Adding periodic images

   /*!
    * @brief Add periodic images to a MappedBoxLevel.
    *
    * This method is a no-op in the case of non-periodic domains.
    *
    * Any periodic image within a certain distance of the domain is
    * added, Those farther out are not added.  The threshold distance
    * is @c threshold_distance.
    *
    * @param[in,out] mapped_box_level MappedBoxLevel subject to the
    * addition of periodic MappedBoxes.
    *
    * @param[in] domain_search_tree Domain description in the reference
    * index space.  This tree must NOT include periodic images.
    *
    * @param[in] threshold_distance
    */
   void
   addPeriodicImages(
      MappedBoxLevel& mapped_box_level,
      const MappedBoxTree& domain_search_tree,
      const IntVector &threshold_distance) const;

   /*!
    * @brief Add periodic images to a MappedBoxLevel and add new
    * relationships to the periodic images.
    *
    * This method is a no-op in the case of non-periodic domains.
    *
    * Any periodic images within a certain distance of the domain is
    * added, but the rest are not added.  The threshold distance is
    * the width of the Connector @c mapped_box_level_to_anchor.
    *
    * This method updates the overlap Connectors between the
    * MappedBoxLevel getting new periodic MappedBoxes and an "anchor"
    * MappedBoxLevel.  New periodic overlap relationships generated
    * are added to the overlap Connector @c
    * mapped_box_level_to_anchor.  If you don't need to have a
    * Connector updated, use addPeriodicImages() instead of this
    * method.
    *
    * @param[in,out] mapped_box_level MappedBoxLevel subject to the
    * addition of periodic MappedBoxes.
    *
    * @param[in,out] mapped_box_level_to_anchor Overlap Connector to
    * be updated with new relationships.
    *
    * @param[in,out] anchor_to_mapped_box_level Overlap Connector to
    * be updated with new relationships.
    *
    * @param[in] domain_search_tree Domain description in the
    * reference index space.  This tree must NOT include periodic
    * images.
    *
    * @param[in] anchor_to_anchor Self overlap Connector for anchor
    * MappedBoxLevel.  Must be a complete overlap Connector with
    * periodic relationships.
    */
   void
   addPeriodicImagesAndRelationships(
      MappedBoxLevel& mapped_box_level,
      Connector& mapped_box_level_to_anchor,
      Connector& anchor_to_mapped_box_level,
      const MappedBoxTree& domain_search_tree,
      const Connector& anchor_to_anchor) const;

   //@}



private:

   /*!
    * @brief Delegated work of computeInternalPartsForMultiblock and
    * computeExternalPartsForMultiblock.
    */
   void computeInternalOrExternalPartsForMultiblock(
      hier::MappedBoxLevel& parts,
      hier::Connector& input_to_parts,
      char internal_or_external,
      const hier::Connector& input_to_reference,
      const hier::IntVector& nesting_width,
      const hier::MultiblockMappedBoxTree& domain) const;

   /*!
    * @brief Call-back function to sort boxes.
    */
   static int
   qsortBoxCompare(
      const void* v,
      const void* w);

   /*!
    * @brief Allocate statics
    *
    * Only called by StartupShutdownManager.
    */
   static void
   initializeCallback();

   /*!
    * @brief Delete statics.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   finalizeCallback();

   static tbox::Pointer<tbox::Timer> t_make_sorting_map;
   static tbox::Pointer<tbox::Timer> t_compute_external_parts;
   static tbox::Pointer<tbox::Timer> t_compute_external_parts_intersection;
   static tbox::Pointer<tbox::Timer> t_compute_internal_parts;
   static tbox::Pointer<tbox::Timer> t_compute_internal_parts_intersection;

   bool d_sanity_check_precond;
   bool d_sanity_check_postcond;

   static tbox::StartupShutdownManager::Handler
   s_initialize_finalize_handler;

};

}
}

#ifdef SAMRAI_INLINE
#include "SAMRAI/hier/MappedBoxLevelConnectorUtils.I"
#endif

#endif  // included_hier_MappedBoxLevelConnectorUtils
