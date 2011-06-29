/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Utilities for working on DLBG edges. 
 *
 ************************************************************************/
#ifndef included_hier_MappedBoxLevelConnectorUtils_C
#define included_hier_MappedBoxLevelConnectorUtils_C

#include "SAMRAI/hier/MappedBoxLevelConnectorUtils.h"
#include "SAMRAI/hier/MappedBoxContainerUtils.h"
#include "SAMRAI/hier/MappedBoxSetSingleBlockIterator.h"
#include "SAMRAI/hier/MappingConnectorAlgorithm.h"
#include "SAMRAI/hier/OverlapConnectorAlgorithm.h"
#include "SAMRAI/hier/PeriodicShiftCatalog.h"
#include "SAMRAI/hier/RealMappedBoxConstIterator.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"
#include "SAMRAI/tbox/TimerManager.h"

#include <limits>
#include <cstdlib>
#include <list>

#ifndef SAMRAI_INLINE
#include "SAMRAI/hier/MappedBoxLevelConnectorUtils.I"
#endif

namespace SAMRAI {
namespace hier {

tbox::Pointer<tbox::Timer> MappedBoxLevelConnectorUtils::t_make_sorting_map;
tbox::Pointer<tbox::Timer> MappedBoxLevelConnectorUtils::t_compute_external_parts;
tbox::Pointer<tbox::Timer> MappedBoxLevelConnectorUtils::t_compute_external_parts_intersection;
tbox::Pointer<tbox::Timer> MappedBoxLevelConnectorUtils::t_compute_internal_parts;
tbox::Pointer<tbox::Timer> MappedBoxLevelConnectorUtils::t_compute_internal_parts_intersection;

tbox::StartupShutdownManager::Handler
MappedBoxLevelConnectorUtils::s_initialize_finalize_handler(
   MappedBoxLevelConnectorUtils::initializeCallback,
   0,
   0,
   MappedBoxLevelConnectorUtils::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);

/*
 ***********************************************************************
 ***********************************************************************
 */
MappedBoxLevelConnectorUtils::MappedBoxLevelConnectorUtils():
   d_sanity_check_precond(false),
   d_sanity_check_postcond(false)
{
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void MappedBoxLevelConnectorUtils::setSanityCheckMethodPreconditions(
   bool do_check)
{
   d_sanity_check_precond = do_check;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void MappedBoxLevelConnectorUtils::setSanityCheckMethodPostconditions(
   bool do_check)
{
   d_sanity_check_postcond = do_check;
}





/*
 ***********************************************************************
 * Given the base and head levels, determine whether the base nests
 * nests in the head.
 ***********************************************************************
 */
bool MappedBoxLevelConnectorUtils::baseNestsInHead(
   bool* locally_nests,
   const MappedBoxLevel& base,
   const MappedBoxLevel& head,
   const IntVector& base_swell,
   const IntVector& head_swell,
   const IntVector& head_nesting_margin,
   const hier::MappedBoxTree* domain) const
{

   tbox::Dimension dim(head.getDim());

   TBOX_DIM_ASSERT_CHECK_ARGS2(head_nesting_margin, base_swell);

   TBOX_ASSERT(head.getMPI() == base.getMPI());

#ifdef DEBUG_CHECK_ASSERTIONS
   const IntVector &zero_vector(IntVector::getZero(dim));
   TBOX_ASSERT(base_swell >= zero_vector);
   TBOX_ASSERT(head_swell >= zero_vector);
   TBOX_ASSERT(head_nesting_margin >= zero_vector);
#endif

   hier::IntVector required_gcw = base_swell;
   if (head.getRefinementRatio() <= base.getRefinementRatio()) {
      const IntVector ratio = base.getRefinementRatio()
         / head.getRefinementRatio();
      required_gcw += (head_swell + head_nesting_margin) * ratio;
   } else if (head.getRefinementRatio() >= base.getRefinementRatio()) {
      const IntVector ratio = head.getRefinementRatio()
         / base.getRefinementRatio();
      required_gcw += IntVector::ceiling((head_swell + head_nesting_margin), ratio);
   } else {
      TBOX_ERROR("MappedBoxLevelConnectorUtils::baseNestsInHead: head index space\n"
         << "must be either a refinement or a coarsening of\n"
         << "base, but not both.");
   }

   Connector base_to_head(
      base,
      head,
      required_gcw);

   OverlapConnectorAlgorithm oca;
   oca.findOverlaps(base_to_head);

   bool rval = baseNestsInHead(
         locally_nests,
         base_to_head,
         base_swell,
         head_swell,
         head_nesting_margin,
         domain);

   return rval;
}


bool MappedBoxLevelConnectorUtils::baseNestsInHeadForMultiblock(
   bool* locally_nests,
   const MappedBoxLevel& base,
   const MappedBoxLevel& head,
   const IntVector& base_swell,
   const IntVector& head_swell,
   const IntVector& head_nesting_margin,
   const hier::MultiblockMappedBoxTree* domain) const
{

   tbox::Dimension dim(head.getDim());

   TBOX_DIM_ASSERT_CHECK_ARGS2(head_nesting_margin, base_swell);

   TBOX_ASSERT(head.getMPI() == base.getMPI());

#ifdef DEBUG_CHECK_ASSERTIONS
   const IntVector &zero_vector(IntVector::getZero(dim));
   TBOX_ASSERT(base_swell >= zero_vector);
   TBOX_ASSERT(head_swell >= zero_vector);
   TBOX_ASSERT(head_nesting_margin >= zero_vector);
#endif

   hier::IntVector required_gcw = base_swell;
   if (head.getRefinementRatio() <= base.getRefinementRatio()) {
      const IntVector ratio = base.getRefinementRatio()
         / head.getRefinementRatio();
      required_gcw += (head_swell + head_nesting_margin) * ratio;
   } else if (head.getRefinementRatio() >= base.getRefinementRatio()) {
      const IntVector ratio = head.getRefinementRatio()
         / base.getRefinementRatio();
      required_gcw += IntVector::ceiling((head_swell + head_nesting_margin), ratio);
   } else {
      TBOX_ERROR("MappedBoxLevelConnectorUtils::baseNestsInHead: head index space\n"
         << "must be either a refinement or a coarsening of\n"
         << "base, but not both.");
   }

   Connector base_to_head(
      base,
      head,
      required_gcw);

   OverlapConnectorAlgorithm oca;
   oca.findOverlaps(base_to_head);

   bool rval = baseNestsInHeadForMultiblock(
         locally_nests,
         base_to_head,
         base_swell,
         head_swell,
         head_nesting_margin,
         domain);

   return rval;
}

/*
 ***********************************************************************
 * Given a Connector, determine the extent to which the base nests in the
 * head.  The Connector is assumed, without verification, to be complete.
 *
 * This method returns true if the base, grown by (non-negative)
 * base_swell nests inside the head by a margin of head_nesting_margin.
 * base_swell should be in the base index space and head_nesting_margin
 * should be in the head index space.  The Connector GCW must be at least
 * the sum of the base_swell and the appropriately converted
 * head_nesting_margin.  base_swell and head_nesting_margin can be
 * interchangable if you do the index space conversion yourself.
 *
 * If the domain is given, disregard non-nesting parts that are outside
 * the domain.
 ***********************************************************************
 */
bool MappedBoxLevelConnectorUtils::baseNestsInHead(
   bool* locally_nests,
   const Connector& connector,
   const IntVector& base_swell,
   const IntVector& head_swell,
   const IntVector& head_nesting_margin,
   const hier::MappedBoxTree* domain) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(
      connector.getBase(), base_swell, head_nesting_margin);
   TBOX_ASSERT(connector.isInitialized());
   const tbox::Dimension& dim(connector.getBase().getDim());
   TBOX_ASSERT(base_swell >= IntVector::getZero(dim));
   TBOX_ASSERT(head_nesting_margin >= IntVector::getZero(dim));

   /*
    * To ensure correct results, connector must be sufficiently wide.
    * It should be at least as wide as the combination of base_swell,
    * head_swell and head_nesting_margin.
    */
   const hier::IntVector required_gcw =
      base_swell
      + (connector.getHeadCoarserFlag() ?
         head_swell * connector.getRatio() :
         IntVector::ceiling(head_swell, connector.getRatio()))
      + (connector.getHeadCoarserFlag() ?
         head_nesting_margin * connector.getRatio() :
         IntVector::ceiling(head_nesting_margin, connector.getRatio()))
      ;
   if (!(connector.getConnectorWidth() >= required_gcw)) {
      TBOX_ERROR("MappedBoxLevelConnectorUtils::baseNestsInHead: connector lacks sufficient\n"
         << "ghost cell width for determining whether its base nests\n"
         << "inside its head.");
   }

   hier::OverlapConnectorAlgorithm oca;

   const MappedBoxLevel& base = connector.getBase();
   const MappedBoxLevel& head = connector.getHead();
   const tbox::ConstPointer<GridGeometry>& grid_geom = base.getGridGeometry();

   tbox::Pointer<hier::MappedBoxTree> refined_domain;
   if (domain != NULL) {
      refined_domain = domain->createRefinedTree(head.getRefinementRatio());
   }

   /*
    * We swell the base then check for the parts outside the head if
    * the head domain is grown by head_swell then shrunken by
    * head_nesting_margin.
    *
    * TODO: We can probably remove the base swelling step by converting
    * the base_swell into head index space and add it to
    * head_nesting_margin.
    */

   MappedBoxLevel swelledbase(dim);
   if (base_swell == IntVector::getZero(dim)) {
      swelledbase = base;
   } else {
      const MappedBoxSet& base_mapped_boxes = base.getMappedBoxes();
      MappedBoxSet swelledbase_mapped_boxes;
      for (MappedBoxSet::const_iterator ni = base_mapped_boxes.begin();
           ni != base_mapped_boxes.end(); ++ni) {
         MappedBox swelledbase_mapped_box = *ni;
         swelledbase_mapped_box.getBox().grow(base_swell);
         swelledbase_mapped_boxes.insert(
            swelledbase_mapped_boxes.end(), swelledbase_mapped_box);
      }
      swelledbase.swapInitialize(swelledbase_mapped_boxes,
         base.getRefinementRatio(),
         grid_geom,
         base.getMPI());
   }

   MappedBoxLevel swelledhead(dim);
   NeighborhoodSet swelledbase_eto_swelledhead;
   if (head_swell == IntVector::getZero(dim)) {
      swelledhead = head;
      swelledbase_eto_swelledhead = connector.getNeighborhoodSets();
   } else {
      const MappedBoxSet& head_mapped_boxes = head.getMappedBoxes();

      MappedBoxSet swelledhead_mapped_boxes;

      for (MappedBoxSet::const_iterator ni = head_mapped_boxes.begin();
           ni != head_mapped_boxes.end(); ++ni) {
         MappedBox swelledhead_mapped_box = *ni;

         swelledhead_mapped_box.getBox().grow(head_swell);
         swelledhead_mapped_boxes.insert(
            swelledhead_mapped_boxes.end(), swelledhead_mapped_box);
      }
      swelledhead.swapInitialize(swelledhead_mapped_boxes,
         head.getRefinementRatio(),
         grid_geom,
         head.getMPI());

      connector.getNeighborhoodSets().growNeighbors(
         swelledbase_eto_swelledhead,
         head_swell);
   }


   Connector swelledbase_to_swelledhead;
   swelledbase_to_swelledhead.swapInitialize(
      swelledbase,
      swelledhead,
      connector.getConnectorWidth() - base_swell,
      swelledbase_eto_swelledhead,
      MappedBoxLevel::DISTRIBUTED);
   if ( d_sanity_check_precond &&
        head_swell == hier::IntVector::getZero(dim) ) {
      /*
       * If head was swelled, it may generate undetected overlaps that
       * cannot be compensated for by shrinking the connector width.
       * The additional overlaps do not matter to the nesting check,
       * so it does not affect our result.  Nevertheless, because they
       * are not detected, don't make this check if head was swelled.
       */
      oca.assertOverlapCorrectness(swelledbase_to_swelledhead);
   }

   MappedBoxLevel external(dim);
   Connector swelledbase_to_external;
   if (domain) {
      computeExternalParts(
         external,
         swelledbase_to_external,
         swelledbase_to_swelledhead,
         -head_nesting_margin,
         *domain);
   } else {
      computeExternalParts(
         external,
         swelledbase_to_external,
         swelledbase_to_swelledhead,
         -head_nesting_margin,
         MappedBoxTree(dim) );
   }
   if (domain) {
      /*
       * If domain is given, do not count external parts that are
       * outside the domain.  In many usages, part of base is outside
       * the domain and we want to ignore those parts.
       */
      hier::MappingConnectorAlgorithm mca;
      std::vector<MappedBox> domain_mapped_boxes;
      domain->getMappedBoxes(domain_mapped_boxes);
      MappedBoxLevel domain_mapped_box_level(
         IntVector::getOne(dim),
         grid_geom,
         connector.getMPI(),
         MappedBoxLevel::GLOBALIZED);
      for (size_t i = 0; i < domain_mapped_boxes.size(); ++i) {
         domain_mapped_box_level.addMappedBox(domain_mapped_boxes[i]);
      }
      Connector external_to_domain(
         external,
         domain_mapped_box_level,
         base_swell);
      oca.findOverlaps(external_to_domain);
      MappedBoxLevel finalexternal(dim);
      Connector external_to_finalexternal;
      computeInternalParts(
         finalexternal,
         external_to_finalexternal,
         external_to_domain,
         IntVector::getZero(dim),
         *domain);
      mca.modify(swelledbase_to_external,
         external_to_finalexternal,
         &external,
         &finalexternal);
   }

   if (locally_nests) {
      *locally_nests = external.getLocalNumberOfBoxes() == 0;
   }
   bool globally_nests = external.getGlobalNumberOfBoxes() == 0;

   return globally_nests;
}



bool MappedBoxLevelConnectorUtils::baseNestsInHeadForMultiblock(
   bool* locally_nests,
   const Connector& connector,
   const IntVector& base_swell,
   const IntVector& head_swell,
   const IntVector& head_nesting_margin,
   const hier::MultiblockMappedBoxTree* domain) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(
      connector.getBase(), base_swell, head_nesting_margin);
   TBOX_ASSERT(connector.isInitialized());
   const tbox::Dimension& dim(connector.getBase().getDim());
   TBOX_ASSERT(base_swell >= IntVector::getZero(dim));
   TBOX_ASSERT(head_nesting_margin >= IntVector::getZero(dim));

   /*
    * To ensure correct results, connector must be sufficiently wide.
    * It should be at least as wide as the combination of base_swell,
    * head_swell and head_nesting_margin.
    */
   const hier::IntVector required_gcw =
      base_swell
      + (connector.getHeadCoarserFlag() ?
         head_swell * connector.getRatio() :
         IntVector::ceiling(head_swell, connector.getRatio()))
      + (connector.getHeadCoarserFlag() ?
         head_nesting_margin * connector.getRatio() :
         IntVector::ceiling(head_nesting_margin, connector.getRatio()))
      ;
   if (!(connector.getConnectorWidth() >= required_gcw)) {
      TBOX_ERROR("MappedBoxLevelConnectorUtils::baseNestsInHead: connector lacks sufficient\n"
         << "ghost cell width for determining whether its base nests\n"
         << "inside its head.");
   }

   hier::OverlapConnectorAlgorithm oca;

   const MappedBoxLevel& base = connector.getBase();
   const MappedBoxLevel& head = connector.getHead();
   const tbox::ConstPointer<GridGeometry>& grid_geom = base.getGridGeometry();

   tbox::Pointer<hier::MappedBoxTree> refined_domain;
   if (domain != NULL) {
      refined_domain = domain->createRefinedTree(head.getRefinementRatio());
   }

   /*
    * We swell the base then check for the parts outside the head if
    * the head domain is grown by head_swell then shrunken by
    * head_nesting_margin.
    *
    * TODO: We can probably remove the base swelling step by converting
    * the base_swell into head index space and add it to
    * head_nesting_margin.
    */

   MappedBoxLevel swelledbase(dim);
   if (base_swell == IntVector::getZero(dim)) {
      swelledbase = base;
   } else {
      const MappedBoxSet& base_mapped_boxes = base.getMappedBoxes();
      MappedBoxSet swelledbase_mapped_boxes;
      for (MappedBoxSet::const_iterator ni = base_mapped_boxes.begin();
           ni != base_mapped_boxes.end(); ++ni) {
         MappedBox swelledbase_mapped_box = *ni;
         swelledbase_mapped_box.getBox().grow(base_swell);
         swelledbase_mapped_boxes.insert(
            swelledbase_mapped_boxes.end(), swelledbase_mapped_box);
      }
      swelledbase.swapInitialize(swelledbase_mapped_boxes,
         base.getRefinementRatio(),
         grid_geom,
         base.getMPI());
   }

   MappedBoxLevel swelledhead(dim);
   NeighborhoodSet swelledbase_eto_swelledhead;
   if (head_swell == IntVector::getZero(dim)) {
      swelledhead = head;
      swelledbase_eto_swelledhead = connector.getNeighborhoodSets();
   } else {
      const MappedBoxSet& head_mapped_boxes = head.getMappedBoxes();

      MappedBoxSet swelledhead_mapped_boxes;

      for (MappedBoxSet::const_iterator ni = head_mapped_boxes.begin();
           ni != head_mapped_boxes.end(); ++ni) {
         MappedBox swelledhead_mapped_box = *ni;

         swelledhead_mapped_box.getBox().grow(head_swell);
         swelledhead_mapped_boxes.insert(
            swelledhead_mapped_boxes.end(), swelledhead_mapped_box);
      }
      swelledhead.swapInitialize(swelledhead_mapped_boxes,
         head.getRefinementRatio(),
         grid_geom,
         head.getMPI());

      connector.getNeighborhoodSets().growNeighbors(
         swelledbase_eto_swelledhead,
         head_swell);
   }


   Connector swelledbase_to_swelledhead;
   swelledbase_to_swelledhead.swapInitialize(
      swelledbase,
      swelledhead,
      connector.getConnectorWidth() - base_swell,
      swelledbase_eto_swelledhead,
      MappedBoxLevel::DISTRIBUTED);
   if ( d_sanity_check_precond &&
        head_swell == hier::IntVector::getZero(dim) ) {
      /*
       * If head was swelled, it may generate undetected overlaps that
       * cannot be compensated for by shrinking the connector width.
       * The additional overlaps do not matter to the nesting check,
       * so it does not affect our result.  Nevertheless, because they
       * are not detected, don't make this check if head was swelled.
       */
      oca.assertOverlapCorrectness(swelledbase_to_swelledhead);
   }

   MappedBoxLevel external(dim);
   Connector swelledbase_to_external;
   if (domain) {
      computeExternalPartsForMultiblock(
         external,
         swelledbase_to_external,
         swelledbase_to_swelledhead,
         -head_nesting_margin,
         *domain);
   } else {
      computeExternalPartsForMultiblock(
         external,
         swelledbase_to_external,
         swelledbase_to_swelledhead,
         -head_nesting_margin,
         MultiblockMappedBoxTree() );
   }
   if (domain) {
      /*
       * If domain is given, do not count external parts that are
       * outside the domain.  In many usages, part of base is outside
       * the domain and we want to ignore those parts.
       */
      hier::MappingConnectorAlgorithm mca;
      std::vector<MappedBox> domain_mapped_boxes;
      domain->getMappedBoxes(domain_mapped_boxes);
      MappedBoxLevel domain_mapped_box_level(
         IntVector::getOne(dim),
         grid_geom,
         connector.getMPI(),
         MappedBoxLevel::GLOBALIZED);
      for (size_t i = 0; i < domain_mapped_boxes.size(); ++i) {
         domain_mapped_box_level.addMappedBox(domain_mapped_boxes[i]);
      }
      Connector external_to_domain(
         external,
         domain_mapped_box_level,
         base_swell);
      oca.findOverlaps(external_to_domain);
      MappedBoxLevel finalexternal(dim);
      Connector external_to_finalexternal;
      computeInternalPartsForMultiblock(
         finalexternal,
         external_to_finalexternal,
         external_to_domain,
         IntVector::getZero(dim),
         *domain);
      mca.modify(swelledbase_to_external,
         external_to_finalexternal,
         &external,
         &finalexternal);
   }

   if (locally_nests) {
      *locally_nests = external.getLocalNumberOfBoxes() == 0;
   }
   bool globally_nests = external.getGlobalNumberOfBoxes() == 0;

   return globally_nests;
}

/*
 ***********************************************************************
 * Make a Connector object for changing the MappedBox indices of a MappedBoxLevel.
 *
 * If sequentialize_global_indices is true, the indices are changed
 * such that they become globally sequential, with processor n
 * starting where processor n-1 ended.  In order to determine what the
 * global indices should be, an allgather communication is used to
 * determine how many mapped_boxes each processor has.  This is a
 * utility function for resetting MappedBox indices to correspond to
 * patch indices while we try to be backward compatible with non-DLBG
 * parts of SAMRAI.
 *
 * If sort_mapped_boxes_by_corner is true, the local MappedBoxes are
 * sorted by their corner indices.  This helps to de-randomize
 * MappedBoxes that may be randomly ordered by non-deterministic
 * algorithms.
 ***********************************************************************
 */
void MappedBoxLevelConnectorUtils::makeSortingMap(
   MappedBoxLevel& sorted_mapped_box_level,
   Connector& output_map,
   const MappedBoxLevel& unsorted_mapped_box_level,
   bool sort_mapped_boxes_by_corner,
   bool sequentialize_global_indices,
   LocalId initial_sequential_index) const
{
   const tbox::Dimension& dim(unsorted_mapped_box_level.getDim());

   if (!sort_mapped_boxes_by_corner && !sequentialize_global_indices) {
      // Make a blank map.
      sorted_mapped_box_level = unsorted_mapped_box_level;
      output_map.initialize(
         unsorted_mapped_box_level,
         sorted_mapped_box_level,
         hier::IntVector::getZero(dim),
         MappedBoxLevel::DISTRIBUTED);
      output_map.setConnectorType(Connector::MAPPING);
      return;
   }

   t_make_sorting_map->start();

   const MappedBoxSet& cur_mapped_boxes =
      unsorted_mapped_box_level.getMappedBoxes();
   int n_mapped_boxes = 
      static_cast<int>(unsorted_mapped_box_level.getLocalNumberOfBoxes());

   LocalId last_index = initial_sequential_index - 1;

   if (sequentialize_global_indices) {

      const int nproc = unsorted_mapped_box_level.getNproc();
      const int rank = unsorted_mapped_box_level.getRank();

      std::vector<int> all_n_mapped_boxes(nproc);
      tbox::SAMRAI_MPI mpi(unsorted_mapped_box_level.getMPI());
      if (mpi.getSize() > 1) {
         mpi.Allgather(&n_mapped_boxes,
            1,
            MPI_INT,
            &all_n_mapped_boxes[0],
            1,
            MPI_INT);
      } else {
         all_n_mapped_boxes[0] = n_mapped_boxes;
      }

      LocalId new_start_index = initial_sequential_index;
      for (int i = 0; i < rank; ++i) {
         new_start_index += all_n_mapped_boxes[i];
      }
      last_index = new_start_index - 1;

   }

   std::vector<MappedBox> real_mapped_box_vector;
   std::vector<MappedBox> periodic_image_mapped_box_vector;
   if (!cur_mapped_boxes.empty()) {
      /*
       * Bypass qsort if we have no mapped_boxes (else there is a memory warning).
       */
      cur_mapped_boxes.separatePeriodicImages(
         real_mapped_box_vector,
         periodic_image_mapped_box_vector);
      if (sort_mapped_boxes_by_corner) {
         qsort((void *)&real_mapped_box_vector[0],
            real_mapped_box_vector.size(),
            sizeof(MappedBox),
            qsortBoxCompare);
      }
   }

   MappedBoxSet new_mapped_boxes;
   NeighborhoodSet edges;

   for (std::vector<MappedBox>::const_iterator ni = real_mapped_box_vector.begin();
        ni != real_mapped_box_vector.end(); ++ni) {

      const MappedBox& cur_mapped_box = *ni;
      const MappedBox new_mapped_box(cur_mapped_box.getBox(),
                                     ++last_index,
                                     cur_mapped_box.getOwnerRank(),
                                     cur_mapped_box.getBlockId(),
                                     cur_mapped_box.getPeriodicId());
      new_mapped_boxes.insert(new_mapped_boxes.end(), new_mapped_box);

      /*
       * Now, add cur_mapped_box's periodic images, but give them cur_mapped_box's
       * new LocalId.  In finding the image mapped_boxes, we use the fact
       * that a real mapped_box's image follows the real mapped_box in a MappedBoxSet.
       */
      MappedBoxSet::const_iterator ini = cur_mapped_boxes.find(cur_mapped_box);
      TBOX_ASSERT(ini != cur_mapped_boxes.end());
      ++ini; // Skip the real mapped_box to look for its image mapped_boxes.
      while (ini != cur_mapped_boxes.end() &&
             ini->getGlobalId() == cur_mapped_box.getGlobalId()) {
         const MappedBox& image_mapped_box = *ini;
         const MappedBox new_image_mapped_box(image_mapped_box.getBox(),
                                              new_mapped_box.getLocalId(),
                                              new_mapped_box.getOwnerRank(),
                                              new_mapped_box.getBlockId(),
                                              image_mapped_box.
                                              getPeriodicId());
         new_mapped_boxes.insert(new_mapped_boxes.end(), new_image_mapped_box);
         ++ini;
      }

      /*
       * Edge for the mapping.  By convention, image mapped_boxes are
       * not explicitly mapped.  Also by convention, we don't create
       * edges unless there is a change.
       */
      if (cur_mapped_box.getLocalId() != new_mapped_box.getLocalId()) {
         edges[cur_mapped_box.getId()].insert(new_mapped_box);
      }
   }

   sorted_mapped_box_level.swapInitialize(
      new_mapped_boxes,
      unsorted_mapped_box_level.getRefinementRatio(),
      unsorted_mapped_box_level.getGridGeometry(),
      unsorted_mapped_box_level.getMPI(),
      MappedBoxLevel::DISTRIBUTED);
   output_map.swapInitialize(
      unsorted_mapped_box_level,
      sorted_mapped_box_level,
      hier::IntVector::getZero(dim),
      edges,
      MappedBoxLevel::DISTRIBUTED);
   output_map.setConnectorType(Connector::MAPPING);

   t_make_sorting_map->stop();
}

/*
 *************************************************************************
 * for use when sorting integers using the C-library qsort
 *************************************************************************
 */
int MappedBoxLevelConnectorUtils::qsortBoxCompare(
   const void* v,
   const void* w)
{
   const hier::MappedBox &mapped_box_v(*(const hier::MappedBox *)v);
   const hier::MappedBox &mapped_box_w(*(const hier::MappedBox *)w);

   if ( mapped_box_v.getBlockId() > mapped_box_w.getBlockId() ) return  1;
   if ( mapped_box_v.getBlockId() < mapped_box_w.getBlockId() ) return -1;

   const tbox::Dimension& dim(mapped_box_v.getDim());

   const hier::IntVector& lowv = mapped_box_v.getBox().lower();
   const hier::IntVector& loww = mapped_box_w.getBox().lower();
   for (int i = 0; i < dim.getValue(); ++i) {
      if (lowv[i] > loww[i]) return 1;
      if (lowv[i] < loww[i]) return -1;
   }

   const hier::IntVector& upv = mapped_box_v.getBox().upper();
   const hier::IntVector& upw = mapped_box_w.getBox().upper();
   for (int i = 0; i < dim.getValue(); ++i) {
      if (upv[i] > upw[i]) return 1;
      if (upv[i] < upw[i]) return -1;
   }

   return 0;
}

/*
 *************************************************************************
 *************************************************************************
 */

void MappedBoxLevelConnectorUtils::computeExternalParts(
   hier::MappedBoxLevel& external,
   hier::Connector& input_to_external,
   const hier::Connector& input_to_reference,
   const hier::IntVector& nesting_width,
   const hier::MappedBoxTree& domain) const
{
   const tbox::Dimension& dim(nesting_width.getDim());

   const BlockId& block_id = domain.getBlockId();

   t_compute_external_parts->start();

   const MappedBoxLevel& input = input_to_reference.getBase();

   const hier::IntVector& zero_vec = IntVector::getZero(dim);
   const hier::IntVector& one_vec = IntVector::getOne(dim);

#ifdef DEBUG_CHECK_ASSERTIONS
   if (!(nesting_width >= zero_vec) &&
       !(nesting_width <= zero_vec)) {
      TBOX_ERROR(
         "MappedBoxLevelConnectorUtils::computeExternalParts: internal error:\n"
         << "nesting_width cannot have mix of positive\n"
         << "and negative values.");
   }
   if (nesting_width >= zero_vec) {
      if (!(input_to_reference.getConnectorWidth() >=
            nesting_width)) {
         TBOX_ERROR(
            "MappedBoxLevelConnectorUtils::computeExternalParts: internal error:\n"
            << "nesting_width "
            << nesting_width << " exceeds\n"
            << "ghost cell width " << input_to_reference.getConnectorWidth()
            << ",\n"
            << "which can lead to possible missed overlaps.");
      }
   }
#endif

   const NeighborhoodSet& input_eto_reference = input_to_reference.getNeighborhoodSets();

   external.initialize(input.getRefinementRatio(),
      input.getGridGeometry(), input.getMPI());

   /*
    *  Keep track of last index so we don't give external an index
    *  used by input.  This is required because if we allow external
    *  to select an index it is not using, it may select one that is
    *  used by input, creating an invalid mapping that prevents
    *  MappingConnectorAlgorithm::modify() from working correctly.
    *
    *  Potential optimization: This loop uses many removeIntersections()
    *  calls.  If it gets too slow, consider precomputing a MappedBoxTree
    *  of all the visible neighbors (or one of the complement of all the
    *  visible neighbors, depending on the sign of nesting_width)
    *  to speed up intersection removals.
    */
   LocalId last_used_index = input.getLastLocalId();

   Connector::NeighborSet reference_nabrs;
   input_eto_reference.getNeighbors(reference_nabrs);
   std::vector<MappedBox> reference_mapped_box_vec;
   for (MappedBoxSetSingleBlockIterator ni(reference_nabrs, block_id);
        ni.isValid();
        ++ni) {
      reference_mapped_box_vec.push_back(*ni);
   }
   reference_nabrs.clear();

   /*
    * Bring reference_mapped_box_vec into refinement ratio of input
    * (for intersection checks).
    */
   if (input_to_reference.getRatio() != hier::IntVector::getOne(dim)) {
      if (input_to_reference.getHeadCoarserFlag()) {
         hier::MappedBoxContainerUtils::refineMappedBoxVectorBoxes(
            reference_mapped_box_vec,
            input_to_reference.getRatio());
      } else {
         hier::MappedBoxContainerUtils::coarsenMappedBoxVectorBoxes(
            reference_mapped_box_vec,
            input_to_reference.getRatio());
      }
   }

   /*
    * Build the tree used for searching for either the internal or external
    * parts (as specified by search_tree_repesents_internal).  The search tree
    * is built differently, depending on the sign of nesting_width.
    */

   hier::MappedBoxTree search_tree(dim);
   bool search_tree_represents_internal;

   if (nesting_width == zero_vec) {
      /*
       * Zero nesting_width.
       * The reference mapped_box_level represents the internl parts.
       */
      search_tree_represents_internal = true;
      search_tree.generateTree(reference_mapped_box_vec);
   } else if (nesting_width >= zero_vec) {
      /*
       * nesting_width is non-negative, grow reference
       * mapped_boxes, which will then represent the internal parts.
       */
      search_tree_represents_internal = true;
      hier::MappedBoxContainerUtils::growMappedBoxVectorBoxes(
         reference_mapped_box_vec,
         nesting_width);
      search_tree.generateTree(reference_mapped_box_vec);
   } else {
      /*
       * nesting_width is non-positive.  The external
       * parts are given by the complex formula above for (A\B).
       */

      search_tree_represents_internal = false;

      hier::BoxList reference_boundary;
      MappedBoxContainerUtils::convertMappedBoxVectorToBoxList(
         reference_mapped_box_vec,
         reference_boundary);
      reference_boundary.grow(one_vec);
      // ... reference_boundary is now (R^1)
      hier::MappedBoxTree reference_mapped_boxes_tree(dim,
                                                    reference_mapped_box_vec);
      t_compute_external_parts_intersection->start();
      reference_boundary.removeIntersections(reference_mapped_boxes_tree);
      // ... reference_boundary is now ( (R^1) \ R )
      if (domain.isInitialized()) {
         if (input.getRefinementRatio() == one_vec) {
            reference_boundary.intersectBoxes(domain);
         } else {
            tbox::Pointer<MappedBoxTree> refined_domain =
               domain.createRefinedTree(input.getRefinementRatio());
            reference_boundary.intersectBoxes(*refined_domain);
         }
      }
      // ... reference_boundary is now ( (R^1) \ R ) <intersection> O )
      t_compute_external_parts_intersection->stop();
      reference_boundary.grow(-nesting_width);
      // ... reference_boundary is now ( (R^1) \ R ) <intersection> O ) - g
      std::vector<MappedBox> reference_boundary_vector;
      MappedBoxContainerUtils::convertBoxListToMappedBoxVector(
         reference_boundary,
         reference_boundary_vector,
         block_id);
      search_tree.generateTree(reference_boundary_vector);

   }

   reference_mapped_box_vec.clear();

   const MappedBoxSet& input_mapped_boxes = input.getMappedBoxes();

   NeighborhoodSet input_eto_external;

   /*
    * For each MappedBox in input_mapped_boxes, compare it to the
    * search tree to compute its external parts.
    */

   for (RealMappedBoxConstIterator ni(input_mapped_boxes); ni.isValid();
        ++ni) {

      const hier::MappedBox& input_mapped_box = *ni;
      const hier::BlockId& input_block_id = input_mapped_box.getBlockId();

      NeighborhoodSet::const_iterator ei =
         input_eto_reference.find(ni->getId());

      if (input_block_id != block_id) {
         external.addMappedBox(input_mapped_box);
      } else if (ei == input_eto_reference.end()) {
         /*
          * Input mapped_box does not overlap reference, so it must be
          * completely external.
          */
         external.addMappedBox(input_mapped_box);
         /*
          * The input_mapped_box should be mapped to itself.
          * We can create such a map, but a missing map
          * means the same thing, so we omit the map
          */
      } else {

         /*
          * external_list will be what remains of input_mapped_box after
          * removing parts we know (by comparing with the reference
          * mapped_box_level) to be internal.
          */
         hier::BoxList external_list(input_mapped_box.getBox());

         t_compute_external_parts_intersection->start();
         if (search_tree_represents_internal) {
            external_list.removeIntersections(search_tree);
         } else {
            external_list.intersectBoxes(search_tree);
         }
         t_compute_external_parts_intersection->stop();

         /*
          * Make external mapped_boxes from external_list and create
          * connectivity with input.
          *
          * We force a neighbor container to be created for input,
          * even if it is replaced with nothing, because
          * MappingConnectorAlgorithm::modify() interprets the presence
          * of a neighbor list as a mapping.  Missing neighbor
          * list mean "no change", which is not what we want.
          */
         if (external_list.size() == 1 &&
             external_list.getFirstItem() == input_mapped_box.getBox()) {
            /*
             * The entire input_mapped_box is external.
             * The input_mapped_box should be mapped to itself.
             * We can create such a map, but a missing map
             * means the same thing, so we omit the map
             */
            external.addMappedBox(input_mapped_box);
         } else {
            hier::Connector::NeighborSet
            & replacements = input_eto_external[input_mapped_box.getId()];
            for (hier::BoxList::Iterator bi(external_list); bi; bi++) {
               const hier::MappedBox
               external_mapped_box((*bi),
                                   ++last_used_index,
                                   input_mapped_box.getOwnerRank(),
                                   input_mapped_box.getBlockId());
               external.addMappedBox(external_mapped_box);
               // Set connectivities between input and external.
               replacements.insert(external_mapped_box);
            }
         }
      }

   }

#ifdef DEBUG_CHECK_ASSERTIONS
   if (external.getMappedBoxes().empty()) {
      /*
       * If there are no external parts, then all in input
       * should be mapped to empty neighbor containers according
       * to the definition of a map in MappingConnectorAlgorithm::modify().
       */
      TBOX_ASSERT(input_eto_external.size() ==
         input.getLocalNumberOfBoxes());
      for (NeighborhoodSet::const_iterator ci = input_eto_external.begin();
           ci != input_eto_external.end(); ++ci) {
         TBOX_ASSERT((*ci).second.empty());
      }
   }
#endif

   /*
    * input_to_external has zero width because a non-zero width means
    * that we should have some edges from an input MappedBox and the
    * external parts of a nearby-input MappedBox, which we don't.
    * Moreover, we cannot compute edges unless we know input<==>input.
    */
   input_to_external.swapInitialize(
      input,
      external,
      zero_vec,
      input_eto_external,
      hier::MappedBoxLevel::DISTRIBUTED);
   input_to_external.setConnectorType(Connector::BASE_GENERATED);

   TBOX_ASSERT(input_to_external.isLocal());

   t_compute_external_parts->stop();
}



/*
*************************************************************************
*************************************************************************
*/
void MappedBoxLevelConnectorUtils::computeExternalPartsForMultiblock(
   hier::MappedBoxLevel& external,
   hier::Connector& input_to_external,
   const hier::Connector& input_to_reference,
   const hier::IntVector& nesting_width,
   const hier::MultiblockMappedBoxTree& domain) const
{
   t_compute_external_parts->start();

   computeInternalOrExternalPartsForMultiblock(
      external,
      input_to_external,
      'e',
      input_to_reference,
      nesting_width,
      domain);

   t_compute_external_parts->stop();
}




/*
 *************************************************************************
 *************************************************************************
 */

void MappedBoxLevelConnectorUtils::computeInternalParts(
   hier::MappedBoxLevel& internal,
   hier::Connector& input_to_internal,
   const hier::Connector& input_to_reference,
   const hier::IntVector& nesting_width,
   const hier::MappedBoxTree& domain) const
{
   const tbox::Dimension& dim(nesting_width.getDim());

   const BlockId& block_id = domain.getBlockId();

   t_compute_internal_parts->start();

   const MappedBoxLevel& input = input_to_reference.getBase();

   const IntVector& zero_vec(IntVector::getZero(dim));
   const IntVector& one_vec(IntVector::getOne(dim));

#ifdef DEBUG_CHECK_ASSERTIONS
   if (!(nesting_width >= zero_vec) &&
       !(nesting_width <= zero_vec)) {
      TBOX_ERROR(
         "MappedBoxLevelConnectorUtils::computeInternalParts: input parameter error:\n"
         << "nesting_width cannot have mix of positive\n"
         << "and negative values.");
   }
   if (nesting_width >= zero_vec) {
      if (!(input_to_reference.getConnectorWidth() >=
            nesting_width)) {
         TBOX_ERROR(
            "MappedBoxLevelConnectorUtils::computeInternalParts: internal error:\n"
            << "nesting_width "
            << nesting_width << " exceeds\n"
            << "ghost cell width " << input_to_reference.getConnectorWidth()
            << ",\n"
            << "which can lead to possible missed overlaps.");
      }
   }
#endif

   const NeighborhoodSet& input_eto_reference = input_to_reference.getNeighborhoodSets();

   internal.initialize(input.getRefinementRatio(),
      input.getGridGeometry(), input.getMPI());

   /*
    * Keep track of last index so we don't give internal an index
    * used by input.  This is required because if we allow internal
    * to select an index it is not using, it may select one that is
    * used by input, creating an invalid mapping that prevents
    * MappingConnectorAlgorithm::modify() from working correctly.
    */
   LocalId last_used_index = input.getLastLocalId();

   Connector::NeighborSet reference_nabrs;
   input_eto_reference.getNeighbors(reference_nabrs);

   /*
    * We eventually put the visible reference neighbors in the
    * reference_mapped_box_vec for the generation of a search tree.
    * However, if the user gave a domain, we have to manipulate these
    * MappedBoxes first, so we put them in a BoxList instead.
    */
   int num_blocks = input.getGridGeometry()->getNumberBlocks();
   BoxList reference_box_list;
   std::vector<std::vector<MappedBox> > reference_mapped_box_vec;
   reference_mapped_box_vec.resize(num_blocks);
   if ( domain.isInitialized() ) {
      for (MappedBoxSetSingleBlockIterator ni(reference_nabrs, block_id);
           ni.isValid();
           ++ni) {
         reference_box_list.appendItem(ni->getBox());
      }
   }
   else {
      for (int b = 0; b < num_blocks; b++) {
         hier::BlockId bid(b);
         for (MappedBoxSetSingleBlockIterator ni(reference_nabrs, bid);
              ni.isValid(); ++ni) {
            reference_mapped_box_vec[b].push_back(*ni);
         }
      }
   }

   reference_nabrs.clear(); // Finished with this.

   /*
    * Bring reference boxes into refinement ratio of input (for
    * intersection checks).  These boxes can be in either
    * reference_box_list or reference_mapped_box_vector.
    */
   if (input_to_reference.getRatio() != hier::IntVector::getOne(dim)) {
      if (input_to_reference.getHeadCoarserFlag()) {
         if ( reference_box_list.size() ) {
            reference_box_list.refine(input_to_reference.getRatio());
         }
         for (int b = 0; b < num_blocks; b++) {
            if ( reference_mapped_box_vec[b].size() ) {
               hier::MappedBoxContainerUtils::refineMappedBoxVectorBoxes(
                  reference_mapped_box_vec[b],
                  input_to_reference.getRatio());
            }
         }
      } else {
         if ( reference_box_list.size() ) {
            reference_box_list.coarsen(input_to_reference.getRatio());
         }
         for (int b = 0; b < num_blocks; b++) {
            if ( reference_mapped_box_vec[b].size() ) {
               hier::MappedBoxContainerUtils::coarsenMappedBoxVectorBoxes(
                  reference_mapped_box_vec[b],
                  input_to_reference.getRatio());
            }
         }
      }
   }

   /*
    * Build the search tree containing either the internal or external
    * parts (as specified by search_tree_repesents_internal).  The
    * search tree is built differently, depending on the sign of
    * nesting_width.
    */

   tbox::Array<hier::MappedBoxTree> search_tree(num_blocks, MappedBoxTree(dim));
   bool search_tree_represents_internal;

   if (nesting_width >= zero_vec) {
      /*
       * nesting_width is non-negative, grow reference
       * MappedBoxes, which will then represent the internal parts.
       */
      search_tree_represents_internal = true;
      if ( domain.isInitialized() ) {
         if (nesting_width != zero_vec) {
            reference_box_list.grow(nesting_width);
         }
         if (input.getRefinementRatio() == one_vec) {
            reference_box_list.intersectBoxes(domain);
         } else {
            tbox::Pointer<MappedBoxTree> refined_domain =
               domain.createRefinedTree(input.getRefinementRatio());
            reference_box_list.intersectBoxes(*refined_domain);
         }
         MappedBoxContainerUtils::convertBoxListToMappedBoxVector(
            reference_box_list,
            reference_mapped_box_vec[block_id.getBlockValue()],
            block_id);
         search_tree[block_id.getBlockValue()].generateTree(
            reference_mapped_box_vec[block_id.getBlockValue()]);
      }
      else {
         for (int b = 0; b < num_blocks; b++) {
            if (nesting_width != zero_vec) {
               hier::MappedBoxContainerUtils::growMappedBoxVectorBoxes(
                  reference_mapped_box_vec[b],
                  nesting_width);
            }
            search_tree[b].generateTree(reference_mapped_box_vec[b]); 
         }
      }
   } else {
      TBOX_ERROR(
         "For lack of need, this method is not currently supporting non-positive g.");
   }

   const MappedBoxSet& input_mapped_boxes = input.getMappedBoxes();

   NeighborhoodSet input_eto_internal;

   /*
    * For each MappedBox in input_mapped_boxes, compare it to the
    * search tree to compute its internal parts.
    */

   for (RealMappedBoxConstIterator ni(input_mapped_boxes); ni.isValid();
        ++ni) {

      const hier::MappedBox& input_mapped_box = *ni;

      NeighborhoodSet::const_iterator ei = input_eto_reference.find(ni->getId());

      if (ei == input_eto_reference.end()) {
         /*
          * Input mapped_box does not overlap reference, so it must be
          * completely external.  Create a blank neighbor set to map it
          * to nothing.
          */
         input_eto_internal[input_mapped_box.getId()];
      } else {

         /*
          * internal_list will be the intersection of input_mapped_box
          * neighbors on the reference mapped_box_level.
          */
         hier::BoxList internal_list(input_mapped_box.getBox());
         int block_num = input_mapped_box.getBlockId().getBlockValue();

         t_compute_internal_parts_intersection->start();
         if (search_tree_represents_internal) {
            internal_list.intersectBoxes(search_tree[block_num]);
         } else {
            internal_list.removeIntersections(search_tree[block_num]);
         }
         t_compute_internal_parts_intersection->stop();

         /*
          * Make internal mapped_boxes from internal_list and create
          * connectivity with input.
          *
          * We force a neighbor container to be created for input,
          * even if it is replaced with nothing, because
          * MappingConnectorAlgorithm::modify() interprets the presence
          * of a neighbor list as a mapping.  Missing neighbor
          * list mean "no change", which is not what we want.
          */
         if (internal_list.size() == 1 &&
             internal_list.getFirstItem() == input_mapped_box.getBox()) {
            /*
             * The entire input_mapped_box is internal.
             * The input_mapped_box should be mapped to itself.
             * We can create such a map, but a missing map
             * means the same thing, so we omit the map
             */
            internal.addMappedBox(input_mapped_box);
         } else {
            hier::Connector::NeighborSet
            & replacements = input_eto_internal[input_mapped_box.getId()];
            for (hier::BoxList::Iterator bi(internal_list); bi; bi++) {
               const hier::MappedBox
               internal_mapped_box((*bi),
                                   ++last_used_index,
                                   input_mapped_box.getOwnerRank(),
                                   input_mapped_box.getBlockId());
               internal.addMappedBox(internal_mapped_box);
               // Set connectivities between input and internal.
               replacements.insert(internal_mapped_box);
            }
         }
      }

   }

#ifdef DEBUG_CHECK_ASSERTIONS
   if (internal.getMappedBoxes().empty()) {
      /*
       * If there are no internal parts, then all in input
       * should be mapped to empty neighbor containers according
       * to the definition of a map in MappingConnectorAlgorithm::modify().
       */
      if (0) {
         tbox::plog << "input/reference/internal:\n"
                    << "input:\n" << input.format("", 2)
                    << "reference:\n" << input_to_reference.getHead().format("", 2)
                    << "input_to_reference:\n" << input_to_reference.format("", 3)
                    << "internal:\n" << internal.format("", 2)
                    << input_eto_internal.size() << ' '
                    << input.getMappedBoxes().size() << std::endl;
      }
      TBOX_ASSERT(input_eto_internal.size() == input.getLocalNumberOfBoxes());
      for (NeighborhoodSet::const_iterator ci = input_eto_internal.begin();
           ci != input_eto_internal.end(); ++ci) {
         TBOX_ASSERT((*ci).second.empty());
      }
   }
#endif

   input_to_internal.swapInitialize(
      input,
      internal,
      zero_vec,
      input_eto_internal,
      hier::MappedBoxLevel::DISTRIBUTED);

   TBOX_ASSERT(input_to_internal.isLocal());

   t_compute_internal_parts->stop();
}


/*
*************************************************************************
*************************************************************************
*/
void MappedBoxLevelConnectorUtils::computeInternalPartsForMultiblock(
   hier::MappedBoxLevel& internal,
   hier::Connector& input_to_internal,
   const hier::Connector& input_to_reference,
   const hier::IntVector& nesting_width,
   const hier::MultiblockMappedBoxTree& domain) const
{
   t_compute_internal_parts->start();

   computeInternalOrExternalPartsForMultiblock(
      internal,
      input_to_internal,
      'i',
      input_to_reference,
      nesting_width,
      domain);

   t_compute_internal_parts->stop();
}




/*
*************************************************************************
* This version of computeInternalParts does not require a domain.
* The domain is taken to be big enough that it does not affect
* the definition of "internal".
*************************************************************************
*/
void MappedBoxLevelConnectorUtils::computeInternalParts(
   hier::MappedBoxLevel& internal,
   hier::Connector& input_to_internal,
   const hier::Connector& input_to_reference,
   const hier::IntVector& nesting_width) const
{
   const tbox::Dimension& dim(nesting_width.getDim());
   const hier::MappedBoxTree dummy_domain(dim);
   computeInternalParts(internal,
                        input_to_internal,
                        input_to_reference,
                        nesting_width,
                        dummy_domain);
   return;
}



/*
*************************************************************************
* Methods computeInternalPartsForMultiblock and
* computeExternalPartsForMultiblock delegates to this method.
*
* Compare an input MappedBoxLevel to a "reference" MappedBoxLevel.
* Identify parts of the input that are internal or external (depending
* on the value of internal_or_external) to the reference
* MappedBoxLevel, and store the in/external parts in a MappedBoxLevel.
* Create the input_to_parts Connector between the input and these
* parts.
*
* For generality, the reference MappedBoxLevel can be grown a
* specified amount (nesting_width) before comparing.  nesting_width
* must be in the index space of the input MappedBoxLevel (not the
* reference MappedBoxLevel, despite the name).  A negative growth
* indicates shrinking the reference layer at its boundary.
*
* As a practical consideration of how this method is used, we do not
* shrink the reference layer where it touches the domain boundary.
* This feature can be disabled by specifying an uninitialized domain
* object.
*
* On return, input_to_parts is set to an appropriate mapping for use
* in MappingConnectorAlgorithm::modify().
*
* This method does not require any communication.
*
* Formula for computing external parts:
*
* Definitions:
* L = input MappedBoxLevel
* R = reference MappedBoxLevel
* g = nesting width (non-negative or non-positive, but not mixed)
* E = parts of L external to R^g (R^g is R grown by g)
* I = parts of L internal to R^g (R^g is R grown by g)
* O = domain (without periodic images).  Universe, if not specified.
* \ = set theory notation.  x \ y means set x with y removed from it.
*
* For non-negative g:
*
* E := L \ { (R^g) <intersection> O }
* I := L <intersection> { (R^g) <intersection> O }
*
* For non-positive g:
*
* E := L <intersection> { ( ( (R^1) \ R ) <intersection> O )^(-g) }
* I := L \ { ( ( (R^1) \ R ) <intersection> O )^(-g) }
*
* A requirement of the computation for negative g is that input must
* nest in R^(1-g).  In other words: L \ (R^(1-g)} = <empty>.  If not
* satisfied, this method may classify some external parts as
* internal.
*
*************************************************************************
*/
void MappedBoxLevelConnectorUtils::computeInternalOrExternalPartsForMultiblock(
   hier::MappedBoxLevel& parts,
   hier::Connector& input_to_parts,
   char internal_or_external,
   const hier::Connector& input_to_reference,
   const hier::IntVector& nesting_width,
   const hier::MultiblockMappedBoxTree& domain) const
{
   const MappedBoxLevel& input = input_to_reference.getBase();

   const tbox::ConstPointer<GridGeometry> &grid_geometry(
      input.getGridGeometry());

   const tbox::Dimension &dim(input.getDim());
   const IntVector& zero_vec(IntVector::getZero(input.getDim()));
   const IntVector& one_vec(IntVector::getOne(dim));

#ifdef DEBUG_CHECK_ASSERTIONS
   const char *caller = internal_or_external == 'i' ?
      "computInternalPartsForMultiblock" : "computeExternalpartsForMultiblock";

   // Sanity check inputs.

   if (!(nesting_width >= zero_vec) && !(nesting_width <= zero_vec)) {
      TBOX_ERROR(
         "MappedBoxLevelConnectorUtils::" << caller << ": error:\n"
         << "nesting_width may not have mix of positive\n"
         << "and negative values.");
   }

   if (nesting_width >= zero_vec) {
      if (!(input_to_reference.getConnectorWidth() >=
            nesting_width)) {
         TBOX_ERROR(
            "MappedBoxLevelConnectorUtils::" << caller << ": error:\n"
            << "nesting_width "
            << nesting_width << " exceeds\n"
            << "ghost cell width " << input_to_reference.getConnectorWidth()
            << ",\n"
            << "which can lead to erroneous results.");
      }
   }
#endif

   parts.initialize(input.getRefinementRatio(),
      input.getGridGeometry(), input.getMPI());


   /*
    * Get the set of neighboring boxes on the reference
    * MappedBoxLevel.  We first store these boxes in a NeighborSet in
    * order to remove duplicate entries.  Then we move them into a
    * vector for box manipulation.
    */
   Connector::NeighborSet reference_nabrs;
   input_to_reference.getNeighborhoodSets().getNeighbors(reference_nabrs);
   std::vector<MappedBox> reference_mapped_box_vec(
      reference_nabrs.begin(), reference_nabrs.end());
   reference_nabrs.clear();

   /*
    * Bring reference_mapped_box_vec into refinement ratio of input
    * (for intersection checks).
    */
   if (input_to_reference.getRatio() != hier::IntVector::getOne(dim)) {

      if (input_to_reference.getHeadCoarserFlag()) {
         hier::MappedBoxContainerUtils::refineMappedBoxVectorBoxes(
            reference_mapped_box_vec,
            input_to_reference.getRatio());
      } else {
         hier::MappedBoxContainerUtils::coarsenMappedBoxVectorBoxes(
            reference_mapped_box_vec,
            input_to_reference.getRatio());
      }

   }


   /*
    * Build a search tree containing either the internal or external
    * parts of the reference MappedBoxLevel, depending on the sign of
    * nesting_width.  If the nesting_width is non-negative, the
    * internal parts are the same as the reference, possibly after
    * growing.  If it is negative, shrinking the reference boxes does
    * not work.  We grow its complement and grow the complement by
    * -nesting_width.  The result represents the external parts of the
    * reference.
    */

   const bool search_tree_represents_internal = nesting_width >= zero_vec;

   hier::MultiblockMappedBoxTree search_tree;

   if (search_tree_represents_internal) {

      if ( !(nesting_width == zero_vec) ) {
         hier::MappedBoxContainerUtils::growMappedBoxVectorBoxes(
            reference_mapped_box_vec,
            nesting_width);
      }

      search_tree.generateTree(grid_geometry,
                               reference_mapped_box_vec);

   } else {

      /*
       * nesting_width is non-positive.  The external parts are given
       * by the grown boundary boxes.
       */

      std::map<BlockId,BoxList> reference_boundary;
      computeBoxesAroundBoundary(
         reference_boundary,
         reference_mapped_box_vec,
         input.getRefinementRatio(),
         grid_geometry );
      // ... reference_boundary is now ( (R^1) \ R )

      if (domain.isInitialized()) {

         if (input.getRefinementRatio() == one_vec) {
            for ( std::map<BlockId,BoxList>::iterator mi=reference_boundary.begin();
                  mi!=reference_boundary.end(); ++mi ) {
               mi->second.intersectBoxes(
                  mi->first,
                  input.getRefinementRatio(),
                  domain);
            }
         } else {
            tbox::Pointer<MultiblockMappedBoxTree> refined_domain =
               domain.createRefinedTree(input.getRefinementRatio());
            for ( std::map<BlockId,BoxList>::iterator mi=reference_boundary.begin();
                  mi!=reference_boundary.end(); ++mi ) {
               mi->second.intersectBoxes(
                  mi->first,
                  input.getRefinementRatio(),
                  *refined_domain);
            }
         }

      }
      // ... reference_boundary is now ( ( (R^1) \ R ) <intersection> O )

      for ( std::map<BlockId,BoxList>::iterator mi=reference_boundary.begin();
            mi!=reference_boundary.end(); ++mi ) {
         mi->second.grow(-nesting_width);
      }
      // ... reference_boundary is now ( ( (R^1) \ R ) <intersection> O )^(-g)

      std::vector<MappedBox> reference_boundary_vector;
      for ( std::map<BlockId,BoxList>::const_iterator mi=reference_boundary.begin();
            mi!=reference_boundary.end(); ++mi ) {
         MappedBoxContainerUtils::convertBoxListToMappedBoxVector(
            mi->second,
            reference_boundary_vector,
            mi->first);
      }

      search_tree.generateTree(grid_geometry, reference_boundary_vector);

   } // search_tree_represents_internal == false

   reference_mapped_box_vec.clear();


   const NeighborhoodSet& input_eto_reference = input_to_reference.getNeighborhoodSets();

   /*
    * Keep track of last index so we don't give parts an index
    * used by input.  This is required because if we allow parts
    * to select an index it is not using, it may select one that is
    * used by input, creating an invalid mapping that prevents
    * MappingConnectorAlgorithm::modify() from working correctly.
    */
   LocalId last_used_index = input.getLastLocalId();

   NeighborhoodSet input_eto_parts;

   const bool compute_overlaps =
      search_tree_represents_internal == (internal_or_external == 'i');

   /*
    * For each MappedBox in input, compare it to the search tree to
    * compute its overlapping (or non-overlapping) parts.
    */

   const MappedBoxSet& input_mapped_boxes = input.getMappedBoxes();

   for (RealMappedBoxConstIterator ni(input_mapped_boxes); ni.isValid();
        ++ni) {

      const hier::MappedBox& input_mapped_box = *ni;

      NeighborhoodSet::const_iterator ei = input_eto_reference.find(ni->getId());

      if (ei == input_eto_reference.end()) {
         /*
          * Absence of a neighbor set in the overlap Connector means
          * the input MappedBox does not overlap the reference
          * MappedBoxLevel.
          */
         if ( compute_overlaps ) {
            /*
             * Trying to get the overlapping parts.  Create empty
             * neighbor list to indicate there is no such parts.
             */
            input_eto_parts[input_mapped_box.getId()];
         }
         else {
            /*
             * Trying to get the non-overlapping parts.
             * Non-overlapping parts is the whole box.
             */
            parts.addMappedBox(input_mapped_box);
         }

      } else {

         hier::BoxList parts_list(input_mapped_box.getBox());
         /*
          * Compute parts of input_mapped_box either overlapping
          * or nor overlapping the search_tree.
          *
          * Note about intersections in singularity neighbor blocks:
          * Cells from multiple singularity neighbor blocks can
          * coincide when transformed into input_mapped_box's block.
          * There is no way to specify that a cell in input_mapped_box
          * intersects in some singularity block neighbors but not
          * others.  By comparing to singularity neighbor blocks, we
          * take the convention that intersection in one singularity
          * block neighbor is considered intersection in all at the
          * same singularity.  When compute_overlaps == true,
          * this can lead to overspecifying parts and
          * underspecifying external parts, and vice versa.
          */

         t_compute_internal_parts_intersection->start();
         if (compute_overlaps) {
            parts_list.intersectBoxes(
               input_mapped_box.getBlockId(),
               input.getRefinementRatio(),
               search_tree,
               true /* Count singularity neighbors */ );
         } else {
            parts_list.removeIntersections(
               input_mapped_box.getBlockId(),
               input.getRefinementRatio(),
               search_tree,
               true /* Count singularity neighbors */ );
         }
         t_compute_internal_parts_intersection->stop();


         /*
          * Make MappedBoxes from parts_list and create
          * Connector from input.
          */
         if (parts_list.size() == 1 &&
             parts_list.getFirstItem() == input_mapped_box.getBox()) {

            /*
             * The entire input_mapped_box is the part we want.
             * The input_mapped_box should be mapped to itself.
             * We can create such a map, but a missing map
             * means the same thing, so we omit the map
             */
            parts.addMappedBox(input_mapped_box);

         } else {

            hier::Connector::NeighborSet
            & replacements = input_eto_parts[input_mapped_box.getId()];
            for (hier::BoxList::Iterator bi(parts_list); bi; bi++) {
               const hier::MappedBox
               parts_mapped_box((*bi),
                                ++last_used_index,
                                input_mapped_box.getOwnerRank(),
                                input_mapped_box.getBlockId());
               parts.addMappedBox(parts_mapped_box);

               // Set connectivities between input and internal.
               replacements.insert(parts_mapped_box);
            }

         } // parts_list

      } // ei != input_eto_reference.end()

   } // Loop through input_mapped_boxes


#ifdef DEBUG_CHECK_ASSERTIONS
   if (parts.getMappedBoxes().empty()) {
      /*
       * If there are no parts, then all in input
       * should be mapped to empty neighbor containers according
       * to the definition of a map in MappingConnectorAlgorithm::modify().
       */
      if (input_eto_parts.size() != input.getLocalNumberOfBoxes()) {
         tbox::perr <<"MappedBoxLevelConnectorUtils::" << caller << ": library error:\n"
                    <<"There are no parts, so input MappedBoxLevel should be completely mapped away.\n"
                    <<"However, not all input MappedBoxes have been mapped.\n"
                    <<"input MappedBoxLevel:\n" << input.format("",2)
                    <<"input_eto_parts:\n" << input_eto_parts.format("",2)
            ;
         TBOX_ERROR("Library error\n");
      }
      for (NeighborhoodSet::const_iterator ci = input_eto_parts.begin();
           ci != input_eto_parts.end(); ++ci) {
         TBOX_ASSERT((*ci).second.empty());
      }
   }
#endif

   /*
    * The output mapping has zero width.  For a mapping, zero width
    * means that no MappedBox is mapped to something outside its
    * extent.
    */
   input_to_parts.swapInitialize(
      input,
      parts,
      zero_vec,
      input_eto_parts,
      hier::MappedBoxLevel::DISTRIBUTED);

   TBOX_ASSERT(input_to_parts.isLocal());

   return;
}



/*
*************************************************************************
Given a MappedBoxSet, compute its boundary as a set of boxes located
just outside it.

Given a set of boxes R, the boundary is computed as (R^1)\R.
R^1 means grown boxes in R by a width of one.
*************************************************************************
*/
void MappedBoxLevelConnectorUtils::computeBoxesAroundBoundary(
   std::map<BlockId,BoxList> &boundary,
   const std::vector<MappedBox>& mapped_boxes,
   const IntVector &refinement_ratio,
   const tbox::ConstPointer<GridGeometry> &grid_geometry,
   const bool simplify_boundary_boxes ) const
{
   const tbox::Dimension &dim(grid_geometry->getDim());
   const IntVector& one_vec(IntVector::getOne(dim));

   boundary.clear();

   for ( std::vector<MappedBox>::const_iterator bi=mapped_boxes.begin();
         bi!=mapped_boxes.end(); ++bi ) {
      boundary[bi->getBlockId()].appendItem(bi->getBox());
   }
   // ... boundary is now R

   for ( std::map<BlockId,BoxList>::iterator mi=boundary.begin();
         mi!=boundary.end(); ++mi ) {
      mi->second.grow(one_vec);
   }
   // ... boundary is now (R^1)

   hier::MultiblockMappedBoxTree reference_mapped_boxes_tree(
      grid_geometry,
      mapped_boxes);

   for ( std::map<BlockId,BoxList>::iterator mi=boundary.begin();
         mi!=boundary.end(); ++mi ) {
      const BlockId &block_id(mi->first);
      BoxList &box_list(mi->second);
      /*
       * Leave boundary boxes in singularity neighbor blocks.
       * These are specially handled in the following if-block.
       */
      const bool include_singularity_neighbors(false);
      box_list.removeIntersections(
         block_id,
         refinement_ratio,
         reference_mapped_boxes_tree,
         include_singularity_neighbors );
   }
   // ... boundary is now ( (R^1) \ R )

   if ( grid_geometry->getNumberOfBlockSingularities() > 0 ) {
      /*
       * The boundary obtained by the formula (R^1)\R can have
       * errors at multiblock singularities.  Fix it here.
       *
       * Boundaries with codimension > 1 (node boundaries in 2D; node
       * and edge boundaries in 3D) passing through a singularity
       * cannot be computed correctly by the formula.
       *
       * What we want to determine is whether the boundary touches the
       * singularity point.
       *
       * - Problem 1: If R touches the singularity in all blocks that
       * touch the singularity, then the boundary does not pass
       * through the singularity.  However, at a reduced connectivity,
       * the formula computes that the boundary does pass through the
       * singularity.
       *
       * - Problem 2: If R touches the singularity in some, but not
       * all, blocks that touch the singularity, then the boundary
       * passes through the singularity point.  However, the formula
       * (R^1)\R may compute that the boundary does not.
       *
       * In both cases, the problem is that the cells representing the
       * codimension > 1 boundaries are removed or left behind when
       * they should not be.  At reduced connectivity singularities,
       * the fix is to always remove them because they do not live in
       * the index space of any block.  At enhanced connectivity, we determine
       * through some box calculus portions of the boundary that are touched
       * by MappedBoxes in all blocks.  The boundary cannot touch these portions
       * so we remove these parts of the boundary.
       */

      std::set<BlockId> blocks_with_mapped_boxes;
      for ( std::vector<MappedBox>::const_iterator bi=mapped_boxes.begin();
            bi!=mapped_boxes.end(); ++bi ) {
         blocks_with_mapped_boxes.insert(bi->getBlockId());
      }


      for ( std::set<BlockId>::const_iterator bi=blocks_with_mapped_boxes.begin();
            bi!=blocks_with_mapped_boxes.end(); ++bi ) {

         const BlockId &block_id(*bi);

         /*
          * Compute a version of singularity boxes for reduced
          * connectivity by removing enhanced connectivity boxes from
          * the singularity box list.  We will remove reduced
          * connectivity singularity boxes from the boundary
          * description, because they do not live in a valid index
          * space.
          */
         BoxList reduced_connectivity_singularity_boxes =
            grid_geometry->getSingularityBoxList(block_id.getBlockValue());

         const tbox::List<GridGeometry::Neighbor> &neighbors(
            grid_geometry->getNeighbors(block_id.getBlockValue()));

         for ( tbox::List<GridGeometry::Neighbor>::Iterator ni(neighbors);
               ni; ni++ ) {
            const GridGeometry::Neighbor &neighbor(*ni);
            if ( neighbor.isSingularity() ) {
               reduced_connectivity_singularity_boxes.removeIntersections(
                  neighbor.getTranslatedDomain());
            }
         }

         if ( ! reduced_connectivity_singularity_boxes.isEmpty() ) {
            if ( refinement_ratio != one_vec ) {
               reduced_connectivity_singularity_boxes.refine(refinement_ratio);
            }
            boundary[block_id].removeIntersections(reduced_connectivity_singularity_boxes);
         }

         /*
          * Intersect singularity_boxes with MappedBoxes from each
          * singularity neighbor.  What remains is where all
          * singularity neighbors have MappedBoxes touching the
          * singularity.  The remains tell us where the boundary does
          * not touch the singularity, overriding what the (R^1)\R
          * formula says.
          */

         BoxList singularity_boxes =
            grid_geometry->getSingularityBoxList(block_id.getBlockValue());
         if ( refinement_ratio != one_vec ) {
            singularity_boxes.refine(refinement_ratio);
         }

         for ( tbox::List<GridGeometry::Neighbor>::Iterator ni(neighbors);
               ni; ni++ ) {
            const GridGeometry::Neighbor &neighbor(*ni);
            const BlockId neighbor_block_id(neighbor.getBlockNumber());

            if ( neighbor.isSingularity() &&
                 reference_mapped_boxes_tree.hasMappedBoxInBlock(neighbor_block_id) ) {

               grid_geometry->translateBoxList(singularity_boxes,
                                               refinement_ratio,
                                               neighbor_block_id,
                                               block_id);

               singularity_boxes.intersectBoxes(
                  reference_mapped_boxes_tree.getSingleBlockMappedBoxTree(neighbor_block_id));

               grid_geometry->translateBoxList(singularity_boxes,
                                               refinement_ratio,
                                               block_id,
                                               neighbor_block_id);
            }
         }

         boundary[block_id].removeIntersections(singularity_boxes);
      }

   } // grid_geometry->getNumberOfBlockSingularities() > 0

   if ( simplify_boundary_boxes ) {
      for ( std::map<BlockId,BoxList>::iterator mi=boundary.begin();
            mi!=boundary.end(); ++mi ) {
         mi->second.simplifyBoxes();
      }
   }

   return;
}


/*
 *************************************************************************
 * Given a mapping from an original MappedBoxLevel to parts to be
 * removed, construct the remainder MappedBoxLevel and the mapping from
 * the original to a remainder.
 *************************************************************************
 */

void MappedBoxLevelConnectorUtils::makeRemainderMap(
   MappedBoxLevel& remainder,
   Connector& orig_to_remainder,
   const Connector& orig_to_rejection) const
{
   TBOX_ASSERT(orig_to_rejection.isLocal());

   const tbox::Dimension& dim(remainder.getDim());

   /*
    * remainder_nodes starts as a copy of orig codes.
    * It will be modified to become the remainder version.
    *
    * orig_to_remainder is the mapping between orig and
    * its properly remainder version.
    */

   const MappedBoxLevel& orig = orig_to_rejection.getBase();
   const MappedBoxSet& orig_nodes = orig.getMappedBoxes();
   const int rank = orig.getRank();
   MappedBoxSet remainder_nodes = orig_nodes;

   /*
    * Track last used index to ensure we use unique indices for new
    * nodes, so that MappingConnectorAlgorithm::modify() works
    * properly.
    */
   LocalId last_used_index = orig.getLastLocalId();

   NeighborhoodSet orig_eto_remainder;

   const NeighborhoodSet& orig_eto_rejection = orig_to_rejection.getNeighborhoodSets();

   for (MappedBoxSet::const_iterator ni = orig_nodes.begin();
        ni != orig_nodes.end(); ++ni) {

      const hier::MappedBox& orig_node = *ni;
      const MappedBoxId mapped_box_id = orig_node.getId();

      NeighborhoodSet::const_iterator ci =
         orig_eto_rejection.find(orig_node.getId());

      if (ci == orig_eto_rejection.end()) {
         /*
          * By the definition of a mapping Connector, no mapping means
          * the entire orig_node is rejected.
          *
          * - Erase rejected node from remainder_nodes
          * - Build connectivities in orig_to_remainder (empty neighbor list).
          */
         remainder_nodes.erase(orig_node);

         TBOX_ASSERT(orig_eto_remainder.find(
               mapped_box_id) == orig_eto_remainder.end());

         orig_eto_remainder[mapped_box_id];
      } else if (ci->second.empty()) {
         /*
          * By the definition of a mapping Connector, empty mapping
          * means entire orig_node remains.
          *
          * No orig<==>remainder mapping is required.
          */
      } else {
         /*
          * The orig_node is partially rejected.
          *
          * - Erase rejected node from remainder_nodes
          * - Remove rejected parts to obtain remainder parts
          * - Add remaining parts to remainder_nodes
          * - Build connectivities in orig_to_remainder
          */

         const hier::Connector::NeighborSet& rejections = (*ci).second;

         remainder_nodes.erase(orig_node);

         hier::BoxList remaining_parts_list(orig_node.getBox());

         for (hier::Connector::NeighborSet::const_iterator
              vi = rejections.begin(); vi != rejections.end(); ++vi) {
            remaining_parts_list.removeIntersections((*vi).getBox());
         }
         /*
          * Coalesce the remaining_parts_list, because it may have unneeded cuts.
          * The coalesce algorithm is O(N^2) or O(N^3), but we expect the
          * length of remaining_parts_list to be very small.
          */
         if (remaining_parts_list.size() > 1) {
            remaining_parts_list.coalesceBoxes();
         }

         /*
          * Create orig_eto_remainder[mapped_box_id] even if remaining_parts_list
          * is empty because its existence defines the required mapping from
          * the orig node to a (possibly empty) container of nesting parts.
          */
         Connector::NeighborSet& remaining_parts_neighbors = orig_eto_remainder[mapped_box_id];
         for (hier::BoxList::Iterator bi(remaining_parts_list);
              bi; bi++) {
            hier::Box new_box = (*bi);
            hier::MappedBox new_node(new_box, ++last_used_index, rank);
            remainder_nodes.insert(remainder_nodes.end(), new_node);
            remaining_parts_neighbors.insert(remaining_parts_neighbors.end(), new_node);
         }
      }

   }

   remainder.initialize(
      remainder_nodes,
      orig.getRefinementRatio(),
      orig.getGridGeometry(),
      orig.getMPI(),
      hier::MappedBoxLevel::DISTRIBUTED);

   orig_to_remainder.initialize(
      orig,
      remainder,
      hier::IntVector::getZero(dim),
      orig_eto_remainder,
      hier::MappedBoxLevel::DISTRIBUTED);
}

/*
 *************************************************************************
 * Add periodic images to a MappedBoxLevel.
 *
 * We add the periodic images by examining real Mapped_boxes in the
 * MappedBoxLevel.  For each real mapped_box, consider all of its
 * possible periodic images and add those that are within the
 * given width of the domain.
 *************************************************************************
 */

void MappedBoxLevelConnectorUtils::addPeriodicImages(
   MappedBoxLevel& mapped_box_level,
   const hier::MappedBoxTree& domain_search_tree,
   const IntVector &threshold_distance) const
{
   const hier::PeriodicShiftCatalog* shift_catalog =
      hier::PeriodicShiftCatalog::getCatalog(mapped_box_level.getDim());

   if (!shift_catalog->isPeriodic()) {
      return; // No-op.
   }

   tbox::Pointer<hier::MappedBoxTree> domain_tree_for_mapped_box_level =
      domain_search_tree.createRefinedTree(mapped_box_level.getRefinementRatio());

   const hier::MappedBoxTree& domain_tree =
      *domain_tree_for_mapped_box_level;

   const hier::IntVector& mapped_box_level_growth = threshold_distance;

   for (hier::RealMappedBoxConstIterator ni(mapped_box_level.getMappedBoxes());
        ni.isValid(); ++ni) {

      const MappedBox& mapped_box = *ni;
      for (int s = 1; s < shift_catalog->getNumberOfShifts(); ++s) {
         PeriodicId id(s);
         const hier::IntVector& try_shift =
            shift_catalog->shiftNumberToShiftDistance(id);
         hier::Box box = mapped_box.getBox();
         box.shift(try_shift);
         box.grow(mapped_box_level_growth);
         if (domain_tree.hasOverlap(box)) {
            mapped_box_level.addPeriodicMappedBox(mapped_box, id);
         }
      }
   }

}

/*
 *************************************************************************
 * Add periodic images to a MappedBoxLevel, and update Connectors that
 * require new edges incident on the additions to the MappedBoxLevel.
 *
 * We add the periodic images by examining real Mapped_boxes in the
 * MappedBoxLevel.  For each real mapped_box, consider all of its
 * possible periodic images and add those that are within the
 * Connector width distance of the domain.  (We are not interested in
 * periodic images so far from the domain that they are never used.)
 *
 * After adding periodic images, we bridge through
 * mapped_box_level<==>anchor<==>anchor so bridge can find the periodic
 * edges.
 *************************************************************************
 */

void MappedBoxLevelConnectorUtils::addPeriodicImagesAndRelationships(
   MappedBoxLevel& mapped_box_level,
   Connector& mapped_box_level_to_anchor,
   Connector& anchor_to_mapped_box_level,
   const hier::MappedBoxTree& domain_search_tree,
   const Connector& anchor_to_anchor) const
{
   hier::OverlapConnectorAlgorithm oca;

   if (d_sanity_check_precond) {
      if (!mapped_box_level_to_anchor.isTransposeOf(anchor_to_mapped_box_level)) {
         TBOX_ERROR(
            "MappedBoxLevelConnectorUtils::addPeriodicImages: non-transposed connector inputs.\n"
            << "mapped_box_level_to_anchor and anchor_to_mapped_box_level\n"
            << "must be mutual transposes.");
      }
      mapped_box_level_to_anchor.assertTransposeCorrectness(
         anchor_to_mapped_box_level);
      if (oca.checkOverlapCorrectness(anchor_to_anchor)) {
         TBOX_ERROR(
            "MappedBoxLevelConnectorUtils::addPeriodicImages: input anchor_to_anchor\n"
            << "Connector failed edge correctness check.");
      }
      if (oca.checkOverlapCorrectness(anchor_to_mapped_box_level, false, true,
             true)) {
         TBOX_ERROR(
            "MappedBoxLevelConnectorUtils::addPeriodicImages: input anchor_to_mapped_box_level\n"
            << "Connector failed edge correctness check.");
      }
      if (oca.checkOverlapCorrectness(mapped_box_level_to_anchor, false, true,
             true)) {
         TBOX_ERROR(
            "MappedBoxLevelConnectorUtils::addPeriodicImages: input mapped_box_level_to_anchor\n"
            << "Connector failed edge correctness check.");
      }
   }
   if (!(anchor_to_anchor.getConnectorWidth() >=
         anchor_to_mapped_box_level.getConnectorWidth())) {
      TBOX_ERROR("MappedBoxLevelConnectorUtils::addPeriodicImages: anchor_to_anchor width\n"
         << anchor_to_anchor.getConnectorWidth() << " is insufficient for\n"
         <<
         "generating periodic edges for anchor_to_mapped_box_level's width of "
         << anchor_to_mapped_box_level.getConnectorWidth() << ".\n");
   }

   const hier::PeriodicShiftCatalog* shift_catalog =
      hier::PeriodicShiftCatalog::getCatalog(anchor_to_anchor.getConnectorWidth(
            ).getDim());

   if (!shift_catalog->isPeriodic()) {
      return; // No-op.
   }

   const MappedBoxLevel& anchor = anchor_to_mapped_box_level.getBase();

   tbox::Pointer<hier::MappedBoxTree> domain_tree_for_mapped_box_level =
      domain_search_tree.createRefinedTree(mapped_box_level.getRefinementRatio());

   {
      /*
       * Add the periodic image mapped_boxes for mapped_box_level.
       *
       * Adding images to a MappedBox without neighbors in the anchor
       * means that this method will not find any edges to the added
       * images.
       */
      mapped_box_level.clearForBoxChanges(false);
      if (0) {
         tbox::perr << "mapped_box_level:\n"
                    << mapped_box_level.format("BEFORE-> ", 3);
      }
      const hier::MappedBoxTree& domain_tree =
         *domain_tree_for_mapped_box_level;

      const hier::IntVector& mapped_box_level_growth =
         mapped_box_level_to_anchor.getConnectorWidth();

      for (hier::RealMappedBoxConstIterator ni(mapped_box_level.getMappedBoxes());
           ni.isValid(); ++ni) {

         const MappedBox& mapped_box = *ni;
         hier::Box grown_box = mapped_box.getBox();
         grown_box.grow(mapped_box_level_growth);
         bool images_added(false);
         for (int s = 1; s < shift_catalog->getNumberOfShifts(); ++s) {
            PeriodicId id(s);
            const hier::IntVector& try_shift =
               shift_catalog->shiftNumberToShiftDistance(id);
            hier::Box box = grown_box;
            box.shift(try_shift);
            if (domain_tree.hasOverlap(box)) {
               mapped_box_level.addPeriodicMappedBox(mapped_box, id);
               images_added = true;
            }
         }
         if (d_sanity_check_precond) {
            if ( images_added &&
                 ( !mapped_box_level_to_anchor.hasNeighborSet(mapped_box.getId()) ||
                   mapped_box_level_to_anchor.getNeighborSet(mapped_box.getId()).empty() ) ) {
               TBOX_WARNING("MappedBoxLevelConnectorUtils::addPeriodicImages: MappedBox " << mapped_box
                            <<"\nhas periodic images in or close to the domain\n"
                            <<"but it does not have any neighbors in the anchor MappedBoxLevel.\n"
                            <<"This will lead to missing neighbors in the output.\n"
                            <<"If post-condition checking is enabled, this will\n"
                            <<"result in an error.\n");
            }
         }

      }
      if (0) {
         tbox::perr << "mapped_box_level:\n" << mapped_box_level.format("AFTER-> ", 3);
      }
   }

   if (0) {
      tbox::plog << "Before bridging for periodic edges:\n"
                 << "anchor_to_anchor:\n" << anchor_to_anchor.format("DBG-> ", 3)
                 << "anchor_to_anchor:\n" << anchor_to_anchor.format("DBG-> ", 3)
                 << "mapped_box_level_to_anchor:\n" << mapped_box_level_to_anchor.format("DBG-> ", 3)
                 << "anchor_to_mapped_box_level:\n" << anchor_to_mapped_box_level.format("DBG-> ", 3);
   }

   Connector new_anchor_to_mapped_box_level;
   Connector new_mapped_box_level_to_anchor;

   hier::IntVector width_limit =
      anchor_to_mapped_box_level.getHeadCoarserFlag() ?
      mapped_box_level_to_anchor.getConnectorWidth() :
      anchor_to_mapped_box_level.getConnectorWidth();

   oca.setSanityCheckMethodPreconditions(d_sanity_check_precond);
   oca.bridge(new_mapped_box_level_to_anchor,
      new_anchor_to_mapped_box_level,
      mapped_box_level_to_anchor,
      anchor_to_anchor,
      anchor_to_anchor,
      anchor_to_mapped_box_level,
      width_limit);
   new_anchor_to_mapped_box_level.eraseEmptyNeighborSets();
   new_mapped_box_level_to_anchor.eraseEmptyNeighborSets();

   if (d_sanity_check_postcond) {
      // Expensive sanity check for consistency.
      size_t err1 = new_anchor_to_mapped_box_level.checkConsistencyWithBase();
      if (err1) {
         tbox::perr << "OverlapConnectorAlgorithm found " << err1
                    << " edge-base consistency errors in\n"
                    <<
         "anchor_to_mapped_box_level after computing periodic images.\n";
      }
      size_t err2 = new_mapped_box_level_to_anchor.checkConsistencyWithBase();
      if (err2) {
         tbox::perr << "OverlapConnectorAlgorithm found " << err2
                    << " edge-base consistency errors in\n"
                    <<
         "mapped_box_level_to_anchor after computing periodic images.\n";
      }
      size_t err3 = new_anchor_to_mapped_box_level.checkConsistencyWithHead();
      if (err3) {
         tbox::perr << "OverlapConnectorAlgorithm found " << err3
                    << " edge-mapped_box consistency errors in\n"
                    <<
         "anchor_to_mapped_box_level after computing periodic images.\n";
      }
      size_t err4 = new_mapped_box_level_to_anchor.checkConsistencyWithHead();
      if (err4) {
         tbox::perr << "OverlapConnectorAlgorithm found " << err4
                    << " edge-mapped_box consistency errors in\n"
                    <<
         "mapped_box_level_to_anchor after computing periodic images.\n";
      }
      if (err1 + err2 + err3 + err4) {
         TBOX_ERROR(
            "OverlapConnectorAlgorithm found consistency errors in\n"
            << "addPeriodicImages\n"
            << "anchor:\n" << anchor.format("ERR-> ", 3)
            << "mapped_box_level:\n" << mapped_box_level.format("ERR-> ", 3)
            << "anchor_to_anchor:\n" << anchor_to_anchor.format("ERR-> ", 3)
            << "new_anchor_to_mapped_box_level:\n" << new_anchor_to_mapped_box_level.format("ERR-> ", 3)
            << "new_mapped_box_level_to_anchor:\n" << new_mapped_box_level_to_anchor.format("ERR-> ", 3));
      }
   }
   if (d_sanity_check_postcond) {
      // Expensive sanity check for correctness.
      size_t err1 = oca.checkOverlapCorrectness(new_anchor_to_mapped_box_level);
      if (err1) {
         tbox::perr << "MappedBoxLevelConnectorUtils::addPeriodicImages found " << err1
                    << " errors\n"
                    << "in new_anchor_to_mapped_box_level after\n"
                    << "computing periodic images.  If you enabled\n"
                    << "precondition checking, this is probably a\n"
                    << "library error.\n";
      }
      size_t err2 = oca.checkOverlapCorrectness(new_mapped_box_level_to_anchor);
      if (err2) {
         tbox::perr << "MappedBoxLevelConnectorUtils::addPeriodicImages found " << err2
                    << " errors\n"
                    << "in new_mapped_box_level_to_anchor after\n"
                    << "computing periodic images.  If you enabled\n"
                    << "precondition checking, this is probably a\n"
                    << "library error.\n";
      }
      if (err1 + err2) {
         TBOX_ERROR(
            "MappedBoxLevelConnectorUtils::addPeriodicImages found edge errors\n"
            << "in output data\n"
            << "anchor:\n" << anchor.format("ERR-> ", 3)
            << "mapped_box_level:\n" << mapped_box_level.format("ERR-> ", 3)
            << "anchor_to_anchor:\n" << anchor_to_anchor.format("ERR-> ", 3)
            << "anchor_to_mapped_box_level:\n" << anchor_to_mapped_box_level.format("ERR-> ", 3)
            << "mapped_box_level_to_anchor:\n" << mapped_box_level_to_anchor.format("ERR-> ", 3)
            << "new_anchor_to_mapped_box_level:\n" << new_anchor_to_mapped_box_level.format("ERR-> ", 3)
            << "new_mapped_box_level_to_anchor:\n" << new_mapped_box_level_to_anchor.format("ERR-> ", 3));
      }
   }

   Connector::swap(new_anchor_to_mapped_box_level, anchor_to_mapped_box_level);
   Connector::swap(new_mapped_box_level_to_anchor, mapped_box_level_to_anchor);
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void MappedBoxLevelConnectorUtils::initializeCallback()
{
   t_make_sorting_map = tbox::TimerManager::getManager()->
      getTimer("hier::MappedBoxLevelConnectorUtils::makeSortingMap()");
   t_compute_external_parts = tbox::TimerManager::getManager()->
      getTimer("hier::MappedBoxLevelConnectorUtils::computeExternalParts()");
   t_compute_external_parts_intersection =
      tbox::TimerManager::getManager()->
      getTimer("mesh::MappedBoxLevelConnectorUtils::computeExternalParts()_intersection");
   t_compute_internal_parts = tbox::TimerManager::getManager()->
      getTimer("hier::MappedBoxLevelConnectorUtils::computeInternalParts()");
   t_compute_internal_parts_intersection =
      tbox::TimerManager::getManager()->
      getTimer("mesh::MappedBoxLevelConnectorUtils::computeInternalParts()_intersection");
}

/*
 ***************************************************************************
 * * Release static timers.  To be called by shutdown registry to make sure  *
 * * memory for timers does not leak.                                        *
 ***************************************************************************
 */

void MappedBoxLevelConnectorUtils::finalizeCallback()
{
   t_make_sorting_map.setNull();
   t_compute_external_parts.setNull();
   t_compute_external_parts_intersection.setNull();
   t_compute_internal_parts.setNull();
   t_compute_internal_parts_intersection.setNull();
}

}
}
#endif
