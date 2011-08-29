/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Algorithms for working with mapping Connectors.
 *
 ************************************************************************/
#ifndef included_hier_MappingConnectorAlgorithm_C
#define included_hier_MappingConnectorAlgorithm_C

#include "SAMRAI/hier/MappingConnectorAlgorithm.h"
#include "SAMRAI/hier/OverlapConnectorAlgorithm.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/hier/RealBoxConstIterator.h"
#include "SAMRAI/tbox/AsyncCommStage.h"
#include "SAMRAI/tbox/AsyncCommPeer.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"
#include "SAMRAI/tbox/TimerManager.h"

#include <algorithm>

namespace SAMRAI {
namespace hier {

const int MappingConnectorAlgorithm::MAPPING_CONNECTOR_ALGORITHM_FIRST_DATA_LENGTH = 1000;

char MappingConnectorAlgorithm::s_print_modify_steps = '\0';

tbox::Pointer<tbox::Timer> MappingConnectorAlgorithm::t_modify;
tbox::Pointer<tbox::Timer> MappingConnectorAlgorithm::t_modify_shortcut;
tbox::Pointer<tbox::Timer> MappingConnectorAlgorithm::t_modify_setup_comm;
tbox::Pointer<tbox::Timer> MappingConnectorAlgorithm::t_modify_remove_and_cache;
tbox::Pointer<tbox::Timer> MappingConnectorAlgorithm::t_modify_discover_and_send;
tbox::Pointer<tbox::Timer> MappingConnectorAlgorithm::t_modify_receive_and_unpack;
tbox::Pointer<tbox::Timer> MappingConnectorAlgorithm::t_modify_MPI_wait;

const std::string MappingConnectorAlgorithm::s_dbgbord;

int MappingConnectorAlgorithm::s_operation_mpi_tag = 0;
/*
 * Do we even need to use different tags each time we modify???
 * Unique tags were used to help debug, but the methods may work
 * with reused tags anyway.
 */

tbox::SAMRAI_MPI MappingConnectorAlgorithm::s_class_mpi(tbox::SAMRAI_MPI::commNull);

tbox::StartupShutdownManager::Handler
MappingConnectorAlgorithm::s_initialize_finalize_handler(
   MappingConnectorAlgorithm::initializeCallback,
   0,
   0,
   MappingConnectorAlgorithm::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);

/*
 ***********************************************************************
 ***********************************************************************
 */

MappingConnectorAlgorithm::MappingConnectorAlgorithm():
   d_sanity_check_inputs(false),
   d_sanity_check_outputs(false),
   d_shortcut_trivial_maps(true)
{
   /*
    * While we figure out how to use multiple communicators in SAMRAI,
    * we are still assuming that all communications use congruent
    * communicators.  This class just makes a duplicate communicator
    * to protect itself from unrelated communications in shared
    * communicators.
    */
   if (s_class_mpi.getCommunicator() == tbox::SAMRAI_MPI::commNull) {
      if (tbox::SAMRAI_MPI::usingMPI()) {
         const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
         s_class_mpi.dupCommunicator(mpi);
      }
   }

}

/*
 ***********************************************************************
 ***********************************************************************
 */

MappingConnectorAlgorithm::~MappingConnectorAlgorithm()
{
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void MappingConnectorAlgorithm::setSanityCheckMethodPreconditions(
   bool do_check)
{
   d_sanity_check_inputs = do_check;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void MappingConnectorAlgorithm::setSanityCheckMethodPostconditions(
   bool do_check)
{
   d_sanity_check_outputs = do_check;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void MappingConnectorAlgorithm::shortcutTrivialMapsInModify(
   bool do_shortcut)
{
   d_shortcut_trivial_maps = do_shortcut;
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void MappingConnectorAlgorithm::modify(
   Connector& anchor_to_mapped,
   Connector& mapped_to_anchor,
   const Connector& old_to_new,
   const Connector& new_to_old,
   BoxLevel* mutable_new,
   BoxLevel* mutable_old) const
{

   /*
    * Ensure that Connectors incident to and from old agree on
    * what the old mapped_box_level is.
    */
   const Connector& anchor_to_old = anchor_to_mapped;
   const Connector& old_to_anchor = mapped_to_anchor;

   const BoxLevel* old = &old_to_new.getBase();

   if (old != &new_to_old.getHead() ||
       old != &anchor_to_mapped.getHead() ||
       old != &mapped_to_anchor.getBase()) {
      TBOX_ERROR("Bad input for MappingConnectorAlgorithm::modify:\n"
         << "Given Connectors to base and head of modify are not incident\n"
         << "from the same old in MappingConnectorAlgorithm::modify:\n"
         << "anchor_to_old is  TO  " << &anchor_to_old.getHead() << "\n"
         << "old_to_new is FROM " << &old_to_new.getBase()
         << "\n"
         << "new_to_old is  TO  " << &new_to_old.getHead()
         << "\n"
         << "old_to_anchor is FROM " << &old_to_anchor.getBase() << "\n"
         );
   }
   /*
    * Ensure that head and base mapped_box_levels in argument agree with
    * head and base in the object.
    */
   if (&anchor_to_old.getBase() != &old_to_anchor.getHead()) {
      TBOX_ERROR("Bad input for MappingConnectorAlgorithm::modify:\n"
         << "Given Connectors to and from base of modify do not refer\n"
         << "to the base of the modify in:\n"
         << "anchor_to_old is FROM " << &anchor_to_old.getBase() << "\n"
         << "old_to_anchor is  TO  " << &old_to_anchor.getHead() << "\n"
         );
   }
   if (&old_to_new.getHead() != &new_to_old.getBase()) {
      TBOX_ERROR("Bad input for MappingConnectorAlgorithm::modify:\n"
         << "Given Connectors to and from head of modify do not refer\n"
         << "to the head of the modify in MappingConnectorAlgorithm::modify:\n"
         << "new_to_old is FROM " << &new_to_old.getBase()
         << "\n"
         << "old_to_new is  TO  " << &old_to_new.getHead()
         << "\n"
         );
   }
   if (!anchor_to_old.isTransposeOf(old_to_anchor)) {
      TBOX_ERROR("Bad input for MappingConnectorAlgorithm::modify:\n"
         << "Given Connectors between base and old of modify\n"
         << "are not transposes of each other.\n"
         << "See MappingConnectorAlgorithm::isTransposeOf().\n"
         );
   }
   if (!new_to_old.isTransposeOf(old_to_new)) {
      TBOX_ERROR("Bad input for MappingConnectorAlgorithm::modify:\n"
         << "Given Connectors between head and old of modify\n"
         << "are not transposes of each other.\n"
         << "See MappingConnectorAlgorithm::isTransposeOf().\n"
         );
   }
   if (anchor_to_old.getParallelState() != BoxLevel::DISTRIBUTED) {
      TBOX_ERROR("Bad input for MappingConnectorAlgorithm::modify:\n"
         << "bridging is currently set up for DISTRIBUTED\n"
         << "mode only.\n");
   }

   if (s_print_modify_steps == 'y') {
      tbox::plog
      << "MappingConnectorAlgorithm::modify: old mapped_box_level:\n"
      << old_to_new.getBase().format(s_dbgbord, 2)
      << "MappingConnectorAlgorithm::modify: new mapped_box_level:\n"
      << new_to_old.getBase().format(s_dbgbord, 2)
      << "MappingConnectorAlgorithm::modify: old_to_new:\n"
      << old_to_new.format(s_dbgbord, 2)
      << "MappingConnectorAlgorithm::modify: new_to_old:\n"
      << new_to_old.format(s_dbgbord, 2);
   }

   if (0) {
      // Expensive input checking.
      const BoxLevel& anchor_mapped_box_level = anchor_to_old.getBase();
      const BoxLevel& old_mapped_box_level = old_to_new.getBase();

      tbox::plog
      << "anchor mapped_box_level:\n" << anchor_mapped_box_level.format(s_dbgbord, 2)
      << "anchor_to_old:\n" << anchor_to_old.format(s_dbgbord, 2)
      << "old mapped_box_level:\n" << old_mapped_box_level.format(s_dbgbord, 2)
      << "old_to_new:\n" << old_to_new.format(s_dbgbord, 2)
      << "new_to_old:\n" << new_to_old.format(s_dbgbord, 2);

      hier::OverlapConnectorAlgorithm oca;
      TBOX_ASSERT(oca.checkOverlapCorrectness(anchor_to_old) == 0);
      TBOX_ASSERT(oca.checkOverlapCorrectness(old_to_anchor) == 0);
      TBOX_ASSERT(old_to_anchor.checkTransposeCorrectness(anchor_to_old,
            true) == 0);
      TBOX_ASSERT(oca.checkOverlapCorrectness(old_to_new, true, false) == 0);
      TBOX_ASSERT(oca.checkOverlapCorrectness(new_to_old, true, false) == 0);
      TBOX_ASSERT(old_to_new.checkTransposeCorrectness(new_to_old,
            true) == 0);
   }

   privateModify(anchor_to_mapped,
      mapped_to_anchor,
      old_to_new,
      new_to_old,
      mutable_new,
      mutable_old);

   if (d_sanity_check_outputs) {
      anchor_to_mapped.assertTransposeCorrectness(mapped_to_anchor);
      mapped_to_anchor.assertTransposeCorrectness(anchor_to_mapped);
   }
}

/*
 *****************************************************************************
 * Version of modify requiring only the forward map
 * and allows only local mappings.
 *
 * This version modifies two transpose Connectors.
 *****************************************************************************
 */

void MappingConnectorAlgorithm::modify(
   Connector& anchor_to_mapped,
   Connector& mapped_to_anchor,
   const Connector& old_to_new,
   BoxLevel* mutable_new,
   BoxLevel* mutable_old) const
{

#ifdef DEBUG_CHECK_ASSERTIONS
   /*
    * Ensure that mapping has only local neighbors.
    */
   const NeighborhoodSet& old_eto_new = old_to_new.getNeighborhoodSets();
   for (NeighborhoodSet::const_iterator ci = old_eto_new.begin();
        ci != old_eto_new.end(); ++ci) {
      const NeighborSet& new_nabrs = (*ci).second;
      for (NeighborSet::const_iterator na = new_nabrs.begin();
           na != new_nabrs.end(); ++na) {
         if ((*na).getOwnerRank() != old_to_new.getRank()) {
            const Box& mapped_box(
               *old_to_new.getBase().getBoxes().
               find(Box(na->getDim(), (*ci).first)));
            TBOX_ERROR("MappingConnectorAlgorithm::modify: this version of modify\n"
               "only allows local mappings.  The local mapped_box\n"
               << mapped_box << " has a non-local map to\n"
               << *na << "\n"
               << "To modify using non-local maps, the\n"
               << "reverse mapping must be provide and\n"
               << "the four-argument version of modify\n"
               << "must be used.");
         }
      }
   }
#endif

   const Connector& anchor_to_old = anchor_to_mapped;
   const Connector& old_to_anchor = mapped_to_anchor;

   const BoxLevel* old = &old_to_new.getBase();

   /*
    * Ensure that Connectors incident to and from old agree on
    * what the old mapped_box_level is.
    */

   if (old != &old_to_new.getBase() ||
       old != &anchor_to_old.getHead()) {
      TBOX_ERROR("Bad input for MappingConnectorAlgorithm::modify:\n"
         << "Given Connectors to base and head of modify are not incident\n"
         << "from the same old in MappingConnectorAlgorithm::modify:\n"
         << "anchor_to_old is  TO   " << &anchor_to_old.getHead() << "\n"
         << "old_to_anchor is FROM  " << &old_to_anchor.getBase() << "\n"
         << "old_to_new is FROM " << &old_to_new.getBase()
         << "\n"
         );
   }
   /*
    * Ensure that head and base mapped_box_levels in argument agree with
    * head and base in the object.
    */
   if (&anchor_to_old.getBase() != &old_to_anchor.getHead()) {
      TBOX_ERROR("Bad input for MappingConnectorAlgorithm::modify:\n"
         << "Given Connectors to and from base of modify do not refer\n"
         << "to the base of the modify in MappingConnectorAlgorithm::modify:\n"
         << "anchor_to_old is FROM " << &anchor_to_old.getBase() << "\n"
         << "old_to_anchor is  TO  " << &old_to_anchor.getHead() << "\n"
         );
   }
   if (!anchor_to_old.isTransposeOf(old_to_anchor)) {
      TBOX_ERROR("Bad input for MappingConnectorAlgorithm::modify:\n"
         << "Given Connectors between base and old of modify\n"
         << "are not transposes of each other.\n"
         << "See MappingConnectorAlgorithm::isTransposeOf().\n"
         );
   }
   if (anchor_to_old.getParallelState() != BoxLevel::DISTRIBUTED) {
      TBOX_ERROR("Bad input for MappingConnectorAlgorithm::modify:\n"
         << "modifying is currently set up for DISTRIBUTED\n"
         << "mode only.\n");
   }

   if (s_print_modify_steps == 'y') {
      const BoxLevel& anchor_mapped_box_level = anchor_to_mapped.getBase();
      const BoxLevel& old_mapped_box_level = old_to_new.getBase();

      tbox::plog
      << "MappingConnectorAlgorithm::modify: anchor mapped_box_level:\n"
      << anchor_mapped_box_level.format(s_dbgbord, 2)
      << "MappingConnectorAlgorithm::modify: anchor_to_old:\n"
      << anchor_to_old.format(s_dbgbord, 2)
      << "MappingConnectorAlgorithm::modify: old mapped_box_level:\n"
      << old_mapped_box_level.format(s_dbgbord, 2)
      << "MappingConnectorAlgorithm::modify: old_to_new:\n"
      << old_to_new.format(s_dbgbord, 2)
      << "MappingConnectorAlgorithm::modify: No new_to_old.\n";
   }

   Connector dummy_new_to_old;
   privateModify(anchor_to_mapped,
      mapped_to_anchor,
      old_to_new,
      dummy_new_to_old,
      mutable_new,
      mutable_old);
}

/*
 *****************************************************************************
 * Version of modify requiring only the forward map
 * and allows only local mappings.
 *
 * This version modifies just one Connector.
 *****************************************************************************
 */

void MappingConnectorAlgorithm::modify(
   Connector& anchor_to_mapped,
   const Connector& old_to_new,
   BoxLevel* mutable_new,
   BoxLevel* mutable_old) const
{
   /*
    * Ensure that connectors incident to and from old agree on
    * what the old mapped_box_level is.
    */
   if (&anchor_to_mapped.getHead() != &old_to_new.getBase()) {
      TBOX_ERROR(
         "Bad input for MappingConnectorAlgorithm::modify:\n"
         << "Input MappingConnectorAlgorithm do not agree on what the old mapped_box_level is.\n"
         << "anchor_to_mapped is  TO  " << &anchor_to_mapped.getHead() << "\n"
         << "old_to_new is FROM " << &old_to_new.getBase()
         << "\n"
         );
   }

   if (anchor_to_mapped.getParallelState() != BoxLevel::DISTRIBUTED) {
      TBOX_ERROR("Bad input for MappingConnectorAlgorithm::modify:\n"
         << "modifying is currently set up for DISTRIBUTED\n"
         << "mode only.\n");
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   /*
    * Ensure that mapping has only local neighbors.
    */
   const NeighborhoodSet& old_eto_new = old_to_new.getNeighborhoodSets();
   for (NeighborhoodSet::const_iterator ci = old_eto_new.begin();
        ci != old_eto_new.end(); ++ci) {
      const NeighborSet& new_nabrs = (*ci).second;
      for (NeighborSet::const_iterator na = new_nabrs.begin();
           na != new_nabrs.end(); ++na) {
         if ((*na).getOwnerRank() != old_to_new.getRank()) {
            const Box& mapped_box(
               *old_to_new.getBase().getBoxes().
               find(Box(na->getDim(), (*ci).first)));
            TBOX_ERROR("MappingConnectorAlgorithm::modify: this version of modify\n"
               "only allows local mappings.  The local mapped_box\n"
               << mapped_box << " has a non-local map to\n"
               << *na << "\n"
               << "To modify using non-local maps, the\n"
               << "reverse mapping must be provide and\n"
               << "the four-argument version of modify\n"
               << "must be used.");
         }
      }
   }
#endif

   /*
    * Generate temporary mapped_to_anchor needed by the modify
    * operation.  This step is the reason anchor_to_mapped
    * may not have non-local neighbors.
    *
    * No need to check that anchor_to_mapped is strictly local,
    * because initializeToLocalTranspose checks that.
    */
   Connector mapped_to_anchor;
   mapped_to_anchor.initializeToLocalTranspose(anchor_to_mapped);

   if (s_print_modify_steps == 'y') {
      const BoxLevel& anchor_mapped_box_level = anchor_to_mapped.getBase();
      const BoxLevel& old_mapped_box_level = old_to_new.getBase();
      const BoxLevel& new_mapped_box_level =
         old_to_new.getHead();

      tbox::plog
      << "MappingConnectorAlgorithm::modify: anchor mapped_box_level:\n"
      << anchor_mapped_box_level.format(s_dbgbord, 2)
      << "MappingConnectorAlgorithm::modify: anchor_to_old:\n"
      << anchor_to_mapped.format(s_dbgbord, 2)
      << "MappingConnectorAlgorithm::modify: old mapped_box_level:\n"
      << old_mapped_box_level.format(s_dbgbord, 2)
      << "MappingConnectorAlgorithm::modify: old_to_new:\n"
      << old_to_new.format(s_dbgbord, 2)
      << "MappingConnectorAlgorithm::modify: new mapped_box_level:\n"
      << new_mapped_box_level.format(s_dbgbord, 2);
   }

   Connector dummy_new_to_old;
   privateModify(anchor_to_mapped,
      mapped_to_anchor,
      old_to_new,
      dummy_new_to_old,
      mutable_new,
      mutable_old);
}

/*
 ***********************************************************************
 * This method modifies overlap Connectors based on the described
 * changes to their base or head BoxLevels.  The change is
 * described as a mapping from the old state to the new.
 *
 * Essential nomenclature:
 * - mapped: The BoxLevel that is being changed.
 * - anchor: A BoxLevel that is NOT being changed.
 *   This method modifies the overlap Connectors between anchor
 *   and mapped.
 * - old: The state of mapped before the change.
 * - new: the state of mapped after the change.
 * - old_to_new: Desription of the change.  The NeighborSet of
 *   an old Box is what the old Box will become.
 *   By convention, the NeighborSet of a un-changing Box
 *   is not required.  However, an empty NeighborSet means that
 *   the old Box will disappear.
 *
 * While modify adjusts the Connector to reflect changes in the mapped
 * BoxLevel, it does NOT manipulate any BoxLevel objects
 * (other than initializing one with another; see in-place changes).
 *
 * The Connector width of the mapping must be at least equal to the
 * ammount that new boxes protrude from their old boxes.  The protusion
 * may generate undetected overlaps between anchor and mapped.  To
 * avoid generating undetected overlaps, the width of anchor<==>mapped
 * are shrunken by the width of the mapping.
 *
 * At this time, new_to_old is only used to determine the remote
 * owners that must be notified of a local mapped_box being mapped.
 * The Connector object contains more info than required, so it can be
 * replaced by something more concise.  However, requiring a transpose
 * Connector lets us use some canned sanity checks.  The use of the
 * transpose Connector is fairly efficient.
 ***********************************************************************
 */

void MappingConnectorAlgorithm::privateModify(
   Connector& anchor_to_mapped,
   Connector& mapped_to_anchor,
   const Connector& old_to_new,
   const Connector& new_to_old,
   BoxLevel* mutable_new,
   BoxLevel* mutable_old) const
{
   t_modify->barrierAndStart();

   if (s_print_modify_steps == 'y') {
      tbox::plog
      << "MappingConnectorAlgorithm::privateModify: anchor_to_old:\n"
      << anchor_to_mapped.format(s_dbgbord, 3)
      << "MappingConnectorAlgorithm::privateModify: mapped_to_anchor:\n"
      << mapped_to_anchor.format(s_dbgbord, 3)
      << "MappingConnectorAlgorithm::privateModify: old_to_new:\n"
      << old_to_new.format(s_dbgbord, 3)
      << "MappingConnectorAlgorithm::privateModify: new_to_old:\n"
      << old_to_new.format(s_dbgbord, 3);
   }

   checkModifyParameters(
      anchor_to_mapped,
      mapped_to_anchor,
      old_to_new,
      new_to_old);

   /*
    * anchor<==>mapped start out as Connectors between
    * old and anchor.  Make copies of Connectors to and from old
    * so we can modify anchor_to_mapped/mapped_to_anchor without
    * losing needed data.
    */
   const Connector anchor_to_old = anchor_to_mapped;
   const Connector old_to_anchor = mapped_to_anchor;

   /*
    * We will modify the mapped BoxLevel to make it the new
    * BoxLevel.  "mapped" is the name for the changing
    * BoxLevel.  Because this BoxLevel will be
    * modified for output, it is synonymous with "new".
    */
   Connector& anchor_to_new = anchor_to_mapped;
   Connector& new_to_anchor = mapped_to_anchor;

   /*
    * Shorthand for the three BoxLevels and their Refinement
    * ratios.
    */
   const BoxLevel& old = old_to_anchor.getBase();
   const BoxLevel& new_mapped_box_level = old_to_new.getHead();
   const BoxLevel& anchor = anchor_to_new.getBase();
   const IntVector& old_ratio = old.getRefinementRatio();
   const IntVector& anchor_ratio = anchor.getRefinementRatio();
   const IntVector& new_ratio = new_mapped_box_level.getRefinementRatio();

   const tbox::Dimension dim(old.getDim());
   const tbox::SAMRAI_MPI& mpi(old.getMPI());

   /*
    * The width of old-->new indicates the maximum amount of box
    * growth caused by the change.  A value of zero means no growth.
    *
    * Growing old Boxes may generate new overlaps that we cannot
    * detect because they lie outside the current anchor<==>mapped
    * widths.  To reflect that we do not see any farther just because
    * the Boxes have grown, we shrink the widths by (nominally)
    * the amount of growth.  To ensure the shrinkage is consistent
    * between transpose pairs of Connectors, it is converted to the
    * coarser index space then converted into the final index space.
    * Calls to Connector::convertHeadWidthToBase() perform the
    * conversions.
    */
   const IntVector shrinkage_in_new_index_space =
      Connector::convertHeadWidthToBase(
         old_ratio,
         new_ratio,
         old_to_new.getConnectorWidth());
   const IntVector shrinkage_in_anchor_index_space =
      Connector::convertHeadWidthToBase(
         anchor_ratio,
         new_ratio,
         shrinkage_in_new_index_space);
   const IntVector anchor_to_new_width =
      anchor_to_old.getConnectorWidth() - shrinkage_in_anchor_index_space;
   if (!(anchor_to_new_width >= hier::IntVector::getZero(dim))) {
      TBOX_ERROR(
         "MappingConnectorAlgorithm::privateModify error:\n"
         << "Mapping connector allows mapped BoxLevel to grow\n"
         << "too much.  The growth may generate new overlaps\n"
         << "that cannot be detected by mapping algorithm, thus\n"
         << "causing output overlap Connectors to be incomplete.\n"
         << "Mapping Connector:\n" << old_to_new.format("", 0)
         << "anchor--->mapped:\n" << anchor_to_new.format("", 0)
         << "Connector width of anchor--->mapped will shrink\n"
         << "by " << shrinkage_in_anchor_index_space << " which\n"
         << "will result in a non-positive width.");
   }
   const IntVector new_to_anchor_width =
      Connector::convertHeadWidthToBase(
         new_mapped_box_level.getRefinementRatio(),
         anchor.getRefinementRatio(),
         anchor_to_new_width);

   bool do_shortcut(false);
   if (d_shortcut_trivial_maps) {
      t_modify_shortcut->start();
      do_shortcut = (old_to_new.getGlobalNumberOfNeighborSets() == 0);
      t_modify_shortcut->stop();
   }

   if (do_shortcut) {

      t_modify_shortcut->start();
      /*
       * Shortcut for trivial mapping.
       *
       * Initialize the output Connectors to connect to the new
       * BoxLevel.  Shrink the ghost widths to those computed
       * above (if needed).
       */

      if (d_sanity_check_inputs) {
         TBOX_ASSERT(old_to_new.getBase() == old_to_new.getHead());
      }

      anchor_to_new.initialize(
         anchor_to_new.getBase(),
         old_to_new.getHead(),
         anchor_to_new.getConnectorWidth(),
         anchor_to_new.getNeighborhoodSets(),
         BoxLevel::DISTRIBUTED);
      new_to_anchor.initialize(
         old_to_new.getHead(),
         anchor_to_new.getBase(),
         new_to_anchor.getConnectorWidth(),
         new_to_anchor.getNeighborhoodSets(),
         BoxLevel::DISTRIBUTED);

      hier::OverlapConnectorAlgorithm oca;
      if (anchor_to_new_width != anchor_to_new.getConnectorWidth()) {
         oca.shrinkConnectorWidth(anchor_to_new, anchor_to_new_width);
      }
      if (new_to_anchor_width != new_to_anchor.getConnectorWidth()) {
         oca.shrinkConnectorWidth(new_to_anchor, new_to_anchor_width);
      }

      t_modify_shortcut->stop();

   } else {

      /*
       * The essential modify algorithm is in this block.
       */

      /*
       * Initialize the output connectors with the correct new
       * BoxLevel.  (As inputs, they were referencing the old
       * BoxLevel.)  At the same time, swap out the neighbor data
       * (into new_eto_anchor and anchor_eto_new) for direct
       * modification.  We'll swap these back in after the mods.
       */
      NeighborhoodSet new_eto_anchor;
      NeighborhoodSet anchor_eto_new;
      anchor_to_new.swapInitialize(
         anchor_to_old.getBase(),
         old_to_new.getHead(),
         anchor_to_new_width,
         anchor_eto_new,
         BoxLevel::DISTRIBUTED);
      new_to_anchor.swapInitialize(
         anchor_to_new.getHead(),
         anchor_to_old.getBase(),
         new_to_anchor_width,
         new_eto_anchor,
         BoxLevel::DISTRIBUTED);

      /*
       * Determine which ranks we have to communicate with.  They are
       * the ones who owns Boxes that the local Boxes will
       * be mapped to.
       */
      std::set<int> incoming_ranks, outgoing_ranks;
      old_to_new.getNeighborhoodSets().getOwners(incoming_ranks);
      old_to_anchor.getNeighborhoodSets().getOwners(incoming_ranks);

      anchor_to_old.getNeighborhoodSets().getOwners(outgoing_ranks);
      if (new_to_old.isInitialized()) {
         new_to_old.getNeighborhoodSets().getOwners(outgoing_ranks);
      }

      // We don't need to communicate locally.
      incoming_ranks.erase(mpi.getRank());
      outgoing_ranks.erase(mpi.getRank());

      /*
       * Object for communicating relationship changes.
       */
      tbox::AsyncCommStage comm_stage;
      tbox::AsyncCommPeer<int> * all_comms(NULL);
      tbox::AsyncCommStage::MemberVec completed;

      /*
       * Set up communication mechanism (and post receives).
       */
      privateModify_setupCommunication(
         all_comms,
         comm_stage,
         completed,
         anchor.getMPI(),
         incoming_ranks,
         outgoing_ranks);

      /*
       * There are three major parts to computing the new neighbors:
       * (1) remove relationships for Boxes being mapped into
       * something else, (2) discover new relationships for the new
       * Boxes and (3) share information with other processes.
       *
       * In steps 1 and 2, the owners of the Boxes being
       * modified determine which relationships to remove and which to
       * add.  Some of this information is kept locally while the rest
       * are to be sent to the owners of the affected anchor and new
       * Boxes.
       *
       * The three parts are done in the 3 steps following.  Note that
       * these steps do not correspond to the parts.
       */

      /*
       * Data for caching relationship removal messages.  This
       * temporary object holds data computed when removing neighbors,
       * to be used when sending out the relationship removal
       * information.
       */
      std::map<int, std::vector<int> > neighbor_removal_mesg;

      /*
       * First step: Remove neighbor data for Boxes that are
       * going away and cache information to be sent out.
       */
      privateModify_removeAndCache(
         neighbor_removal_mesg,
         anchor_eto_new,
         new_eto_anchor,
         old_to_new);

      /*
       * Second step: Discover overlaps for new Boxes.  Send
       * messages with information about what relationships to remove
       * or create.
       */
      privateModify_discoverAndSend(
         neighbor_removal_mesg,
         anchor_to_new,
         new_to_anchor,
         anchor_eto_new,
         new_eto_anchor,
         incoming_ranks,
         outgoing_ranks,
         all_comms,
         completed,
         old_to_anchor,
         anchor_to_old,
         old_to_new);

      /*
       * Third step: Receive and unpack messages from incoming_ranks.
       * These message contain information about what relationships to
       * remove or add.
       */
      privateModify_receiveAndUnpack(
         anchor_eto_new,
         new_eto_anchor,
         outgoing_ranks,
         all_comms,
         comm_stage,
         completed,
         dim,
         mpi);

      delete[] all_comms;

      TBOX_ASSERT(anchor_to_new.getNeighborhoodSets().empty());
      TBOX_ASSERT(new_to_anchor.getNeighborhoodSets().empty());

      /*
       * We have finished manually changing the NeighborhoodSets
       * anchor_eto_new and new_eto_anchor.  Swap them back into the
       * Connectors.
       */
      anchor_to_new.swapInitialize(
         anchor_to_new.getBase(),
         anchor_to_new.getHead(),
         anchor_to_new.getConnectorWidth(),
         anchor_eto_new,
         BoxLevel::DISTRIBUTED);
      new_to_anchor.swapInitialize(
         new_to_anchor.getBase(),
         new_to_anchor.getHead(),
         new_to_anchor.getConnectorWidth(),
         new_eto_anchor,
         BoxLevel::DISTRIBUTED);

      if (!anchor_to_new.ratioIsExact()) {
         TBOX_WARNING("MappingConnectorAlgorithm::privateModify: generated\n"
            << "overlap Connectors with non-integer ratio between\n"
            << "the base and head.  The results are not guaranteed\n"
            << "to be complete overlap Connectors.");
      }

      TBOX_ASSERT(anchor_eto_new.empty());
      TBOX_ASSERT(new_eto_anchor.empty());

   }

   /*
    * Optional in-place changes:
    *
    * Note that the old and new BoxLevels gotten from Connectors
    * are const, so this method cannot modify them.  Only by specifying
    * the mutable BoxLevels, can this method modify them.
    *
    * If users provide mutable object to initialize to the new
    * BoxLevel, this method initializes it to the new
    * BoxLevel and uses it in the output Connectors.
    */
   if (mutable_new == &old_to_new.getBase() &&
       mutable_old == &old_to_new.getHead()) {
      /*
       * Since mutable_new is old and mutable_old is new, shortcut
       * two assignments by swapping.
       */
      BoxLevel::swap(*mutable_new, *mutable_old);
      new_to_anchor.initialize(
         *mutable_new,
         new_to_anchor.getHead(),
         new_to_anchor.getConnectorWidth(),
         new_to_anchor.getNeighborhoodSets(),
         BoxLevel::DISTRIBUTED);
      anchor_to_new.initialize(
         anchor_to_new.getBase(),
         *mutable_new,
         anchor_to_new.getConnectorWidth(),
         anchor_to_new.getNeighborhoodSets(),
         BoxLevel::DISTRIBUTED);
   } else {
      if (mutable_new != NULL) {
         *mutable_new = old_to_new.getHead();
         new_to_anchor.initialize(
            *mutable_new,
            new_to_anchor.getHead(),
            new_to_anchor.getConnectorWidth(),
            new_to_anchor.getNeighborhoodSets(),
            BoxLevel::DISTRIBUTED);
         anchor_to_new.initialize(
            anchor_to_new.getBase(),
            *mutable_new,
            anchor_to_new.getConnectorWidth(),
            anchor_to_new.getNeighborhoodSets(),
            BoxLevel::DISTRIBUTED);
      }
      if (mutable_old != NULL) {
         *mutable_old = old_to_new.getBase();
      }
   }

   t_modify->barrierAndStop();
}

/*
 ***********************************************************************
 * Receive messages and unpack info sent from other processes.
 ***********************************************************************
 */
void MappingConnectorAlgorithm::privateModify_setupCommunication(
   tbox::AsyncCommPeer<int> *& all_comms,
   tbox::AsyncCommStage& comm_stage,
   tbox::AsyncCommStage::MemberVec& completed,
   const tbox::SAMRAI_MPI& mpi,
   const std::set<int>& incoming_ranks,
   const std::set<int>& outgoing_ranks) const
{
   t_modify_setup_comm->start();

   /*
    * Set up communication mechanism (and post receives).  We lump all
    * communication objects into one array, all_comms.  all_comms is
    * ordered with the incoming first and the outgoing afterward.
    */
   comm_stage.setCommunicationWaitTimer(t_modify_MPI_wait);
   const int n_comm = static_cast<int>(
         outgoing_ranks.size() + incoming_ranks.size());
   all_comms = new tbox::AsyncCommPeer<int>[n_comm];

   completed.clear();

   const int tag0 = ++s_operation_mpi_tag;
   const int tag1 = ++s_operation_mpi_tag;

   std::set<int>::const_iterator owneri;
   size_t peer_idx = 0;
   for (owneri = outgoing_ranks.begin(); owneri != outgoing_ranks.end(); ++owneri,
        ++peer_idx) {
      const int peer_rank = *owneri;
      tbox::AsyncCommPeer<int>& incoming_comm = all_comms[peer_idx];
      incoming_comm.initialize(&comm_stage);
      incoming_comm.setPeerRank(peer_rank);
      incoming_comm.setMPI(mpi);
      incoming_comm.setMPITag(tag0, tag1);
      incoming_comm.limitFirstDataLength(MAPPING_CONNECTOR_ALGORITHM_FIRST_DATA_LENGTH);
      incoming_comm.beginRecv();
      if (s_print_modify_steps == 'y') tbox::plog << "Post recv to "
                                                  << incoming_comm.getPeerRank()
                                                  << std::endl;
      if (incoming_comm.isDone()) {
         completed.insert(completed.end(), &incoming_comm);
      }
   }

   for (owneri = incoming_ranks.begin(); owneri != incoming_ranks.end(); ++owneri,
        ++peer_idx) {
      const int peer_rank = *owneri;
      tbox::AsyncCommPeer<int>& outgoing_comm = all_comms[peer_idx];
      outgoing_comm.initialize(&comm_stage);
      outgoing_comm.setPeerRank(peer_rank);
      outgoing_comm.setMPI(mpi);
      outgoing_comm.setMPITag(tag0, tag1);
      outgoing_comm.limitFirstDataLength(MAPPING_CONNECTOR_ALGORITHM_FIRST_DATA_LENGTH);
   }

   t_modify_setup_comm->stop();
}

/*
 ***********************************************************************
 * Receive messages and unpack info sent from other processes.
 ***********************************************************************
 */
void MappingConnectorAlgorithm::privateModify_receiveAndUnpack(
   NeighborhoodSet& anchor_eto_new,
   NeighborhoodSet& new_eto_anchor,
   std::set<int>& outgoing_ranks,
   tbox::AsyncCommPeer<int> all_comms[],
   tbox::AsyncCommStage& comm_stage,
   tbox::AsyncCommStage::MemberVec& completed,
   const tbox::Dimension& dim,
   const tbox::SAMRAI_MPI& mpi) const
{
   t_modify_receive_and_unpack->start();

   const int rank(mpi.getRank());

   do {
      for (unsigned int i = 0; i < completed.size(); ++i) {

         tbox::AsyncCommPeer<int>* peer =
            dynamic_cast<tbox::AsyncCommPeer<int> *>(completed[i]);
         TBOX_ASSERT(completed[i] != NULL);
         TBOX_ASSERT(peer != NULL);

         if ((unsigned int)(peer - all_comms) < outgoing_ranks.size()) {
            // Received from this peer.
            const int sender = peer->getPeerRank();
            const int* ptr = peer->getRecvData();

            const int mapped_box_com_buffer_size = Box::commBufferSize(
                  dim);
            /*
             * Content of send_mesg:
             * - neighbor-removal section cached in neighbor_removal_mesg.
             * - offset to the reference section (see below)
             * - number of anchor mapped_boxes for which neighbors are found
             * - number of new mapped_boxes for which neighbors are found
             * - index of base/head mapped_box
             * - number of neighbors found for base/head mapped_box.
             *   - owner and local indices of neighbors found.
             *     Boxes of found neighbors are given in the
             *     reference section of the message.
             * - reference section: all the mapped_boxes referenced as
             *   neighbors (accumulated in referenced_anchor_nabrs
             *   and referenced_new_nabrs).
             *   - number of referenced base neighbors
             *   - number of referenced head neighbors
             *   - referenced base neighbors
             *   - referenced head neighbors
             */

            // Unpack neighbor-removal section.
            const int num_removed_mapped_boxes = *(ptr++);
            for (int ii = 0; ii < num_removed_mapped_boxes; ++ii) {
               const LocalId id_gone(*(ptr++));
               const BlockId block_id_gone(*(ptr++));
               const int number_affected = *(ptr++);
               const Box mapped_box_gone(dim, id_gone, sender, block_id_gone);
               if (s_print_modify_steps == 'y') tbox::plog << "Box "
                                                           << mapped_box_gone
                                                           << " removed, affecting "
                                                           << number_affected
                                                           << " mapped_boxes."
                                                           << std::endl;
//TODO: Is BoxId usage in this method correct regarding block id?
               for (int iii = 0; iii < number_affected; ++iii) {
                  const LocalId id_affected(*(ptr++));
                  const BlockId block_id_affected(*(ptr++));
                  NeighborhoodSet::iterator ci =
                     anchor_eto_new.find(BoxId(id_affected, rank, block_id_affected));
                  NeighborSet& nabrs = (*ci).second;
                  if (s_print_modify_steps == 'y') tbox::plog
                     << " Removing " << mapped_box_gone
                     << " from nabr list for " << id_affected
                     << std::endl;
                  NeighborSet::iterator ni = nabrs.find(mapped_box_gone);
                  TBOX_ASSERT(ni != nabrs.end());
                  nabrs.erase(ni);
               }
               TBOX_ASSERT(ptr != peer->getRecvData() + peer->getRecvSize());
            }

            // Get the referenced neighbor Boxes.
            NeighborSet referenced_base_nabrs;
            NeighborSet referenced_head_nabrs;
            const int offset = *(ptr++);
            const int* ref_mapped_box_ptr = peer->getRecvData() + offset;
            const int n_reference_base_mapped_boxes = *(ref_mapped_box_ptr++);
            const int n_reference_head_mapped_boxes = *(ref_mapped_box_ptr++);
            for (int ii = 0; ii < n_reference_base_mapped_boxes; ++ii) {
               Box mapped_box(dim);
               mapped_box.getFromIntBuffer(ref_mapped_box_ptr);
               referenced_base_nabrs.insert(
                  referenced_base_nabrs.end(), mapped_box);
               ref_mapped_box_ptr += mapped_box_com_buffer_size;
            }
            for (int ii = 0; ii < n_reference_head_mapped_boxes; ++ii) {
               Box mapped_box(dim);
               mapped_box.getFromIntBuffer(ref_mapped_box_ptr);
               referenced_head_nabrs.insert(
                  referenced_head_nabrs.end(), mapped_box);
               ref_mapped_box_ptr += mapped_box_com_buffer_size;
            }
            TBOX_ASSERT(
               ref_mapped_box_ptr == peer->getRecvData() + peer->getRecvSize());

            /*
             * Unpack neighbor data for base mapped_boxes and head
             * mapped_box_levels.  The neighbor info given includes
             * only owner and local index.  Refer to reference data to
             * get the box info.
             */
            const int n_base_mapped_boxes = *(ptr++);
            const int n_head_mapped_boxes = *(ptr++);
            for (int ii = 0; ii < n_base_mapped_boxes; ++ii) {
               LocalId local_id(*(ptr++));
               BlockId block_id(*(ptr++));
               const BoxId gid(local_id, rank, block_id);
               NeighborSet& nabrs = anchor_eto_new[gid];
               const int n_nabrs_found = *(ptr++);
               for (int j = 0; j < n_nabrs_found; ++j) {
                  Box tmp_nabr(dim, LocalId(ptr[1]), ptr[0], BlockId(ptr[2]));
                  ptr += 3;
                  NeighborSet::const_iterator na =
                     referenced_head_nabrs.find(tmp_nabr);
                  TBOX_ASSERT(na != referenced_head_nabrs.end());
                  const Box& nabr = *na;
                  nabrs.insert(nabr);
               }
            }
            for (int ii = 0; ii < n_head_mapped_boxes; ++ii) {
               LocalId local_id(*(ptr++));
               BlockId block_id(*(ptr++));
               const BoxId gid(local_id, rank, block_id);
               NeighborSet& nabrs = new_eto_anchor[gid];
               const int n_nabrs_found = *(ptr++);
               for (int j = 0; j < n_nabrs_found; ++j) {
                  Box tmp_nabr(dim, LocalId(ptr[1]), ptr[0], BlockId(ptr[2]));
                  ptr += 3;
                  NeighborSet::const_iterator na =
                     referenced_base_nabrs.find(tmp_nabr);
                  TBOX_ASSERT(na != referenced_base_nabrs.end());
                  const Box& nabr = *na;
                  nabrs.insert(nabr);
               }
            }
         } else {
            // Sent to this peer.  No need to do anything.
         }
      }

      completed.clear();
      comm_stage.advanceSome(completed);

   } while (completed.size() > 0);

   t_modify_receive_and_unpack->stop();
}

/*
 ***********************************************************************
 * Discover overlaps with new Boxes and send outgoing messages.
 * We interlace discovery with sends rather than wait until all
 * discovery is completed before sending.  This makes sure sends are
 * started asap to maximize communication and computation.
 ***********************************************************************
 */
void MappingConnectorAlgorithm::privateModify_discoverAndSend(
   std::map<int, std::vector<int> >& neighbor_removal_mesg,
   Connector& anchor_to_new,
   Connector& new_to_anchor,
   NeighborhoodSet& anchor_eto_new,
   NeighborhoodSet& new_eto_anchor,
   std::set<int>& incoming_ranks,
   std::set<int>& outgoing_ranks,
   tbox::AsyncCommPeer<int> all_comms[],
   tbox::AsyncCommStage::MemberVec& completed,
   const Connector& old_to_anchor,
   const Connector& anchor_to_old,
   const Connector& old_to_new) const
{
   t_modify_discover_and_send->start();

   const BoxLevel& old(old_to_new.getBase());
   const BoxLevel& new_mapped_box_level(old_to_new.getHead());
   const IntVector& anchor_to_new_width(anchor_to_new.getConnectorWidth());
   const IntVector& new_to_anchor_width(new_to_anchor.getConnectorWidth());
   const NeighborhoodSet& old_eto_anchor(old_to_anchor.getNeighborhoodSets());
   const NeighborhoodSet& old_eto_new(old_to_new.getNeighborhoodSets());
   const tbox::ConstPointer<GridGeometry>& grid_geometry(old.getGridGeometry());

   const tbox::Dimension& dim(old.getDim());
   const int rank = old.getMPI().getRank();
   const int nproc = old.getMPI().getSize();

   /*
    * visible_anchor_nabrs, visible_new_nabrs are the neighbors that
    * are seen by the local process.  For communication efficiency, we
    * use a looping construct that assumes that they are sorted by
    * owners first.  Note the comparator BoxOwnerFirst used to
    * achieve this ordering.
    */
   BoxSet visible_anchor_nabrs, visible_new_nabrs;
   InvertedNeighborhoodSet anchor_eto_old, new_eto_old;
   for (NeighborhoodSet::const_iterator ei = old_eto_anchor.begin();
        ei != old_eto_anchor.end(); ++ei) {
      const BoxId& old_gid = (*ei).first;
      const NeighborSet& anchor_nabrs = (*ei).second;
      for (NeighborSet::const_iterator na = anchor_nabrs.begin();
           na != anchor_nabrs.end(); ++na) {
         visible_anchor_nabrs.insert(*na);
         if (old_eto_new.find(old_gid) != old_eto_new.end()) {
            /*
             * anchor_eto_old is an InvertedNeighborhoodSet mapping visible anchor
             * Boxes to local old Boxes that are changing (excludes
             * old Boxes that do not change).
             */
            anchor_eto_old[*na].insert(old_gid);
         }
      }
   }
   for (NeighborhoodSet::const_iterator ei = old_eto_new.begin();
        ei != old_eto_new.end(); ++ei) {
      const BoxId& old_gid = (*ei).first;
      const NeighborSet& new_nabrs = (*ei).second;
      for (NeighborSet::const_iterator na = new_nabrs.begin();
           na != new_nabrs.end(); ++na) {
         visible_new_nabrs.insert(visible_new_nabrs.end(), *na);
         new_eto_old[*na].insert(old_gid);
      }
   }

   /*
    * Compute relationships.  Relationships are either locally
    * stored or packed into a message for sending.
    */

   /*
    * anchor_ni and new_ni point to the base/head mapped_box whose
    * neighbors are being sought.
    *
    * new_ni is used only if we are trying to compute new_to_anchor.
    * If we are not interested in that connector, then new_ni plays no
    * role.
    */
   NeighborSet::iterator anchor_ni;
   NeighborSet::iterator new_ni;
   /*
    * Local process can find neighbors for the owners of mapped_boxes
    * in visible_anchor_nabrs and visible_new_nabrs.  As an
    * optimization measure, start loop on the first owner with higher
    * rank than the local rank.  This avoid the higher-end ranks from
    * having to wait at the beginning and the lower-end ranks from
    * having to wait at the end.  After the highest rank owner has
    * been handled, start at the beginning and do the remaining ranks.
    * If local rank is highest of all owners of the visible
    * mapped_boxes, start at the begining.
    */

   const Box start_loop_here(dim, LocalId::getZero(), rank + 1);
   anchor_ni = visible_anchor_nabrs.lower_bound(start_loop_here);
   new_ni = visible_new_nabrs.lower_bound(start_loop_here);

   if (anchor_ni == visible_anchor_nabrs.end() &&
       new_ni == visible_new_nabrs.end()) {
      /*
       * There are no visible mapped_boxes owned by rank higher than
       * local process.  So loop from the beginning.
       */
      anchor_ni = visible_anchor_nabrs.begin();
      new_ni = visible_new_nabrs.begin();
   }

   /*
    * Set peer_idx to index the first incoming_ranks peer in all_comms.  It
    * will be incremented to correpond to the peer whose neighbors are
    * being searched for.
    */
   int peer_idx = static_cast<int>(outgoing_ranks.size());

   if (s_print_modify_steps == 'y') {
      tbox::plog << "visible_anchor_nabrs:" << std::endl;
      for (NeighborSet::const_iterator na = visible_anchor_nabrs.begin();
           na != visible_anchor_nabrs.end(); ++na) {
         tbox::plog << "  " << *na << std::endl;
      }
      tbox::plog << "visible_new_nabrs:" << std::endl;
      for (NeighborSet::const_iterator na = visible_new_nabrs.begin();
           na != visible_new_nabrs.end(); ++na) {
         tbox::plog << "  " << *na << std::endl;
      }
   }

   // Owners that we have sent messages to.
   std::set<int> owners_sent_to;

   /*
    * Loop until all visible anchor/new neighbors have their
    * new/anchor neighbors searched for.
    */
   while ((anchor_ni != visible_anchor_nabrs.end()) ||
          (new_ni != visible_new_nabrs.end())) {

      /*
       * curr_owner is the owner whose neighbors is currently being
       * searched for.  It should be the owner of the next anchor or
       * new Box in our cyclic-type looping.
       */
      int curr_owner = nproc; // Start with invalid value.
      if (anchor_ni != visible_anchor_nabrs.end() &&
          curr_owner > (*anchor_ni).getOwnerRank()) {
         curr_owner = (*anchor_ni).getOwnerRank();
      }
      if (new_ni != visible_new_nabrs.end() &&
          curr_owner > (*new_ni).getOwnerRank()) {
         curr_owner = (*new_ni).getOwnerRank();
      }

      std::vector<int> send_mesg;
      /*
       * Set up send_message to contain info discovered
       * locally and should be sent to curr_owner.
       *
       * Content of send_mesg:
       * - neighbor-removal section cached in neighbor_removal_mesg.
       * - offset to the reference section (see below)
       * - number of anchor mapped_boxes for which neighbors are found
       * - number of new mapped_boxes for which neighbors are found
       *   - index of base/head mapped_box
       *   - number of neighbors found for base/head mapped_box.
       *     - owner and local indices of neighbors found.
       *       Boxes of found neighbors are given in the
       *       reference section of the message.
       * - reference section: all the mapped_boxes referenced as
       *   neighbors (accumulated in referenced_anchor_nabrs
       *   and referenced_new_nabrs).
       *   - number of referenced base neighbors
       *   - number of referenced head neighbors
       *   - referenced base neighbors
       *   - referenced head neighbors
       *
       * The purpose of factoring out info on the neighbors referenced
       * is to reduce redundant data that can eat up lots of memory
       * and message passing bandwidth when there are lots of mapped_boxes
       * with the same neighbors.
       */

      /*
       * The first section of the send_mesg is the remote neighbor-removal
       * section (computed above).
       */
      if (curr_owner != rank) {
         // swap( send_mesg, neighbor_removal_mesg[curr_owner] );
         // We could use swap, but assign instead to leave data for debugging.
         send_mesg = neighbor_removal_mesg[curr_owner];
         if (send_mesg.empty()) {
            // No neighbor-removal data found for curr_owner.
            send_mesg.insert(send_mesg.end(), 0);
         }
      }

      // Indices of certain positions in send_mesg.
      const size_t idx_offset_to_ref = send_mesg.size();
      const size_t idx_num_anchor_mapped_boxes = idx_offset_to_ref + 1;
      const size_t idx_num_new_mapped_boxes = idx_offset_to_ref + 2;
      send_mesg.insert(send_mesg.end(), 3, 0);

      // Mapped_boxes referenced in the message, used when adding ref section.
      BoxSet referenced_anchor_nabrs; // Referenced neighbors in anchor.
      BoxSet referenced_new_nabrs; // Referenced neighbors in new.

      /*
       * Find locally visible new neighbors for all anchor
       * Boxes owned by curr_owner.
       */
      while (anchor_ni != visible_anchor_nabrs.end() &&
             (*anchor_ni).getOwnerRank() == curr_owner) {
         const Box& anchor_mapped_box = *anchor_ni;
         if (s_print_modify_steps == 'y')
            tbox::plog << "Finding neighbors for anchor_mapped_box "
                       << anchor_mapped_box << std::endl;
         Box compare_box = anchor_mapped_box;
         compare_box.grow(anchor_to_new_width);
         if (anchor_to_old.getHeadCoarserFlag()) compare_box.coarsen(
               anchor_to_old.getRatio());
         else if (old_to_anchor.getHeadCoarserFlag()) compare_box.refine(
               old_to_anchor.getRatio());

         BlockId compare_box_block_id(anchor_mapped_box.getBlockId());
         Box transformed_compare_box(compare_box);

         std::vector<Box> found_nabrs;
         const BoxIdSet& old_indices = anchor_eto_old[anchor_mapped_box];

         for (BoxIdSet::const_iterator na = old_indices.begin();
              na != old_indices.end(); ++na) {
            const NeighborSet& new_nabrs = old_eto_new.find(*na)->second;
            for (NeighborSet::const_iterator naa = new_nabrs.begin();
                 naa != new_nabrs.end(); ++naa) {
               const Box& new_nabr(*naa);
               if (compare_box_block_id != new_nabr.getBlockId()) {
                  // Re-transform compare_box and note its new BlockId.
                  transformed_compare_box = compare_box;
                  compare_box_block_id = new_nabr.getBlockId();
                  if (compare_box_block_id != anchor_mapped_box.getBlockId()) {
                     grid_geometry->transformBox(transformed_compare_box,
                        old.getRefinementRatio(),
                        new_nabr.getBlockId(),
                        anchor_mapped_box.getBlockId());
                  }
               }
               if (transformed_compare_box.intersects(new_nabr)) {
                  found_nabrs.insert(found_nabrs.end(), *naa);
               }
            }
         }
         if (s_print_modify_steps == 'y') {
            tbox::plog << "Found " << found_nabrs.size() << " neighbors :"
                       << std::endl;
            for (std::vector<Box>::const_iterator naa = found_nabrs.begin();
                 naa != found_nabrs.end(); ++naa) {
               tbox::plog << "  " << *naa << std::endl;
            }
         }
         if (!found_nabrs.empty()) {
            if (anchor_mapped_box.getOwnerRank() != rank) {
               // Pack up info for sending.
               ++send_mesg[idx_num_anchor_mapped_boxes];
               int subsize = 3 + 3 * static_cast<int>(found_nabrs.size());
               send_mesg.insert(send_mesg.end(), subsize, -1);
               int* submesg = &send_mesg[send_mesg.size() - subsize];
               *(submesg++) = anchor_mapped_box.getLocalId().getValue();
               *(submesg++) = anchor_mapped_box.getBlockId().getBlockValue();
               *(submesg++) = static_cast<int>(found_nabrs.size());
               for (std::vector<Box>::const_iterator na =
                       found_nabrs.begin();
                    na != found_nabrs.end(); ++na) {
                  const Box& nabr = *na;
                  referenced_new_nabrs.insert(nabr);
                  submesg[0] = nabr.getOwnerRank();
                  submesg[1] = nabr.getLocalId().getValue();
                  submesg[2] = nabr.getBlockId().getBlockValue();
                  submesg += 3;
               }
            } else {
               /*
                * Save neighbor info locally.
                *
                * To improve communication time, we should really send
                * the head neighbors before doing anything locally.
                */
               NeighborSet& local_nabrs =
                  anchor_eto_new[anchor_mapped_box.getId()];
               for (std::vector<Box>::const_iterator na =
                       found_nabrs.begin();
                    na != found_nabrs.end(); ++na) {
                  const Box& nabr = *na;
                  local_nabrs.insert(nabr);
               }
            }
         }
         if (s_print_modify_steps == 'y')
            tbox::plog << "Erasing visible base nabr " << (*anchor_ni)
                       << std::endl;
         visible_anchor_nabrs.erase(anchor_ni++);
         if (s_print_modify_steps == 'y') {
            if (anchor_ni == visible_anchor_nabrs.end())
               tbox::plog << "Next base nabr: end" << std::endl;
            else
               tbox::plog << "Next base nabr: " << *anchor_ni << std::endl;
         }

      }

      /*
       * Find locally visible anchor neighbors for all new
       * Boxes owned by curr_owner.
       */
      while (new_ni != visible_new_nabrs.end() &&
             (*new_ni).getOwnerRank() == curr_owner) {
         const Box& new_mapped_box = *new_ni;
         if (s_print_modify_steps == 'y')
            tbox::plog << "Finding neighbors for new_mapped_box "
                       << new_mapped_box << std::endl;
         Box compare_box = new_mapped_box;
         compare_box.grow(new_to_anchor_width);
         if (anchor_to_new.getHeadCoarserFlag()) compare_box.refine(
               anchor_to_new.getRatio());
         else if (new_to_anchor.getHeadCoarserFlag()) compare_box.coarsen(
               new_to_anchor.getRatio());

         BlockId compare_box_block_id(new_mapped_box.getBlockId());
         Box transformed_compare_box(compare_box);

         std::vector<Box> found_nabrs;
         const BoxIdSet& old_indices = new_eto_old[new_mapped_box];

         for (BoxIdSet::const_iterator na = old_indices.begin();
              na != old_indices.end(); ++na) {
            NeighborhoodSet::const_iterator anchor_nabrs_i = old_eto_anchor.find(*na);
            if (anchor_nabrs_i != old_eto_anchor.end()) {
               /*
                * There are anchor Boxes with relationships to
                * the old Box identified by *na.
                */
               const NeighborSet& anchor_nabrs = anchor_nabrs_i->second;
               for (NeighborSet::const_iterator naa = anchor_nabrs.begin();
                    naa != anchor_nabrs.end(); ++naa) {
                  const Box& anchor_nabr(*naa);
                  if (compare_box_block_id != anchor_nabr.getBlockId()) {
                     // Re-transform compare_box and note its new BlockId.
                     transformed_compare_box = compare_box;
                     compare_box_block_id = anchor_nabr.getBlockId();
                     if (compare_box_block_id != new_mapped_box.getBlockId()) {
                        grid_geometry->transformBox(transformed_compare_box,
                           new_mapped_box_level.getRefinementRatio(),
                           anchor_nabr.getBlockId(),
                           new_mapped_box.getBlockId());
                     }
                  }
                  if (transformed_compare_box.intersects(anchor_nabr)) {
                     found_nabrs.insert(found_nabrs.end(), *naa);
                  }
               }
            }
         }

         if (s_print_modify_steps == 'y') {
            tbox::plog << "Found " << found_nabrs.size() << " neighbors:"
                       << std::endl;
            for (std::vector<Box>::const_iterator naa = found_nabrs.begin();
                 naa != found_nabrs.end(); ++naa) {
               tbox::plog << "  " << *naa << std::endl;
            }
         }

         if (!found_nabrs.empty()) {
            if (new_mapped_box.getOwnerRank() != rank) {
               // Pack up info for sending.
               ++send_mesg[idx_num_new_mapped_boxes];
               int subsize = 3 + 3 * static_cast<int>(found_nabrs.size());
               send_mesg.insert(send_mesg.end(), subsize, -1);
               int* submesg = &send_mesg[send_mesg.size() - subsize];
               *(submesg++) = new_mapped_box.getLocalId().getValue();
               *(submesg++) = new_mapped_box.getBlockId().getBlockValue();
               *(submesg++) = static_cast<int>(found_nabrs.size());
               for (std::vector<Box>::const_iterator na =
                       found_nabrs.begin();
                    na != found_nabrs.end(); ++na) {
                  const Box& nabr(*na);
                  referenced_anchor_nabrs.insert(nabr);
                  submesg[0] = nabr.getOwnerRank();
                  submesg[1] = nabr.getLocalId().getValue();
                  submesg[2] = nabr.getBlockId().getBlockValue();
                  submesg += 3;
               }
            } else {
               /*
                * Save neighbor info locally.
                */
               NeighborSet& local_nabrs =
                  new_eto_anchor[new_mapped_box.getId()];
               for (std::vector<Box>::const_iterator na =
                       found_nabrs.begin();
                    na != found_nabrs.end(); ++na) {
                  const Box& nabr = *na;
                  local_nabrs.insert(nabr);
               }
            }
         }
         if (s_print_modify_steps == 'y')
            tbox::plog << "Erasing visible head nabr " << (*new_ni)
                       << std::endl;
         visible_new_nabrs.erase(new_ni++);
         if (s_print_modify_steps == 'y') {
            if (new_ni == visible_new_nabrs.end())
               tbox::plog << "Next head nabr: end" << std::endl;
            else
               tbox::plog << "Next head nabr: " << *new_ni << std::endl;
         }
      }

      if (curr_owner != rank) {
         /*
          * Send discoveries to the curr_owner.
          */

         /*
          * Fill the messages's reference section with neighbors
          * that have been referenced.
          */
         const int offset = send_mesg[idx_offset_to_ref] =
               static_cast<int>(send_mesg.size());
         const int n_referenced_nabrs =
            static_cast<int>(
               referenced_new_nabrs.size() + referenced_anchor_nabrs.size());
         const int reference_section_size =
            2 + n_referenced_nabrs * Box::commBufferSize(dim);
         send_mesg.insert(send_mesg.end(),
            reference_section_size,
            -1);
         int* ptr = &send_mesg[offset];
         *(ptr++) = static_cast<int>(referenced_anchor_nabrs.size());
         *(ptr++) = static_cast<int>(referenced_new_nabrs.size());
         for (BoxSet::const_iterator ni = referenced_anchor_nabrs.begin();
              ni != referenced_anchor_nabrs.end(); ++ni) {
            const Box& mapped_box = *ni;
            mapped_box.putToIntBuffer(ptr);
            ptr += Box::commBufferSize(dim);
         }
         for (BoxSet::const_iterator ni = referenced_new_nabrs.begin();
              ni != referenced_new_nabrs.end(); ++ni) {
            const Box& mapped_box = *ni;
            mapped_box.putToIntBuffer(ptr);
            ptr += Box::commBufferSize(dim);
         }

         TBOX_ASSERT(ptr == &send_mesg[send_mesg.size() - 1] + 1);

         /*
          * Find the communication object by increasing peer_idx
          * (cyclically) until it corresponds to curr_owner.
          */
         while (all_comms[peer_idx].getPeerRank() != curr_owner) {
            ++peer_idx;
            if (peer_idx == static_cast<int>(outgoing_ranks.size()
                                             + incoming_ranks.size()))
               peer_idx -= static_cast<int>(incoming_ranks.size());
         }
         /*
          * Send message.
          */
         tbox::AsyncCommPeer<int>& outgoing_comm = all_comms[peer_idx];
         // Make sure we picked the correct peer.
         TBOX_ASSERT(outgoing_comm.getPeerRank() == curr_owner);

         outgoing_comm.beginSend(&send_mesg[0],
            static_cast<int>(send_mesg.size()));
         if (outgoing_comm.isDone()) {
            completed.insert(completed.end(), &outgoing_comm);
         }

         TBOX_ASSERT(owners_sent_to.find(curr_owner) == owners_sent_to.end());

         owners_sent_to.insert(curr_owner);

      }

      /*
       * If we come to the end of visible mapped_boxes, go back and
       * work on the mapped_boxes owned by processors with lower rank
       * than the local rank.  (This is part of the optimization to
       * reduce communication time.)
       */
      if (anchor_ni == visible_anchor_nabrs.end() &&
          new_ni == visible_new_nabrs.end()) {
         /*
          * There are no mapped_boxes that are owned by rank higher
          * than local process and that we want to find neighbors for.
          * So loop from the beginning.
          */
         anchor_ni = visible_anchor_nabrs.begin();
         new_ni = visible_new_nabrs.begin();
      }

   }

   t_modify_discover_and_send->stop();
}

/*
 ***********************************************************************
 * Remove relationships made obsolete by mapping.  Cache outgoing
 * information in message buffers.
 ***********************************************************************
 */
void MappingConnectorAlgorithm::privateModify_removeAndCache(
   std::map<int, std::vector<int> >& neighbor_removal_mesg,
   NeighborhoodSet& anchor_eto_new,
   NeighborhoodSet& new_eto_anchor,
   const Connector& old_to_new) const
{
   t_modify_remove_and_cache->start();

   const tbox::Dimension& dim(old_to_new.getBase().getDim());
   const tbox::SAMRAI_MPI& mpi(old_to_new.getBase().getMPI());
   const int rank(mpi.getRank());

   const NeighborhoodSet& old_eto_new = old_to_new.getNeighborhoodSets();

   /*
    * Remove relationships with old mapped_boxes (because
    * they are going away). These are Boxes mapped by
    * old_to_new.
    *
    * Erase local old Boxes from new_eto_anchor.
    *
    * If the old mapped_boxes have neighbors in the anchor
    * BoxLevel, some relationships from anchor_eto_old should be
    * erased also.  For each neighbor from a remote anchor mapped_box to a
    * local old mapped_box, add data to mesg_to_owners saying what
    * Box is disappearing and what anchor Box should
    * no longer reference it.
    */
   for (NeighborhoodSet::const_iterator iold = old_eto_new.begin();
        iold != old_eto_new.end(); ++iold) {

      const BoxId& old_gid_gone = iold->first;
      const Box old_mapped_box_gone(dim, old_gid_gone);

      NeighborhoodSet::iterator inew =
         new_eto_anchor.find(old_gid_gone);

      if (inew != new_eto_anchor.end()) {
         // old_gid_gone exists in new_to_anchor.  Remove it.

         const NeighborSet& affected_anchor_nabrs = inew->second;

         if (s_print_modify_steps == 'y') tbox::plog << "Box "
                                                     << old_mapped_box_gone
                                                     << " is gone."
                                                     << std::endl;

         for (NeighborSet::const_iterator ianchor = affected_anchor_nabrs.begin();
              ianchor != affected_anchor_nabrs.end(); /* incremented in loop */) {

            if (s_print_modify_steps == 'y') tbox::plog << "  Box "
                                                        << *ianchor
                                                        << " is affected."
                                                        << std::endl;

            const int anchor_nabr_owner = ianchor->getOwnerRank();
            if (anchor_nabr_owner == rank) {
               // Erase local relationship from anchor to old_gid_gone.
               do {

                  if (s_print_modify_steps == 'y') tbox::plog
                     << "  Fixing affected mapped_box " << *ianchor
                     << std::endl;

                  NeighborhoodSet::iterator ck = anchor_eto_new.find(ianchor->getId());
                  TBOX_ASSERT(ck != anchor_eto_new.end());

                  NeighborSet& nabrs_nabrs = ck->second;

                  if (s_print_modify_steps == 'y') {
                     tbox::plog << nabrs_nabrs.format("XX-> ") << std::endl;
                  }

                  NeighborSet::iterator nb = nabrs_nabrs.find(old_mapped_box_gone);
                  if (nb != nabrs_nabrs.end()) {
                     if (s_print_modify_steps == 'y') tbox::plog
                        << "    Removing neighbor " << *nb
                        << " from list for " << *ianchor << std::endl;
                     nabrs_nabrs.erase(nb);
                  }

                  ++ianchor;

                  // Skip past periodic image Boxes.
                  while (ianchor != affected_anchor_nabrs.end() &&
                         ianchor->isPeriodicImage()) {
                     ++ianchor;
                  }

               } while (ianchor != affected_anchor_nabrs.end() &&
                        ianchor->getOwnerRank() == rank);
            } else {
               // Tell owner of nabr to erase references to old_gid_gone.
               std::vector<int>& mesg = neighbor_removal_mesg[anchor_nabr_owner];
               // mesg[0] is the counter for how many mapped_boxes are removed.
               if (mesg.empty()) {
                  mesg.insert(mesg.end(), 1);
               } else { ++mesg[0];
               }
               mesg.insert(mesg.end(), old_gid_gone.getLocalId().getValue());
               mesg.insert(mesg.end(), old_gid_gone.getBlockId().getBlockValue());
               int i_count = static_cast<int>(mesg.size());
               mesg.insert(mesg.end(), 0);
               do {
                  mesg.insert(mesg.end(), ianchor->getLocalId().getValue());
                  mesg.insert(mesg.end(), ianchor->getBlockId().getBlockValue());
                  ++mesg[i_count];
                  if (s_print_modify_steps == 'y') tbox::plog
                     << "    Request change " << mesg[i_count]
                     << " to neighbors fo " << *ianchor << std::endl;
                  ++ianchor;
               } while (ianchor != affected_anchor_nabrs.end() &&
                        ianchor->getOwnerRank() == anchor_nabr_owner);
            }
         }

         /*
          * Erase relationships from old_mapped_box_gone to anchor
          * mapped_box_level.
          */
         new_eto_anchor.erase(inew);

      }

   }

   t_modify_remove_and_cache->stop();
}

/*
 ***********************************************************************
 * Do some standard (not expensive) checks on the arguments of modify.
 ***********************************************************************
 */

void MappingConnectorAlgorithm::checkModifyParameters(
   const Connector& anchor_to_mapped,
   const Connector& mapped_to_anchor,
   const Connector& old_to_new,
   const Connector& new_to_old) const
{
   const BoxLevel& old = mapped_to_anchor.getBase();

   /*
    * Ensure that Connectors incident to and from the old agree on
    * what the old is.
    */
   if (&old != &old_to_new.getBase() ||
       (new_to_old.isInitialized() && &old !=
        &new_to_old.getHead()) ||
       &old != &anchor_to_mapped.getHead()) {
      TBOX_ERROR("Bad input for MappingConnectorAlgorithm::modify:\n"
         << "Given Connectors to anchor and new of modify are not incident\n"
         << "from the same old in MappingConnectorAlgorithm::modify:\n"
         << "anchor_to_mapped is  TO  " << &anchor_to_mapped.getHead() << "\n"
         << "old_to_new is FROM " << &old_to_new.getBase()
         << "\n"
         << "new_to_old is  TO  " << &new_to_old.getHead()
         << "\n"
         << "mapped_to_anchor is FROM " << &mapped_to_anchor.getBase() << "\n"
         );
   }
   /*
    * Ensure that new and anchor mapped_box_levels in argument agree with
    * new and anchor in the object.
    */
   if (&mapped_to_anchor.getHead() != &anchor_to_mapped.getBase()) {
      TBOX_ERROR("Bad input for MappingConnectorAlgorithm::modify:\n"
         << "Given Connectors to and from anchor of modify do not refer\n"
         << "to the anchor of the modify in MappingConnectorAlgorithm::modify:\n"
         << "anchor_to_mapped is FROM " << &anchor_to_mapped.getBase() << "\n"
         << "mapped_to_anchor is  TO  " << &mapped_to_anchor.getHead() << "\n"
         << "anchor of modify is    " << &anchor_to_mapped.getBase() << "\n"
         );
   }
   if (new_to_old.isInitialized() && &old_to_new.getHead() !=
       &new_to_old.getBase()) {
      TBOX_ERROR("Bad input for MappingConnectorAlgorithm::modify:\n"
         << "Given Connectors to and from new of modify do not refer\n"
         << "to the new of the modify in MappingConnectorAlgorithm::modify:\n"
         << "new_to_old is FROM " << &new_to_old.getBase()
         << "\n"
         << "old_to_new is  TO  " << &old_to_new.getHead()
         << "\n"
         << "new of modify is    " << &anchor_to_mapped.getHead() << "\n"
         );
   }
   if (!anchor_to_mapped.isTransposeOf(mapped_to_anchor)) {
      TBOX_ERROR("Bad input for MappingConnectorAlgorithm::modify:\n"
         << "Given Connectors between anchor and mapped of modify\n"
         << "are not transposes of each other.\n"
         << "See Connector::isTransposeOf().\n"
         );
   }
   if (new_to_old.isInitialized() &&
       !new_to_old.isTransposeOf(old_to_new)) {
      TBOX_ERROR("Bad input for MappingConnectorAlgorithm::modify:\n"
         << "Given Connectors between new and old of modify\n"
         << "are not transposes of each other.\n"
         << "See Connector::isTransposeOf().\n"
         );
   }
   if (anchor_to_mapped.getParallelState() != BoxLevel::DISTRIBUTED) {
      TBOX_ERROR("Bad input for MappingConnectorAlgorithm::modify:\n"
         << "bridging is currently set up for DISTRIBUTED\n"
         << "mode only.\n");
   }
   if (mapped_to_anchor.getParallelState() != BoxLevel::DISTRIBUTED) {
      TBOX_ERROR("Bad input for MappingConnectorAlgorithm::modify:\n"
         << "bridging is currently set up for DISTRIBUTED\n"
         << "mode only.\n");
   }

   // Expensive sanity checks:
   if (d_sanity_check_inputs) {
      anchor_to_mapped.assertTransposeCorrectness(mapped_to_anchor);
      mapped_to_anchor.assertTransposeCorrectness(anchor_to_mapped);
      if (new_to_old.isInitialized()) {
         /*
          * Not sure if the following are valid checks for modify operation.
          * Modify *may* have different restrictions on the mapping Connector.
          */
         new_to_old.assertTransposeCorrectness(old_to_new);
         old_to_new.assertTransposeCorrectness(new_to_old);
      }
      size_t nerrs = findMappingErrors(old_to_new);
      if (nerrs != 0) {
         TBOX_ERROR("MappingConnectorUtil::privateModify: found errors in\n"
            << "mapping Connector.");
      }
   }
}

/*
 ***********************************************************************
 * Run findMappingErrors and assert that no errors are found.
 ***********************************************************************
 */

void MappingConnectorAlgorithm::assertMappingValidity(
   const Connector& connector,
   char is_local_map) const
{
   size_t nerr = findMappingErrors(connector, is_local_map);
   if (nerr != 0) {
      tbox::perr << "MappingConnectorAlgorithm::assertMappingValidity found\n"
                 << nerr << " errors.\n"
                 << "mapping connector:\n" << connector.format("MAP: ", 2)
                 << "pre-map:\n" << connector.getBase().format("PRE: ", 2)
                 << "post-map:\n" << connector.getHead().format("POST: ", 2)
                 << std::endl;
      TBOX_ERROR("MappingConnectorAlgorithm::assertMappingValidity exiting due\n"
         << "to above errors.");
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */

size_t MappingConnectorAlgorithm::findMappingErrors(
   const Connector& connector,
   char is_local_map) const
{
   const tbox::SAMRAI_MPI& mpi(connector.getMPI());
   const NeighborhoodSet& neighborhoods = connector.getNeighborhoodSets();

   // Need to know whether this is a local map.
   if (is_local_map == '\0') {
      if (mpi.getSize() > 1) {
         for (NeighborhoodSet::const_iterator ei = neighborhoods.begin();
              ei != neighborhoods.end(); ++ei) {
            const NeighborSet& nabrs = (*ei).second;
            for (NeighborSet::const_iterator ni = nabrs.begin();
                 ni != nabrs.end(); ++ni) {
               if ((*ni).getOwnerRank() != connector.getRank()) {
                  is_local_map = 'n';
                  break;
               }
            }
            if (is_local_map == 'n') {
               break;
            }
         }
         if (is_local_map != 'n') {
            is_local_map = 'y';
         }
         int tmpi = is_local_map == 'y' ? 0 : 1;
         int tmpj; // For some reason, MPI_IN_PLACE is undeclared!
         mpi.Allreduce(&tmpi, &tmpj, 1, MPI_INT, MPI_MAX);
         if (tmpj > 0) is_local_map = 'n';
      } else {
         is_local_map = 'y';
      }
   }

   /*
    * If not a local map, we need a globalized copy of the head.
    */
   const BoxLevel& new_mapped_box_level =
      is_local_map == 'y' ? connector.getHead() :
      connector.getHead().getGlobalizedVersion();

   int error_count = 0;

   /*
    * Find old Boxes that changed or disappeared on
    * the new BoxLevel.  There should be a mapping for each
    * Box that changed or disappeared.
    */
   const BoxSet& old_mapped_boxes = connector.getBase().getBoxes();
   for (RealBoxConstIterator ni(old_mapped_boxes); ni.isValid();
        ++ni) {
      const Box& old_mapped_box = *ni;
      if (!new_mapped_box_level.hasBox(old_mapped_box)) {
         // old_mapped_box disappeared.  Require a mapping for old_mapped_box.
         if (!connector.hasNeighborSet(old_mapped_box.getId())) {
            ++error_count;
            tbox::perr << "MappingConnectorAlgorithm::findMappingError ("
                       << error_count
                       << "): old mapped_box " << old_mapped_box
                       << " disappeared without being mapped." << std::endl;
         }
      } else {
         const Box& new_mapped_box =
            *(new_mapped_box_level.getBoxStrict(old_mapped_box));
         if (!new_mapped_box.isSpatiallyEqual(old_mapped_box)) {
            // old_mapped_box has changed its box.  A mapping must exist for it.
            if (!connector.hasNeighborSet(old_mapped_box.getId())) {
               ++error_count;
               tbox::perr << "MappingConnectorAlgorithm::findMappingError ("
                          << error_count
                          << "): old mapped_box " << old_mapped_box
                          << " changed to " << new_mapped_box
                          << " without being mapped." << std::endl;
            }
         }
      }
   }

   /*
    * All mappings should point from a old Box to a new
    * set of Boxes.
    */
   for (NeighborhoodSet::const_iterator ei = neighborhoods.begin();
        ei != neighborhoods.end(); ++ei) {

      const BoxId& gid = (*ei).first;
      const NeighborSet& nabrs = (*ei).second;

      if (!connector.getBase().hasBox(gid)) {
         // Mapping does not go from a old mapped_box.
         ++error_count;
         tbox::perr << "MappingConnectorAlgorithm::findMappingError ("
                    << error_count
                    << "): mapping given for nonexistent index " << gid
                    << std::endl;
      } else {
         const Box& old_mapped_box =
            *(connector.getBase().getBoxStrict(gid));

         Box grown_box(old_mapped_box);
         grown_box.grow(connector.getConnectorWidth());

         for (NeighborSet::const_iterator ni = nabrs.begin();
              ni != nabrs.end(); ++ni) {
            const Box& nabr = *ni;

            if (!grown_box.contains(nabr)) {
               ++error_count;
               tbox::perr << "MappingConnectorAlgorithm::findMappingError ("
                          << error_count
                          << "): old mapped_box " << old_mapped_box
                          << " grown by " << connector.getConnectorWidth()
                          << " to " << grown_box << " does not contain neighbor "
                          << nabr << std::endl;
            }

            if (!new_mapped_box_level.hasBox(nabr)) {
               ++error_count;
               tbox::perr << "MappingConnectorAlgorithm::findMappingError ("
                          << error_count
                          << "): old mapped_box " << old_mapped_box
                          << " mapped to nonexistent new mapped_box "
                          << nabr << std::endl;
            } else {
               const Box& head_mapped_box =
                  *(new_mapped_box_level.getBoxStrict(nabr.getId()));
               if (!nabr.isSpatiallyEqual(head_mapped_box)) {
                  ++error_count;
                  tbox::perr << "MappingConnectorAlgorithm::findMappingError ("
                             << error_count
                             << "): old mapped_box " << old_mapped_box
                             << " mapped to neighbor " << nabr
                             << " inconsistent new mapped_box "
                             << head_mapped_box << std::endl;
               }
            }

         }
      }
   }

   return error_count;
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void MappingConnectorAlgorithm::initializeCallback()
{
   /*
    * The first constructor:
    * - gets timers from the TimerManager.
    * - sets up their deallocation.
    * - sets up debugging flags.
    */

   if (s_print_modify_steps == 0) {
      if (tbox::InputManager::inputDatabaseExists()) {
         s_print_modify_steps = 'n';
         tbox::Pointer<tbox::Database> idb =
            tbox::InputManager::getInputDatabase();
         if (idb->isDatabase("MappingConnectorAlgorithm")) {
            tbox::Pointer<tbox::Database> ocu_db =
               idb->getDatabase("MappingConnectorAlgorithm");
            s_print_modify_steps =
               ocu_db->getCharWithDefault("print_modify_steps",
                  s_print_modify_steps);
         }
      }
   }

   t_modify = tbox::TimerManager::getManager()->
      getTimer("hier::MappingConnectorAlgorithm::privateModify()");
   t_modify_shortcut = tbox::TimerManager::getManager()->
      getTimer("hier::MappingConnectorAlgorithm::privateModify()_shortcut");
   t_modify_setup_comm = tbox::TimerManager::getManager()->
      getTimer("hier::MappingConnectorAlgorithm::privateModify_setupCommunication()");
   t_modify_remove_and_cache = tbox::TimerManager::getManager()->
      getTimer("hier::MappingConnectorAlgorithm::privateModify_removeAndCache()");
   t_modify_discover_and_send = tbox::TimerManager::getManager()->
      getTimer("hier::MappingConnectorAlgorithm::privateModify_discoverAndSend()");
   t_modify_receive_and_unpack = tbox::TimerManager::getManager()->
      getTimer("hier::MappingConnectorAlgorithm::privateModify_receiveAndUnpack()");
   t_modify_MPI_wait = tbox::TimerManager::getManager()->
      getTimer("hier::MappingConnectorAlgorithm::privateModify()_MPI_wait");
}

/*
 ***************************************************************************
 * Free statics.  To be called by shutdown registry to make sure
 * memory for statics do not leak.
 ***************************************************************************
 */

void MappingConnectorAlgorithm::finalizeCallback()
{
   t_modify.setNull();
   t_modify_setup_comm.setNull();
   t_modify_remove_and_cache.setNull();
   t_modify_discover_and_send.setNull();
   t_modify_receive_and_unpack.setNull();
   t_modify_MPI_wait.setNull();

   if (s_class_mpi.getCommunicator() != tbox::SAMRAI_MPI::commNull) {
      s_class_mpi.freeCommunicator();
   }
}

}
}
#endif
