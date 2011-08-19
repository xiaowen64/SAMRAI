/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Algorithms to work with maping Connectors. 
 *
 ************************************************************************/
#ifndef included_hier_MappingConnectorAlgorithm
#define included_hier_MappingConnectorAlgorithm

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Connector.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/MappedBoxLevel.h"
#include "SAMRAI/tbox/AsyncCommPeer.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"

#include <map>
#include <string>

namespace SAMRAI {
namespace hier {

/*!
 * @brief Algorithms for using Connectors representing changes
 * to a MappedBoxLevel.
 *
 * MappingConnectorAlgorithm objects check and apply mappings.
 */
class MappingConnectorAlgorithm
{

public:
   /*!
    * @brief Constructor.
    *
    * The default constructor creates an uninitialized object in
    * distributed state.
    *
    * @see initialize()
    * @see swapInitialize()
    */
   explicit MappingConnectorAlgorithm();

   /*!
    * @brief Deallocate internal data.
    */
   virtual ~MappingConnectorAlgorithm(
      void);

   /*!
    * @brief Set whether to run expensive sanity checks on input
    * parameters when at the beginning of certain methods.
    *
    * The checks are expensive and meant mainly for debugging.
    *
    * @param[in] do_check
    */
   void
   setSanityCheckMethodPreconditions(
      bool do_check);

   /*!
    * @brief Set whether to run expensive sanity checks on outputs
    * before returning from certain methods.
    *
    * The checks are expensive and meant mainly for debugging.
    *
    * @param[in] do_check
    */
   void
   setSanityCheckMethodPostconditions(
      bool do_check);

   /*!
    * @brief Most general version for modifying Connectors using
    * mapping Connectors.  Modification is the changing of existing
    * Connectors when mapped_boxes in a MappedBoxLevel changes according to specified
    * mapping connectors.
    *
    * The change is represented by a the mapper @c old_to_new
    * and its transpose @c new_to_old.  The Connectors to be
    * modified are @c anchor_to_mapped and @c mapped_to_anchor, which
    * on input, go between a anchor (not mapped) MappedBoxLevel and the old
    * MappedBoxLevel.  On output, these Connectors will go from the anchor mapped_box_level
    * to the new mapped_box_level.
    *
    * @code
    * Input:
    *
    *                                    (anchor)
    *                                   ^ /
    *                                  / /
    *              mapped_to_anchor-> / /
    *                                / / <--anchor_to_mapped
    *                               / v
    * mapped mapped_box_level:    (old) ---------> (new)
    *                                   <---------
    *
    *
    * Output:
    *
    *                                    (anchor)
    *                                          \ ^
    *                                           \ \
    *                                            \ \ <-mapped_to_anchor
    *                        anchor_to_mapped --> \ \
    *                                              v \
    *                                              (new)
    *
    * @endcode
    *
    * The MappedBoxLevels before and after the mapping are represented
    * by the base and head of the mapping.  No MappedBoxLevel is modified,
    * other than the "mutable" MappedBoxLevels in the argument.
    *
    * An important constraint in the old_to_new Connectors is
    * that this method cannot handle multiple maps at once.  For
    * example, it cannot map mapped_box J to mapped_box K and at the
    * same time map mapped_box I to mapped_box J.  Box J in the
    * old mapped_box_level and mapped_box J on the new
    * mapped_box_level are considered entirely different mapped_boxes.
    *
    * After modifying, the output Connectors that had referenced old
    * MappedBoxLevels will be reset to reference the new
    * MappedBoxedLevel.  This is the end of the modify operation.
    *
    * The following "in-place" modification is provided for users'
    * convenience.  Often times, the "new" level is a temporary object
    * and users often reset output equivalent to:
    *
    * @li old = new
    * @li reset Connectors using old in place of new.
    *
    * The modify methods support these optional steps as follows: If
    * mutable versions of some MappedBoxLevel are given, the output
    * Connectors can be reset to reference these versions instead.
    *
    * If mutable_new points to the old MappedBoxLevel and mutable_old
    * points to the new, then do an in-place switch as follows:
    * @li Swap the mutable_old and mutable_new MappedBoxLevels.
    * @li Use mutable_old (which actually the new MappedBoxLevel after
    *    the swap) as the mapped MappedBoxLevel in the output Connectors.
    * Otherwise:
    * @li If mutable_new is non-NULL, set it equal to new and use
    *    it as the mapped MappedBoxLevel in the output Connectors.
    * @li If mutable_old is non-NULL, set it equal to the old MappedBoxLevel.
    *
    * @param[in,out] anchor_to_mapped Connector to be modified.  On input, this
    *   points to the MappedBoxLevel being mapped.
    * @param[in,out] mapped_to_anchor Reverse (transpose) of anchor_to_mapped.
    *   points to the MappedBoxLevel being mapped.
    * @param[in] old_to_new Mapping from the old MappedBoxLevel to the
    *   new MappedBoxLevel.
    *   The width of old-->new should indicate the maximum amount of box
    *   growth caused by the change.  A value of zero means no growth.
    * @param[in] new_to_old Reverse (transpose) of old_to_new.
    * @param[in,out] mutable_new See comments.
    * @param[in,out] mutable_old See comments.
    */
   void
   modify(
      Connector& anchor_to_mapped,
      Connector& mapped_to_anchor,
      const Connector& old_to_new,
      const Connector& new_to_old,
      MappedBoxLevel* mutable_new = NULL,
      MappedBoxLevel* mutable_old = NULL) const;

   /*!
    * @brief Version of modify requiring only the forward map
    * and allows only local mappings.
    *
    * This version does not require the reverse mapping, but the
    * mapping must be totally local, e.g.  old_to_new must
    * contain no remote neighbor.
    *
    * @code
    * Input:
    *
    *                                    (anchor)
    *                                   ^ /
    *                                  / /
    *              mapped_to_anchor-> / /
    *                                / / <--anchor_to_mapped
    *                               / v
    * mapped mapped_box_level:    (old) ---------> (new)
    *
    *
    * Output:
    *
    *                                    (anchor)
    *                                          \ ^
    *                                           \ \
    *                                            \ \ <-mapped_to_anchor
    *                        anchor_to_mapped --> \ \
    *                                              v \
    *                                              (new)
    *
    * @endcode
    *
    * The modify methods support these optional steps as follows: If
    * mutable versions of some MappedBoxLevel are given, the output
    * Connectors can be reset to reference these versions instead.
    *
    * If mutable_new points to the old MappedBoxLevel and mutable_old
    * points to the new, then do an in-place switch as follows:
    * @li Swap the mutable_old and mutable_new MappedBoxLevels.
    * @li Use mutable_old (which actually the new MappedBoxLevel after
    *    the swap) as the mapped MappedBoxLevel in the output Connectors.
    * Otherwise:
    * @li If mutable_new is non-NULL, set it equal to new and use
    *    it as the mapped MappedBoxLevel in the output Connectors.
    * @li If mutable_old is non-NULL, set it equal to the old MappedBoxLevel.
    *
    * @param[in,out] anchor_to_mapped Connector to be modified.  On input, this
    *   points to the MappedBoxLevel being mapped.
    * @param[in,out] mapped_to_anchor Reverse (transpose) of anchor_to_mapped.
    *   points to the MappedBoxLevel being mapped.
    * @param[in] old_to_new Mapping from the old MappedBoxLevel to the
    *   new MappedBoxLevel.
    *   The width of old-->new should indicate the maximum amount of box
    *   growth caused by the change.  A value of zero means no growth.
    * @param[in,out] mutable_new See comments.
    * @param[in,out] mutable_old See comments.
    */
   void
   modify(
      Connector& anchor_to_mapped,
      Connector& mapped_to_anchor,
      const Connector& old_to_new,
      MappedBoxLevel* mutable_new = NULL,
      MappedBoxLevel* mutable_old = NULL) const;

   /*!
    * @brief Version of modify requiring only the forward map
    * and does not set the Connector from the mapped MappedBoxLevel back
    * to the anchor MappedBoxLevel.
    *
    * This version does not require the reverse mapping,
    * but the forward mapping must be totally local, e.g.
    * old_to_new must contain no remote neighbor.
    *
    * This version does not update the Connector from
    * the mapped mapped_box_level back to the anchor mapped_box_level.
    * (This is implicit in the fact that mapped_to_anchor
    * is absent from the interface.)
    *
    * @code
    * Input:
    *
    *                                    (anchor)
    *                                     /
    *                                    /
    *                                   /
    *                                  / <--anchor_to_mapped
    *                                 v
    * mapped mapped_box_level:    (old) ---------> (new)
    *
    *
    * Output:
    *
    *                                    (anchor)
    *                                          \
    *                                           \
    *                                            \
    *                        anchor_to_mapped --> \
    *                                              v
    *                                              (new)
    *
    * @endcode
    *
    *
    * If mutable_new points to the old MappedBoxLevel and mutable_old
    * points to the new, then do an in-place switch as follows:
    * @li Swap the mutable_old and mutable_new MappedBoxLevels.
    * @li Use mutable_old (which actually the new MappedBoxLevel after
    *    the swap) as the mapped MappedBoxLevel in the output Connectors.
    * Otherwise:
    * @li If mutable_new is non-NULL, set it equal to new and use
    *    it as the mapped MappedBoxLevel in the output Connectors.
    * @li If mutable_old is non-NULL, set it equal to the old MappedBoxLevel.
    *
    * @param[in,out] anchor_to_mapped Connector to be modified.  On input, this
    *   points to the MappedBoxLevel being mapped.
    * @param[in,out] mapped_to_anchor Reverse (transpose) of anchor_to_mapped.
    *   points to the MappedBoxLevel being mapped.
    * @param[in] new_to_old Reverse (transpose) of old_to_new.
    * @param[in,out] mutable_new See comments.
    * @param[in,out] mutable_old See comments.
    */
   void
   modify(
      Connector& anchor_to_mapped,
      const Connector& old_to_new,
      MappedBoxLevel* mutable_new = NULL,
      MappedBoxLevel* mutable_old = NULL) const;

   /*
    * @brief Set whether to check for trivial mappings in modify
    * operations.
    *
    * Many common maps describe small (possibly zero) changes to a
    * bigger MappedBoxLevel.  If @c do_shortcut is set, we check whether
    * the mapping is globally empty.  (This costs a global reduction
    * unless it is already cached in the mapping Connector.)  If
    * empty, the shortcut improves performance.
    *
    * @param[in] do_shortcut
    */
   void
   shortcutTrivialMapsInModify(
      bool do_shortcut);

   // TODO:  Create an enum to replace the char for clarity.
   /*!
    * @brief Check if the Connector has a valid mapping.
    *
    * This function can be called prior to calling modify().  It
    * is intended to check the Connector that is the @c old_to_new
    * argument in modify() to determine if it has a valid mapping.
    * In other words, it checks to see if the mapping can be used
    * in modify() without logic errors.  It does no other checks.
    *
    * @param[in] connector
    * @param[in] is_local_map 'y' means assume the mapping is local.  'n'
    * means the mapping is not local.  '\0' means find out whether the
    * map is local or not (communication required) and act
    * accordingly.
    *
    * @return number of errors found.
    */
   size_t
   findMappingErrors(
      const Connector& connector,
      char is_local_map = '\0') const;

   /*!
    * @brief Run findMappingErrors and abort if any errors are found.
    *
    * @param[in] connector
    * @param[in] is_local_map 'y' means assume the mapping is local.  'n'
    * means the mapping is not local.  '\0' means find out whether the
    * map is local or not (communication required) and act
    * accordingly.
    */
   void
   assertMappingValidity(
      const Connector& connector,
      char is_local_map = '\0') const;

   //@}

private:
   static const int MAPPING_CONNECTOR_ALGORITHM_FIRST_DATA_LENGTH;

   // Internal shorthand.
   typedef Connector::NeighborSet NeighborSet;

   /*!
    * @brief BoxIdSet is a clarifying typedef.
    */
   typedef std::set<BoxId> BoxIdSet;

   /*!
    * @brief Mapping from a (potentially remote) Box to a
    * set of BoxIds, representing an inverted information
    * from a NeighborhoodSet.
    */
   typedef std::map<Box, BoxIdSet, Box::id_less> InvertedNeighborhoodSet;

   /*!
    * @brief Most general version of method to modify existing
    * Connectors objects by using another to map the head mapped_boxes.
    *
    * This version does no checking of the inputs.  The three
    * public versions do input checking and setting up temporaries
    * (where needed), then call this function.
    *
    * If new_to_old is uninitialized, treat it as a dummy
    * and assume that all mappings are local.
    *
    * If mutable_new points to the old MappedBoxLevel and mutable_new
    * points to the new, then do an in-place switch as follows:
    * -# Swap the old and new MappedBoxLevels.
    * -# Use mutable_old (which actually the new MappedBoxLevel after
    *    the swap) as the mapped MappedBoxLevel in the output Connectors.
    * Otherwise:
    * -# If mutable_new is non-NULL, set it equal to new and use
    *    it as the mapped MappedBoxLevel in the output Connectors.
    * -# If mutable_old is non-NULL, set it equal to the old MappedBoxLevel.
    */
   void
   privateModify(
      Connector& anchor_to_mapped,
      Connector& mapped_to_anchor,
      const Connector& old_to_new,
      const Connector& new_to_old,
      MappedBoxLevel* mutable_new,
      MappedBoxLevel* mutable_old) const;

   /*!
    * @brief Set up communication objects for use in privateModify.
    */
   void
   privateModify_setupCommunication(
      tbox::AsyncCommPeer<int> *& all_comms,
      tbox::AsyncCommStage& comm_stage,
      tbox::AsyncCommStage::MemberVec& completed,
      const tbox::SAMRAI_MPI &mpi,
      const std::set<int>& incoming_ranks,
      const std::set<int>& outgoing_ranks) const;

   /*!
    * @brief Relationship removal part of modify algorithm, caching
    * outgoing information in message buffers.
    */
   void
   privateModify_removeAndCache(
      std::map<int, std::vector<int> >& neighbor_removal_mesg,
      NeighborhoodSet& anchor_eto_new,
      NeighborhoodSet& new_eto_anchor,
      const Connector& old_to_new) const;

   /*!
    * @brief Remove relationships made obsolete by mapping and send
    * outgoing information.
    */
   void
   privateModify_discoverAndSend(
      std::map<int, std::vector<int> >& neighbor_removal_mesg,
      Connector& anchor_to_new,
      Connector& new_to_anchor,
      NeighborhoodSet& anchor_eto_new,
      NeighborhoodSet& new_eto_anchor,
      std::set<int>& incoming_ranks,
      std::set<int>& outoing_ranks,
      tbox::AsyncCommPeer<int> all_comms[],
      tbox::AsyncCommStage::MemberVec& completed,
      const Connector& old_to_anchor,
      const Connector& anchor_to_old,
      const Connector& old_to_new) const;

   /*!
    * @brief Receive messages and unpack info sent from other processes.
    */
   void
   privateModify_receiveAndUnpack(
      NeighborhoodSet& anchor_eto_new,
      NeighborhoodSet& new_eto_anchor,
      std::set<int>& outgoing_ranks,
      tbox::AsyncCommPeer<int> all_comms[],
      tbox::AsyncCommStage& comm_stage,
      tbox::AsyncCommStage::MemberVec& completed,
      const tbox::Dimension& dim,
      const tbox::SAMRAI_MPI &mpi) const;

   /*
    * @brief Performn checks on the arguments of modify.
    */
   void
   checkModifyParameters(
      const Connector& anchor_to_mapped,
      const Connector& mapped_to_anchor,
      const Connector& old_to_new,
      const Connector& new_to_old) const;

   /*!
    * @brief Set up things for the entire class.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   initializeCallback();

   /*!
    * Free statics.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   finalizeCallback();

   /*
    * @brief Border for debugging output.
    */
   static const std::string s_dbgbord;

   /*!
    * @brief Private communicator object shared by all objects in class.
    * Protects internal communication from mixing with external.
    *
    * For communication, we usually use the SAMRAI_MPI of the
    * Connectors, and this object is mainly for debugging.  If we
    * suspect interference from unrelated communication calls, we
    * resort to this exclusive SAMRAI_MPI object to rule out that
    * possibility.
    */
   static tbox::SAMRAI_MPI s_class_mpi;
   /*!
    * @brief Tag to use (and increment) at begining of operations that
    * require nearest-neighbor communication, to aid in eliminating
    * mixing of messages from different internal operations.
    */
   static int s_operation_mpi_tag;

   // Extra checks independent of optimization/debug.
   static char s_print_modify_steps;

   static tbox::Pointer<tbox::Timer> t_modify;
   static tbox::Pointer<tbox::Timer> t_modify_shortcut;
   static tbox::Pointer<tbox::Timer> t_modify_setup_comm;
   static tbox::Pointer<tbox::Timer> t_modify_remove_and_cache;
   static tbox::Pointer<tbox::Timer> t_modify_discover_and_send;
   static tbox::Pointer<tbox::Timer> t_modify_receive_and_unpack;
   static tbox::Pointer<tbox::Timer> t_modify_MPI_wait;

   bool d_sanity_check_inputs;
   bool d_sanity_check_outputs;
   bool d_shortcut_trivial_maps;

   static tbox::StartupShutdownManager::Handler
   s_initialize_finalize_handler;

};

}
}

#endif // included_hier_MappingConnectorAlgorithm
