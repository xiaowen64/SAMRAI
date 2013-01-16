/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2013 Lawrence Livermore National Security, LLC
 * Description:   Algorithms to work with MapingConnectors.
 *
 ************************************************************************/
#ifndef included_hier_MappingConnectorAlgorithm
#define included_hier_MappingConnectorAlgorithm

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/BaseConnectorAlgorithm.h"
#include "SAMRAI/hier/MappingConnector.h"

#include <map>
#include <string>

namespace SAMRAI {
namespace hier {

/*!
 * @brief Algorithms for using MappingConnectors representing changes
 * to a BoxLevel.
 *
 * MappingConnectorAlgorithm objects check and apply mappings.
 */
class MappingConnectorAlgorithm : public BaseConnectorAlgorithm
{

public:
   /*!
    * @brief Constructor.
    *
    * The default constructor creates an uninitialized object in
    * distributed state.
    */
   MappingConnectorAlgorithm();

   /*!
    * @brief Deallocate internal data.
    */
   virtual ~MappingConnectorAlgorithm();

   /*!
    * @brief Read extra debugging flag from input database.
    */
   void
   getFromInput();

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
      bool do_check)
   {
      d_sanity_check_inputs = do_check;
   }

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
      bool do_check)
   {
      d_sanity_check_outputs = do_check;
   }

   /*!
    * @brief Most general version for modifying Connectors using
    * MappingConnectors.  Modification is the changing of existing
    * Connectors when boxes in a BoxLevel changes according to specified
    * MappingConnectors.
    *
    * The change is represented by a the mapper @c old_to_new
    * and its transpose @c new_to_old.  The Connectors to be
    * modified are @c anchor_to_mapped and @c mapped_to_anchor, which
    * on input, go between a anchor (not mapped) BoxLevel and the old
    * BoxLevel.  On output, these Connectors will go from the anchor box_level
    * to the new box_level.
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
    * mapped box_level:           (old) ---------> (new)
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
    * The BoxLevels before and after the mapping are represented
    * by the base and head of the mapping.  No BoxLevel is modified,
    * other than the "mutable" BoxLevels in the argument.
    *
    * An important constraint in the old_to_new MappingConnectors is
    * that this method cannot handle multiple maps at once.  For
    * example, it cannot map box J to box K and at the
    * same time map box I to box J.  Box J in the
    * old box_level and box J on the new
    * box_level are considered entirely different boxes.
    *
    * After modifying, the output Connectors that had referenced old
    * BoxLevels will be reset to reference the new
    * BoxLevel.  This is the end of the modify operation.
    *
    * The following "in-place" modification is provided for users'
    * convenience.  Often times, the "new" level is a temporary object
    * and users often reset output equivalent to:
    *
    * @li old = new
    * @li reset Connectors using old in place of new.
    *
    * The modify methods support these optional steps as follows: If
    * mutable versions of some BoxLevel are given, the output
    * Connectors can be reset to reference these versions instead.
    *
    * If mutable_new points to the old BoxLevel and mutable_old
    * points to the new, then do an in-place switch as follows:
    * @li Swap the mutable_old and mutable_new BoxLevels.
    * @li Use mutable_old (which actually the new BoxLevel after
    *    the swap) as the mapped BoxLevel in the output Connectors.
    * Otherwise:
    * @li If mutable_new is non-NULL, set it equal to new and use
    *    it as the mapped BoxLevel in the output Connectors.
    * @li If mutable_old is non-NULL, set it equal to the old BoxLevel.
    *
    * @param[in,out] anchor_to_mapped Connector to be modified.  On input, this
    *   points to the BoxLevel being mapped.
    * @param[in,out] mapped_to_anchor Reverse (transpose) of anchor_to_mapped.
    *   points to the BoxLevel being mapped.
    * @param[in] old_to_new Mapping from the old BoxLevel to the
    *   new BoxLevel.
    *   The width of old-->new should indicate the maximum amount of box
    *   growth caused by the change.  A value of zero means no growth.
    * @param[in] new_to_old Reverse (transpose) of old_to_new.
    * @param[in,out] mutable_new See comments.
    * @param[in,out] mutable_old See comments.
    *
    * @pre (&old_to_new.getBase() == &new_to_old.getHead()) &&
    *      (&old_to_new.getBase() == &anchor_to_mapped.getHead()) &&
    *      (&old_to_new.getBase() == &mapped_to_anchor.getBase())
    * @pre &anchor_to_mapped.getBase() == &mapped_to_anchor.getHead()
    * @pre &old_to_new.getHead() == &new_to_old.getBase()
    * @pre anchor_to_mapped.isTransposeOf(mapped_to_anchor)
    * @pre new_to_old.isTransposeOf(old_to_new)
    * @pre anchor_to_mapped.getParallelState() == BoxLevel::DISTRIBUTED
    */
   void
   modify(
      Connector& anchor_to_mapped,
      Connector& mapped_to_anchor,
      const MappingConnector& old_to_new,
      const MappingConnector& new_to_old,
      BoxLevel* mutable_new = 0,
      BoxLevel* mutable_old = 0) const;

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
    * mapped box_level:           (old) ---------> (new)
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
    * mutable versions of some BoxLevel are given, the output
    * Connectors can be reset to reference these versions instead.
    *
    * If mutable_new points to the old BoxLevel and mutable_old
    * points to the new, then do an in-place switch as follows:
    * @li Swap the mutable_old and mutable_new BoxLevels.
    * @li Use mutable_old (which actually the new BoxLevel after
    *    the swap) as the mapped BoxLevel in the output Connectors.
    * Otherwise:
    * @li If mutable_new is non-NULL, set it equal to new and use
    *    it as the mapped BoxLevel in the output Connectors.
    * @li If mutable_old is non-NULL, set it equal to the old BoxLevel.
    *
    * @param[in,out] anchor_to_mapped Connector to be modified.  On input, this
    *   points to the BoxLevel being mapped.
    * @param[in,out] mapped_to_anchor Reverse (transpose) of anchor_to_mapped.
    *   points to the BoxLevel being mapped.
    * @param[in] old_to_new Mapping from the old BoxLevel to the
    *   new BoxLevel.
    *   The width of old-->new should indicate the maximum amount of box
    *   growth caused by the change.  A value of zero means no growth.
    * @param[in,out] mutable_new See comments.
    * @param[in,out] mutable_old See comments.
    *
    * @pre (&old_to_new.getBase() == &old_to_new.getBase()) &&
    *      (&old_to_new.getBase() == &anchor_to_mapped.getHead())
    * @pre &anchor_to_mapped.getBase() == &mapped_to_anchor.getHead()
    * @pre anchor_to_mapped.isTransposeOf(mapped_to_anchor)
    * @pre anchor_to_mapped.getParallelState() == BoxLevel::DISTRIBUTED
    */
   void
   modify(
      Connector& anchor_to_mapped,
      Connector& mapped_to_anchor,
      const MappingConnector& old_to_new,
      BoxLevel* mutable_new = 0,
      BoxLevel* mutable_old = 0) const;

   /*!
    * @brief Version of modify requiring only the forward map
    * and does not set the Connector from the mapped BoxLevel back
    * to the anchor BoxLevel.
    *
    * This version does not require the reverse mapping,
    * but the forward mapping must be totally local, e.g.
    * old_to_new must contain no remote neighbor.
    *
    * This version does not update the Connector from
    * the mapped box_level back to the anchor box_level.
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
    * mapped box_level:           (old) ---------> (new)
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
    * If mutable_new points to the old BoxLevel and mutable_old
    * points to the new, then do an in-place switch as follows:
    * @li Swap the mutable_old and mutable_new BoxLevels.
    * @li Use mutable_old (which actually the new BoxLevel after
    *    the swap) as the mapped BoxLevel in the output Connectors.
    * Otherwise:
    * @li If mutable_new is non-NULL, set it equal to new and use
    *    it as the mapped BoxLevel in the output Connectors.
    * @li If mutable_old is non-NULL, set it equal to the old BoxLevel.
    *
    * @param[in,out] anchor_to_mapped Connector to be modified.  On input, this
    *   points to the BoxLevel being mapped.
    * @param[in,out] old_to_new Reverse (transpose) of anchor_to_mapped.
    *   points to the BoxLevel being mapped.
    * @param[in,out] mutable_new See comments.
    * @param[in,out] mutable_old See comments.
    *
    * @pre &anchor_to_mapped.getHead() == &old_to_new.getBase()
    * @pre anchor_to_mapped.getParallelState() == BoxLevel::DISTRIBUTED
    */
   void
   modify(
      Connector& anchor_to_mapped,
      const MappingConnector& old_to_new,
      BoxLevel* mutable_new = 0,
      BoxLevel* mutable_old = 0) const;

   /*!
    * @brief Get the name of this object.
    */
   const std::string
   getObjectName() const
   {
      return "MappingConnectorAlgorithm";
   }

private:
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
    * Connectors objects by using MappingConnectors to map the head boxes.
    *
    * This version does no checking of the inputs.  The three
    * public versions do input checking and setting up temporaries
    * (where needed), then call this function.
    *
    * If new_to_old is uninitialized, treat it as a dummy
    * and assume that all mappings are local.
    *
    * If mutable_new points to the old BoxLevel and mutable_new
    * points to the new, then do an in-place switch as follows:
    * -# Swap the old and new BoxLevels.
    * -# Use mutable_old (which actually the new BoxLevel after
    *    the swap) as the mapped BoxLevel in the output Connectors.
    * Otherwise:
    * -# If mutable_new is non-NULL, set it equal to new and use
    *    it as the mapped BoxLevel in the output Connectors.
    * -# If mutable_old is non-NULL, set it equal to the old BoxLevel.
    */
   void
   privateModify(
      Connector& anchor_to_mapped,
      Connector& mapped_to_anchor,
      const MappingConnector& old_to_new,
      const MappingConnector* new_to_old,
      BoxLevel* mutable_new,
      BoxLevel* mutable_old) const;

   /*
    * @brief Perform checks on the arguments of modify.
    */
   void
   privateModify_checkParameters(
      const Connector& anchor_to_mapped,
      const Connector& mapped_to_anchor,
      const MappingConnector& old_to_new,
      const MappingConnector* new_to_old) const;

   /*!
    * @brief Relationship removal part of modify algorithm, caching
    * outgoing information in message buffers.
    */
   void
   privateModify_removeAndCache(
      std::map<int, std::vector<int> >& neighbor_removal_mesg,
      Connector& anchor_to_new,
      Connector* new_to_anchor,
      const MappingConnector& old_to_new) const;

   /*!
    * @brief Remove relationships made obsolete by mapping and send
    * outgoing information.
    */
   void
   privateModify_discoverAndSend(
      std::map<int, std::vector<int> >& neighbor_removal_mesg,
      Connector& anchor_to_new,
      Connector* new_to_anchor,
      std::set<int>& incoming_ranks,
      std::set<int>& outoing_ranks,
      tbox::AsyncCommPeer<int> all_comms[],
      BoxContainer& visible_new_nabrs,
      BoxContainer& visible_anchor_nabrs,
      InvertedNeighborhoodSet& anchor_eto_old,
      InvertedNeighborhoodSet& new_eto_old,
      const Connector& old_to_anchor,
      const Connector& anchor_to_old,
      const MappingConnector& old_to_new) const;

   /*!
    * @brief Find overlap and save in mapping connector or pack
    * into send message, used in privateModify().
    */
   void
   privateModify_findOverlapsForOneProcess(
      const int owner_rank,
      BoxContainer& visible_base_nabrs,
      BoxContainer::iterator& base_ni,
      std::vector<int>& send_mesg,
      const int remote_box_counter_index,
      Connector& mapped_connector,
      BoxContainer& referenced_head_nabrs,
      const Connector& unmapped_connector,
      const Connector& unmapped_connector_transpose,
      const Connector& mapping_connector,
      InvertedNeighborhoodSet& inverted_nbrhd,
      const IntVector& refinement_ratio) const;

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

   // Extra checks independent of optimization/debug.
   static char s_print_steps;

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
    *
    * Using unique tags was important to the early development of this
    * class but may not be necessary anymore.
    */
   static int s_operation_mpi_tag;

   static boost::shared_ptr<tbox::Timer> t_modify;
   static boost::shared_ptr<tbox::Timer> t_modify_setup_comm;
   static boost::shared_ptr<tbox::Timer> t_modify_remove_and_cache;
   static boost::shared_ptr<tbox::Timer> t_modify_discover_and_send;
   static boost::shared_ptr<tbox::Timer> t_modify_receive_and_unpack;
   static boost::shared_ptr<tbox::Timer> t_modify_MPI_wait;
   static boost::shared_ptr<tbox::Timer> t_modify_misc;

   bool d_sanity_check_inputs;
   bool d_sanity_check_outputs;

   static tbox::StartupShutdownManager::Handler
      s_initialize_finalize_handler;

};

}
}

#endif // included_hier_MappingConnectorAlgorithm
