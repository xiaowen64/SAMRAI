/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Set of distributed box-graph relationships from one MappedBoxLevel to another. 
 *
 ************************************************************************/
#ifndef included_hier_Connector
#define included_hier_Connector

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/MappedBox.h"
#include "SAMRAI/hier/MappedBoxLevel.h"
#include "SAMRAI/hier/NeighborhoodSet.h"
#include "SAMRAI/tbox/Timer.h"

#include <vector>
#include <string>

namespace SAMRAI {
namespace hier {

class MappedBoxLevelHandle;
class NeighborhoodSet;
class SAMRAI_MPI;

/*!
 * @brief A container which holds relationship connections between two
 * MappedBoxLevels.
 *
 * Connectors have a notion of a "base" and a "head", representing a
 * directional relationship between two MappedBoxLevels.  The relationships
 * are a collection of MappedBoxes in the head, pointed to by a MappedBox in
 * the base.   The association between base and head relationships is
 * 1 .. 0-many.  That is, one MappedBox in the base can be related to zero or
 * more MappedBoxes (called its NeighborSet) in the head.
 *
 * @par Usage
 * Connections in a Connector can have three possible relationships:
 *
 * # A MappedBox in the base has no related NeighborSet.  In this case,
 *   the MappedBox in the base will exist as-is in the head.
 * # A MappedBox in the base has NeighborSet which is empty.  In this case, the
 *   MappedBox from the base will not exist in the head MappedBoxLevel.
 * # A MappedBox in the base has a corresponding NeighborSet which in
 *   non-empty.  In this case, the NeighborSet contains the set of MappedBoxes
 *   to which the MappedBox in the base is related.
 */

class Connector:public tbox::DescribedClass
{
public:
   /*!
    * @brief NeighborsSet is a clarifying typedef.
    */
   typedef MappedBoxSet NeighborSet;


   /// TODO:  Possible refactor?  Since Connectors do not imply relationship
   // meanings, why is this even defined?  The "getConnectorType function
   // is never called; The individual types are never used except in
   // conjunction with the "setConnectorType" function.  This seems to be
   // a useless enumeration producing unused code.  SGH.

   /*!
    * @brief Types of Connectors
    *
    * The types describe the meaning of the relationships in a Connector.
    *
    * @b COMPLETE_OVERLAP: The relationships represent overlaps, and every
    * overlap is represented by an relationship, including overlaps with
    * periodic images.
    *
    * @b COMPLETE_OVERLAP_NO_PERIODIC: The relationships represent overlaps,
    * and every overlap is represented by an relationship.  Overlaps with
    * periodic images are omitted.
    *
    * @b INCOMPLETE_OVERLAP: The relationships represent overlaps, but not
    * all overlaps are represented.
    *
    * @b BASE_GENERATED: The head is generated from the base.  Each
    * head MappedBox comes from a base MappedBox and there is an relationship
    * from the base MappedBox to the head MappedBox.
    *
    * @b MAPPING: relationships indicate a mapping relationship.  Applying
    * the map would change Connectors incident to the base into
    * Connectors incident to the head.
    *
    * @b UNKNOWN: Meaning of relationships are unknown.
    *
    * See setConnectorType(), getConnectorType().
    *
    * The Connector types are not exclusive.  For example, a mapping
    * Connector may also be used as an overlap Connector.
    */
   enum ConnectorType {
      COMPLETE_OVERLAP = 1,
      COMPLETE_OVERLAP_NO_PERIODIC = 2,
      INCOMPLETE_OVERLAP = 3,
      BASE_GENERATED = 4,
      MAPPING = 5,
      UNKNOWN = 6
   };

   /*!
    * @brief Creates an uninitialized Connector object in the
    * distributed state.
    *
    * @see initialize()
    * @see swapInitialize()
    */
   Connector();

   /*!
    * @brief Copy constructor.
    *
    * @param[in] other
    */
   Connector(
      const Connector& other);

   /*!
    * @brief Creates an initialized Connector with the given relationships.
    *
    * @param[in] base_mapped_box_level
    * @param[in] head_mapped_box_level
    * @param[in] base_width
    * @param[in] relationships
    * @param[in] parallel_state
    */
   explicit Connector(
      const MappedBoxLevel& base_mapped_box_level,
      const MappedBoxLevel& head_mapped_box_level,
      const IntVector& base_width,
      const NeighborhoodSet& relationships,
      const MappedBoxLevel::ParallelState parallel_state = MappedBoxLevel::DISTRIBUTED);

   /*!
    * @brief Initialize a Connector with no defined relationships.
    *
    * The Connector's relationships are initialized to a dummy state.
    *
    * @param[in] base_mapped_box_level
    * @param[in] head_mapped_box_level
    * @param[in] base_width
    * @param[in] parallel_state
    */
   explicit Connector(
      const MappedBoxLevel& base_mapped_box_level,
      const MappedBoxLevel& head_mapped_box_level,
      const IntVector& base_width,
      const MappedBoxLevel::ParallelState parallel_state = MappedBoxLevel::DISTRIBUTED);

   /*!
    * @brief Destructor.
    */
   virtual ~Connector();

   /*!
    * @brief Initializes the Connector.
    *
    * The Connector is initialized with the provided information.
    *
    * Consistency with the base and head requires that the relationships be
    * from an existing MappedBox in the base to an existing MappedBox
    * in the head.
    *
    * @par Assertions
    * With assertion checking turned on, this function will check for
    * base consistency only.  Consistency with the head is not checked
    * explicitly.
    *
    * @param[in] base
    * @param[in] head
    * @param[in] base_width The Connector width associated with the meaning of the
    *   relationships, specified in base refinement ratio.
    * @param[in] parallel_state Either DISTRIBUTED or GLOBALIZED.
    *   If state is GLOBALIZED, base must be in GLOBALIZED mode.
    * @param[in] relationships The input relationship data.  For DISTRIBUTED state, we
    *   disregard relationships from remote MappedBoxes.
    *
    * @see initializePrivate()
    * @see checkConsistencyWithBase()
    * @see checkConsistencyWithHead()
    *
    */
   void
   initialize(
      const MappedBoxLevel& base,
      const MappedBoxLevel& head,
      const IntVector& base_width,
      const NeighborhoodSet& relationships,
      const MappedBoxLevel::ParallelState parallel_state = MappedBoxLevel::DISTRIBUTED);

   /*!
    * @brief Initializes the Connector without any relationships.
    *
    * @param[in] base
    * @param[in] head
    * @param[in] base_width The Connector width, specified in the base
    *   refinement ratio, associated with the meaning of the relationships to be
    *   added to the Connector.
    * @param[in] parallel_state Either DISTRIBUTED or GLOBALIZED.
    *   If state is GLOBALIZED, base must be in GLOBALIZED mode.
    *
    * @see swapInitialize()
    * @see initializePrivate()
    *
    */
   void
   initialize(
      const MappedBoxLevel& base,
      const MappedBoxLevel& head,
      const IntVector& base_width,
      const MappedBoxLevel::ParallelState parallel_state = MappedBoxLevel::DISTRIBUTED);

   /*!
    * @brief Set data defining the relationship set.
    *
    * @param[in] base
    * @param[in] head
    * @param[in] base_width
    * @param[out] relationships
    * @param[in] parallel_state
    *
    * POST-CONDITION:  @c relationships contains the NeighborhoodSet.
    *
    * @see initializePrivate()
    *
    */
   void
   swapInitialize(
      const MappedBoxLevel& base,
      const MappedBoxLevel& head,
      const IntVector& base_width,
      NeighborhoodSet& relationships,
      const MappedBoxLevel::ParallelState parallel_state = MappedBoxLevel::DISTRIBUTED);

   /*!
    * @brief Clear the Connector, putting it into an uninitialized state.
    */
   void
   clear();

   /*!
    * @brief Returns true if the object has been initialized
    */
   bool
   isInitialized() const;

   /*!
    * @brief Return relationships from local base MappedBoxes.
    */
   const NeighborhoodSet&
   getNeighborhoodSets() const;

   /*!
    * @brief Return the globalized relationship data.
    *
    * @par Assertions
    * Throws an unrecoverable assertion if not in GLOBALIZED mode.
    */
   const NeighborhoodSet&
   getGlobalNeighborhoodSets() const;

   /*!
    * @brief Return true if a neighbor set exists for the specified
    * MappedBoxId.
    *
    * @param[in] mapped_box_id
    */
   bool
   hasNeighborSet(
      const MappedBoxId& mapped_box_id) const;

   /*!
    * @brief Return the neighbor set for the specified MappedBoxId.
    *
    * @param[in] mapped_box_id
    */
   const NeighborSet&
   getNeighborSet(
      const MappedBoxId& mapped_box_id) const;

   //@{
   /*!
    * @name Algorithms for changing individual MappedBox's neighbor data
    */

   /*!
    * @brief Insert additional neighbors for the specified MappedBox.
    *
    * @param[in] neighbors
    * @param[in] mapped_box_id
    */
   void
   insertNeighbors(
      const NeighborSet& neighbors,
      const MappedBoxId& mapped_box_id);

   /*!
    * @brief Erase neighbor of the specified MappedBoxId.
    *
    * @note Assertions
    * It is an error to to specify a non-existent MappedBoxId.
    *
    * @param[in] neighbor
    * @param[in] mapped_box_id
    */
   void
   eraseNeighbor(
      const MappedBox& neighbor,
      const MappedBoxId& mapped_box_id);

   /*!
    * @brief Set the neighbors for the specified MappedBoxId to the
    * given set by swapping the sets.
    *
    * If no neighbor set exists for the specified MappedBoxId, an empty
    * set is first created for swapping.
    *
    * An assertion failure will occur if the MappedBoxId has a non-zero
    * PeriodicId.
    *
    * @param[out] neighbors
    * @param[in] mapped_box_id
    */
   void
   swapNeighbors(
      NeighborSet& neighbors,
      const MappedBoxId& mapped_box_id);

   /*!
    * @brief Remove empty sets of neighbors.
    */
   void
   eraseEmptyNeighborSets();

   //@}

   /*!
    * @brief Return a reference to the base MappedBoxLevel.
    */
   const MappedBoxLevel&
   getBase() const;

   /*!
    * @brief Return a reference to the head MappedBoxLevel.
    */
   const MappedBoxLevel&
   getHead() const;

   /*!
    * @brief Get the refinement ratio between the base and head
    * MappedBoxLevels.
    *
    * The ratio is the same regardless of which is the coarser of the two.
    * Use getHeadCoarserFlag() to determine which is coarser.  If the ratio
    * cannot be represented by an IntVector, truncated.  @see ratioIsExact().
    */
   const IntVector&
   getRatio() const;

   /*!
    * @brief Whether the ratio given by getRatio() is exact.
    *
    * The ratio is exact if it can be represented by an IntVector.
    * @see getRatio().
    */
   bool ratioIsExact() const;

   /*!
    * @brief Return true if head MappedBoxLevel is coarser than base
    * MappedBoxLevel.
    */
   bool
   getHeadCoarserFlag() const;

   /*!
    * @brief Return true if the Connector contains only relationships to local
    * MappedBoxes.
    *
    * The check only applies to neighbors of local base MappedBoxes,
    * so it is possible for the results to be different on different
    * processors.
    */
   bool
   isLocal() const;

   /*!
    * @brief Initialize to the transpose of a given Connector object,
    * assuming that all relationships are local (no remote neighbors).
    *
    * If any remote neighbor is found an unrecoverable assertion is
    * thrown.
    *
    * Non-periodic relationships in @c connector are simply reversed to get the
    * transpose relationship.  For each periodic relationships in @c connector, we create a
    * periodic relationship incident from @c connector's unshifted head neighbor to
    * @c connectors's shifted base neighbor.  This is because all relationships must be
    * incident from a real (unshifted) MappedBox.
    *
    * @param[in] connector
    */
   void
   initializeToLocalTranspose(
      const Connector& connector);

   /*!
    * @brief Assignment operator
    */
   const Connector&
   operator= (
      const Connector& rhs);

   //  TODO:  need to find out what the use case is for this, especially
   //  considering the caution statement.
   /*!
    * @brief Equality operator checks relationship data, Connector width and
    * equality of base and head MappedBox pointers.
    *
    * @par CAUTION
    * Equality here means just the local parts are equal.
    * This means that one processor may see the equality differently
    * from another.
    *
    * The cost for the comparison is on the order of the local relationship
    * count.  However, an object may be compared to itself, an
    * efficient operation that always returns true.  When comparing
    * Connector objects, if you expect equality to hold, using the
    * same objects would improve performance.
    */
   bool
   operator== (
      const Connector& rhs) const;

   /*!
    * @brief Inequality operator checks the same data that equality
    * operator checks.
    *
    * @see operator==( const Connector &rhs );
    */
   bool
   operator != (
      const Connector& rhs) const;

   /*!
    * @brief Swap the contents of two Connector objects.
    */
   static void
   swap(
      Connector& a,
      Connector& b);

   /*!
    * @brief Set the parallel distribution state.
    *
    * Before a Connector can be in a GLOBALIZED state, The base
    * MappedBoxLevel given in initialize() must already be in
    * GLOBALIZED mode.  The base MappedBoxLevel should remain in
    * GLOBALIZED mode for compatibility with the Connector.
    *
    * This method is not necessarily trivial.  More memory is required
    * to store additional relationships.
    *
    * For serial (one processor) runs, there is no difference between
    * the parallel states (except for the names), and there is no real
    * cost for switching parallel states.
    *
    * @param[in] parallel_state
    */
   void
   setParallelState(
      const MappedBoxLevel::ParallelState parallel_state);

   /*!
    * @brief Return the current parallel state.
    */
   MappedBoxLevel::ParallelState
   getParallelState() const;

   /*!
    * @brief Returns the MPI communication object, which is always
    * that of the base MappedBoxLevel.
    */
   const tbox::SAMRAI_MPI&
   getMPI() const;

   /*!
    * @return Processor rank in the internal MPI communicator.
    *
    * This is redundant information.  You can obtain it from the MPI
    * Communicator.
    */
   int
   getRank() const;

   /*!
    * @return Number of processes for the internal MPI communicator.
    *
    * This is redundant information.  You can obtain it from the MPI
    * Communicator.
    */
   int
   getNproc() const;

   /*!
    * @brief Return the Connector width associated with the relationships.
    *
    * For overlap Connectors, an relationship exists between a base and head
    * MappedBoxes if the base mapped_box, grown by this width,
    * overlaps the head mapped_box.  For mapping Connectors, the width
    * the amount that a pre-map box must grow to nest the post-map
    * boxes.
    */
   const IntVector&
   getConnectorWidth() const;

   //@{
   /*!
    * @name For outputs, error checking and debugging.
    */


   /*
    * @brief output data
    *
    * @param[in,out] os The output stream
    * @param[in] border
    * @param[in] detail_depth
    */
   void
   recursivePrint(
      std::ostream& os,
      const std::string& border,
      int detail_depth = 0) const;

   /*!
    * @brief Print out statistics on the relationships.
    *
    * Requires communication, so all processors must call this.
    *
    * @param[in,out] os The output stream
    * @param[in] border
    */
   void
   printNeighborStats(
      std::ostream& co,
      const std::string& border) const;

   /*!
    * @brief Return true if two Connector objects are
    * transposes of each other.
    *
    * Each Connector represents a set of directed relationships incident from its base
    * to its head.  The transpose represent the relationships in the opposite
    * direction.  In order for two Connector objects to be transpose of each
    * other, their Connector widths and base refinement ratios must be
    * such that an relationship in one set also appears in the other set.
    * A transpose set must have
    * @li base and head MappedBoxLevels reversed from the untransposed set.
    * @li the same Connector width, although it is described in the index
    *     space of a different base MappedBoxLevel.
    *
    * @param[in] other
    */
   bool
   isTransposeOf(
      const Connector& other) const;

   /*!
    * @brief Given the Connector width in the head index space, convert
    * it to the base index space.
    *
    * This method is useful for computing Connector widths for
    * transpose Connectors.  It handles negative refinement ratios. By
    * SAMRAI convention, a refinement ratio of -N is interpreted as
    * 1/N.)
    *
    * This method is static because (1) it has nothing to do with an
    * existing Connector object, and (2) it is often used to compute a
    * Connector's initializing data.
    *
    * @param[in] base_refinement_ratio
    * @param[in] head_refinement_ratio
    * @param[in] head_gcw The connector width in the head index space.
    *
    * @return A copy of the connector width converted to the base index
    * space.
    */
   static IntVector
   convertHeadWidthToBase(
      const IntVector& base_refinement_ratio,
      const IntVector& head_refinement_ratio,
      const IntVector& head_gcw);


   // TODO: refactor use of size_t as return type.  This could be
   // problematic.

   /*!
    * @brief Check for consistency between the relationship data and
    * base MappedBoxes, and return the number of consistency errors.
    *
    * Consistency stipulates that each neighbor list must correspond to
    * a base mapped box.
    *
    * Relationship consistency errors should be treated as fatal
    * because many operations assume consistency.
    *
    * @return number of inconsistencies found.
    */
   size_t
   checkConsistencyWithBase() const;

   /*!
    * @brief Run checkConsistencyWithBase().  If any inconsistency is
    * found, write out diagnostic information and throw an
    * unrecoverable assertion.
    */
   void
   assertConsistencyWithBase() const;

   /*!
    * @brief Check that the neighbors exist in the head
    * MappedBoxLevel.
    *
    * Consistency stipulates that each neighbor referenced in the
    * Connector must exist in the head.
    *
    * Relationship consistency errors should be treated as fatal
    * because many operations assume consistency.
    *
    * @return number of inconsistencies found.
    */

   size_t
   checkConsistencyWithHead() const;

   /*!
    * @brief Run checkConsistencyWithBase().  If any inconsistency is
    * found, write out diagnostic information and throw an
    * unrecoverable assertion.
    */
   void
   assertConsistencyWithHead() const;

   /*!
    * @brief Check that MappedBoxes referenced by the given
    * NeighborhoodSet match those in the given MappedBoxLevel.
    *
    * This method is static so users can check data without having to
    * put it in a Connector object.  Connectors prohibit initializing
    * with inconsistent data.
    *
    * @param[in] relationships
    *
    * @param[in] head_mapped_box_level
    *
    * @return number of inconsistencies found.
    */
   static size_t
   checkConsistencyWithHead(
      const NeighborhoodSet& relationships,
      const MappedBoxLevel& head_mapped_box_level);

   /*!
    * @brief Compute the differences between two relationship sets.
    *
    * Given Connectors @c left_connector and @c right_connector,
    * compute the relationships that are in @c left_connector but not in
    * @c right_connector.
    *
    * @param[out] left_minus_right
    * @param[in] left_connector
    * @param[in] right_connector
    */
   static void
   computeNeighborhoodDifferences(
      Connector& left_minus_right,
      const Connector& left_connector,
      const Connector& right_connector);

   /*!
    * @brief Check that the relationships are a correct transpose of another
    * Connector and return the number of erroneous relationships.
    *
    * For every relationship in this Connector, there should be a corresponding relationship
    * in the transpose Connector.  Any missing or extra relationship constitutes
    * an error.
    *
    * Errors found are written to perr.
    *
    * @param[in] transpose
    * @param[in] ignore_periodic_relationships
    *
    * @return Number of errors in assuming that @c transpose is a
    * transpose of @c *this.
    */
   size_t
   checkTransposeCorrectness(
      const Connector& transpose,
      const bool ignore_periodic_relationships = false) const;

   /*!
    * @brief Run checkTransposeCorrectness.  If any errors are found,
    * print out diagnostic information and throw an unrecoverable
    * assertion.
    *
    * @param[in] transpose
    * @param[in] ignore_periodic_relationships
    */
   void
   assertTransposeCorrectness(
      const Connector& transpose,
      const bool ignore_periodic_relationships = false) const;

   //@}

   /*!
    * @brief Set the Connector type.
    *
    * @param[in] connector_type
    */
   void
   setConnectorType(
      ConnectorType connector_type);

   /*!
    * @brief Return the Connector type.
    */
   ConnectorType
   getConnectorType() const;

   // TODO:  refactor size_t
   /*!
    * @brief Return local number of neighbor sets.
    */
   size_t
   getLocalNumberOfNeighborSets() const;

   // TODO: refactor size_t
   /*!
    * @brief Return local number of relationships.
    */
   size_t
   getLocalNumberOfRelationships() const;

   /*!
    * @brief Return global number of neighbor sets.
    *
    * This requires a global sum reduction, if the global size has not
    * been computed and cached.  When communication is required, all
    * processors must call this method.  To ensure that no
    * communication is needed, call cacheGlobalReducedData() first.
    */
   int
   getGlobalNumberOfNeighborSets() const;

   /*!
    * @brief Return global number of relationships.
    *
    * This requires a global sum reduction if the global size has not
    * been computed and cached.  When communication is required, all
    * processors must call this method.  To ensure that no
    * communication is needed, call cacheGlobalReducedData() first.
    */
   int
   getGlobalNumberOfRelationships() const;

   /*!
    * @brief If global reduced data (global number of relationships,
    * etc.) has not been updated, compute and cache them
    * (communication required).
    *
    * After this method is called, data requiring global reduction can
    * be accessed without further communications, until the object
    * changes.
    *
    * Sets d_global_data_up_to_date;
    */
   void
   cacheGlobalReducedData() const;

   /*!
    * @brief A class for outputting Connector.
    *
    * This class simplifies the insertion of a Connector into a stream
    * while letting the user control how the Connector is formatted
    * for output.
    *
    * Each Outputter is a light-weight object constructed with a
    * Connector and output parameters.  The Outputter is capable of
    * outputting its Connector, formatted according to the parameters.
    *
    * To use, @see Connector::format().
    */
   class Outputter {
      friend std::ostream& operator << ( std::ostream& s, const Outputter& f);
   private:
      friend class Connector;
      /*!
       * @brief Construct the Outputter with a Connector and the
       * parameters needed to output the Connector to a stream.
       */
      Outputter( const Connector &connector,
                 const std::string& border,
                 int detail_depth = 0);
      void operator=( const Outputter &r ); // Unimplemented private.
      const Connector &d_conn;
      const std::string d_border;
      const int d_detail_depth;
   };

   /*!
    * @brief Return an Outputter object that is formatted
    * insertion into output streams.
    *
    * Usage example:
    * @code
    *    cout << "my connector:\n"
    *         << connector.format("  ", 2) << endl;
    * @endcode
    *
    * @param[in] border
    * @param[in] detail_depth
    */
   Outputter format( const std::string& border=std::string(),
                     int detail_depth = 0 ) const;

private:
   /*
    * Static integer constant descibing class's version number.
    */
   static const int HIER_CONNECTOR_VERSION;

   enum { BAD_INT = (1 << (8 * sizeof(int) - 2)) };

   /*
    * @brief Create a copy of a DISTRIBUTED Connector and
    * change its state to GLOBALIZED.
    *
    * The returned object should be deleted to prevent memory leaks.
    */
   Connector *
   makeGlobalizedCopy(
      const Connector& other) const;

   /*!
    * @brief Get and store info on remote MappedBoxes.
    *
    * This requires global communication (all gather).
    * Call acquireRemoteNeighborhoods_pack to pack up messages.
    * Do an all-gather.  Call acquireRemoteNeighborhoods_unpack
    * to unpack data from other processors.
    */
   void
   acquireRemoteNeighborhoods();

   //! @brief Pack local MappedBoxes into an integer array.
   void
   acquireRemoteNeighborhoods_pack(
      std::vector<int>& send_mesg,
      int offset) const;

   //! @brief Unpack MappedBoxes from an integer array into internal storage.
   void
   acquireRemoteNeighborhoods_unpack(
      const std::vector<int>& recv_mesg,
      const std::vector<int>& proc_offset);

   /*! @brief Encapsulates functionality common to all initialization
    *  functions.
    */
   void
   initializePrivate(
      const MappedBoxLevel& base,
      const MappedBoxLevel& head,
      const IntVector& base_width,
      const IntVector& baseRefinementRatio,
      const IntVector& headRefinementRatio,
      const MappedBoxLevel::ParallelState parallel_state = MappedBoxLevel::DISTRIBUTED);

   /*!
    * @brief Set up things for the entire class.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   initializeCallback();

   /*!
    * Free static timers.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   finalizeCallback();

   //@{ @name Private utilities.

   //@}

   /*!
    * @brief Handle for access to the base MappedBoxLevel.
    *
    * We don't use a pointer to the MappedBoxLevel, because it would
    * become dangling when the MappedBoxLevel goes out of scope.
    */
   tbox::Pointer<MappedBoxLevelHandle> d_base_handle;

   /*!
    * @brief Handle for access to the base MappedBoxLevel.
    *
    * We don't use a pointer to the MappedBoxLevel, because it would
    * become dangling when the MappedBoxLevel goes out of scope.
    */
   tbox::Pointer<MappedBoxLevelHandle> d_head_handle;

   /*!
    * @brief Connector width for the base MappedBoxLevel.
    *
    * This is the amount of growth applied to a mapped_box in the base MappedBoxLevel
    * before checking if the mapped_box overlaps a mapped_box in the head MappedBoxLevel.
    */
   IntVector d_base_width;

   /*!
    * @brief Refinement ratio between base and head.
    *
    * If d_head_coarser is false, the head is not coarser than
    * the base and this is the refinement ratio from base to head.
    * If d_head_coarser is true, this is the coarsen ratio
    * from base to head.
    *
    * This is redundant information.  You can compute it
    * from the base and head MappedBoxLevels.
    */
   IntVector d_ratio;

   /*!
    * @brief Whether the ratio between the base and head
    * MappedBoxLevel refinement ratios are exactly as given by
    * d_ratio.  It can only be exact if it can be represented as an
    * IntVector.
    */
   bool d_ratio_is_exact;

   /*!
    * @brief Whether the base MappedBoxLevel is at a finer index space.
    *
    * When this is true, d_ratio is the refinement ratio going
    * from the head to the base.
    *
    * This is redundant information.  You can compute it
    * from the base and head MappedBoxLevels.
    */
   bool d_head_coarser;

   /*!
    * @brief Neighbor data for local MappedBoxes.
    */
   NeighborhoodSet d_relationships;

   /*!
    * @brief Neighbor data for global MappedBoxes in GLOBALIZED mode.
    */
   NeighborhoodSet d_global_relationships;

   /*!
    * @brief State flag.
    *
    * Modified by setParallelState().
    */
   MappedBoxLevel::ParallelState d_parallel_state;

   /*!
    * @brief Number of NeighborSets in d_relationships globally.
    */
   mutable int d_global_number_of_neighbor_sets;

   /*!
    * @brief Number of relationships in d_relationships globally.
    */
   mutable int d_global_number_of_relationships;

   /*!
    * @brief Whether globally reduced data is up to date or needs
    * recomputing using cacheGlobalReducedData().
    */
   mutable bool d_global_data_up_to_date;

   /*!
    * @brief Process rank (id),
    * for convenience and data management use after MPI_Finalize.
    */
   int d_rank;
   /*!
    * @brief Number of processes,
    * for convenience and data management use after MPI_Finalize.
    *
    * If d_nproc == BAD_INT, the object is in uninitialized state.
    */
   int d_nproc;

   ConnectorType d_connector_type;

   static tbox::Pointer<tbox::Timer> t_initialize;
   static tbox::Pointer<tbox::Timer> t_acquire_remote_relationships;

   static tbox::StartupShutdownManager::Handler
   s_initialize_finalize_handler;

};

}
}

#ifdef SAMRAI_INLINE
#include "SAMRAI/hier/Connector.I"
#endif

#endif // included_hier_Connector
