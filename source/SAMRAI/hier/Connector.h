/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Set of distributed box-graph relationships from one BoxLevel to another.
 *
 ************************************************************************/
#ifndef included_hier_Connector
#define included_hier_Connector

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/BoxLevel.h"
#include "SAMRAI/hier/NeighborhoodSet.h"
#include "SAMRAI/tbox/Timer.h"

#include <set>
#include <string>
#include <vector>

namespace SAMRAI {
namespace hier {

class BoxLevelHandle;

/*!
 * @brief A container which holds relationship connections between two
 * BoxLevels.
 *
 * Connectors have a notion of a "base" and a "head", representing a
 * directional relationship between two BoxLevels.  The relationships
 * are a collection of Boxes in the head, pointed to by a Box in
 * the base.   The association between base and head relationships is
 * 1 .. 0-many.  That is, one Box in the base can be related to zero or
 * more Boxes (called its NeighborSet) in the head.
 *
 * @par Usage
 * Connections in a Connector can have three possible relationships:
 *
 * # A Box in the base has no related NeighborSet.  In this case,
 *   the Box in the base will exist as-is in the head.
 * # A Box in the base has NeighborSet which is empty.  In this case, the
 *   Box from the base will not exist in the head BoxLevel.
 * # A Box in the base has a corresponding NeighborSet which in
 *   non-empty.  In this case, the NeighborSet contains the set of Boxes
 *   to which the Box in the base is related.
 */

class Connector:public tbox::DescribedClass
{
public:
   /*!
    * @brief NeighborsSet is a clarifying typedef.
    */
   typedef BoxSet NeighborSet;

   /*!
    * @brief Type of the iterator over neighborhoods.
    */
   typedef NeighborhoodSet::const_iterator ConstNeighborhoodIterator;

   /*!
    * @brief Type of the iterator over neighbors in a neighborhood.
    */
   typedef NeighborSet::const_iterator ConstNeighborIterator;

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
    * head Box comes from a base Box and there is an relationship
    * from the base Box to the head Box.
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
      const BoxLevel& base_mapped_box_level,
      const BoxLevel& head_mapped_box_level,
      const IntVector& base_width,
      const BoxLevel::ParallelState parallel_state = BoxLevel::DISTRIBUTED);

   /*!
    * @brief Destructor.
    */
   virtual ~Connector();

   /*!
    * @brief Clear the Connector, putting it into an uninitialized state.
    */
   void
   clear();

   /*!
    * @brief Clear the Connector's neighborhood relations.
    */
   void
   clearNeighborhoods();

   /*!
    * @brief Returns true if the object has been finalized
    */
   bool
   isFinalized() const;

   /*!
    * @brief Iterator pointing to the first neighborhood.
    */
   ConstNeighborhoodIterator
   begin() const;

   /*!
    * @brief Iterator pointing one past the last neighborhood.
    */
   ConstNeighborhoodIterator
   end() const;

   /*!
    * @brief Iterator pointing to the first neighbor in nbrhd.
    *
    * @param nbrhd The neighborhood whose neighbors are to be iterated.
    */
   ConstNeighborIterator
   begin(ConstNeighborhoodIterator& nbrhd) const;

   /*!
    * @brief Iterator pointing one past the last neighbor in nbrhd.
    *
    * @param nbrhd The neighborhood whose neighbors are to be iterated.
    */
   ConstNeighborIterator
   end(ConstNeighborhoodIterator& nbrhd) const;

   /*!
    * @brief Returns an Iterator pointing to the neighborhood of box_id--
    * localized version.
    *
    * @param[in] box_id
    */
   ConstNeighborhoodIterator
   findLocal(const BoxId& box_id) const;

   /*!
    * @brief Returns an Iterator pointing to the neighborhood of box_id--
    * globalized version.
    *
    * @param[in] box_id
    */
   ConstNeighborhoodIterator
   find(const BoxId& box_id) const;

   /*!
    * @brief Returns true if the local neighborhoods of this and other are the
    * same.
    *
    * @param[in] other
    */
   bool
   localNeighborhoodsEqual(const Connector& other) const;

   /*!
    * @brief Returns true if the neighborhood of the supplied BoxId of this
    * and other are the same.
    *
    * @param[in] box_id
    * @param[in] other
    */
   bool
   neighborhoodEqual(const BoxId& box_id, const Connector& other) const;

   /*!
    * @brief Return true if a neighbor set exists for the specified
    * BoxId.
    *
    * @param[in] mapped_box_id
    */
   bool
   hasNeighborSet(
      const BoxId& mapped_box_id) const;

   /*!
    * @brief Return true if the supplied box is in the neighborhood of the
    * supplied BoxId.
    *
    * @param[in] box_id
    * @param[in] neighbor
    */
   bool
   hasLocalNeighbor(
      const BoxId& box_id,
      const Box& neighbor) const;

   /*!
    * @brief Return the neighbor set for the specified BoxId.
    *
    * @param[in] mapped_box_id
    */
   void
   getNeighborBoxes(
      const BoxId& mapped_box_id,
      BoxList& nbr_boxes) const;

   /*!
    * @brief Return all neighbors for all neighborhoods.
    *
    * @param[out] neighbors
    */
   void
   getLocalNeighbors(
      NeighborSet& neighbors) const;

   /*!
    * @brief Return all neighbors for all neighborhoods.
    *
    * @param[out] neighbors
    */
   void
   getLocalNeighbors(
      BoxList& neighbors) const;

   /*!
    * @brief Return all neighbors for all neighborhoods segragated by BlockId.
    *
    * @param[out] neighbors
    */
   void
   getLocalNeighbors(
      std::map<BlockId, BoxList>& neighbors) const;

   /*!
    * @brief Returns the number of neighbors in the neighborhood with the
    * supplied BoxId.
    *
    * @param[in] box_id
    */
   int
   numLocalNeighbors(const BoxId& box_id) const;

   /*!
    * @brief Returns the number of empty neighborhoods in the Connector.
    *
    * @return The number of empty neighborhoods in the Connector.
    */
   int
   numLocalEmptyNeighborhoods() const;

   /*!
    * @brief Places the ranks of the processors owning neighbors into owners.
    *
    * @param[out] owners
    */
   void
   getLocalOwners(std::set<int>& owners) const;

   //@{
   /*!
    * @name Algorithms for changing individual Box's neighbor data
    */

   /*!
    * @brief Insert additional neighbors for the specified Box.
    *
    * @param[in] neighbors
    * @param[in] mapped_box_id
    */
   void
   insertNeighbors(
      const NeighborSet& neighbors,
      const BoxId& mapped_box_id);

   /*!
    * @brief Erase neighbor of the specified BoxId.
    *
    * @note Assertions
    * It is an error to to specify a non-existent BoxId.
    *
    * @param[in] neighbor
    * @param[in] mapped_box_id
    */
   void
   eraseNeighbor(
      const Box& neighbor,
      const BoxId& mapped_box_id);

   /*!
    * @brief Adds a neighbor of the specified BoxId.
    *
    * @param[in] neighbor
    * @param[in] box_id
    */
   void
   insertLocalNeighbor(
      const Box& neighbor,
      const BoxId& box_id);

   /*!
    * @brief Erases the neighborhood of the specified BoxId.
    *
    * @param[in] box_id
    */
   void
   eraseLocalNeighborhood(
      const BoxId& box_id);

   /*!
    * @brief Remove all the periodic relationships in the Connector.
    */
   void
   removePeriodicRelationships();

   /*!
    * @brief Remove all the periodic neighbors in all local neighborhoods.
    */
   void
   removePeriodicLocalNeighbors();

   /*!
    * @brief Check for any neighborhood roots which are periodic.
    *
    * @return true if any neighborhood root is periodic
    */
   bool
   hasPeriodicLocalNeighborhoodRoots() const;

   /*!
    * @brief Make and empty set of neighbors of the box with the supplied
    * box_id.
    *
    * @param[in] box_id
    */
   void
   makeEmptyLocalNeighborhood(
      const BoxId& box_id);

   /*!
    * @brief Remove empty sets of neighbors.
    */
   void
   eraseEmptyNeighborSets();

   /*!
    * @brief Returns true is the neighborhood of the supplied BoxId is empty.
    *
    * @param[in] box_id
    */
   bool
   isEmptyNeighborhood(
      const BoxId& box_id) const;

   /*!
    * @brief Coarsen all neighbors of this connector by ratio placing result
    * into coarser.
    *
    * @param[out] coarser
    * @param[in] ratio
    */
   void
   coarsenLocalNeighbors(
      Connector& coarser,
      const IntVector& ratio) const;

   /*!
    * @brief Refine all neighbors of this connector by ratio placing result
    * into finer.
    *
    * @param[out] finer
    * @param[in] ratio
    */
   void
   refineLocalNeighbors(
      Connector& finer,
      const IntVector& ratio) const;

   /*!
    * @brief Grow all neighbors of this connector by growth placing result
    * into grown.
    *
    * @param[out] grown
    * @param[in] growth
    */
   void
   growLocalNeighbors(
      Connector& grown,
      const IntVector& growth) const;

   //@}

   /*!
    * @brief Enforces implicit class invariants and removes non-local neighbors
    * from relationships.
    *
    * To be called after modifying a Connector's context through setBase,
    * setHead, or setWidth methods.
    */
   void
   finalizeContext();

   /*!
    * @brief Change the Connector base to new_base.  If finalize_context is
    * true then this is the last atomic change being made and finalizeContext
    * should be called.
    *
    * @param new_level
    * @param finalize_context
    */
   void
   setBase(
      const BoxLevel& new_base,
      bool finalize_context = false);

   /*!
    * @brief Return a reference to the base BoxLevel.
    */
   const BoxLevel&
   getBase() const;

   /*!
    * @brief Change the Connector head to new_head.  If finalize_context is
    * true then this is the last atomic change being made and finalizeContext
    * should be called.
    *
    * @param new_head
    * @param finalize_context
    */
   void
   setHead(
      const BoxLevel& new_head,
      bool finalize_context = false);

   /*!
    * @brief Return a reference to the head BoxLevel.
    */
   const BoxLevel&
   getHead() const;

   /*!
    * @brief Get the refinement ratio between the base and head
    * BoxLevels.
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
   bool
   ratioIsExact() const;

   /*!
    * @brief Return true if head BoxLevel is coarser than base
    * BoxLevel.
    */
   bool
   getHeadCoarserFlag() const;

   /*!
    * @brief Return true if the Connector contains only relationships to local
    * Boxes.
    *
    * The check only applies to neighbors of local base Boxes,
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
    * incident from a real (unshifted) Box.
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
   operator = (
      const Connector& rhs);

   //  TODO:  need to find out what the use case is for this, especially
   //  considering the caution statement.
   /*!
    * @brief Equality operator checks relationship data, Connector width and
    * equality of base and head Box pointers.
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
   operator == (
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
    * @brief Set the parallel distribution state.
    *
    * Before a Connector can be in a GLOBALIZED state, The base
    * BoxLevel given in initialize() must already be in
    * GLOBALIZED mode.  The base BoxLevel should remain in
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
      const BoxLevel::ParallelState parallel_state);

   /*!
    * @brief Return the current parallel state.
    */
   BoxLevel::ParallelState
   getParallelState() const;

   /*!
    * @brief Returns the MPI communication object, which is always
    * that of the base BoxLevel.
    */
   const tbox::SAMRAI_MPI&
   getMPI() const;

   /*!
    * @brief Change the Connector width to new_width.  If finalize_context is
    * true then this is the last atomic change being made and finalizeContext
    * should be called.
    *
    * @param new_width
    * @param finalize_context
    */
   void
   setWidth(
      const IntVector& new_width,
      bool finalize_context = false);

   /*!
    * @brief Return the Connector width associated with the relationships.
    *
    * For overlap Connectors, an relationship exists between a base and head
    * Boxes if the base mapped_box, grown by this width,
    * overlaps the head mapped_box.  For mapping Connectors, the width
    * the amount that a pre-map box must grow to nest the post-map
    * boxes.
    */
   const IntVector&
   getConnectorWidth() const;

   /*!
    * @brief Shrink the width of the connector modifying the proximity
    * relationships as needed.
    *
    * @param[in] new_width
    */
   void
   shrinkWidth(
      const IntVector& new_width);

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
    * @brief Return true if two Connector objects are
    * transposes of each other.
    *
    * Each Connector represents a set of directed relationships incident from its base
    * to its head.  The transpose represent the relationships in the opposite
    * direction.  In order for two Connector objects to be transpose of each
    * other, their Connector widths and base refinement ratios must be
    * such that an relationship in one set also appears in the other set.
    * A transpose set must have
    * @li base and head BoxLevels reversed from the untransposed set.
    * @li the same Connector width, although it is described in the index
    *     space of a different base BoxLevel.
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
    * @brief Check for consistency between the relationship data and base
    * mapped boxes, and return the number of consistency errors.
    *
    * Consistency stipulates that each neighbor list must correspond to
    * a base mapped box.
    *
    * relationship consistency errors should be treated as fatal because many
    * operations assume consistency.
    */
   size_t
   checkConsistencyWithBase() const;

   /*!
    * @brief Run checkConsistencyWithBase().
    *
    * If any inconsistency is
    * found, write out diagnostic information and throw an
    * unrecoverable assertion is found.
    */
   void
   assertConsistencyWithBase() const;

   /*!
    * @brief Check that the neighbors specified by the relationships exist in
    * the head BoxLevel.
    *
    * If the head is not GLOBALIZED, a temporary copy is made and
    * globalized for checking, triggering communication.
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
    * @brief Write the neighborhoods to a database.
    *
    * @param[in] database
    */
   void
   putNeighborhoodsToDatabase(tbox::Database& database) const;

   /*!
    * @brief Read the neighborhoods from a database.
    *
    * @param[in] database
    */
   void
   getNeighborhoodsFromDatabase(tbox::Database& database);

   /*!
    *
    * @brief Computes refinement ratio between head and base, whether that
    * ratio is exact and whether the head is coarser than the base.
    *
    * @param[in] baseRefinementRatio
    * @param[in] headRefinementRatio
    * @param[out] ratio
    * @param[out[ head_coarser
    * @param[out] ratio_exact
    */
   static void
   computeRatioInfo(
      const IntVector& baseRefinementRatio,
      const IntVector& headRefinementRatio,
      IntVector& ratio,
      bool& head_coarser,
      bool& ratio_is_exact);

   /*!
    * @brief Writes the neighborhoods to tbox::perr.
    *
    * @param[in] border Left border of the output.
    */
   void
   writeNeighborhoodsToErrorStream(
      const std::string& border) const;

   /*!
    * @brief Writes the requested neighborhood to tbox::perr.
    *
    * @param[in[ box_id
    * @param[in] border Left border of the output.
    */
   void
   writeNeighborhoodToErrorStream(
      const BoxId& box_id,
      const std::string& border) const;

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
   class Outputter
   {
      friend std::ostream&
      operator << (
         std::ostream& s,
         const Outputter& f);
private:
      friend class Connector;
      /*!
       * @brief Construct the Outputter with a Connector and the
       * parameters needed to output the Connector to a stream.
       */
      Outputter(
         const Connector& connector,
         const std::string& border,
         int detail_depth = 0);
      void
      operator = (
         const Outputter& r);               // Unimplemented private.
      const Connector& d_conn;
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
   Outputter
   format(
      const std::string& border = std::string(),
      int detail_depth = 0) const;

private:
   /*
    * Static integer constant descibing class's version number.
    */
   static const int HIER_CONNECTOR_VERSION;

   enum { BAD_INT = (1 << (8 * sizeof(int) - 2)) };

   /*!
    * @brief Return the globalized relationship data.
    *
    * @par Assertions
    * Throws an unrecoverable assertion if not in GLOBALIZED mode.
    */
   const NeighborhoodSet&
   getGlobalNeighborhoodSets() const;

   /*!
    * @brief Return the neighbor set for the specified BoxId.
    *
    * @param[in] mapped_box_id
    */
   const NeighborSet&
   getNeighborSet(
      const BoxId& mapped_box_id) const;

   /*!
    * @brief Create a copy of a DISTRIBUTED Connector and
    * change its state to GLOBALIZED.
    *
    * The returned object should be deleted to prevent memory leaks.
    */
   Connector *
   makeGlobalizedCopy(
      const Connector& other) const;

   /*!
    * @brief Get and store info on remote Boxes.
    *
    * This requires global communication (all gather).
    * Call acquireRemoteNeighborhoods_pack to pack up messages.
    * Do an all-gather.  Call acquireRemoteNeighborhoods_unpack
    * to unpack data from other processors.
    */
   void
   acquireRemoteNeighborhoods();

   //! @brief Pack local Boxes into an integer array.
   void
   acquireRemoteNeighborhoods_pack(
      std::vector<int>& send_mesg,
      int offset) const;

   //! @brief Unpack Boxes from an integer array into internal storage.
   void
   acquireRemoteNeighborhoods_unpack(
      const std::vector<int>& recv_mesg,
      const std::vector<int>& proc_offset);

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
    * @brief Handle for access to the base BoxLevel.
    *
    * We don't use a pointer to the BoxLevel, because it would
    * become dangling when the BoxLevel goes out of scope.
    */
   tbox::Pointer<BoxLevelHandle> d_base_handle;

   /*!
    * @brief Handle for access to the base BoxLevel.
    *
    * We don't use a pointer to the BoxLevel, because it would
    * become dangling when the BoxLevel goes out of scope.
    */
   tbox::Pointer<BoxLevelHandle> d_head_handle;

   /*!
    * @brief Connector width for the base BoxLevel.
    *
    * This is the amount of growth applied to a mapped_box in the base BoxLevel
    * before checking if the mapped_box overlaps a mapped_box in the head BoxLevel.
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
    * from the base and head BoxLevels.
    */
   IntVector d_ratio;

   /*!
    * @brief Whether the ratio between the base and head
    * BoxLevel refinement ratios are exactly as given by
    * d_ratio.  It can only be exact if it can be represented as an
    * IntVector.
    */
   bool d_ratio_is_exact;

   /*!
    * @brief Whether the base BoxLevel is at a finer index space.
    *
    * When this is true, d_ratio is the refinement ratio going
    * from the head to the base.
    *
    * This is redundant information.  You can compute it
    * from the base and head BoxLevels.
    */
   bool d_head_coarser;

   /*!
    * @brief Neighbor data for local Boxes.
    */
   NeighborhoodSet d_relationships;

   /*!
    * @brief Neighbor data for global Boxes in GLOBALIZED mode.
    */
   NeighborhoodSet d_global_relationships;

   /*!
    * @brief SAMRAI_MPI object.
    *
    * This is a copy of the getBase().getMPI().  We maintain a copy to
    * allow continued limited functionality should the base detaches
    * itself.
    */
   tbox::SAMRAI_MPI d_mpi;

   /*!
    * @brief State flag.
    *
    * Modified by setParallelState().
    */
   BoxLevel::ParallelState d_parallel_state;

   /*!
    * @brief true when Container's context has been finalized--base, head and
    * width are all defined and consistent/valid.
    */
   bool d_finalized;

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

   ConnectorType d_connector_type;

   static tbox::Pointer<tbox::Timer> t_acquire_remote_relationships;
   static tbox::Pointer<tbox::Timer> t_cache_global_reduced_data;

   static tbox::StartupShutdownManager::Handler
      s_initialize_finalize_handler;

};

}
}

#ifdef SAMRAI_INLINE
#include "SAMRAI/hier/Connector.I"
#endif

#endif // included_hier_Connector
