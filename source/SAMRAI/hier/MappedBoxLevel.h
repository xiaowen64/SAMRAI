/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Set of Boxes in the same "level". 
 *
 ************************************************************************/
#ifndef included_hier_MappedBoxLevel
#define included_hier_MappedBoxLevel

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/GridGeometry.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/MappedBoxLevelHandle.h"
#include "SAMRAI/hier/MappedBoxSet.h"
#include "SAMRAI/hier/PersistentOverlapConnectors.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/ConstPointer.h"
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/tbox/DescribedClass.h"
#include "SAMRAI/tbox/Timer.h"

#include <vector>

namespace SAMRAI {
namespace hier {

/*!
 * @brief A distributed set of Box objects which reside in the
 * same index space.
 *
   *
   * TODO: Are we eliminating DLBG terminology?
   *
 * This class is a part of the distributed layered box graph (DLBG) for
 * managing SAMR meshes in parallel.  A MappedBoxLevel is a set of 
 * boxes in the same index space. Relationships (e.g., neighbor adjacency)
 * among boxes is contained in a Connector object. Also, each MappedBoxLevel 
 * has an refinement ratio vector describing the relationship of the 
 * index space to that of a reference level in a patch hierarchy (typically
 * the coarsest level or level zero).  
 *
 * Like a PatchLevel, a MappedBoxLevel is a parallel object.  The
 * Boxes of a MappedBoxLevel may be distributed across all the
 * processors in an MPI communicator and can be in one
 * of two parallel states:
 *
 * - @b DISTRIBUTED: Each MPI process knows only the Boxes in the set
 * that are "owned" by that process.  This is analogous to a PatchLevel 
 * which owns only the Patches that reside on a process.
 *
 * - @b GLOBALIZED: All processes know all Boxes in the set.
 * This is analogous to PatchLevel BoxList state when it is globalized
 * (@see PatchLevel::getBoxes()).
 *
 * @par Performance notes
 * <li> The parallel state is changed by calling setParallelState().  Going
 * from DISTRIBUTED to GLOBALIZED state is an expensive operation requiring
 * all-to-all communication.  Using this state can incur a significant 
 * performance penalty.
 *
 * <li> The GLOBALIZED state requires more memory.
 *
 * <li> Transitioning from GLOBALIZED state to DISTRIBUTED state is
 * cheap.
 *
 * @note
 * The general attributes of a MappedBoxLevel are
 * <li> the set of Box objects with unique BoxIds,
 * <li> the refinement ratio defining their index space, and
 * <li> the parallel state.
 *
 * Box object uniqueness is based on the Box equality operator,
 * which compares owner MPI ranks and local indices.  Therefore, 
 * a valid MappedBoxLevel does not contain two Boxes with the same
 * owner and index.
 */
class MappedBoxLevel:public tbox::DescribedClass
{

public:
   /*!
    * @brief Names of parallel states.
    */
   /*
    * TODO: We should only have one enum for this. Currently, we have
    * this one and a similar one in Connector.h.
    */
   enum ParallelState { DISTRIBUTED, GLOBALIZED };

   /*!
    * @brief Construct uninitialized object.
    *
    * Uninitialized objects can be initialized by calling initialize()
    * or swapInitialize().
    *
    * @see initialize()
    * @see swapInitialize()
    *
    * @param[in] dim
    */
   explicit MappedBoxLevel(
      const tbox::Dimension& dim);

   /*!
    * @brief Copy constructor.
    *
    * New object has the same parallel state as original.
    *
    * Persistent Connectors are not duplicated.  This
    * decision was based on expected usage, which is that
    * copies are either for short term usage or meant to
    * be changed in some way and will invalidate Connectors.
    *
    * @param[in] rhs
    */
   MappedBoxLevel(
      const MappedBoxLevel& rhs);

   /*!
    * @brief Construct initialized and populated object.
    *
    * @see initialize()
    *
    * @param[in] mapped_boxes
    * @param[in] ratio
    * @param[in] grid_geom
    * @param[in] mpi
    * @param[in] parallel_state
    */
   explicit MappedBoxLevel(
      const MappedBoxSet& mapped_boxes,
      const IntVector& ratio,
      const tbox::ConstPointer<GridGeometry> &grid_geom,
      const tbox::SAMRAI_MPI& mpi = tbox::SAMRAI_MPI::getSAMRAIWorld(),
      const ParallelState parallel_state = DISTRIBUTED);

   /*!
    * @brief Constructs an empty, initialized object.
    *
    * @see addMappedBox()
    * @see addBox()
    * @see initialize()
    *
    * @param[in] ratio
    * @param[in] grid_geom
    * @param[in] mpi
    * @param[in] parallel_state
    */
   explicit MappedBoxLevel(
      const IntVector& ratio,
      const tbox::ConstPointer<GridGeometry> &grid_geom,
      const tbox::SAMRAI_MPI& mpi = tbox::SAMRAI_MPI::getSAMRAIWorld(),
      const ParallelState parallel_state = DISTRIBUTED);

   /*!
    * @brief Destructor.
    *
    * Deallocate internal data.
    */
   ~MappedBoxLevel(
      void);

   //@{
   //! @name Initialization and clearing methods

   /*!
    * @brief Initialize the MappedBoxLevel
    *
    * If @c parallel_state is GLOBALIZED, mapped_boxes must be
    * the global set of Boxes.
    *
    * If @c parallel_state is DISTRIBUTED, mapped_boxes should contain
    * all local Boxes.  Non-local Boxes are ignored.
    *
    * Once the object is initialized, you can further modify it by
    * adding and removing boxes.
    *
    * Initializing is a modifying operation, causing the
    * PersistentOverlapConnectors to be cleared.
    *
    * The content and state of the object before calling this function
    * is discarded.
    *
    * @see initializePrivate()
    * @see getPersistentOverlapConnectors().
    *
    * @param[in] mapped_boxes
    * @param[in] ratio
    * @param[in] grid_geom
    * @param[in] mpi
    * @param[in] parallel_state
    */
   void
   initialize(
      const MappedBoxSet& mapped_boxes,
      const IntVector& ratio,
      const tbox::ConstPointer<GridGeometry> &grid_geom,
      const tbox::SAMRAI_MPI& mpi = tbox::SAMRAI_MPI::getSAMRAIWorld(),
      const ParallelState parallel_state = DISTRIBUTED);

  /*!
    * @brief Initialize the MappedBoxLevel without and Boxes
    *
    * The content and state of the object before calling this function
    * is discarded.
    *
    * @see addMappedBox()
    * @see addBox()
    * @see initialize(const MappedBoxSet&, const IntVector&, const tbox::SAMRAI_MPI&, const ParallelState)
    *
    * @param[in] ratio
    * @param[in] grid_geom
    * @param[in] mpi
    * @param[in] parallel_state
    */
   void
   initialize(
      const IntVector& ratio,
      const tbox::ConstPointer<GridGeometry> &grid_geom,
      const tbox::SAMRAI_MPI& mpi = tbox::SAMRAI_MPI::getSAMRAIWorld(),
      const ParallelState parallel_state = DISTRIBUTED);

   /*!
    * @brief Initialize the MappedBoxLevel.
    *
    * Similar to initialize(const MappedBoxSet&, const IntVector&, const tbox::SAMRAI_MPI&, const ParallelState), except that the @c mapped_boxes are mutable.
    *
    * The state of the object before calling this function is
    * discarded.  The Box content before calling this function
    * is returned via the @c mapped_boxes argument.
    *
    * @see initializePrivate()
    *
    * @param[in,out] mapped_boxes On input, this should contain the
    * Boxes to place in the MappedBoxLevel.  On output, it
    * contains the Boxes that were in the MappedBoxLevel before 
    * the call.
    *
    * @param[in] mapped_boxes
    * @param[in] ratio
    * @param[in] grid_geom
    * @param[in] mpi
    * @param[in] parallel_state
    */
   void
   swapInitialize(
      MappedBoxSet& mapped_boxes,
      const IntVector& ratio,
      const tbox::ConstPointer<GridGeometry> &grid_geom,
      const tbox::SAMRAI_MPI& mpi = tbox::SAMRAI_MPI::getSAMRAIWorld(),
      const ParallelState parallel_state = DISTRIBUTED);

   /*!
    * @brief Returns True if the object has been initialized.
    */
   bool
   isInitialized() const;

   /*!
    * @brief Clear the internal state of the MappedBoxLevel.
    *
    * The MappedBoxLevel will be in an uninitialized state
    * after a call to this method.
    */
   void
   clear();

   /*!
    * @brief Clear the globalized version and the persistent overlap
    * connectors for data consistency.
    *
    * Most of the time, this method is automatically called by methods
    * that know when some data is stale and needs to be cleared.  
    * For example, adding a box makes the global number of
    * boxes stale.  However, sometimes it is necessary to call this
    * method manually.  For example, when only some processes add
    * boxes while others do not, resulting in some processes not
    * knowing that the global number of boxes is inconsistent.
    *
    * @param[in] isInvalid A flag indicating that boxes have been (or will
    *            be) removed, thus invalidating the handle.
    */
   void clearForBoxChanges( bool isInvalid = true );


   //@{

   //! @name Parallelism

   /*!
    * @brief Set the parallel state.
    *
    * This method is potentially expensive.
    * Acquiring remote Box information (when going
    * to GLOBALIZED mode) triggers all-gather communication.
    * More memory is required to store additional Boxes.
    *
    * Data not used by the new state gets deallocated.
    *
    * @param[in] parallel_state
    */
   void
   setParallelState(
      const ParallelState parallel_state);

   /*!
    * @brief Returns the ParallelState of the object.
    */
   ParallelState
   getParallelState() const;

   /*!
    * @brief If global reduced data (global Box count, global
    * cell count and global bounding box) have not been updated,
    * compute and cache them (communication required).
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
    * @brief Return the globalized version of the MappedBoxLevel,
    * creating it if needed.
    *
    * If the MappedBoxLevel is in globalized state, return @c *this.
    * If not, create and cache a globalized version (if necessary) and
    * return that.
    *
    * The cached version remains until it is removed by
    * deallocateGlobalizedVersion() or a method that can potentially
    * change the Boxes is called.  Note that globalizing and
    * globalized data is not scalable.  Use only when necessary.
    *
    * Obviously, when the globalized version must be created (when the
    * MappedBoxLevel is in DISTRIBUTED state and there is no cached
    * version yet), all processes must make this call at the same
    * point.
    */
   const MappedBoxLevel&
   getGlobalizedVersion() const;

   /*!
    * @brief Deallocate the internal globalized version of the
    * MappedBoxLevel, if there is any.
    */
   void
   deallocateGlobalizedVersion() const;

   /*!
    * @brief Returns the SAMRAI_MPI communicator over which the Boxes
    * are distributed.
    */
   const tbox::SAMRAI_MPI&
   getMPI() const;

   /*
    * TODO: Why are these methods here? Shouldn't we get this information
    * from the SAMRAI_MPI object (method above).
    */
   /*!
    * @brief Return the processor rank for the internal MPI communicator.
    */
   int
   getRank() const;

   /*!
    * @brief Return the number of processes for the internal MPI communicator.
    */
   int
   getNproc() const;

   //@}


   /*!
    * @brief Assignment operator duplicates all internal data,
    * including parallel mode.
    *
    * Assignment is a modifying operation, causing the
    * PersistentOverlapConnectors to be cleared.
    *
    * Persistent Connectors are not duplicated.  This
    * decision was based on expected usage, which is 
    * that copies are either for short term usage or meant to
    * be changed in some way (thus invalidating current
    * Connectors anyway).
    *
    * @see getPersistentOverlapConnectors()
    *
    * @param[in] rhs
    */
   MappedBoxLevel&
   operator = (
      const MappedBoxLevel& rhs);

   /*!
    * @brief Swap the contents of two MappedBoxLevel objects.
    *
    * Swapping is a modifying operation, so the
    * PersistentOverlapConnectorss of the operands are cleared.
    *
    * Persistent Connectors are not swapped.  This decision
    * was based on expected usage, which is that
    * copies are either for short term usage or meant to be
    * changed in some way (thus invalidating current
    * Connectors anyway).
    *
    * @param[in,out] level_a
    * @param[in,out] level_b
    */
   static void
   swap(
      MappedBoxLevel& level_a,
      MappedBoxLevel& level_b);

   //@}

   /*!
    * @brief Equality comparison.
    *
    * All data required to initialize the object is compared, except
    * for the parallel state.  Thus equality here means just the local
    * parts are equal.  @b BEWARE!  This means that one processor may
    * see the equality differently from another.
    *
    * The cost for the comparison is on the order of the local
    * Box count.  An object may be compared to itself, an
    * efficient operation that always returns true.
    *
    * @param[in] rhs
    */
   bool
   operator == (
      const MappedBoxLevel& rhs) const;

   /*!
    * @brief Inequality comparison.
    *
    * All data required to initialize the object is compared, except
    * for the parallel state.  Thus equality here means just the local
    * parts are equal.  @b BEWARE!  This means that one processor may
    * see the inequality differently from another.
    *
    * The cost for the comparison is on the order of the local
    * Box count.  However, an object may be compared to itself,
    * an efficient operation that always returns false.
    *
    * @param[in] rhs
    */
   bool
   operator != (
      const MappedBoxLevel& rhs) const;


   //@{
   /*!
    * @name Accessors
    */

   /*!
    * @brief Returns the container of local Boxes.
    *
    * @par Important
    * The MappedBoxSet returned contains periodic image
    * Boxes (if any).  To iterate through real Boxes only, see
    * RealMappedBoxConstIterator.
    *
    * You cannot directly modify the MappedBoxSet because it may
    * invalidate other internal data.  Use other methods for modifying
    * the MappedBoxSet.
    *
    * @see getGlobalNumberOfBoxes()
    * @see getLocalNumberOfBoxes()
    *
    */
   const MappedBoxSet&
   getMappedBoxes() const;

   /*!
    * @brief Returns the container of global Boxes.
    *
    * @par Assertions
    * Throws an unrecoverable assertion if not in GLOBALIZED mode.
    */
   const MappedBoxSet&
   getGlobalMappedBoxes() const;

   /*!
    * @brief Returns the container of global Boxes.
    *
    * @par Assertions
    * Throws an unrecoverable assertion if not in GLOBALIZED mode.
    */
   void
   getGlobalBoxes(BoxList &global_boxes) const;

   /*
    * TODO: Why are the following two methods here?  Returning local id 
    * information like this is dangerous in that it seems to imply that 
    * they can be used as an integer range, or a count (last - first + 1).  
    * Since the first method is not used and the second is used in a few 
    * places, wouldn't it be better to just use the previous method?
    */
   /*!
    * @brief Returns the first LocalId, or one with a value of -1 if
    * no local Box exists.
    */
   LocalId
   getFirstLocalId() const;

   /*!
    * @brief Returns the last LocalId, or one with a value of -1 if no
    * local Box exists.
    */
   LocalId
   getLastLocalId() const;

   /*!
    * @brief Get const access to MappedBoxLevel's refinement ratio
    * (with respect to a reference level).
    */
   const IntVector&
   getRefinementRatio() const;

   /*!
    * @brief Return local number of boxes.
    *
    * Periodic image Boxes are excluded.
    */
   size_t
   getLocalNumberOfBoxes() const;

   /*!
    * @brief Return number of boxes local to the given rank.
    *
    * Periodic image Boxes are excluded.
    *
    * Object must be in GLOBALIZED mode to use this method.
    *
    * @param[in] rank
    */
   size_t
   getLocalNumberOfBoxes(
      int rank) const;

   /*!
    * @brief Return global number of Boxes.
    *
    * This requires a global reduction, if the global-reduced data has
    * not been computed and cached.  When communication is required,
    * all processors must call this method.  To ensure that no
    * communication is needed, call cacheGlobalReducedData() first.
    *
    * Periodic image Boxes are excluded.
    */
   int
   getGlobalNumberOfBoxes() const;

   /*!
    * @brief Return maximum number of Boxes over all processes.
    *
    * This requires a global reduction, if the global-reduced data has
    * not been computed and cached.  When communication is required,
    * all processors must call this method.  To ensure that no
    * communication is needed, call cacheGlobalReducedData() first.
    *
    * Periodic image Boxes are excluded.
    */
   int
   getMaxNumberOfBoxes() const;

   /*!
    * @brief Return maximum number of Boxes over all processes.
    *
    * This requires a global reduction, if the global-reduced data has
    * not been computed and cached.  When communication is required,
    * all processors must call this method.  To ensure that no
    * communication is needed, call cacheGlobalReducedData() first.
    *
    * Periodic image Boxes are excluded.
    */
   int
   getMinNumberOfBoxes() const;

   /*!
    * @brief Return local number of cells.
    *
    * Cells in periodic image Boxes are excluded.
    */
   size_t
   getLocalNumberOfCells() const;

   /*!
    * @brief Return maximum number of cells over all processes.
    *
    * This requires a global reduction, if the global-reduced data has
    * not been computed and cached.  When communication is required,
    * all processors must call this method.  To ensure that no
    * communication is needed, call cacheGlobalReducedData() first.
    *
    * Periodic image Boxes are excluded.
    */
   int
   getMaxNumberOfCells() const;

   /*!
    * @brief Return maximum number of cells over all processes.
    *
    * This requires a global reduction, if the global-reduced data has
    * not been computed and cached.  When communication is required,
    * all processors must call this method.  To ensure that no
    * communication is needed, call cacheGlobalReducedData() first.
    *
    * Periodic image Boxes are excluded.
    */
   int
   getMinNumberOfCells() const;

   /*!
    * @brief Return number of cells local to the given rank.
    *
    * Cells in periodic image Boxes are excluded.
    *
    * Object must be in GLOBALIZED mode to use this method.
    *
    * @param[in] rank
    */
   size_t
   getLocalNumberOfCells(
      int rank) const;

   /*!
    * @brief Return global number of cells.
    *
    * This requires a global reduction if the global-reduced data has
    * not been computed and cached.  When communication is required,
    * all processors must call this method.  To ensure that no
    * communication is needed, call cacheGlobalReducedData() first.
    *
    * Cells in periodic image Boxes are excluded.
    */
   int
   getGlobalNumberOfCells() const;

   /*!
    * @brief Return bounding box for local Boxes in a block.
    */
   const Box&
   getLocalBoundingBox( int block_number ) const;

   /*!
    * @brief Return bounding box for global Boxes in a block.
    *
    * This requires a global reduction if the global bounding box has
    * not been computed and cached.  When communication is required,
    * all processors must call this method.  To ensure that no
    * communication is needed, call cacheGlobalReducedData() first.
    */
   const Box&
   getGlobalBoundingBox( int block_number ) const;

   /*!
    * @brief Return size of the largest local Box in a block.
    */
   const IntVector&
   getLocalMaxBoxSize( int block_number ) const;

   /*!
    * @brief Return size of the smallest local Box in a block.
    */
   const IntVector&
   getLocalMinBoxSize( int block_number ) const;

   /*!
    * @brief Return size of the largest Box globally in a block.
    *
    * This requires a global reduction if the global bounding box has
    * not been computed and cached.  When communication is required,
    * all processors must call this method.  To ensure that no
    * communication is needed, call cacheGlobalReducedData() first.
    */
   const IntVector&
   getGlobalMaxBoxSize( int block_number ) const;

   /*!
    * @brief Return size of the smallest Box globally in a block.
    *
    * This requires a global reduction if the global bounding box has
    * not been computed and cached.  When communication is required,
    * all processors must call this method.  To ensure that no
    * communication is needed, call cacheGlobalReducedData() first.
    */
   const IntVector&
   getGlobalMinBoxSize( int block_number ) const;

   /*!
    * @brief Return the dimension of this object.
    *
    * If object has never been initialized, return
    * tbox::Dimension::getInvalidDimension().
    */
   const tbox::Dimension&
   getDim() const;

   /*!
    * @brief Return the GridGeometry associated with this object.
    *
    * If object has never been initialized, return NULL pointer.
    */
   const tbox::ConstPointer<GridGeometry>&
   getGridGeometry() const;

   //@}


   //@{

   //! @name Individual Box methods.

   /*
    * TODO: Why the vacant index thing?  The comments say that the
    * box will be assigned an unused index.  Isn't this a "vacant" index?
    *
    * TODO: Why does the first method require "distributed" state and the ones
    * after require GLOBALIZED state for a remote box?  These comments
    * are inconsistent and confusing.  I think it would be best to have the
    * same pre/post conditions apply to all similar methods (e.g., "add box" 
    * methods), then describe them once rather than repeat (potentially 
    * inconsistently for each method).
    */


   /*!
    * @brief Create new local Box from given Box and add it to this
    * level.
    *
    * The new Box will be assigned an unused local index.  To be
    * efficient, no communication will be used.  Therefore, the state
    * must be distributed.
    *
    * The new Box will have a periodic shift number
    * corresponding to zero-shift.
    *
    * It is faster not to request a vacant index when adding a box.
    *
    * @param[in] box
    * @param[in] block_id 
    * @param[in] use_vacant_index
    *
    * @return iterator to the new Box
    */
   MappedBoxSet::iterator
   addBox(
      const Box& box,
      const BlockId& block_id,
      const bool use_vacant_index = true);

   /*!
    * @brief Add a Box to this level.
    *
    * Adding a remote Box is allowed if the object is in
    * GLOBALIZED mode.
    *
    * @par CAUTION
    * To be efficient, no checks are made to make sure the
    * MappedBoxLevel representation is consistent across all
    * processors.  Setting inconsistent data leads potentially
    * elusive bugs.
    *
    * @par Errors
    * It is an error to add a periodic image of a Box that is
    * not a part of the MappedBoxLevel.
    *
    * It is an error to add any Box that already exists.
    *
    * FIXME: Should we prevent this operation if persistent overlap
    * Connectors are attached to this object?
    *
    * @param[in] mapped_box
    */
   void
   addMappedBox(
      const Box& mapped_box);

   /*!
    * @brief Insert given periodic image of an existing Box.
    *
    * Adding a remote Box is allowed if the object is in
    * GLOBALIZED mode.
    *
    * Unlike adding a regular Box, it is OK to add a periodic
    * image Box that already exists.  However, that is a no-op.
    *
    * @par CAUTION
    * To be efficient, no checks are made to make sure the
    * MappedBoxLevel representation is consistent across all
    * processors.  Setting inconsistent data leads potentially elusive
    * bugs.
    *
    * @par Errors
    * It is an error to add a periodic image of a Box that does
    * not exist.
    *
    * FIXME: Should we prevent this operation if persistent overlap
    * Connectors are attached to this object?
    *
    * @param[in] existing_mapped_box  An existing Box for reference.
    *      This Box must be in the MappedBoxLevel.  The Box added
    *      is an image of the reference Box but shifted to another
    *      position.
    * @param[in] shift_number The valid shift number for the Box being
    *      added.  The shift amount is taken from the PeriodicShiftCatalog.
    */
   void
   addPeriodicMappedBox(
      const Box& existing_mapped_box,
      const PeriodicId& shift_number);

   /*!
    * @brief Erase the existing Box specified by its iterator.
    *
    * The given iterator @em MUST be a valid iterator pointing to a
    * Box currently in this object.  After erasing, the iterator
    * is advanced to the next valid Box (or the end of its
    * MappedBoxSet).
    *
    * Erasing a Box also erases all of its periodic images.
    *
    * FIXME: Should we prevent this operation if the object has
    * persistent overlap Connectors?
    *
    * @param[in] ibox The iterator of the Box to erase.
    */
   void
   eraseMappedBox(
      MappedBoxSet::iterator& ibox);

   /*!
    * @brief Erase the Box matching the one given.
    *
    * The given Box @em MUST match a Box currently in this
    * object.  Matching means that the BoxId's match
    * (disregarding the Boxes).
    *
    * Erasing a Box also erases all of its periodic images.
    *
    * FIXME: Should we prevent this operation if the object has
    * persistent overlap Connectors?
    *
    * @param[in] mapped_box
    */
   void
   eraseMappedBox(
      const Box& mapped_box);

   /*!
    * @brief Find the Box matching the one given.
    *
    * Only the BoxId matters in matching, so the actual Box can
    * be anything.
    *
    * If @c mapped_box is not a local Box, the state must be
    * GLOBALIZED.
    *
    * @param[in] mapped_box
    *
    * @return Iterator to the mapped_box, or @c
    * getMappedBoxes(owner).end() if mapped_box does not exist in set.
    */
   MappedBoxSet::const_iterator
   getMappedBox(
      const Box& mapped_box) const;

   /*!
    * @brief Find the Box specified by the given BoxId and
    * periodic shift.
    *
    * If @c mapped_box is not a local Box, the state must be
    * GLOBALIZED.
    * @param[in] mapped_box_id
    *
    * @return Iterator to the mapped_box, or @c
    * getMappedBoxes(owner).end() if mapped_box does not exist in set.
    */
   MappedBoxSet::const_iterator
   getMappedBox(
      const BoxId& mapped_box_id) const;

   /*
    * TODO: What is different about these "strict" methods compared to
    * the preceding ones?  I can't tell from the comments.
    */

   /*!
    * @brief Find the Box matching the one given.
    *
    * Only the BoxId matters in matching, so the actual Box can
    * be anything.
    *
    * If @c mapped_box is not owned by the local process, the state
    * must be GLOBALIZED.
    *
    * You cannot directly modify the MappedBoxSet because it may
    * invalidate other internal data.  Use other methods for modifying
    * the MappedBoxSet.
    *
    * @par Assertions
    * Throws an unrecoverable assertion if the Box does not
    * exist.
    *
    * @param[in] mapped_box
    *
    * @return Iterator to the mapped_box.
    */
   MappedBoxSet::const_iterator
   getMappedBoxStrict(
      const Box& mapped_box) const;

   /*!
    * @brief Find the Box specified by the given BoxId and
    * periodic shift.
    *
    * You cannot directly modify the MappedBoxSet because it may
    * invalidate other internal data.  Use other methods for modifying
    * the MappedBoxSet.
    *
    * @par Assertions
    * Throw an unrecoverable assertion if the Box does not exist.
    *
    * @param[in] mapped_box_id
    *
    * @return Iterator to the mapped_box.
    */
   MappedBoxSet::const_iterator
   getMappedBoxStrict(
      const BoxId& mapped_box_id) const;

   /*!
    * @brief Returns true when the object has a Box specified by the
    * BoxId.
    *
    * @param[in] mapped_box_id
    */
   bool
   hasMappedBox(
      const BoxId& mapped_box_id) const;

   /*!
    * @brief Returns true when the object has a Box consistent with all
    * of the arguments
    *
    * @param[in] global_id
    * @param[in] block_id
    * @param[in] periodic_id
    */
   bool
   hasMappedBox(
      const GlobalId& global_id,
      const BlockId& block_id,
      const PeriodicId &periodic_id) const;

   /*!
    * @brief Returns true when the object has a Box matching the
    * one parameter.
    *
    * @param[in] mapped_box
    */
   bool
   hasMappedBox(
      const Box& mapped_box) const;

   //@}

   //@{
   /*!
    * @name IO support.
    */

   /*!
    * @brief Write the MappedBoxLevel to a database.
    *
    * Write only local parts regardless of parallel state (to avoid
    * writing tons of repetitive data).
    *
    * @par Assertions
    * Check that database is a non-null Pointer.
    *
    * @param[in,out] database
    */
   void
   putToDatabase(
      tbox::Database& database) const;

   /*!
    * @brief Read the MappedBoxLevel from a database.
    *
    * Put the MappedBoxLevel in the DISTRIBUTED parallel state and
    * read only local parts.
    *
    * If the MappedBoxLevel is initialized, use its SAMRAI_MPI object
    * and require its refinement ratio to match that in the database.
    * If the MappedBoxLevel is uninitialized, it will be initialized
    * to use tbox::SAMRAI_MPI::getSAMRAIWorld() for the SAMRAI_MPI
    * object.  Note that these behaviors have not been extensively
    * discussed by the SAMRAI developers and may be subject to change.
    *
    * @par Assertions
    * Check that database is a non-null Pointer.
    *
    * @param[in,out] database
    */
   void
   getFromDatabase(
      tbox::Database& database,
      const tbox::ConstPointer<GridGeometry> &grid_geom);

   //@}


   /*!
    * @brief Get the collection of overlap Connectors dedicated to
    * provide overlap neighbors for this MappedBoxLevel.
    *
    * The PersistentOverlapConnectors provides overlap neighbors for
    * this MappedBoxLevel.  Its role is to create and manage
    * persistent overlap Connectors based at this MappedBoxLevel and
    * persisting until the MappedBoxLevel changes (so they should not
    * be set up until the MappedBoxLevel is in its final state).  This
    * is the mechanism by which code that can efficiently generate the
    * overlap Connectors (usually the code that generated the
    * MappedBoxLevel) provides overlap data to code using the
    * MappedBoxLevel.  The PersistentOverlapConnectors are guaranteed
    * to be correct, so any changes to the MappedBoxLevel will cause
    * current Connectors to be deallocated.
    *
    * @see PersistentOverlapConnectors for instructions on creating
    * the Connectors.
    */
   PersistentOverlapConnectors&
   getPersistentOverlapConnectors() const;

   /*
    * TODO: The following method is "not for general use" and indeed
    * is only used in two Connector classes.  Would it be better to 
    * make the method private and make this class a friend of those?
    */

   /*!
    * @brief Get the handle with which Connectors
    * reference the MappedBoxLevel instead of referencing the
    * MappedBoxLevel itself.  Not for general use.
    *
    * Connectors referencing their base and head MappedBoxLevels should
    * reference their handles instead of the MappedBoxLevels themselves.
    * As long as the MappedBoxLevel does not change in a way
    * that can invalidate Connector data, you can access
    * the MappedBoxLevel from the MappedBoxLevelHandle.
    *
    * If the MappedBoxLevel go out of scope before the
    * Connector disconnects, this tbox::Pointer object will
    * stay around until all Connectors have disconnected.
    *
    * Operations that can invalidate Connector data are those
    * that remove information from the MappedBoxLevel.  These
    * are:
    *
    * @li initialize()
    * @li swapInitialize()
    * @li swap()
    * @li clear()
    * @li operator=() (assignment) (Exception: assigning to
    *     self is a no-op, which does not invalidate Connector
    *     data.
    * @li eraseMappedBox() (Note that adding a Box
    *     does not invalidate Connector data.)
    * @li going out of scope
    *
    * @see MappedBoxLevelHandle.
    *
    * @return A tbox::Pointer to the MappedBoxLevelHandle
    */
   const tbox::Pointer<MappedBoxLevelHandle>&
   getMappedBoxLevelHandle() const;

   //@{

   /*!
    * @name Methods for outputs, error checking and debugging.
    */

   /*!
    * @brief Print Box info from this level
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
    * @brief Print out statistics on the Boxes.
    *
    * @param[in,out] os The output stream
    * @param[in] border
    */
   void
   printMappedBoxStats(
      std::ostream& os,
      const std::string& border) const;

   /*!
    * @brief A class for outputting MappedBoxLevel.
    *
    * To use, see MappedBoxLevel::format().
    *
    * This class simplifies the insertion of a MappedBoxLevel into a
    * stream while letting the user control how the MappedBoxLevel is
    * formatted for output.
    *
    * Each Outputter is a light-weight object constructed with a
    * MappedBoxLevel and output parameters.  The Outputter is capable
    * of outputting its MappedBoxLevel, formatted according to the
    * parameters.
    */
   class Outputter {
      friend std::ostream& operator << ( std::ostream& s, const Outputter& f);
   private:
      friend class MappedBoxLevel;
      /*!
       * @brief Construct the Outputter with a MappedBoxLevel and the
       * parameters needed to output the MappedBoxLevel to a stream.
       *
       * @param[in] mapped_box_level
       * @param[in] border
       * @param[in] detail_depth
       */
      Outputter( const MappedBoxLevel &mapped_box_level,
                 const std::string& border,
                 int detail_depth = 0);
      void operator=( const Outputter &r ); // Unimplemented private.
      const MappedBoxLevel &d_level;
      const std::string d_border;
      const int d_detail_depth;
   };

   /*!
    * @brief Return a object that can format the MappedBoxLevel for
    * inserting into output streams.
    *
    * Usage example:
    * @code
    *    std::cout << "my mapped_box_level:\n"
    *         << mapped_box_level.format("  ", 2) << std::endl;
    * @endcode
    *
    * @param[in] border
    * @param[in] detail_depth
    */
   Outputter format( const std::string& border=std::string(),
                     int detail_depth = 0 ) const;

   //@}

   /*!
    * @brief Allows std::vector to allocate objects with
    * uninitialized dimensions.
    */
   friend class std::vector<MappedBoxLevel>;

   /*!
    * @brief Set up things for the entire class.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   initializeCallback();

   /*!
    * @brief Free static timers.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   finalizeCallback();

private:
   /*
    * Static integer constant describing class's version number.
    */
   static const int HIER_MAPPED_BOX_LEVEL_VERSION;

   /*
    * Static integer constant describing the number of statistics in this
    * class.
    */
   static const int MAPPED_BOX_LEVEL_NUMBER_OF_STATS;

   /*
    * TODO: This same enum is defined in the Connector header.
    * If there is a common use for this, should it be defined in a
    * common location? Also, it seems to be used similarly to the
    * BAD_INTEGER #define in BergerRigoutsosNode.C.  Is the intent
    * really the same?
    */
   enum { BAD_INT = (1 << (8 * sizeof(int) - 2)) };

   /*!
    * @brief Construct uninitialized object.
    *
    * Constructor creates an uninitialized object in distributed state.
    *
    * Private to limit where an uninitialized object can
    * be created.
    */
   MappedBoxLevel();

   /*
    * TODO: The comments for the following method use the phrase
    * "local redundant data" three times, but I still don't know
    * what that is!
    */
   /*!
    * @brief Recompute local redundant data.
    *
    * Local redundant data is usually updated immediately after their
    * dependencies change.  On certain occasions, we recompute all
    * the local redundant data.
    */
   void
   computeLocalRedundantData();

   //@{

   /*!
    * @brief Get and store info on remote Boxes.
    *
    * This requires global communication (all-gather).
    * Call acquireRemoteMappedBoxes_pack to pack up messages.
    * Do an all-gather.  Call acquireRemoteMappedBoxes_unpack
    * to unpack data from other processors.
    */
   void
   acquireRemoteMappedBoxes();

   //! @brief Pack local Boxes into an integer array.
   void
   acquireRemoteMappedBoxes_pack(
      std::vector<int>& send_mesg) const;

   /*!
    * @brief Unpack Boxes from an integer array into internal
    * storage.
    */
   void
   acquireRemoteMappedBoxes_unpack(
      const std::vector<int>& recv_mesg,
      std::vector<int>& proc_offset);

   /*!
    * @brief Get and store info on remote Boxes for multiple
    * MappedBoxLevel objects.
    *
    * This method combines communication for the multiple
    * mapped_box_levels to increase message passing efficiency.
    *
    * Note: This method is stateless (could be static).
    */
   void acquireRemoteMappedBoxes(
      const int num_sets,
      MappedBoxLevel * multiple_mapped_box_level[]);
   //@}

   /*!
    * @brief Deallocate persistent overlap Connectors, if there are any.
    */
   void
   clearPersistentOverlapConnectors();

   /*!
    * @brief Detach this object from the handle it has been using.
    *
    * Postcondition: Objects that cached the handle would no longer
    * be able to access this MappedBoxLevel by the handle.
    */
   void
   detachMyHandle() const;

   /*!
    * @brief Encapsulates functionality common to all initialization
    * functions.
    */
   void
   initializePrivate(
      const IntVector& ratio,
      const tbox::ConstPointer<GridGeometry> &grid_geom,
      const tbox::SAMRAI_MPI& mpi = tbox::SAMRAI_MPI::getSAMRAIWorld(),
      const ParallelState parallel_state = DISTRIBUTED);

   /*!
    * @brief MappedBoxLevel is a parallel object,
    * and this describes its MPI object.
    */
   tbox::SAMRAI_MPI d_mpi;

   /*!
    * @brief Locally-stored Boxes.
    *
    * This is always the container of local Boxes, regardless of
    * parallel mode.
    */
   MappedBoxSet d_mapped_boxes;

   /*!
    * @brief Locally-stored global Boxes (for GLOBALIZED mode).
    *
    * In DISTRIBUTED mode, this is empty.
    */
   MappedBoxSet d_global_mapped_boxes;

   /*
    * TODO: I certainly hope we are not using tests on whether the
    * ratio vector is zero to check whether we have an initialized object.
    */
   /*!
    * @brief Refinement ratio from a reference such as level 0.
    *
    * If d_ratio(0) == 0, the object is in uninitialized state.
    */
   IntVector d_ratio;

   /*!
    * @brief Local cell count, excluding periodic images.
    *
    * Unlike d_global_number_of_cells, this parameter is always current.
    */
   size_t d_local_number_of_cells;

   /*!
    * @brief Global cell count, excluding periodic images.
    *
    * This is mutable because it depends on the Boxes and may be
    * saved by a const object if computed.
    *
    * A value < 0 means it has not been computed.
    */
   mutable int d_global_number_of_cells;

   /*!
    * @brief Local Box count, excluding periodic images.
    *
    * Unlike d_global_number_of_mapped_boxes, this parameter is always current.
    */
   size_t d_local_number_of_mapped_boxes;

   /*!
    * @brief Global box count, excluding periodic images.
    *
    * This is mutable because it depends on the Boxes and may be
    * saved by a const object if computed.
    *
    * A value < 0 means it has not been computed.
    */
   mutable int d_global_number_of_mapped_boxes;

   //! @brief Global max box count on any proc, excluding periodic images.
   mutable int d_max_number_of_mapped_boxes;
   //! @brief Global min box count on any proc, excluding periodic images.
   mutable int d_min_number_of_mapped_boxes;
   //! @brief Global max cell count on any proc, excluding periodic images.
   mutable int d_max_number_of_cells;
   //! @brief Global min cell count on any proc, excluding periodic images.
   mutable int d_min_number_of_cells;

   //! @brief Max size of largest local box, one for each block.
   std::vector<IntVector> d_local_max_box_size;
   //! @brief Max size of largest box globally, one for each block.
   mutable std::vector<IntVector> d_global_max_box_size;
   //! @brief Min size of largest local box, one for each block.
   std::vector<IntVector> d_local_min_box_size;
   //! @brief Min size of largest box globally, one for each block.
   mutable std::vector<IntVector> d_global_min_box_size;

   /*!
    * @brief Bounding box of local Boxes, excluding periodic images.
    * One for each block.
    *
    * This is mutable because it depends on the Boxes and may be
    * saved by a const object if computed.
    */
   mutable std::vector<Box> d_local_bounding_box;

   /*!
    * @brief Whether d_local_bounding_box is up to date (or needs
    * recomputing.
    */
   mutable bool d_local_bounding_box_up_to_date;

   /*!
    * @brief Bounding box of global Boxes, excluding periodic images.
    * One for each block.
    *
    * This is mutable because it depends on the Boxes and may be
    * saved by a const object if computed.
    */
   mutable std::vector<Box> d_global_bounding_box;

   /*!
    * @brief Whether globally reduced data is up to date or needs
    * recomputing using cacheGlobalReducedData().
    */
   mutable bool d_global_data_up_to_date;

   /*!
    * @brief State flag.
    *
    * Modified by setParallelState().
    */
   ParallelState d_parallel_state;

   /*!
    * @brief A globalized version of the MappedBoxLevel.
    *
    * Initialized by getGlobalizedVersion().  Deallocated by
    * deallocateGlobalizedVersion().
    *
    * Like other redundant data, this is automatically removed if any
    * method that can potentially change the MappedBoxLevel is called.
    *
    * This is mutable because it is redundant data and gets
    * automatically set as needed.
    */
   mutable MappedBoxLevel const *d_globalized_version;

   /*!
    * @brief Connectors managed by this MappedBoxLevel,
    * providing overlap neighbor data across multiple
    * scopes.
    *
    * This is mutable so it can be allocated as needed (by
    * getPersistentOverlapConnectors()).  We can make it non-mutable
    * by always allocating the PersistentOverlapConnectors in the
    * constructor, but most MappedBoxLevel won't need it at all.
    */
   mutable PersistentOverlapConnectors * d_persistent_overlap_connectors;

   /*!
    * @brief A Handle for Connectors to reference this
    * MappedBoxLevel, used to help prevent invalid Connector
    * data.
    *
    * Connectors reference the handle instead of the
    * MappedBoxLevel directly.  When the MappedBoxLevel
    * changes in a way that can invalidate Connector data,
    * it detaches its handle from itself.  A detached handle
    * tells Connectors that the MappedBoxLevel has changed
    * in a way that can invalidate their data.
    *
    * Note: The automatic detaching mechanism prevents some
    * logic errors.  It cannot prevent incorrect Connector
    * data because correctness depends on the Connector's
    * intended usage.
    */
   mutable tbox::Pointer<MappedBoxLevelHandle> d_handle;

   static tbox::Pointer<tbox::Timer> t_acquire_remote_mapped_boxes;

   /*!
    * @brief Process rank (id),
    * for convenience and data management use after MPI_Finalize.
    */
   int d_rank;

   /*!
    * @brief Number of processes,
    * for convenience and data management use after MPI_Finalize.
    */
   int d_nproc;

   /*!
    * @brief Pointer to the GridGeometry associated with this object.
    */
   tbox::ConstPointer<GridGeometry> d_grid_geometry;

   /*!
    * @brief A LocalId object with value of -1.
    */
   static const LocalId s_negative_one_local_id;

   static tbox::StartupShutdownManager::Handler
   s_initialize_finalize_handler;

};

}
}

#ifdef SAMRAI_INLINE
#include "SAMRAI/hier/MappedBoxLevel.I"
#endif

#endif  // included_hier_MappedBoxLevel
