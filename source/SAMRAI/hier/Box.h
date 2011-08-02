/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Box representing a portion of the AMR index space 
 *
 ************************************************************************/

#ifndef included_hier_Box
#define included_hier_Box

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/MappedBoxId.h"
#include "SAMRAI/hier/Transformation.h"
#include "SAMRAI/tbox/DatabaseBox.h"

#include <iostream>

namespace SAMRAI {
namespace hier {

class BoxIterator;

/**
 * Class Box represents a n-dimensional box in the AMR index
 * space.  It is defined by lower and upper bounds given by index objects.
 * The box semantics assumes that the box is cell-centered.  A cell-centered
 * convention implies that the index set covered by the box includes both
 * the lower and upper bounds.
 *
 * @see hier::BoxIterator
 * @see hier::Index
 */

class Box
{
public:
   /**
    * Creates an ``empty'' box.
    */
   explicit Box(
      const tbox::Dimension& dim);

   /**
    * Create a box describing the index space between lower and upper.  The
    * box is assumed to be cell centered and include all elements between lower
    * and upper, including the end points.
    */
   explicit Box(
      const Index& lower,
      const Index& upper);

   /**
    * The copy constructor copies the index space of the argument box.
    */
   Box(
      const Box& box);

   /**
    * Construct a Box from a DatabaseBox.
    */
   explicit Box(
      const tbox::DatabaseBox& box);

   /*!
    * @brief Initializing constructor.
    *
    * @param[in] box
    *
    * @param[in] local_id
    *
    * @param[in] owner_rank
    *
    * @param[in] block_id
    *
    * @param[in] periodic_id The periodic shift number.  If
    * periodic_id is non-zero, specify the Box in the position shifted
    * according to the @c periodic_id.  The default argument for @c
    * periodic_id corresponds to the zero-shift.
    */
   explicit Box(
      const hier::Box& box,
      const LocalId &local_id,
      const int owner_rank,
      const BlockId &block_id = BlockId::zero(),
      const PeriodicId &periodic_id = PeriodicId::zero());

   /*!
    * @brief Constructor with undefined box.
    *
    * The box can be initialized using any of the initialize()
    * methods or by assignment.
    *
    * @param[in] dim
    *
    * @param[in] local_id
    *
    * @param[in] owner_rank
    *
    * @param[in] block_id
    *
    * @param[in] periodic_id
    */
   /*
    * TODO: Constructors initializing boxes are only used to construct
    * temporary objects for finding other MappedBoxes in a
    * stl::set<MappedBox>.  We need another way to do it and get rid
    * of these constructors.
    */
   explicit Box(
      const tbox::Dimension& dim,
      const LocalId &local_id,
      const int owner_rank,
      const BlockId &block_id = BlockId::zero(),
      const PeriodicId &periodic_id = PeriodicId::zero());

   /*!
    * @brief Constructor with undefined box.
    *
    * The box can be initialized using any of the initialize()
    * methods or by assignment.
    *
    * @param[in] dim
    *
    * @param[in] global_id
    *
    * @param[in] block_id
    *
    * @param[in] periodic_id
    */
   /*
    * TODO: Constructors initializing boxes are only used to construct
    * temporary objects for finding other MappedBoxes in a
    * stl::set<MappedBox>.  We need another way to do it and get rid
    * of these constructors.
    */
   explicit Box(
      const tbox::Dimension& dim,
      const GlobalId& id,
      const BlockId &block_id = BlockId::zero(),
      const PeriodicId &periodic_id = PeriodicId::zero());


   /*!
    * @brief Constructor with undefined box and a MappedBoxId.
    *
    * The box can be initialized using any of the initialize()
    * methods or by assignment.
    *
    * @param[in] dim
    *
    * @param[in] mapped_box_id
    *
    * @param[in] periodic_id
    */
   /*
    * TODO: Constructors initializing boxes are only used to construct
    * temporary objects for finding other MappedBoxes in a
    * stl::set<MappedBox>.  We need another way to do it and get rid
    * of these constructors.
    */
   explicit Box(
      const tbox::Dimension& dim,
      const MappedBoxId& mapped_box_id);

   /*!
    * @brief "Copy" constructor allowing change in PeriodicId.
    *
    * @param[in] other Make a copy (but not an exact copy) of this
    * MappedBox.
    *
    * @param[in] periodic_id Periodic shift number to use instead of
    * the shift in @c other.  The box will be set to the real box
    * shifted to the position specified by this value.
    *
    * @param[in] refinement_ratio The index space where the MappedBox
    * lives.
    *
    * @see initialize( const MappedBox&, const int, const IntVector&);
    */
   explicit Box(
      const Box& other,
      const PeriodicId &periodic_id,
      const IntVector& refinement_ratio);

   /*!
    * @brief Set all the attributes of the MappedBox.
    *
    * @param[in] box
    *
    * @param[in] local_id
    *
    * @param[in] owner_rank
    *
    * @param[in] block_id
    *
    * @param[in] periodic_id The periodic shift number.  If
    * this is not zero, specify @c box in the shifted position.  The
    * default argument for @c periodic_id corresponds to the
    * zero-shift.
    */
   void
   initialize(
      const hier::Box& box,
      const LocalId &local_id,
      const int owner_rank,
      const BlockId &block_id = BlockId::zero(),
      const PeriodicId &periodic_id = PeriodicId::zero());

   /*!
    * @brief Set all the attributes identical to that of a reference
    * MappedBox, but with a different PeriodicId.
    *
    * @param[in] other Initialize to this MappedBox, but with the shift
    * given by @c periodic_id.
    *
    * @param[in] periodic_id PeriodicId number to use instead of the
    * shift in @c other.  The box will be set to the real box shifted
    * to the position specified by this value.
    *
    * @param[in] refinement_ratio The index space where the MappedBox
    * lives.
    */
   void
   initialize(
      const Box& other,
      const PeriodicId &periodic_id,
      const IntVector& refinement_ratio);

   /**
    * The destructor for Box.
    */
   ~Box();

   //! @brief Get the MappedBoxId.
   MappedBoxId&
   getId();

   //! @brief Get the MappedBoxId.
   const MappedBoxId&
   getId() const;

   //! @brief Get the Block.
   const BlockId &getBlockId() const;

   //! @brief Get the LocalId.
   const LocalId &getLocalId() const;

   //! @brief Get the GlobalId.
   const GlobalId& getGlobalId() const;

   //! @brief Get the owner rank.
   int getOwnerRank() const;

   /*!
    * @brief Get the periodic shift number.
    *
    * @see PeriodicShiftCatalog.
    */
   const PeriodicId &getPeriodicId() const;

   //! @brief Whether the MappedBox is a periodic image.
   bool isPeriodicImage() const;

   bool isIdEqual(const Box& other) const;

   struct id_equal {
      bool operator() (const Box& b1, const Box& b2) const
      {
         return b1.isIdEqual(b2);
      }
   };

   struct id_less {
      bool operator() (const Box& b1, const Box& b2) const
      {
         return b1.getId() < b2.getId();
      }
   };

   /*!
    * @brief Give number of ints required for putting a MappedBox in
    * message passing buffer.
    *
    * This number is independent of instance but dependent on
    * dimension.
    *
    * @see putToIntBuffer(), getFromIntBuffer().
    */
   static int
   commBufferSize(
      const tbox::Dimension& dim);

   /*!
    * @brief Put self into a int buffer.
    *
    * This is the opposite of getFromIntBuffer().  Number of ints
    * written is given by commBufferSize().
    */
   void
   putToIntBuffer(
      int* buffer) const;

   /*!
    * @brief Set attributes according to data in int buffer.
    *
    * This is the opposite of putToIntBuffer().  Number of ints read
    * is given by commBufferSize().
    */
   void
   getFromIntBuffer(
      const int* buffer);

   /**
    * The assignment operator copies the index space of the argument box.
    *
    * An assignment to an uninitialized box is allowed but assigning
    * from an uninitialized box will result in an assert.
    */
   Box&
   operator = (
      const Box& box);

   /**
    * Return a non-const lower index of the box.
    */
   Index&
   lower();

   /**
    * Return a non-const upper index of the box.
    */
   Index&
   upper();

   /**
    * Return a const lower index of the box.
    */
   const Index&
   lower() const;

   /**
    * Return a const upper index of the box.
    */
   const Index&
   upper() const;

   /**
    * Return the i'th component (non-const) of the lower index of the box.
    */
   int&
   lower(
      const int i);

   /**
    * Return the i'th component (non-const) of the upper index of the box.
    */
   int&
   upper(
      const int i);

   /**
    * Return the i'th component (const) of the lower index of the box.
    */
   const int&
   lower(
      const int i) const;

   /**
    * Return the i'th component (const) of the upper index of the box.
    */
   const int&
   upper(
      const int i) const;

   /**
    * Set the index space represented by the box to empty.
    */
   void
   setEmpty();

   /**
    * @brief Return whether the box is ``empty''.
    *
    * isEmpty() is preferred to match "is" standard syntax for
    * boolean methods.
    *
    * @see isEmpty()
    */
   bool
   empty() const;

   /**
    * @brief Return whether the box is ``empty''.
    *
    * A box is empty if any of the lower bounds is greater than the
    * corresponding upper bound.  An empty box has a size of zero.
    */
   bool
   isEmpty() const;

   /**
    * Return the number of cells (an integer) represented by the box in
    * the given coordinate direction.
    */
   int
   numberCells(
      const int i) const;

   /**
    * Return the number of cells (a vector of integers) represented by
    * the box in every coordinate direction.
    */
   IntVector
   numberCells() const;

   /**
    * Calculate the number of indices represented by the box.  If the box
    * is empty, then the number of index points within the box is zero.
    */
   int
   size() const;

   /**
    *  Return the dimension of the box that is longest.
    */
   int
   longestDimension() const;

   /**
    * Given an index, calculate the offset of the index into the box.
    * This function assumes column-major (e.g., Fortran) ordering of
    * the indices within the box.  This operation is a convenience
    * function for later array indexing operations.
    */
   int
   offset(
      const Index& p) const;

   /**
    * Given an offset, calculate the index of the offset into the box.
    * This function assumes column-major (e.g., Fortran) ordering of
    * the indices within the box.  This operation is a convenience
    * function for later array indexing operations.
    */
   Index
   index(
      const int offset) const;

   /**
    * Check whether an index lies within the bounds of the box.
    */
   bool
   contains(
      const Index& p) const;

   /**
    * Check whether a given box lies within the bounds of the box.
    *
    * If @c b is empty, always return true.
    */
   bool
   contains(
      const Box& b) const;

   /**
    * Check whether two boxes represent the same portion of index space.
    */
   bool isSpatiallyEqual(const Box& box) const;
   struct box_equality {
      bool operator() (const Box& b1, const Box& b2) const
      {
         return b1.isSpatiallyEqual(b2);
      }
   };


   /**
    * Calculate the intersection of the index spaces of two boxes.  The
    * intersection with an empty box always yields an empty box.
    */
   Box
   operator * (
      const Box& box) const;

   /**
    * Preallocated box intersection.
    * Calculate the intersection of the index spaces of two boxes.  The
    * intersection with an empty box always yields an empty box.
    */
   void
   intersect(
      const Box& other,
      Box &result) const;

   /**
    * Return true if two boxes have a non-empty intersection.
    * Otherwise, return false.
    */
   bool
   intersects(
      const Box& box) const;

   /**
    * Calculate the bounding box of two boxes.  Note that this is not
    * the union of the two boxes (since union is not closed over boxes),
    * but rather the smallest box that contains both boxes.
    */
   Box
   operator + (
      const Box& box) const;

   /**
    * Increase the bounding box to include the argument box.
    */
   Box&
   operator += (
      const Box& box);

   /**
    * Return true if this box can be coalesced with the argument box,
    * and set this box to the union of the boxes.  Otherwise, return false
    * and leave boxes as is.  Two boxes may be coalesced if their
    * union is a box (recall that index set union is not closed over boxes).
    * If either box is empty, then the return value is true and this box
    * becomes the union of the two.
    */
   bool
   coalesceWith(
      const Box& box);

   /**
    * Grow a box by the specified ghost cell width.  The lower bound is
    * decremented by the width, and the upper bound is incremented by the
    * width.  All dimensions are grown by the corresponding component in
    * the IntVector; ghost cell widths may be different in each dimension.
    * Negative ghost cell widths will shrink the box.
    */
   void
   grow(
      const IntVector& ghosts);

   /**
    * Grow a box by the specified ghost cell width in the given coordinate
    * direction in index space.  The lower bound is decremented by the
    * width, and the upper bound is incremented by the width.  Note that
    * negative ghost cell widths will shrink the box.
    */
   void
   grow(
      const int direction,
      const int ghosts);

   /**
    * Similar to grow() functions. However, box is only grown in lower
    * directions (i.e., only lower index is changed).
    */
   void
   growLower(
      const IntVector& ghosts);

   /**
    * Similar to grow() functions. However, box is only grown in lower
    * bound of given direction in index space.
    */
   void
   growLower(
      const int direction,
      const int ghosts);

   /**
    * Similar to grow() function. However, box is only grown in upper
    * directions (i.e., only upper index is changed).
    */
   void
   growUpper(
      const IntVector& ghosts);

   /**
    * Similar to grow() functions. However, box is only grown in upper
    * bound of given direction in index space.
    */
   void
   growUpper(
      const int direction,
      const int ghosts);

   /**
    * Similar to growUpper() and growLower() functions. However, box is
    * lengthened (never shortened).  The sign of @c ghosts refer to whether
    * the box is lengthened in the upper or lower side.
    */
   void
   lengthen(
      const int direction,
      const int ghosts);

   /**
    * Similar to growUpper() and growLower() functions. However, box is
    * shortened (never lengthened).  The sign of @c ghosts refer to whether
    * the box is shortened in the upper or lower side.
    */
   void
   shorten(
      const int direction,
      const int ghosts);

   /**
    * Shift a box by the specified amount (a vector of integers).
    * The new box is located at (lower+offset, upper+offset).
    */
   void
   shift(
      const IntVector& offset);

   /**
    * Similar to shift() function above, but shift occurs only in specified
    * direction in index space.  The new box is located at (lower+offset,
    * upper+offset) in that direction.
    */
   void
   shift(
      const int direction,
      const int offset);

   /**
    * Rotate 90 degrees around origin.
    */
   void
   rotate(
      const Transformation::RotationIdentifier rotation_ident);

   /**
    * Refine the index space of a box by specified vector ratio.  Each
    * component of the box is multiplied by the refinement ratio,
    * then @c (ratio-1) is added to the upper corner.
    */
   void
   refine(
      const IntVector& ratio);

   /**
    * Coarsen the index space of a box by specified vector ratio.  Each
    * component is divided by the specified coarsening ratio and rounded
    * (if necessary) such that the coarsened box contains the cells that
    * are the parents of the refined box.  In other words, refining a
    * coarsened box will always yield a box that is equal to or larger
    * than the original box.
    */
   void
   coarsen(
      const IntVector& ratio);

   /**
    * This assignment operator constructs a Box given a DatabaseBox.
    */
   Box&
   operator = (
      const tbox::DatabaseBox& box);

   /**
    * Sets a Box from a tbox::DatabaseBox and returns a reference to
    * the Box.
    */
   Box&
   Box_from_DatabaseBox(
      const tbox::DatabaseBox& box);

   /**
    * Sets a Box from a DatabaseBox.
    */
   void
   set_Box_from_DatabaseBox(
      const tbox::DatabaseBox& box);

   /**
    * Returns a tbox::DatabaseBox generated from a Box.
    */
   tbox::DatabaseBox
   DatabaseBox_from_Box() const;

   /**
    * Type conversion from Box to Box
    */
   operator tbox::DatabaseBox ();

   /**
    * Type conversion from Box to Box
    */
   operator tbox::DatabaseBox () const;

   /**
    * Utility function to grow a box by the specified vector ghost cell
    * width.  A new box is returned and the argument is not changed.
    */
   static Box
   grow(
      const Box& box,
      const IntVector& ghosts);

   /**
    * Utility function to shift a box by the specified offset.  A new
    * box is returned and the argument is not changed.
    */
   static Box
   shift(
      const Box& box,
      const IntVector& offset);

   /**
    * Utility function to refine the index space of a box by the specified
    * refinement ratio.  A new box is returned and the argument is not changed.
    */
   static Box
   refine(
      const Box& box,
      const IntVector& ratio);

   /**
    * Utility function to coarsen the index space of a box by the specified
    * coarsening ratio.  A new box is returned and the argument is not changed.
    */
   static Box
   coarsen(
      const Box& box,
      const IntVector& ratio);

   /**
    * Return the dimension of this object.
    */
   const tbox::Dimension&
   getDim() const;

   /**
    * Read the box description in the form [L,U], where L and U are the
    * lower and upper bounds of the box.
    */
   friend std::istream&
   operator >> (
      std::istream& s,
      Box& box);

   /**
    * Output the box description in the form [L,U], where L and U are the
    * lower and upper bounds of the box.
    */
   friend std::ostream&
   operator << (
      std::ostream& s,
      const Box& box);

   /*!
    * @brief Return an empty Box of the specified dimension.
    *
    * Can be used to avoid object creation overheads.
    */
   static const Box&
   getEmptyBox(
      const tbox::Dimension& dim);

   /**
    * Returns a Box that represents the maximum allowed index extents
    * for a given dimension.   The "universe" that can be represented.
    */
   static const Box&
   getUniverse(
      const tbox::Dimension& dim);

   /**
    * A box iterator iterates over the elements of a box.  This class is
    * defined elsewhere, and the typedef is used to point to that class.
    */
   typedef BoxIterator Iterator;

   template<class>
   friend class ::SAMRAI::pdat::ArrayData;
   friend class BoxIterator;
   friend class std::vector<Box>;

#ifdef BOX_TELEMETRY
   // These are to optionally track the cumulative number of Boxes constructed,
   // the cumulative number of Box assignments, and the high water mark of
   // Boxes in existance at any given time.
   static int s_cumulative_constructed_ct;

   static int s_cumulative_assigned_ct;

   static int s_active_ct;

   static int s_high_water;
#endif

private:
   /**
    * The default constructor creates an uninitialized box.
    *
    * This should never be invoked, it will cause assertions
    */
   Box();

   static int
   coarsen(
      const int index,
      const int ratio);

   static bool
   coalesceIntervals(
      const int* lo1,
      const int* hi1,
      const int* lo2,
      const int* hi2,
      const int dim);

   void
   rotateAboutAxis(
      const int axis,
      const int num_rotations);

   Index d_lo;
   Index d_hi;
   MappedBoxId d_id;

   /*!
    * @brief Initialize static objects and register shutdown routine.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   initializeCallback();

   /*!
    * @brief Method registered with ShutdownRegister to cleanup statics.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   finalizeCallback();

   /*
    * Array of empty boxes for each dimension.  Preallocated
    * as a performance enhancement to avoid constructing
    * them in multiple places.
    */

   static Box* s_emptys[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

   /*
    * Box that represents the maximum allowed index extents,
    * the "universe" that can be represented.
    */
   static Box* s_universes[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

   static tbox::StartupShutdownManager::Handler
   s_initialize_finalize_handler;
};

/**
 * Class BoxIterator is an iterator that provides methods for
 * stepping through the index space associated with a box.  The indices
 * are enumerated in column-major (e.g., Fortran) order.  The iterator
 * should be used as follows:
 * \verbatim
 * Box box;
 * ...
 * for (Box::Iterator b(box); b; b++) {
 *    // use index b of the box
 * }
 * \endverbatim
 * Note that the box iterator may not compile to efficient code, depending
 * on your compiler.  Many compilers are not smart enough to optimize the
 * looping constructs and indexing operations.
 *
 * @see hier::Index
 * @see hier::Box
 */

class BoxIterator
{
public:
   /**
    * Default constructor for the box iterator.  The iterator must
    * be initialized before it can be used to iterate over a box.
    *
    * @see initialize().
    */
   BoxIterator();

   /**
    * Constructor for the box iterator.  The iterator will enumerate the
    * indices in the argument box.
    */
   BoxIterator(
      const Box& box);

   /**
    * Copy constructor for the box iterator.
    */
   BoxIterator(
      const BoxIterator& iterator);

   /**
    * Assignment operator for the box iterator.
    */
   BoxIterator&
   operator = (
      const BoxIterator& iterator);

   /**
    * Destructor for the box iterator.
    */
   ~BoxIterator();

   /**
    * @brief Initializer for the box iterator.
    *
    * The iterator will enumerate the indices in the box.
    */
   void
   initialize(
      const Box& box);

   /**
    * Return the current index in the box.  This operation is undefined
    * if the iterator is past the last Index in the box.
    */
   const Index&
   operator * () const;

   /**
    * Return the current index in the box.  This operation is undefined
    * if the iterator is past the last Index in the box.
    */
   const Index&
   operator () () const;

   /**
    * Return true if the iterator points to a valid index in the box.
    */
   operator bool () const;

#ifndef LACKS_BOOL_VOID_RESOLUTION
   /**
    * Return a non-NULL if the iterator points to a valid index in the box.
    */
   operator const void
   * () const;
#endif

   /**
    * Return whether the iterator points to a valid item in the box.  This
    * operator mimics the !p operation applied to a pointer p.
    */
   bool
   operator ! () const;

   /**
    * Increment the iterator to point to the next index in the box.
    */
   void
   operator ++ (
      int);

   /**
    * Test two iterators for equality (same index value).
    */
   bool
   operator == (
      const BoxIterator& iterator) const;

   /**
    * Test two iterators for inequality (different index values).
    */
   bool
   operator != (
      const BoxIterator& iterator) const;

private:
   Index d_index;
   Box d_box;

};

}
}

#ifdef SAMRAI_INLINE
#include "SAMRAI/hier/Box.I"
#endif

#endif
