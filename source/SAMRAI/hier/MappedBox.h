/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Box in the distributed box graph.
 *
 ************************************************************************/
#ifndef included_hier_MappedBox
#define included_hier_MappedBox

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/MappedBoxId.h"

#include <iostream>

namespace SAMRAI {
namespace hier {

/*!
 * @brief A MappedBox represents a Box that is mapped to (owned by) a
 * process.
 *
 * The MappedBox has a MappedBoxId for the purpose of identifying it.
 * The MappedBoxId also has information about the MappedBox's BlockId,
 * its PeriodicId and its GlobalId.
 *
 * A MappedBox describes Box that is "owned" by a certain process.
 * The MappedBox is identified by its MappedBoxId identifier, which
 * includes the owner's rank and its LocalId.  The LocalId is
 * typically assigned to the Box by the process that owns it.
 *
 * Ownership typically means that only the owner should change the
 * object.  Other processes may work with a copy of the object but
 * should not change it.  (Any changes made by the owner should
 * eventually be communicated to other processes for re-establishing
 * consistency.  Such coordinations are the responsibility of
 * algorithms manipulating the MappedBoxes.)
 *
 * Each MappedBox may be real or a periodic image of a real MappedBox.
 * A MappedBox that is a periodic image has the same GlobalId as its
 * real MappedBox but a different periodic shift number.  A real
 * MappedBox has a periodic shift number of zero.  Periodic images
 * should set the Box portion of the MappedBox at the shifted
 * postition rather than the real position.  See PeriodicShiftCatalog.
 *
 * Comparison operators are implemented for sorting MappedBoxes and
 * instantiating STL containers of MappedBoxes.  Comparisons are
 * strictly by the MappedBoxId.  (The Box is ignored by comparisons.)
 */
class MappedBox
{

public:
   /*!
    * @brief Non-initializing constructor.
    *
    * The object can be initialized using any of the initialize()
    * methods or by assignment.
    *
    * @param[in] dim
    */
   explicit MappedBox(
      const tbox::Dimension& dim);

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
   explicit MappedBox(
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
   explicit MappedBox(
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
   explicit MappedBox(
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
   explicit MappedBox(
      const tbox::Dimension& dim,
      const MappedBoxId& mapped_box_id);

   /*!
    * @brief Copy constructor making an exact copy.
    *
    * @param[in] r
    */
   MappedBox(
      const MappedBox& r);

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
   explicit MappedBox(
      const MappedBox& other,
      const PeriodicId &periodic_id,
      const IntVector& refinement_ratio);

   /*!
    * @brief Destructor.
    */
   virtual ~MappedBox();



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
      const MappedBox& other,
      const PeriodicId &periodic_id,
      const IntVector& refinement_ratio);


   //@{

   //! @name Accessors

   //! @brief Get the Box.
   hier::Box& getBox();

   //! @brief Get the Box.
   const hier::Box&
   getBox() const;

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

   /*!
    * @brief Return the dimension of this object.
    */
   const tbox::Dimension&
   getDim() const;

   //@}



   //@{

   //! @name Comparison operators

   /*!
    * @brief Equality operator.
    *
    * All comparison operators use only the MappedBoxId.  The
    * comparison yield the same results that comparing just the
    * MappedBoxIds of the MappedBoxes yield.  The box does not affect
    * comparison results.  This allows finding a specific box by
    * knowing just its id.
    */
   bool operator == (
      const MappedBox& r) const;

   /*!
    * @brief Inequality operator.
    *
    * See note on comparison for operator==(const MappedBox&).
    */
   bool operator != (
      const MappedBox& r) const;

   /*!
    * @brief Less-than operator.
    *
    * See note on comparison for operator==(const MappedBox&).
    */
   bool operator < (
      const MappedBox& r) const;

   /*!
    * @brief Greater-than operator.
    *
    * See note on comparison for operator==(const MappedBox&).
    */
   bool operator > (
      const MappedBox& r) const;

   /*!
    * @brief Less-than-or-equal-to operator.
    *
    * See note on comparison for operator==(const MappedBox&).
    */
   bool operator <= (
      const MappedBox& r) const;

   /*!
    * @brief Greater-than-or-equal-to operator.
    *
    * See note on comparison for operator==(const MappedBox&).
    */
   bool operator >= (
      const MappedBox& r) const;

   //@}


   //@{

   //! @name Support for message passing

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

   //@}


   //! @brief Stream-insert operator.
   friend std::ostream&
   operator << (
      std::ostream& co,
      const MappedBox& r);

   /*!
    * @brief Allows container to construct objects without dimension.
    */
   friend class std::vector<MappedBox>;

private:
   /**
    * Needed by the stupid STL vector class.
    */
   MappedBox();

   /*!
    * @brief MappedBox identifier.
    *
    * This state variable is used in sorting a set of MappedBox objects.
    */
   MappedBoxId d_id;

   Box d_box;
};

}
}

#ifdef SAMRAI_INLINE
#include "SAMRAI/hier/MappedBox.I"
#endif

#endif  // included_hier_MappedBox
