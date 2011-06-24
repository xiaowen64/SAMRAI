/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Identifier for a MappedBox.
 *
 ************************************************************************/

#ifndef included_hier_MappedBoxId
#define included_hier_MappedBoxId

#include <iostream>

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/GlobalId.h"
#include "SAMRAI/hier/BlockId.h"
#include "SAMRAI/hier/PeriodicId.h"

namespace SAMRAI {
namespace hier {

/*!
 * @brief Identifier for a MappedBox, consisting of a GlobalId,
 * a BlockId and a PeriodicId.
 *
 * MappedBoxes are identified by their BlockId, GlobalId and
 * PeriodicId.  SAMRAI does not support multiblock periodic domains so
 * either the BlockId must be zero or the PeriodicId must be zero, or
 * both.
 *
 * A MappedBox and all its periodic images have the same GlobalId but
 * different PeriodicId.
 *
 * Comparison operators are provided to define sorted ordering of
 * objects.  The BlockId, GlobalId and PeriodicId are used for all
 * comparisons.  The BlockIds are compared first, then GlobalIds,
 * followed by the PeriodicId.
 */
class MappedBoxId
{

public:

   /*!
    * @brief Default constructor uses the default constructors for the
    * BlockId, GlobalId and PeriodicId.
    *
    * The object can be changed using initialize() or by assignment.
    */
   MappedBoxId();

   /*!
    * @brief Initializing constructor.
    *
    * @param[in] local_id
    *
    * @param[in] owner_rank
    *
    * @param[in] block_id
    *
    * @param[in] periodic_id
    */
   explicit MappedBoxId(
      const LocalId &local_id,
      const int owner_rank,
      const BlockId &block_id = BlockId::zero(),
      const PeriodicId &periodic_id = PeriodicId::zero());

   /*!
    * @brief Initializing constructor.
    *
    * @param[in] global_id
    *
    * @param[in] block_id
    *
    * @param[in] periodic_id
    */
   explicit MappedBoxId(
      const GlobalId& id,
      const BlockId &block_id = BlockId::zero(),
      const PeriodicId &periodic_id = PeriodicId::zero());

   /*!
    * @brief Copy constructor.
    *
    * @param[in] other
    */
   MappedBoxId(
      const MappedBoxId& other);

   /*!
    * @brief Destructor.
    */
   virtual ~MappedBoxId();

   /*!
    * @brief Set all the attributes to given values.
    *
    * @param[in] local_id
    *
    * @param[in] owner_rank
    *
    * @param[in] block_id
    *
    * @param[in] periodic_id
    */
   void
   initialize(
      const LocalId &local_id,
      const int owner_rank,
      const BlockId &block_id = BlockId::zero(),
      const PeriodicId &periodic_id = PeriodicId::zero());

   /*!
    * @brief Access the GlobalId.
    */
   const GlobalId&
   getGlobalId() const;

   /*!
    * @brief Access the owner rank.
    */
   int
   getOwnerRank() const;

   /*!
    * @brief Access the LocalId.
    */
   const LocalId
   &getLocalId() const;

   /*!
    * @brief Access the BlockId.
    */
   const BlockId &
   getBlockId() const;

   /*!
    * @brief Access the PeriodicId.
    */
   const PeriodicId &
   getPeriodicId() const;

   /*!
    * @brief Whether the PeriodicId refers to a periodic
    * image.
    */
   bool
   isPeriodicImage() const;


   //@{

   //! @name Comparison operators

   /*!
    * @brief Equality operator.
    *
    * All comparison operators use the BlockId, GlobalId and
    * PeriodicId.
    */
   bool
   operator == (
      const MappedBoxId& r) const;

   /*!
    * @brief Inequality operator.
    *
    * See note on comparison for operator==(const MappedBoxId&);
    */
   bool
   operator != (
      const MappedBoxId& r) const;

   /*!
    * @brief Less-than operator.
    *
    * Compare the owner ranks first; if they compare equal, compare
    * the BlockIds next; if they compare equal, compare the LocalId
    * next; if they compare equal, compare the PeriodicIds.
    */
   bool
   operator < (
      const MappedBoxId& r) const;

   /*!
    * @brief Greater-than operator.
    *
    * Compare the owner ranks first; if they compare equal, compare
    * the BlockIds next; if they compare equal, compare the LocalId
    * next; if they compare equal, compare the PeriodicIds.
    */
   bool
   operator > (
      const MappedBoxId& r) const;

   /*!
    * @brief Less-than-or-equal-to operator.
    */
   bool
   operator <= (
      const MappedBoxId& r) const;

   /*!
    * @brief Greater-than-or-equal-to operator.
    */
   bool
   operator >= (
      const MappedBoxId& r) const;

   //@}


   //@{

   //! @name Support for message passing

   /*!
    * @brief Give number of ints required for putting a MappedBoxId in
    * message passing buffer.
    *
    * @see putToIntBuffer(), getFromIntBuffer().
    */
   static int
   commBufferSize();

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


   /*!
    * @brief Format and insert the object into a stream.
    */
   friend std::ostream&
   operator << (
      std::ostream& co,
      const MappedBoxId& r);

private:
   GlobalId d_global_id;

   BlockId d_block_id;

   PeriodicId d_periodic_id;
};

}
}

#ifdef SAMRAI_INLINE
#include "SAMRAI/hier/MappedBoxId.I"
#endif

#endif  // included_hier_MappedBoxId
