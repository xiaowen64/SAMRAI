/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Special iterator for MappedBoxSet.
 *
 ************************************************************************/
#ifndef included_hier_MappedBoxSetSingleBlockIterator
#define included_hier_MappedBoxSetSingleBlockIterator

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/MappedBoxSet.h"
#include "SAMRAI/hier/BlockId.h"

namespace SAMRAI {
namespace hier {

/*!
 * @brief MappedBoxSet iterator picking items with a specified
 * BlockId.
 *
 * This iterator runs through all MappedBoxes in a MappedBoxSet that
 * has the given BlockId.  The iterator runs through the MappedBoxes
 * in the order they appear in the MappedBoxSet, skipping over
 * MappedBoxes that do not have the specified owner rank.
 */
class MappedBoxSetSingleBlockIterator {

public:

   /*!
    * @brief Constructor
    *
    * @param [i] container
    * @param [i] block_id
    */
   MappedBoxSetSingleBlockIterator(
      const MappedBoxSet &container,
      const BlockId &block_id);

   //! @brief Destructor
   ~MappedBoxSetSingleBlockIterator();

   /*!
    * @brief Assignment operator.
    */
   MappedBoxSetSingleBlockIterator&
   operator = (
      const MappedBoxSetSingleBlockIterator& r);

   /*!
    * @brief Dereference operator mimicking a pointer dereference.
    */
   const MappedBox &operator * () const;

   /*!
    * @brief Dereference operator mimicking a pointer dereference.
    */
   const MappedBox *operator -> () const;

   /*!
    * @brief Equality comparison.
    */
   bool operator == (const MappedBoxSetSingleBlockIterator& r) const;

   /*!
    * @brief Inequality comparison.
    */
   bool operator != (const MappedBoxSetSingleBlockIterator& r) const;

   /*!
    * @brief Whether the iterator can be dereferenced.  When the
    * iterator reaches its end, this returns false.
    */
   bool isValid() const;

   /*!
    * @brief Pre-increment iterator.
    *
    * Pre-increment increment the iterator and returns the incremented
    * state.
    */
   MappedBoxSetSingleBlockIterator &operator ++ ();

   /*!
    * @brief Post-increment iterator.
    *
    * Post-increment saves the iterator, increment it and returns the
    * saved iterator.
    */
   MappedBoxSetSingleBlockIterator operator ++ (int);

   /*!
    * @brief Returns the number of MappedBoxSets being iterated through.
    */
   int count() const;

private:

   /*!
    * @brief MappedBoxSet being iterated through.
    */
   const MappedBoxSet* d_mapped_boxes;

   /*!
    * @brief The BlockId.
    */
   BlockId d_block_id;

   /*!
    * @brief The iterator.
    */
   MappedBoxSet::const_iterator d_iter;

};


}
}

#ifdef SAMRAI_INLINE
// #include "SAMRAI/hier/MappedBoxSetSingleBlockIterator.I"
#endif

#endif  // included_hier_MappedBoxSetSingleBlockIterator
