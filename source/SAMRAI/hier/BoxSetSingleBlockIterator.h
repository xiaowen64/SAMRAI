/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Special iterator for BoxSet.
 *
 ************************************************************************/
#ifndef included_hier_BoxSetSingleBlockIterator
#define included_hier_BoxSetSingleBlockIterator

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/BlockId.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BoxContainerConstIterator.h"

namespace SAMRAI {
namespace hier {


/*!
 * @brief BoxContainer iterator picking items with a specified
 * BlockId.
 *
 * This iterator runs through all Boxes in a BoxContainer that
 * has the given BlockId.  The iterator runs through the Boxes
 * in the order they appear in the BoxContainer, skipping over
 * Boxes that do not have the specified owner rank.
 */
class BoxSetSingleBlockIterator
{

public:
   /*!
    * @brief Constructor
    *
    * @param [i] container
    * @param [i] block_id
    */
   BoxSetSingleBlockIterator(
      const BoxContainer& container,
      const BlockId& block_id);

   //! @brief Destructor
   ~BoxSetSingleBlockIterator();

   /*!
    * @brief Assignment operator.
    */
   BoxSetSingleBlockIterator&
   operator = (
      const BoxSetSingleBlockIterator& r);

   /*!
    * @brief Dereference operator mimicking a pointer dereference.
    */
   const Box&
   operator * () const;

   /*!
    * @brief Dereference operator mimicking a pointer dereference.
    */
   const Box *
   operator -> () const;

   /*!
    * @brief Equality comparison.
    */
   bool
   operator == (
      const BoxSetSingleBlockIterator& r) const;

   /*!
    * @brief Inequality comparison.
    */
   bool
   operator != (
      const BoxSetSingleBlockIterator& r) const;

   /*!
    * @brief Whether the iterator can be dereferenced.  When the
    * iterator reaches its end, this returns false.
    */
   bool
   isValid() const;

   /*!
    * @brief Pre-increment iterator.
    *
    * Pre-increment increment the iterator and returns the incremented
    * state.
    */
   BoxSetSingleBlockIterator&
   operator ++ ();

   /*!
    * @brief Post-increment iterator.
    *
    * Post-increment saves the iterator, increment it and returns the
    * saved iterator.
    */
   BoxSetSingleBlockIterator
   operator ++ (
      int);

   /*!
    * @brief Returns the number of BoxContainers being iterated through.
    */
   int
   count() const;

private:
   /*!
    * @brief BoxContainer being iterated through.
    */
   const BoxContainer* d_mapped_boxes;

   /*!
    * @brief The BlockId.
    */
   BlockId d_block_id;

   /*!
    * @brief The iterator.
    */
   BoxContainer::ConstIterator d_iter;

};

}
}

#ifdef SAMRAI_INLINE
// #include "SAMRAI/hier/BoxSetSingleBlockIterator.I"
#endif

#endif  // included_hier_BoxSetSingleBlockIterator
