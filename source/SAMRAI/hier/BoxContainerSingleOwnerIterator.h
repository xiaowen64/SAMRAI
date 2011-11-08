/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Special iterator for BoxContainer.
 *
 ************************************************************************/
#ifndef included_hier_BoxContainerSingleOwnerIterator
#define included_hier_BoxContainerSingleOwnerIterator

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BoxContainerConstIterator.h"

namespace SAMRAI {
namespace hier {

/*!
 * @brief BoxContainer iterator picking items with a specified
 * owner rank.
 *
 * This iterator runs through all Boxes in a BoxContainer that
 * has the given owner rank.  The iterator runs through the
 * Boxes in the order they appear in the BoxContainer, skipping
 * over Boxes that do not have the specified owner rank.
 */
class BoxContainerSingleOwnerIterator
{

public:
   /*!
    * @brief Constructor
    *
    * @param [i] container
    * @param [i] owner_rank
    */
   BoxContainerSingleOwnerIterator(
      const BoxContainer& container,
      const int& owner_rank);

   //! @brief Destructor
   ~BoxContainerSingleOwnerIterator();

   /*!
    * @brief Assignment operator.
    */
   BoxContainerSingleOwnerIterator&
   operator = (
      const BoxContainerSingleOwnerIterator& r);

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
      const BoxContainerSingleOwnerIterator& r) const;

   /*!
    * @brief Inequality comparison.
    */
   bool
   operator != (
      const BoxContainerSingleOwnerIterator& r) const;

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
   BoxContainerSingleOwnerIterator&
   operator ++ ();

   /*!
    * @brief Post-increment iterator.
    *
    * Post-increment saves the iterator, increment it and returns the
    * saved iterator.
    */
   BoxContainerSingleOwnerIterator
   operator ++ (
      int);

private:
   /*!
    * @brief BoxContainer being iterated through.
    */
   const BoxContainer* d_mapped_boxes;

   /*!
    * @brief The owner_rank.
    */
   int d_owner_rank;

   /*!
    * @brief The iterator.
    */
   BoxContainer::ConstIterator d_iter;

};

}
}

#ifdef SAMRAI_INLINE
// #include "SAMRAI/hier/BoxContainerSingleOwnerIterator.I"
#endif

#endif  // included_hier_BoxContainerSingleOwnerIterator
