/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Special iterator for BoxSet.
 *
 ************************************************************************/
#ifndef included_hier_BoxSetSingleOwnerIterator
#define included_hier_BoxSetSingleOwnerIterator

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BoxContainerConstIterator.h"

namespace SAMRAI {
namespace hier {

/*!
 * @brief BoxSet iterator picking items with a specified
 * owner rank.
 *
 * This iterator runs through all Boxes in a BoxSet that
 * has the given owner rank.  The iterator runs through the
 * Boxes in the order they appear in the BoxSet, skipping
 * over Boxes that do not have the specified owner rank.
 */
class BoxSetSingleOwnerIterator
{

public:
   /*!
    * @brief Constructor
    *
    * @param [i] container
    * @param [i] owner_rank
    */
   BoxSetSingleOwnerIterator(
      const BoxContainer& container,
      const int& owner_rank);

   //! @brief Destructor
   ~BoxSetSingleOwnerIterator();

   /*!
    * @brief Assignment operator.
    */
   BoxSetSingleOwnerIterator&
   operator = (
      const BoxSetSingleOwnerIterator& r);

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
      const BoxSetSingleOwnerIterator& r) const;

   /*!
    * @brief Inequality comparison.
    */
   bool
   operator != (
      const BoxSetSingleOwnerIterator& r) const;

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
   BoxSetSingleOwnerIterator&
   operator ++ ();

   /*!
    * @brief Post-increment iterator.
    *
    * Post-increment saves the iterator, increment it and returns the
    * saved iterator.
    */
   BoxSetSingleOwnerIterator
   operator ++ (
      int);

private:
   /*!
    * @brief BoxSet being iterated through.
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
// #include "SAMRAI/hier/BoxSetSingleOwnerIterator.I"
#endif

#endif  // included_hier_BoxSetSingleOwnerIterator
