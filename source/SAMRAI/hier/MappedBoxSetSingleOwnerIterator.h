/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Special iterator for MappedBoxSet.
 *
 ************************************************************************/
#ifndef included_hier_MappedBoxSetSingleOwnerIterator
#define included_hier_MappedBoxSetSingleOwnerIterator

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/MappedBoxSet.h"

namespace SAMRAI {
namespace hier {

/*!
 * @brief MappedBoxSet iterator picking items with a specified
 * owner rank.
 *
 * This iterator runs through all MappedBoxes in a MappedBoxSet that
 * has the given owner rank.  The iterator runs through the
 * MappedBoxes in the order they appear in the MappedBoxSet, skipping
 * over MappedBoxes that do not have the specified owner rank.
 */
class MappedBoxSetSingleOwnerIterator {

public:

   /*!
    * @brief Constructor
    *
    * @param [i] container
    * @param [i] owner_rank
    */
   MappedBoxSetSingleOwnerIterator(
      const MappedBoxSet &container,
      const int &owner_rank);

   //! @brief Destructor
   ~MappedBoxSetSingleOwnerIterator();

   /*!
    * @brief Assignment operator.
    */
   MappedBoxSetSingleOwnerIterator&
   operator = (
      const MappedBoxSetSingleOwnerIterator& r);

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
   bool operator == (const MappedBoxSetSingleOwnerIterator& r) const;

   /*!
    * @brief Inequality comparison.
    */
   bool operator != (const MappedBoxSetSingleOwnerIterator& r) const;

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
   MappedBoxSetSingleOwnerIterator &operator ++ ();

   /*!
    * @brief Post-increment iterator.
    *
    * Post-increment saves the iterator, increment it and returns the
    * saved iterator.
    */
   MappedBoxSetSingleOwnerIterator operator ++ (int);

private:

   /*!
    * @brief MappedBoxSet being iterated through.
    */
   const MappedBoxSet* d_mapped_boxes;

   /*!
    * @brief The owner_rank.
    */
   int d_owner_rank;

   /*!
    * @brief The iterator.
    */
   MappedBoxSet::const_iterator d_iter;

};


}
}

#ifdef SAMRAI_INLINE
// #include "SAMRAI/hier/MappedBoxSetSingleOwnerIterator.I"
#endif

#endif  // included_hier_MappedBoxSetSingleOwnerIterator
