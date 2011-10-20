/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   A const iterator over containers of boxes.
 *
 ************************************************************************/

#ifndef included_hier_BoxContainerConstIterator
#define included_hier_BoxContainerConstIterator

#include <list>
#include <set>

namespace SAMRAI {
namespace hier {

class Box;
class BoxContainer;
class BoxContainerIterator;

/*!
 * @brief A immutable iterator over the boxes in a BoxContainer.
 *
 * If iterating over an ordered BoxContainer, then iteration will follow
 * the sequence of the BoxId-based ordering of the container.  If iterating
 * over an unordered BoxContainer, the sequence of the iteration will be
 * based on how the members of the container were added.
 *
 * @see BoxContainer
 * @see BoxContainerIterator
 */
class BoxContainerConstIterator
{
   friend class BoxContainer;
   friend class BoxContainerIterator;

public:

   /*!
    * @brief Constructor for the BoxContainerConstIterator.
    *
    * The iterator will point to the beginning or the end of the argument
    * container, depending on the from_start argument
    *
    * @param[in] container The container whose members are iterated.
    * @param[in] from_start true if iteration starts at beginning of container.
    */
   explicit BoxContainerConstIterator(
      const BoxContainer& container,
      bool from_start = true);

   /*!
    * @brief Copy constructor.
    *
    * @param[in] other
    */
   BoxContainerConstIterator(
      const BoxContainerConstIterator& other);

   /*!
    * @brief Copy constructor to copy mutable iterator to an immutable
    * iterator.
    */
   BoxContainerConstIterator(
      const BoxContainerIterator& other);

   /*!
    * @brief Assignment operator.
    *
    * @param[in] rhs
    */
   BoxContainerConstIterator&
   operator = (
      const BoxContainerConstIterator& rhs);

   /*!
    * @brief The destructor releases all storage.
    */
   ~BoxContainerConstIterator();

   /*!
    * @brief Get box corresponding to iterator's position in container.
    *
    * @return An immutable reference to the current Box in the iteration.
    */
   const Box&
   operator * () const;

   /*!
    * @brief Get box corresponding to iterator's position in container.
    *
    * @return An immutable reference to the current Box in the iteration.
    */
   const Box&
   operator () () const;

   /*!
    * @brief Get pointer to box at iterator's position in container.
    *
    * @return Const pointer to the current box.
    */
   const Box*
   operator -> () const;

   /*!
    * @brief Post-increment iterator to point to next box in the container.
    *
    * @return Iterator at the position in the container before the increment.
    */
   BoxContainerConstIterator
   operator ++ (
      int);

   /*!
    * @brief Pre-increment iterator to point to next box in the container.
    *
    * @return Reference to iterator at the position in the container after
    * the increment.
    */
   const BoxContainerConstIterator&
   operator ++ ();

   /*!
    * @brief Post-decrement iterator to point to next box in the container.
    *
    * @return Iterator at the position in the container before the decrement.
    */
   BoxContainerConstIterator
   operator -- (
      int);

   /*!
    * @brief Pre-decrement iterator to point to next box in the container.
    *
    * @return Reference to iterator at the position in the container after
    * the decrement.
    */
   const BoxContainerConstIterator&
   operator -- ();

   /*!
    * @brief Equality operator.
    *
    * @return true if both iterators point to the same box.
    *
    * @param[in] other
    */
   bool
   operator == (
      const BoxContainerConstIterator& other) const;

   /*!
    * @brief Inequality operator.
    *
    * @return true if both iterators point to different boxes.
    *
    * @param[in] other
    */
   bool
   operator != (
      const BoxContainerConstIterator& other) const;

private:
   /*
    * Default constructor just to be clear that there is none.
    */
   BoxContainerConstIterator();

   /*
    * Underlying iterator to be used when unordered.
    */
   std::list<Box>::const_iterator d_list_iter;

   /*
    * Underlying iterator to be used when ordered.
    */
   std::set<Box*>::const_iterator d_set_iter;

   bool d_ordered;
};

}
}

#ifdef SAMRAI_INLINE
#include "SAMRAI/hier/BoxContainerConstIterator.I"
#endif

#endif
