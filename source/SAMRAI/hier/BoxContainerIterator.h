/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   An iterator over containers of boxes.
 *
 ************************************************************************/

#ifndef included_hier_BoxContainerIterator
#define included_hier_BoxContainerIterator

#include <list>
#include <set>

namespace SAMRAI {
namespace hier {

class Box;
class BoxContainer;
class BoxContainerConstIterator;

/*!
 * @brief A mutable iterator over the boxes in a BoxContainer.
 *
 * If iterating over an ordered BoxContainer, then iteration will follow
 * the sequence of the BoxId-based ordering of the container.  If iterating
 * over an unordered BoxContainer, the sequence of the iteration will be
 * based on how the members of the container were added.
 *
 * @see BoxContainer
 * @see BoxContainerConstIterator
 */
class BoxContainerIterator
{
   friend class BoxContainer;
   friend class BoxContainerConstIterator;

public:

   /*!
    * @brief Constructor for the BoxContainerIterator.
    *
    * The iterator will point to the beginning or the end of the argument
    * container, depending on the from_start argument
    *
    * @param[in] container The container whose members are iterated.
    * @param[in] from_start true if iteration starts at beginning of container.
    */
   explicit BoxContainerIterator(
      BoxContainer& container,
      bool from_start = true);

   /*!
    * @brief Copy constructor.
    *
    * @param[in] other
    */
   BoxContainerIterator(
      const BoxContainerIterator& other);

   /*!
    * @brief Assignment operator.
    *
    * @param[in] rhs
    */
   BoxContainerIterator&
   operator = (
      const BoxContainerIterator& rhs);

   /*!
    * @brief The destructor releases all storage.
    */
   ~BoxContainerIterator();

   /*!
    * @brief Get box corresponding to iterator's position in container.
    *
    * @return A mutable reference to the current Box in the iteration.
    */
   Box&
   operator * () const;

   /*!
    * @brief Get box corresponding to iterator's position in container.
    *
    * @return A mutable reference to the current Box in the iteration.
    */
   Box&
   operator () () const;

   /*!
    * @brief Get pointer to box at iterator's position in container.
    *
    * @return Pointer to the current box.
    */ 
   Box*
   operator -> () const;

   /*!
    * @brief Post-increment iterator to point to next box in the container.
    *
    * @return Iterator at the position in the container before the increment.
    */
   BoxContainerIterator 
   operator ++ (
      int);

   /*!
    * @brief Pre-increment iterator to point to next box in the container.
    *
    * @return Reference to iterator at the position in the container after
    * the increment.
    */
   const BoxContainerIterator&
   operator ++ ();

   /*!
    * @brief Post-decrement iterator to point to next box in the container.
    *
    * @return Iterator at the position in the container before the decrement.
    */
   BoxContainerIterator
   operator -- (
      int);


   /*!
    * @brief Pre-decrement iterator to point to next box in the container.
    *
    * @return Reference to iterator at the position in the container after
    * the decrement.
    */
   const BoxContainerIterator&
   operator -- ();


   /*!
    * @brief Equality operators
    *
    * @return true if both iterators point to the same box.
    *
    * @param[in] other
    */
   bool
   operator == (
      const BoxContainerIterator& other) const;

   bool
   operator == (
      const BoxContainerConstIterator& other) const;

   /*!
    * @brief Inequality operators.
    *
    * @return true if both iterators point to different boxes.
    *
    * @param[in] other
    */
   bool
   operator != (
      const BoxContainerIterator& other) const;

   bool
   operator != (
      const BoxContainerConstIterator& other) const;

private:

   /*!
    * @brief Default constructor is defined but accessible only by friends.
    */
   BoxContainerIterator();

   /*
    * Underlying iterator to be used when unordered.
    */
   std::list<Box>::iterator d_list_iter;

   /*
    * Underlying iterator to be used when ordered.
    */
   std::set<Box*>::iterator d_set_iter;

   bool d_ordered;

};

}
}

#ifdef SAMRAI_INLINE
#include "SAMRAI/hier/BoxContainerIterator.I"
#endif

#endif
