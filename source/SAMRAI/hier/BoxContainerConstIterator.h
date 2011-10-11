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

/**
 * An immutable iterator over the boxes in a BoxContainer or the underlying
 * boxes in a BoxContainer.
 *
 * @see hier::BoxContainer
 */
class BoxContainerConstIterator
{
   friend class BoxContainer;
   friend class BoxContainerIterator;

public:
   // Constructors.

   /*!
    * @brief Constructor for the BoxContainerConstIterator.
    *
    * The iterator will enumerate the boxes in the argument container.
    *
    * @param[in] container The container whose members are iterated.
    * @param[in] from_start true if iteration starts at front of container.
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

   // Destructor.

   /*!
    * @brief The destructor releases all storage.
    */
   ~BoxContainerConstIterator();

   // Operators.

   /*!
    * @brief Extract box corresponding to iterator's position in container.
    *
    * @return An immutable reference to the current Box in the iteration.
    */
   const Box&
   operator * () const;

   /*!
    * @brief Extract box corresponding to iterator's position in container.
    *
    * @return An immutable reference to the current Box in the iteration.
    */
   const Box&
   operator () () const;

   const Box*
   operator -> () const;
   /*!
    * @brief Determine if iterator points to a valid position in container.
    *
    * @return true if iterator points to a valid position in container.
    */
//   operator bool () const;

   /*!
    * @brief Determine if iterator points to an invalid position in container.
    *
    * @return true if iterator points to an invalid position in container.
    */
//   bool
//   operator ! () const;

   /*!
    * @brief Post-increment iterator to point to next box in the container.
    */
   BoxContainerConstIterator
   operator ++ (
      int);

   /*!
    * @brief Pre-increment iterator to point to next box in the container.
    */
   const BoxContainerConstIterator&
   operator ++ ();


   BoxContainerConstIterator
   operator -- (
      int);

   const BoxContainerConstIterator&
   operator -- ();
   /*!
    * @brief Determine if two iterators are equivalent.
    *
    * @return true if both iterators point to the same box.
    *
    * @param[in] other
    */

   bool
   operator == (
      const BoxContainerConstIterator& other) const;

   /*!
    * @brief Determine if two iterators are not equivalent.
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
    * Underlying iterator for a BoxContainer.  This is a wrapper.
    */
   std::list<Box>::const_iterator d_list_iter;
   std::set<Box*>::const_iterator d_set_iter;

   bool d_ordered;
};

}
}

#ifdef SAMRAI_INLINE
#include "SAMRAI/hier/BoxContainerConstIterator.I"
#endif

#endif
