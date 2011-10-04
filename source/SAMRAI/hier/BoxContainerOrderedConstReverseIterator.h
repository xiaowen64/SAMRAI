/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   An iterator over containers of boxes.
 *
 ************************************************************************/

#ifndef included_hier_BoxContainerOrderedConstReverseIterator
#define included_hier_BoxContainerOrderedConstReverseIterator


#include <set>

namespace SAMRAI {
namespace hier {


class Box;
class BoxContainer;

/**
 * A mutable iterator over the boxes in a BoxContainer or the underlying boxes
 * in a BoxContainer.
 *
 * @see hier::BoxContainer
 */
class BoxContainerOrderedConstReverseIterator
{
   friend class BoxContainer;

public:
   // Constructors.

   /*!
    * @brief Constructor for the BoxContainerOrderedConstReverseIterator.
    *
    * The iterator will enumerate the boxes in the argument container.
    *
    * @param[in] container The container whose members are iterated.
    * @param[in] from_back true if iteration starts at back of container.
    */
   explicit BoxContainerOrderedConstReverseIterator(
      const BoxContainer& container,
      bool from_back = true);

   /*!
    * @brief Copy constructor.
    *
    * @param[in] other
    */
   BoxContainerOrderedConstReverseIterator(
      const BoxContainerOrderedConstReverseIterator& other);

   /*!
    * @brief Assignment operator.
    *
    * @param[in] rhs
    */
   BoxContainerOrderedConstReverseIterator&
   operator = (
      const BoxContainerOrderedConstReverseIterator& rhs);

   // Destructor.

   /*!
    * @brief The destructor releases all storage.
    */
   ~BoxContainerOrderedConstReverseIterator();

   // Operators.

   /*!
    * @brief Extract box corresponding to iterator's position in container.
    *
    * @return A mutable reference to the current Box in the iteration.
    */
   const Box&
   operator * () const;

   const Box*
   operator -> () const;

   /*!
    * @brief Extract box corresponding to iterator's position in container.
    *
    * @return A mutable reference to the current Box in the iteration.
    */
   const Box&
   operator () () const;

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
    * @brief Post-increment iterator to point to previous box in the container.
    */
   BoxContainerOrderedConstReverseIterator
   operator ++ (
      int);

   /*!
    * @brief Pre-increment iterator to point to previous box in the container.
    */
   const BoxContainerOrderedConstReverseIterator&
   operator ++ ();

   /*!
    * @brief Determine if two iterators are equivalent.
    *
    * @return true if both iterators point to the same box.
    *
    * @param[in] other
    */
   bool
   operator == (
      const BoxContainerOrderedConstReverseIterator& other) const;

   /*!
    * @brief Determine if two iterators are not equivalent.
    *
    * @return true if both iterators point to different boxes.
    *
    * @param[in] other
    */
   bool
   operator != (
      const BoxContainerOrderedConstReverseIterator& other) const;

private:
   /*
    * Default constructor just to be clear that there is none.
    */
   BoxContainerOrderedConstReverseIterator();

   /*
    * The BoxContainer being iterated over.
    */
   const BoxContainer* d_box_container;

   /*
    * Underlying iterator for a BoxContainer.  This is a wrapper.
    */
   std::set<const Box*>::const_reverse_iterator d_set_iter;

};

}
}

#ifdef SAMRAI_INLINE
#include "SAMRAI/hier/BoxContainerOrderedConstReverseIterator.I"
#endif

#endif
