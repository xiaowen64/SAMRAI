/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   An iterator over containers of boxes.
 *
 ************************************************************************/

#ifndef included_hier_BoxContainerSetIterator
#define included_hier_BoxContainerSetIterator

#include <set>

namespace SAMRAI {
namespace hier {

class Box;
class BoxContainer;
class BoxContainerSetConstIterator;

/**
 * A mutable iterator over the boxes in a BoxContainer or the underlying boxes
 * in a BoxContainer.
 *
 * @see hier::BoxContainer
 */
class BoxContainerSetIterator
{
   friend class BoxContainer;
   friend class BoxContainerSetConstIterator;

public:
   // Constructors.

   /*!
    * @brief Constructor for the BoxContainerSetIterator.
    *
    * The iterator will enumerate the boxes in the argument container.
    *
    * @param[in] container The container whose members are iterated.
    * @param[in] from_start true if iteration starts at front of container.
    */
   explicit BoxContainerSetIterator(
      BoxContainer& container,
      bool from_start = true);

   explicit BoxContainerSetIterator(
      const BoxContainer& container,
      bool from_start = true);

   /*!
    * @brief Copy constructor.
    *
    * @param[in] other
    */
   BoxContainerSetIterator(
      const BoxContainerSetIterator& other);

   /*!
    * @brief Assignment operator.
    *
    * @param[in] rhs
    */
   BoxContainerSetIterator&
   operator = (
      const BoxContainerSetIterator& rhs);

   // Destructor.

   /*!
    * @brief The destructor releases all storage.
    */
   ~BoxContainerSetIterator();

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
   operator bool () const;

   /*!
    * @brief Determine if iterator points to an invalid position in container.
    *
    * @return true if iterator points to an invalid position in container.
    */
   bool
   operator ! () const;

   /*!
    * @brief Post-increment iterator to point to next box in the container.
    */
   BoxContainerSetIterator
   operator ++ (
      int);

   /*!
    * @brief Pre-increment iterator to point to next box in the container.
    */
   const BoxContainerSetIterator&
   operator ++ ();

   /*!
    * @brief Post-decrement iterator to point to previous box in the container.
    */
   BoxContainerSetIterator
   operator -- (
      int);

   /*!
    * @brief Pre-decrement iterator to point to previous box in the container.
    */
   const BoxContainerSetIterator&
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
      const BoxContainerSetIterator& other) const;

   /*!
    * @brief Determine if two iterators are not equivalent.
    *
    * @return true if both iterators point to different boxes.
    *
    * @param[in] other
    */
   bool
   operator != (
      const BoxContainerSetIterator& other) const;

private:
   /*
    * Default constructor just to be clear that there is none.
    */
   BoxContainerSetIterator();

   /*
    * The BoxContainer being iterated over.
    */
   const BoxContainer* d_box_container;

   /*
    * Underlying iterator for a BoxContainer.  This is a wrapper.
    */
   std::set<const Box*>::iterator d_set_iter;

};

}
}

#ifdef SAMRAI_INLINE
#include "SAMRAI/hier/BoxContainerSetIterator.I"
#endif

#endif
