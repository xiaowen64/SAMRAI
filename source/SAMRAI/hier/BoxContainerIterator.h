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

namespace SAMRAI {
namespace hier {

class Box;
class BoxContainer;

/**
 * A mutable iterator over the boxes in a BoxContainer or the underlying boxes
 * in a MappedBoxContainer.
 *
 * @see hier::BoxContainer
 */
class BoxContainerIterator
{
   friend class BoxContainer;

public:
   // Constructors.

   /*!
    * @brief Constructor for the BoxContainerIterator.
    *
    * The iterator will enumerate the boxes in the argument container.
    *
    * @param[in] container The container whose members are iterated.
    * @param[in] from_start true if iteration starts at front of container.
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

   // Destructor.

   /*!
    * @brief The destructor releases all storage.
    */
   ~BoxContainerIterator();

   // Operators.

   /*!
    * @brief Extract box corresponding to iterator's position in container.
    *
    * @return A mutable reference to the current Box in the iteration.
    */
   Box&
   operator * () const;

   /*!
    * @brief Extract box corresponding to iterator's position in container.
    *
    * @return A mutable reference to the current Box in the iteration.
    */
   Box&
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
   void
   operator ++ (
      int);

   /*!
    * @brief Pre-increment iterator to point to next box in the container.
    */
   const BoxContainerIterator&
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
      const BoxContainerIterator& other) const;

   /*!
    * @brief Determine if two iterators are not equivalent.
    *
    * @return true if both iterators point to different boxes.
    *
    * @param[in] other
    */
   bool
   operator != (
      const BoxContainerIterator& other) const;

private:
   /*
    * Default constructor just to be clear that there is none.
    */
   BoxContainerIterator();

   /*
    * The BoxContainer being iterated over.
    */
   BoxContainer* d_box_container;

   /*
    * Underlying iterator for a BoxContainer.  This is a wrapper.
    */
   std::list<Box>::iterator d_list_iter;

   // Add pointer to MappedBoxContainer and underlying iterator for
   // MappedBoxContainer here.
};

}
}

#ifdef SAMRAI_INLINE
#include "SAMRAI/hier/BoxContainerIterator.I"
#endif

#endif
