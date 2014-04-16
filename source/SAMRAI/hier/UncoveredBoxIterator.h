/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2014 Lawrence Livermore National Security, LLC
 * Description:   A container of boxes with basic domain calculus operations
 *
 ************************************************************************/

#ifndef included_hier_UncoveredBoxIterator
#define included_hier_UncoveredBoxIterator

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/BoxContainer.h"

#include "boost/shared_ptr.hpp"
#include <utility>
#include <vector>

namespace SAMRAI {
namespace hier {

class Patch;
class PatchHierarchy;

/*!
 * @brief An iterator over the uncovered Boxes in a hierarhcy.
 *
 * This iterator points to a pair consisting of the current uncovered Box and
 * the Patch with which it is associated.  Note that in the case of overlapping
 * Patches in which the overlap is uncovered, the overlapping region will
 * appear multiple times in the iteration.  The number of appearances is equal
 * to the number of Patches that overlap that region.
 */
class UncoveredBoxIterator
{
friend class PatchHierarchy;

public:
   /*!
    * @brief Copy constructor.
    *
    * @param[in] other The iterator being copied.
    */
   UncoveredBoxIterator(
      const UncoveredBoxIterator& other);

   /*!
    * Destructor.
    */
   ~UncoveredBoxIterator();

   /*!
    * @brief Assignment operator.
    *
    * @param[in] rhs The right hand side of the assignment.
    */
   UncoveredBoxIterator&
   operator = (
      const UncoveredBoxIterator& rhs);

   /*!
    * @brief Dereference operator mimicking a pointer dereference.
    */
   const std::pair<boost::shared_ptr<Patch>, Box>&
   operator * () const;

   /*!
    * @brief Dereference operator mimicking a pointer dereference.
    */
   const std::pair<boost::shared_ptr<Patch>, Box>*
   operator -> () const;

   /*!
    * @brief Equality comparison.
    *
    * @param[in] rhs The right hand side of the comparison.
    */
   bool
   operator == (
      const UncoveredBoxIterator& rhs) const;

   /*!
    * @brief Inequality comparison.
    *
    * @param[in] rhs The right hand side of the comparison.
    */
   bool
   operator != (
      const UncoveredBoxIterator& rhs) const;

   /*!
    * @brief Pre-increment iterator.
    *
    * Pre-increment increment the iterator and returns the incremented
    * state.
    */
   UncoveredBoxIterator&
   operator ++ ();

   /*!
    * @brief Post-increment iterator.
    *
    * Post-increment saves the iterator, increments it and returns the
    * saved iterator.
    */
   UncoveredBoxIterator
   operator ++ (
      int);

private:
   /*!
    * @brief Unimplemented default constructor.
    */
   UncoveredBoxIterator();

   /*!
    * @brief Constructor.
    *
    * @param[in] hierarchy The hierarchy over which the iteration will occur.
    * @param[in] begin If true the iteration starts from the beginning.
    */
   UncoveredBoxIterator(
      const PatchHierarchy* hierarchy,
      bool begin);

   /*!
    * @brief Private method to do work common to both pre and post increments.
    */
   void
   incrementIterator();

   /*!
    * @brief Private method to find next finest uncovered boxes.
    */
   void
   findNextFinestUncoveredBoxes();

   /*!
    * @brief Find the patch overlapped by the current uncovered box at set
    * d_item appropriatly.
    */
   void
   findOverlappedPatch();

   /* The PatchHierarchy on which this iterator operates. */
   const PatchHierarchy* d_hierarchy;

   /* The current level in the PatchHierarchy. */
   int d_level_num;

   /* The uncovered Boxes for the current level. */
   BoxContainer d_uncovered_boxes;

   /* The iterator over d_uncovered_boxes. */
   BoxContainer::const_iterator d_uncovered_boxes_itr;

   /* The iterator at the end of d_uncovered_boxes. */
   BoxContainer::const_iterator d_uncovered_boxes_itr_end;

   /* All Boxes uncovered or not in the current level of the PatchHierarchy. */
   BoxContainer d_level_boxes;

   /* The current item in the iteration. */
   std::pair<boost::shared_ptr<Patch>, Box>* d_item;

   /* The number of the finest level in the hierarchy. */
   int d_finest_level_num;

   /* The level boxes that the current uncovered box originally came from.
    * Since overlapping patches are allowed, a given uncovered box may be
    * associated with multiple patches and may therefore have multiple
    * overlapping level boxes.  Much depends on how the original level boxes
    * are cut up to form the uncovered boxes.
    */
   BoxContainer d_overlapping_level_boxes;

   /* The originating level box corresponding to the patch that the iterator is
    * currently pointing to.
    */
   BoxContainer::const_iterator d_cur_overlapping_level_box;

   /* The iterator at the end of d_overlapping_level_boxes. */
   BoxContainer::const_iterator d_end_overlapping_level_boxes;
};

}
}

#endif
