/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Iterator over real MappedBoxes in a MappedBoxSet. 
 *
 ************************************************************************/
#ifndef included_hier_RealMappedBoxConstIterator
#define included_hier_RealMappedBoxConstIterator

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/MappedBoxSet.h"

namespace SAMRAI {
namespace hier {

/*
 * TODO: Do we really need a separate class for this?  Couldn't we just
 *       add an argument to the MappedBoxSet iterator construction that
 *       (e.g., an enum with values: All (default value), RealOnly, 
 *       PeriodicImagesOnly, etc.  and extend for Multiblock stuff)?
 *       Then, which boxes in the set are selected for iteration would be 
 *       controlled internally.
 */
/*!
 * @brief Iterator through real MappedBoxes (not periodic images) in a
 * const MappedBoxSet.
 *
 * RealMappedBoxConstIterator is an iterator that provides methods for
 * stepping through a MappedBoxSet, skipping periodic images.
 *
 * Example usage:
 * @verbatim
 *  MappedBoxSet mapped_boxes;
 *  // fill in mapped_boxes
 *  for ( RealMappedBoxConstIterator ni(mapped_boxes); ni.isValid(); ++ni ) {
 *    TBOX_ASSERT( ! ni->isPeriodicImage() );
 *  }
 * @endverbatim
 */
class RealMappedBoxConstIterator
{

public:

   /*!
    * @brief Construct the iterator for the given MappedBoxSet.
    *
    * The iterator will iterate through the items in mapped_boxes.
    *
    * @param[in] mapped_boxes
    */
   explicit RealMappedBoxConstIterator(
      const MappedBoxSet& mapped_boxes);

   /*!
    * @brief Destructor.
    */
   virtual ~RealMappedBoxConstIterator(
      void);

   /*!
    * @brief Assignment operator.
    */
   RealMappedBoxConstIterator&
   operator = (
      const RealMappedBoxConstIterator& r);

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
      const RealMappedBoxConstIterator& r) const;

   /*!
    * @brief Inequality comparison.
    */
   bool
   operator != (
      const RealMappedBoxConstIterator& r) const;

   /*!
    * @brief Pre-increment iterator.
    *
    * Pre-increment increment the iterator and returns the incremented
    * state.
    */
   RealMappedBoxConstIterator
   &operator ++ ();

   /*!
    * @brief Post-increment iterator.
    *
    * Post-increment saves the iterator, increment it and returns the
    * saved iterator.
    */
   RealMappedBoxConstIterator
   operator ++ (
      int);

   /*!
    * @brief Whether the iterator can be dereferenced.  When the
    * iterator reaches its end, this returns false.
    */
   bool
   isValid() const;

private:
   /*!
    * @brief MappedBoxSet being iterated through.
    */
   const MappedBoxSet* d_mapped_boxes;

   /*!
    * @brief The iterator.
    */
   MappedBoxSet::const_iterator d_ni;

};

}
}

#endif  // included_hier_RealMappedBoxConstIterator
