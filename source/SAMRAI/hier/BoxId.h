/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Identifier for a Box.
 *
 ************************************************************************/

#ifndef included_hier_BoxId
#define included_hier_BoxId

#include <iostream>

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/GlobalId.h"
#include "SAMRAI/hier/PeriodicId.h"

namespace SAMRAI {
namespace hier {

/*!
 * @brief Identifier for a Box, consisting of a GlobalId and a PeriodicId
 *
 * Boxes are identified by their GlobalId and PeriodicId.
 * A Box and all its periodic images have the same GlobalId but
 * different PeriodicId.
 *
 * Comparison operators are provided to define sorted ordering of
 * objects.  The GlobalId and PeriodicId are used for all
 * comparisons.  The GlobalIds are compared first, followed by
 * the PeriodicIds.
 */
class BoxId
{

public:
   /*!
    * @brief Default constructor uses the default constructors for the
    * GlobalId and PeriodicId.
    *
    * The object can be changed using initialize() or by assignment.
    */
   BoxId();

   /*!
    * @brief Initializing constructor.
    *
    * @param[in] local_id
    *
    * @param[in] owner_rank
    *
    * @param[in] periodic_id
    */
   explicit BoxId(
      const LocalId& local_id,
      const int owner_rank,
      const PeriodicId& periodic_id = PeriodicId::zero());

   /*!
    * @brief Initializing constructor.
    *
    * @param[in] global_id
    *
    * @param[in] periodic_id
    */
   explicit BoxId(
      const GlobalId& id,
      const PeriodicId& periodic_id = PeriodicId::zero());

   /*!
    * @brief Copy constructor.
    *
    * @param[in] other
    */
   BoxId(
      const BoxId& other);

   /*!
    * @brief Destructor.
    */
   virtual ~BoxId();

   /*!
    * @brief Set all the attributes to given values.
    *
    * @param[in] local_id
    *
    * @param[in] owner_rank
    *
    * @param[in] periodic_id
    */
   void
   initialize(
      const LocalId& local_id,
      const int owner_rank,
      const PeriodicId& periodic_id = PeriodicId::zero());

   /*!
    * @brief Access the GlobalId.
    */
   const GlobalId&
   getGlobalId() const;

   /*!
    * @brief Access the owner rank.
    */
   int
   getOwnerRank() const;

   /*!
    * @brief Access the LocalId.
    */
   const LocalId
   &
   getLocalId() const;

   /*!
    * @brief Access the PeriodicId.
    */
   const PeriodicId&
   getPeriodicId() const;

   /*!
    * @brief Whether the PeriodicId refers to a periodic
    * image.
    */
   bool
   isPeriodicImage() const;

   /*!
    * @brief Whether the BoxId is valid--meaning it has a valid
    * GlobalId and PeriodicId.
    */ 
   bool
   isValid() const;

   //@{

   //! @name Comparison operators

   /*!
    * @brief Equality operator.
    *
    * All comparison operators use the GlobalId and PeriodicId.
    */
   bool
   operator == (
      const BoxId& r) const;

   /*!
    * @brief Inequality operator.
    *
    * See note on comparison for operator==(const BoxId&);
    */
   bool
   operator != (
      const BoxId& r) const;

   /*!
    * @brief Less-than operator.
    *
    * Compare the owner ranks first; if they compare equal, compare the
    * LocalIds next; if they compare equal, compare the PeriodicIds.
    */
   bool
   operator < (
      const BoxId& r) const;

   /*!
    * @brief Greater-than operator.
    *
    * Compare the owner ranks first; if they compare equal, compare the
    * LocalIds next; if they compare equal, compare the PeriodicIds.
    */
   bool
   operator > (
      const BoxId& r) const;

   /*!
    * @brief Less-than-or-equal-to operator.
    */
   bool
   operator <= (
      const BoxId& r) const;

   /*!
    * @brief Greater-than-or-equal-to operator.
    */
   bool
   operator >= (
      const BoxId& r) const;

   //@}

   //@{

   //! @name Support for message passing

   /*!
    * @brief Give number of ints required for putting a BoxId in
    * message passing buffer.
    *
    * @see putToIntBuffer(), getFromIntBuffer().
    */
   static int
   commBufferSize();

   /*!
    * @brief Put self into a int buffer.
    *
    * This is the opposite of getFromIntBuffer().  Number of ints
    * written is given by commBufferSize().
    */
   void
   putToIntBuffer(
      int * buffer) const;

   /*!
    * @brief Set attributes according to data in int buffer.
    *
    * This is the opposite of putToIntBuffer().  Number of ints read
    * is given by commBufferSize().
    */
   void
   getFromIntBuffer(
      const int * buffer);

   //@}

   /*!
    * @brief Format and insert the object into a stream.
    */
   friend std::ostream&
   operator << (
      std::ostream& co,
      const BoxId& r);

private:
   GlobalId d_global_id;

   PeriodicId d_periodic_id;
};

}
}

#ifdef SAMRAI_INLINE
#include "SAMRAI/hier/BoxId.I"
#endif

#endif  // included_hier_BoxId
