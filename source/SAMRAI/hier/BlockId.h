/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Block identifier in multiblock domain.
 *
 ************************************************************************/

#ifndef included_hier_BlockId
#define included_hier_BlockId

#include "SAMRAI/SAMRAI_config.h"

#include <iostream>

namespace SAMRAI {
namespace hier {

/*!
 * @brief Generic identifier for identifying the block id.
 *
 * Comparison operators are provided to define sorted ordering of
 * objects.
 */
class BlockId
{

public:
   /*!
    * @brief Default constructor sets the value to invalid.
    */
   BlockId();

   /*!
    * @brief Copy constructor.
    */
   BlockId(
      const BlockId& other);

   /*!
    * @brief Construct from a numerical value.
    *
    * This method is explicit to prevent automatic conversion.
    */
   explicit BlockId(
      const int& value);

   /*!
    * @brief Default constructor.
    */
   ~BlockId();

   /*!
    * @brief Assignment operator.
    *
    * @param[in] rhs
    *
    * @return @c *this
    */
   BlockId&
   operator = (
      const BlockId& rhs);

   /*!
    * @brief Set to an int value.
    *
    * @param[in] rhs
    */
   void
   setId(
      const int& rhs);

   /*!
    * @brief Whether the value is valid.
    */
   bool
   isValid() const;

   /*!
    * @brief Access the numerical value.
    */
   const int&
   getBlockValue() const;

   /*!
    * @brief Get the BlockId with a numerical value of zero.
    */
   static const BlockId&
   zero();

   /*!
    * @brief Get the designated invalid value for this class.
    */
   static const BlockId&
   invalidId();

   //@{

   //! @name Comparison with another BlockId.

   /*!
    * @brief Equality operator.
    *
    * All comparison operators compare the numerical value.
    *
    * @param[in] rhs
    */
   bool
   operator == (
      const BlockId& rhs) const;

   /*!
    * @brief Inequality operator.
    *
    * See note on comparison for operator==(const BlockId&);
    *
    * @param[in] rhs
    */
   bool
   operator != (
      const BlockId& rhs) const;

   /*!
    * @brief Less-than operator.
    *
    * See note on comparison for operator==(const BlockId&);
    *
    * @param[in] rhs
    */
   bool
   operator < (
      const BlockId& rhs) const;

   /*!
    * @brief Greater-than operator.
    *
    * See note on comparison for operator==(const BlockId&);
    *
    * @param[in] rhs
    */
   bool
   operator > (
      const BlockId& rhs) const;

   /*!
    * @brief Less-than-or-equal-to operator.
    *
    * See note on comparison for operator==(const BlockId&);
    *
    * @param[in] rhs
    */
   bool
   operator <= (
      const BlockId& rhs) const;

   /*!
    * @brief Greater-thanor-equal-to operator.
    *
    * See note on comparison for operator==(const BlockId&);
    *
    * @param[in] rhs
    */
   bool
   operator >= (
      const BlockId& rhs) const;

   //@}

   //@{

   //! @name Comparison with an integer.

   /*!
    * @brief Equality operator.
    *
    * All comparison operators compare the numerical value.
    *
    * @param[in] rhs
    */
   bool
   operator == (
      const int& rhs) const;

   /*!
    * @brief Inequality operator.
    *
    * See note on comparison for operator==(const BlockId&);
    *
    * @param[in] rhs
    */
   bool
   operator != (
      const int& rhs) const;

   /*!
    * @brief Less-than operator.
    *
    * See note on comparison for operator==(const BlockId&);
    *
    * @param[in] rhs
    */
   bool
   operator < (
      const int& rhs) const;

   /*!
    * @brief Greater-than operator.
    *
    * See note on comparison for operator==(const BlockId&);
    *
    * @param[in] rhs
    */
   bool
   operator > (
      const int& rhs) const;

   /*!
    * @brief Less-than-or-equal-to operator.
    *
    * See note on comparison for operator==(const BlockId&);
    *
    * @param[in] rhs
    */
   bool
   operator <= (
      const int& rhs) const;

   /*!
    * @brief Greater-thanor-equal-to operator.
    *
    * See note on comparison for operator==(const BlockId&);
    *
    * @param[in] rhs
    */
   bool
   operator >= (
      const int& rhs) const;

   //@}

   /*!
    * @brief Format and insert object into a stream.
    */
   friend std::ostream&
   operator << (
      std::ostream& co,
      const BlockId& r);

private:
   /*!
    * @brief Numerical value of the identifier.
    */
   int d_value;

   /*!
    * @brief BlockId with a numerical value of zero.
    */
   static const BlockId s_zero_id;

   /*!
    * @brief Definition of invalid BlockId.
    */
   static const BlockId s_invalid_id;

};

}
}

#ifdef SAMRAI_INLINE
#include "SAMRAI/hier/BlockId.I"
#endif

#endif  // included_hier_BlockId
