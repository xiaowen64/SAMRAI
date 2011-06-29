/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Generic identifier used on a single process.
 *
 ************************************************************************/

#ifndef included_hier_LocalId
#define included_hier_LocalId

#include "SAMRAI/SAMRAI_config.h"

#include <iostream>

namespace SAMRAI {
namespace hier {

/*!
 * @brief Generic identifier for identifying things on the local
 * process.
 *
 * The LocalId can be combined with a process rank, as is done in
 * GlobalId, to create global identifiers.
 *
 * Comparison operators are provided to define the sorted ordering of
 * objects.
 */
class LocalId {

public:

   /*!
    * @brief Default constructor.
    */
   LocalId();

   /*!
    * @brief Copy constructor.
    */
   LocalId(const LocalId &other);

   /*!
    * @brief Construct from a numerical value.
    *
    * This method is explicit to prevent automatic conversion.
    */
   explicit LocalId(const int &value);

   /*!
    * @brief Default constructor.
    */
   ~LocalId();

   /*!
    * @brief Assignment operator.
    *
    * @param[in] rhs
    *
    * @return @c *this
    */
   LocalId &operator = ( const LocalId& rhs);

   /*!
    * @brief Assignment operator.
    *
    * @param[in] rhs
    *
    * @return @c *this
    */
   LocalId &operator = ( const int &rhs);

   /*!
    * @brief Access the numerical value.
    */
   int &getValue();

   /*!
    * @brief Access the numerical value.
    */
   const int &getValue() const;

   /*!
    * @brief Get the LocalId with a numerical value of zero.
    */
   static const LocalId &getZero();

   /*!
    * @brief Get the designated invalid value for this class.
    */
   static const LocalId &getInvalidId();


   //@{

   //! @name Numerical operations.

   /*!
    * @brief Pre-increment iterator.
    *
    * Pre-increment increments the value and returns the incremented
    * state.
    */
   LocalId operator ++ ();

   /*!
    * @brief Post-increment iterator.
    *
    * Post-increment saves the value, increment it and returns an
    * object with the saved value.
    */
   LocalId operator ++ (int);


   /*!
    * @brief Addition.
    *
    * @param[in] rhs
    */
   LocalId operator+(const LocalId &rhs) const;

   /*!
    * @brief Subtraction.
    *
    * @param[in] rhs
    */
   LocalId operator-(const LocalId &rhs) const;

   /*!
    * @brief Multiplication.
    *
    * @param[in] rhs
    */
   LocalId operator*(const LocalId &rhs) const;

   /*!
    * @brief Division.
    *
    * @param[in] rhs
    */
   LocalId operator/(const LocalId &rhs) const;

   /*!
    * @brief Modulus.
    *
    * @param[in] rhs
    */
   LocalId operator%(const LocalId &rhs) const;

   /*!
    * @brief Addition and assignment.
    *
    * @param[in] rhs
    */
   LocalId &operator+=(const LocalId &rhs);

   /*!
    * @brief Subtraction and assignment.
    *
    * @param[in] rhs
    */
   LocalId &operator-=(const LocalId &rhs);


   /*!
    * @brief Integer addition.
    *
    * @param[in] rhs
    */
   LocalId operator+(const int &rhs) const;

   /*!
    * @brief Integer subtraction.
    *
    * @param[in] rhs
    */
   LocalId operator-(const int &rhs) const;

   /*!
    * @brief Integer multiplication.
    *
    * @param[in] rhs
    */
   LocalId operator*(const int &rhs) const;

   /*!
    * @brief Integer division.
    *
    * @param[in] rhs
    */
   LocalId operator/(const int &rhs) const;

   /*!
    * @brief Integer modulus.
    *
    * @param[in] rhs
    */
   LocalId operator%(const int &rhs) const;

   /*!
    * @brief Integer addition and assignment.
    *
    * @param[in] rhs
    */
   LocalId &operator+=(const int &rhs);

   /*!
    * @brief Integer subtraction and assignment.
    *
    * @param[in] rhs
    */
   LocalId &operator-=(const int &rhs);


   //@}


   //@{

   //! @name Comparison with another LocalId.

   /*!
    * @brief Equality operator.
    *
    * All comparison operators compare the numerical value.
    *
    * @param[in] rhs
    */
   bool operator == ( const LocalId& rhs) const;

   /*!
    * @brief Inequality operator.
    *
    * See note on comparison for operator==(const LocalId&);
    *
    * @param[in] rhs
    */
   bool operator != ( const LocalId& rhs) const;

   /*!
    * @brief Less-than operator.
    *
    * See note on comparison for operator==(const LocalId&);
    *
    * @param[in] rhs
    */
   bool operator < ( const LocalId& rhs) const;

   /*!
    * @brief Greater-than operator.
    *
    * See note on comparison for operator==(const LocalId&);
    *
    * @param[in] rhs
    */
   bool operator > ( const LocalId& rhs) const;

   /*!
    * @brief Less-than-or-equal-to operator.
    *
    * See note on comparison for operator==(const LocalId&);
    *
    * @param[in] rhs
    */
   bool operator <= ( const LocalId& rhs) const;

   /*!
    * @brief Greater-thanor-equal-to operator.
    *
    * See note on comparison for operator==(const LocalId&);
    *
    * @param[in] rhs
    */
   bool operator >= ( const LocalId& rhs) const;

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
   bool operator == ( const int &rhs) const;

   /*!
    * @brief Inequality operator.
    *
    * See note on comparison for operator==(const LocalId&);
    *
    * @param[in] rhs
    */
   bool operator != ( const int &rhs) const;

   /*!
    * @brief Less-than operator.
    *
    * See note on comparison for operator==(const LocalId&);
    *
    * @param[in] rhs
    */
   bool operator < ( const int &rhs) const;

   /*!
    * @brief Greater-than operator.
    *
    * See note on comparison for operator==(const LocalId&);
    *
    * @param[in] rhs
    */
   bool operator > ( const int &rhs) const;

   /*!
    * @brief Less-than-or-equal-to operator.
    *
    * See note on comparison for operator==(const LocalId&);
    *
    * @param[in] rhs
    */
   bool operator <= ( const int &rhs) const;

   /*!
    * @brief Greater-thanor-equal-to operator.
    *
    * See note on comparison for operator==(const LocalId&);
    *
    * @param[in] rhs
    */
   bool operator >= ( const int &rhs) const;

   //@}


   /*!
    * @brief Format and insert object into a stream.
    */
   friend std::ostream&
   operator << (
      std::ostream& co,
      const LocalId& r);

private:

   /*!
    * @brief Numerical value of the identifier.
    */
   int d_value;

   /*!
    * @brief LocalId with a numerical value of zero.
    */
   static const LocalId s_zero_id;

   /*!
    * @brief Definition of invalid LocalId.
    */
   static const LocalId s_invalid_id;

};

}
}

#ifdef SAMRAI_INLINE
#include "SAMRAI/hier/LocalId.I"
#endif

#endif  // included_hier_LocalId
