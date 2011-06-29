/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Generic identifier used on a single process.
 *
 ************************************************************************/

#include "SAMRAI/tbox/MathUtilities.h"

namespace SAMRAI {
namespace hier {


/*
********************************************************************************
********************************************************************************
*/
SAMRAI_INLINE_KEYWORD
LocalId::LocalId()
   : d_value( getInvalidId().d_value ) {
   return;
}


/*
********************************************************************************
********************************************************************************
*/
SAMRAI_INLINE_KEYWORD
LocalId::LocalId(const LocalId &other)
   : d_value( other.d_value ) {
   return;
}


/*
********************************************************************************
********************************************************************************
*/
SAMRAI_INLINE_KEYWORD
LocalId::LocalId(const int &value)
   : d_value( value ) {
   return;
}


/*
********************************************************************************
********************************************************************************
*/
SAMRAI_INLINE_KEYWORD
LocalId::~LocalId()
{
#ifdef DEBUG_CHECK_ASSERTIONS
   d_value = s_invalid_id.d_value;
#endif
   return;
}


/*
********************************************************************************
* Assignment operator.
********************************************************************************
*/
SAMRAI_INLINE_KEYWORD
LocalId &LocalId::operator = ( const LocalId &rhs)
{
   d_value = rhs.d_value;
   return *this;
}

/*
********************************************************************************
* Assignment operator.
********************************************************************************
*/
SAMRAI_INLINE_KEYWORD
LocalId &LocalId::operator = ( const int &rhs)
{
   d_value = rhs;
   return *this;
}


/*
********************************************************************************
********************************************************************************
*/
SAMRAI_INLINE_KEYWORD
int &LocalId::getValue()
{
   return d_value;
}


/*
********************************************************************************
********************************************************************************
*/
SAMRAI_INLINE_KEYWORD
const int &LocalId::getValue() const
{
   return d_value;
}




/*
****************************************************************************
* Pre-increment operator.
****************************************************************************
 */

SAMRAI_INLINE_KEYWORD
LocalId LocalId::operator ++ ()
{
   ++d_value;
   return *this;
}

/*
****************************************************************************
* Post-increment operator.
****************************************************************************
*/

SAMRAI_INLINE_KEYWORD
LocalId LocalId::operator ++ (int)
{
   int saved = d_value;
   ++d_value;
   return LocalId(saved);
}




/*
****************************************************************************
* Addition
****************************************************************************
 */

SAMRAI_INLINE_KEYWORD
LocalId LocalId::operator + (const LocalId &rhs) const
{
   return LocalId(d_value + rhs.d_value);
}

/*
****************************************************************************
* Subtraction
****************************************************************************
*/

SAMRAI_INLINE_KEYWORD
LocalId LocalId::operator - (const LocalId &rhs) const
{
   return LocalId(d_value - rhs.d_value);
}

/*
****************************************************************************
* Multiplication
****************************************************************************
 */

SAMRAI_INLINE_KEYWORD
LocalId LocalId::operator * (const LocalId &rhs) const
{
   return LocalId(d_value * rhs.d_value);
}

/*
****************************************************************************
* Division
****************************************************************************
*/

SAMRAI_INLINE_KEYWORD
LocalId LocalId::operator / (const LocalId &rhs) const
{
   return LocalId(d_value / rhs.d_value);
}

/*
****************************************************************************
* Modulus
****************************************************************************
*/

SAMRAI_INLINE_KEYWORD
LocalId LocalId::operator % (const LocalId &rhs) const
{
   return LocalId(d_value % rhs.d_value);
}


/*
********************************************************************************
* Addition and assignment.
********************************************************************************
*/
SAMRAI_INLINE_KEYWORD
LocalId &LocalId::operator += ( const LocalId &rhs)
{
   d_value += rhs.d_value;
   return *this;
}

/*
********************************************************************************
* Subtraction and assignment.
********************************************************************************
*/
SAMRAI_INLINE_KEYWORD
LocalId &LocalId::operator -= ( const LocalId &rhs)
{
   d_value -= rhs.d_value;
   return *this;
}




/*
****************************************************************************
* Integer addition
****************************************************************************
 */

SAMRAI_INLINE_KEYWORD
LocalId LocalId::operator + (const int &rhs) const
{
   return LocalId(d_value + rhs);
}

/*
****************************************************************************
* Integer subtraction
****************************************************************************
*/

SAMRAI_INLINE_KEYWORD
LocalId LocalId::operator - (const int &rhs) const
{
   return LocalId(d_value - rhs);
}

/*
****************************************************************************
* Integer multiplication
****************************************************************************
 */

SAMRAI_INLINE_KEYWORD
LocalId LocalId::operator * (const int &rhs) const
{
   return LocalId(d_value * rhs);
}

/*
****************************************************************************
* Integer division
****************************************************************************
*/

SAMRAI_INLINE_KEYWORD
LocalId LocalId::operator / (const int &rhs) const
{
   return LocalId(d_value / rhs);
}

/*
****************************************************************************
* Integer modulus
****************************************************************************
*/

SAMRAI_INLINE_KEYWORD
LocalId LocalId::operator % (const int &rhs) const
{
   return LocalId(d_value % rhs);
}


/*
********************************************************************************
* Integer addition and assignment.
********************************************************************************
*/
SAMRAI_INLINE_KEYWORD
LocalId &LocalId::operator += ( const int &rhs)
{
   d_value += rhs;
   return *this;
}

/*
********************************************************************************
* Integer subtraction and assignment.
********************************************************************************
*/
SAMRAI_INLINE_KEYWORD
LocalId &LocalId::operator -= ( const int &rhs)
{
   d_value -= rhs;
   return *this;
}





/*
********************************************************************************
* Equality operator.
********************************************************************************
*/
SAMRAI_INLINE_KEYWORD
bool LocalId::operator == ( const LocalId &rhs) const
{
   return d_value == rhs.d_value;
}

/*
********************************************************************************
* Equality operator.
********************************************************************************
*/
SAMRAI_INLINE_KEYWORD
bool LocalId::operator == ( const int &rhs) const
{
   return d_value == rhs;
}

/*
********************************************************************************
* Inequality operator.
********************************************************************************
*/
SAMRAI_INLINE_KEYWORD
bool LocalId::operator != ( const LocalId &rhs) const
{
   return d_value != rhs.d_value;
}

/*
********************************************************************************
* Inequality operator.
********************************************************************************
*/
SAMRAI_INLINE_KEYWORD
bool LocalId::operator != ( const int &rhs) const
{
   return d_value != rhs;
}

/*
********************************************************************************
* Less-than operator.
********************************************************************************
*/
SAMRAI_INLINE_KEYWORD
bool LocalId::operator < ( const LocalId &rhs) const
{
   return d_value < rhs.d_value;
}

/*
********************************************************************************
* Less-than operator.
********************************************************************************
*/
SAMRAI_INLINE_KEYWORD
bool LocalId::operator < ( const int &rhs) const
{
   return d_value < rhs;
}

/*
********************************************************************************
* Greater-than operator.
********************************************************************************
*/
SAMRAI_INLINE_KEYWORD
bool LocalId::operator > ( const LocalId &rhs) const
{
   return d_value > rhs.d_value;
}

/*
********************************************************************************
* Greater-than operator.
********************************************************************************
*/
SAMRAI_INLINE_KEYWORD
bool LocalId::operator > ( const int &rhs) const
{
   return d_value > rhs;
}

/*
********************************************************************************
* Less-than-or-equal-to operator.
********************************************************************************
*/
SAMRAI_INLINE_KEYWORD
bool LocalId::operator <= ( const LocalId &rhs) const
{
   return d_value <= rhs.d_value;
}

/*
********************************************************************************
* Less-than-or-equal-to operator.
********************************************************************************
*/
SAMRAI_INLINE_KEYWORD
bool LocalId::operator <= ( const int &rhs) const
{
   return d_value <= rhs;
}

/*
********************************************************************************
* Greater-thanor-equal-to operator.
********************************************************************************
*/
SAMRAI_INLINE_KEYWORD
bool LocalId::operator >= ( const LocalId &rhs) const
{
   return d_value >= rhs.d_value;
}

/*
********************************************************************************
* Greater-thanor-equal-to operator.
********************************************************************************
*/
SAMRAI_INLINE_KEYWORD
bool LocalId::operator >= ( const int &rhs) const
{
   return d_value >= rhs;
}



/*
********************************************************************************
********************************************************************************
*/
SAMRAI_INLINE_KEYWORD
const LocalId &LocalId::getZero()
{
   return s_zero_id;
}

/*
********************************************************************************
********************************************************************************
*/
SAMRAI_INLINE_KEYWORD
const LocalId &LocalId::getInvalidId()
{
   return s_invalid_id;
}

}
}
