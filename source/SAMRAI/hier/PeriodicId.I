/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Periodic shift identifier in periodic domain.
 *
 ************************************************************************/

namespace SAMRAI {
namespace hier {


/*
********************************************************************************
********************************************************************************
*/
SAMRAI_INLINE_KEYWORD
PeriodicId::PeriodicId()
   : d_value( invalidId().d_value ) {
   return;
}


/*
********************************************************************************
********************************************************************************
*/
SAMRAI_INLINE_KEYWORD
PeriodicId::PeriodicId(const PeriodicId &other)
   : d_value( other.d_value ) {
   return;
}


/*
********************************************************************************
********************************************************************************
*/
SAMRAI_INLINE_KEYWORD
PeriodicId::PeriodicId(const int &value)
   : d_value( value ) {
   return;
}


/*
********************************************************************************
********************************************************************************
*/
SAMRAI_INLINE_KEYWORD
PeriodicId::~PeriodicId()
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
PeriodicId &PeriodicId::operator = ( const PeriodicId &rhs)
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
PeriodicId &PeriodicId::operator = ( const int &rhs)
{
   d_value = rhs;
   return *this;
}


/*
********************************************************************************
********************************************************************************
*/
SAMRAI_INLINE_KEYWORD
const int &PeriodicId::getPeriodicValue() const
{
   return d_value;
}


/*
********************************************************************************
* Equality operator.
********************************************************************************
*/
SAMRAI_INLINE_KEYWORD
bool PeriodicId::operator == ( const PeriodicId &rhs) const
{
   return d_value == rhs.d_value;
}

/*
********************************************************************************
* Inequality operator.
********************************************************************************
*/
SAMRAI_INLINE_KEYWORD
bool PeriodicId::operator != ( const PeriodicId &rhs) const
{
   return d_value != rhs.d_value;
}

/*
********************************************************************************
* Less-than operator.
********************************************************************************
*/
SAMRAI_INLINE_KEYWORD
bool PeriodicId::operator < ( const PeriodicId &rhs) const
{
   return d_value < rhs.d_value;
}

/*
********************************************************************************
* Greater-than operator.
********************************************************************************
*/
SAMRAI_INLINE_KEYWORD
bool PeriodicId::operator > ( const PeriodicId &rhs) const
{
   return d_value > rhs.d_value;
}

/*
********************************************************************************
* Less-than-or-equal-to operator.
********************************************************************************
*/
SAMRAI_INLINE_KEYWORD
bool PeriodicId::operator <= ( const PeriodicId &rhs) const
{
   return d_value <= rhs.d_value;
}

/*
********************************************************************************
* Greater-thanor-equal-to operator.
********************************************************************************
*/
SAMRAI_INLINE_KEYWORD
bool PeriodicId::operator >= ( const PeriodicId &rhs) const
{
   return d_value >= rhs.d_value;
}

/*
********************************************************************************
********************************************************************************
*/
SAMRAI_INLINE_KEYWORD
const PeriodicId &PeriodicId::zero()
{
   return s_zero_id;
}

/*
********************************************************************************
********************************************************************************
*/
SAMRAI_INLINE_KEYWORD
const PeriodicId &PeriodicId::invalidId()
{
   return s_invalid_id;
}

/*
********************************************************************************
********************************************************************************
*/
SAMRAI_INLINE_KEYWORD
bool PeriodicId::isValid() const
{
   return d_value >= 0;
}

}
}
